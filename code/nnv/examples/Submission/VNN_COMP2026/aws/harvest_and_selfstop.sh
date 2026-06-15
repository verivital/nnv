#!/bin/bash
# Runs ON the AWS box, launched right after overnight_sweep.sh. Waits for the
# sweep to finish, tars all results onto the (persistent) root EBS, writes a
# DONE sentinel, then self-stops the instance (shutdown behavior = stop, so this
# STOPS not terminates -- EBS + license-anchor ENI both survive). No IAM and no
# local machine needed: cost-saving is guaranteed once the sweep ends.
set -u
LOG=~/sweep_logs/harvest.log
exec >>"$LOG" 2>&1
echo "[harvest] started $(date -u +%FT%TZ)"

# --- Phase 1: confirm the sweep actually started (up to 25 min for git fetch +
#     manifest regen + MATLAB cold start). If MATLAB never appears, do NOT
#     self-stop -- leave the box up for inspection.
seen=0
for i in $(seq 1 150); do            # 150 * 10s = 25 min
  n=$(pgrep -c MATLAB 2>/dev/null || echo 0)
  if [ "$n" -gt 0 ]; then seen=1; echo "[harvest] sweep up: $n MATLAB procs at $(date -u +%T)"; break; fi
  sleep 10
done
if [ "$seen" -eq 0 ]; then
  echo "[harvest] ABORT: no MATLAB proc within 25 min -- leaving box UP for inspection"
  echo "harvest-abort-no-matlab $(date -u +%FT%TZ)" > ~/SWEEP_HARVEST_ABORTED
  exit 1
fi

# --- Phase 2: wait for completion. Done = MATLAB procs == 0 for 3 consecutive
#     checks (guards against a brief gap between workers). Hard cap 14h.
zero=0
for i in $(seq 1 5040); do           # 5040 * 10s = 14h
  n=$(pgrep -c MATLAB 2>/dev/null || echo 0)
  if [ "$n" -eq 0 ]; then zero=$((zero+1)); else zero=0; fi
  [ "$zero" -ge 3 ] && { echo "[harvest] all workers exited at $(date -u +%T)"; break; }
  sleep 10
done

# --- Phase 3: harvest. Collect every results CSV the sweep produced plus all
#     worker logs into a single timestamped tarball on the root EBS.
TS=$(date -u +%Y%m%dT%H%M%SZ)
OUT=~/sweep_results_$TS
mkdir -p "$OUT/csv"
# run_all_benchmarks writes results_*.csv in its cwd; collect from the whole tree
find ~/nnv -name 'results_*.csv' -newermt '-20 hours' -exec cp -t "$OUT/csv" {} + 2>/dev/null
cp ~/sweep_logs/*.log "$OUT/" 2>/dev/null
# quick at-a-glance summary
{
  echo "sweep finished $(date -u +%FT%TZ)"
  echo "commit: $(cd ~/nnv && git log --oneline -1 2>/dev/null)"
  echo "DONE markers: $(grep -lhE 'JOB-.*-DONE|GROUP-.*-DONE' ~/sweep_logs/*.log 2>/dev/null | wc -l) / $(ls ~/sweep_logs/*.m 2>/dev/null | wc -l) jobs"
  echo "result CSVs collected: $(ls "$OUT/csv" 2>/dev/null | wc -l)"
  echo "total result rows: $(cat "$OUT"/csv/results_*.csv 2>/dev/null | grep -vc '^$')"
  echo "--- verdict tally ---"
  cat "$OUT"/csv/results_*.csv 2>/dev/null | grep -v '^subfolder' | cut -d, -f6 | sort | uniq -c | sort -rn
} > "$OUT/SUMMARY.txt"
tar czf ~/sweep_results_$TS.tgz -C ~ "sweep_results_$TS" 2>/dev/null
cp "$OUT/SUMMARY.txt" ~/SWEEP_DONE
echo "[harvest] tarball: ~/sweep_results_$TS.tgz"
cat ~/SWEEP_DONE

# --- Phase 4: self-stop (shutdown behavior = stop). 60s grace so this log flushes.
echo "[harvest] self-stopping in 60s at $(date -u +%T)"
sleep 60
sudo shutdown -h now
