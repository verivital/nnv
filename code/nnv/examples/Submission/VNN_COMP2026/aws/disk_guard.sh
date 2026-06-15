#!/bin/bash
# Disk guard for the overnight sweep. Decompressed .onnx files (each has a .gz
# sibling, so they're reclaimable) accumulate as jobs run and can fill the disk.
# When free space drops below the threshold, delete the decompressed .onnx ONLY
# for benchmarks whose job is already DONE (they won't be re-read), keeping the
# .gz so they could be regenerated. Safe: never touches a still-running
# benchmark's models. Runs until the sweep finishes (SWEEP_DONE) or MATLAB exits.
THRESH_G=15
while true; do
  [ -f ~/SWEEP_DONE ] && { echo "$(date -u +%T) SWEEP_DONE -> guard exit" >> ~/sweep_logs/diskguard.log; break; }
  [ "$(pgrep -c MATLAB)" = "0" ] && { echo "$(date -u +%T) no MATLAB -> guard exit" >> ~/sweep_logs/diskguard.log; break; }
  free_g=$(df -BG --output=avail / | tail -1 | tr -dc '0-9')
  echo "$(date -u +%T) free=${free_g}G" >> ~/sweep_logs/diskguard.log
  if [ "${free_g:-99}" -lt "$THRESH_G" ]; then
    reclaimed=0
    for f in ~/sweep_logs/job*.log; do
      grep -qE 'JOB-.*-DONE' "$f" 2>/dev/null || continue
      base=$(basename "$f" .log)                          # job12_2026_cora_2024
      yr=$(echo "$base" | grep -oE '_20[0-9][0-9]_' | tr -d _ | head -1)
      bench=$(echo "$base" | sed -E 's/^job[0-9]+_20[0-9][0-9]_//')
      d="$HOME/vnncomp${yr}_benchmarks/benchmarks/$bench"
      [ -d "$d" ] || continue
      while IFS= read -r onnx; do
        [ -f "${onnx}.gz" ] && { rm -f "$onnx" && reclaimed=$((reclaimed+1)); }
      done < <(find "$d" -name '*.onnx' 2>/dev/null)
    done
    echo "$(date -u +%T) LOW DISK (${free_g}G<${THRESH_G}G) -> reclaimed $reclaimed .onnx from done benchmarks; now $(df -BG --output=avail / | tail -1 | tr -dc '0-9')G" >> ~/sweep_logs/diskguard.log
  fi
  sleep 420
done
