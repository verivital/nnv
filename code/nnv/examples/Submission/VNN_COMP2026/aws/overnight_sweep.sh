#!/bin/bash
# Overnight FULL 2026 sweep on the m5.16xlarge: 11 concurrent matlab -batch
# workers, each running run_all_benchmarks (per-instance parfeval timeout 300 s;
# vggnet16 gets its full budgets) over a folder group. CSVs are timestamped so
# they never collide; combine afterwards. Usage: bash overnight_sweep.sh
#
# Runs run_all_benchmarks from VNN_COMP2026. Until the 2026 migration branch merges
# to master, default to it; after merge override with NNV_BRANCH=master.
set -u
NNV_BRANCH="${NNV_BRANCH:-master}"
cd ~/nnv || exit 1
git remote get-url ttj >/dev/null 2>&1 || git remote add ttj https://github.com/ttj/nnv.git
git fetch -q ttj "$NNV_BRANCH" || { echo "FATAL: could not fetch $NNV_BRANCH"; exit 1; }
git checkout -q -f "ttj/$NNV_BRANCH" || { echo "FATAL: checkout failed"; exit 1; }
echo "sweep on $(git log --oneline -1)"

# PARSE GUARD: never launch onto an unparseable driver again. An orphan line that
# is just a quote followed by comma/paren = a string split by a literal newline.
RAB=code/nnv/examples/Submission/VNN_COMP2026/run_all_benchmarks.m
if grep -qE "^[[:space:]]*'[,)]" "$RAB"; then
  echo "FATAL: run_all_benchmarks.m still has an orphan-quote (parse error) line -- aborting"; exit 1
fi
echo "parse guard passed"

# regenerate ALL manifest-category .nnv.mat (incl. vit post-#345; nn4sys-mscn + test
# post-error-fix). Handles both .onnx and .onnx.gz, and 1.0/2.0 version dirs.
for b in lsnc_relu traffic_signs_recognition_2023 cgan_2023 cgan2026 soundnessbench soundnessbench_2026 vit_2023 nn4sys test; do
  for f in ~/vnncomp2026_benchmarks/benchmarks/$b/*/onnx/*.onnx ~/vnncomp2026_benchmarks/benchmarks/$b/*/onnx/*.onnx.gz; do
    [ -f "$f" ] || continue
    onnx="$f"; case "$f" in *.gz) gunzip -kf "$f"; onnx="${f%.gz}";; esac
    if [ "$b" = "nn4sys" ]; then
      # only mscn routes via the manifest; mscn_2048d_dual is corrupt upstream (skip)
      case "$onnx" in *mscn_2048d_dual*) continue;; *mscn*) ;; *) continue;; esac
    fi
    [ -f "${onnx%.onnx}.nnv.mat" ] || python3 ~/nnv/code/nnv/tools/onnx2nnv_python/onnx2nnv.py "$onnx" "${onnx%.onnx}.nnv.mat" > /dev/null 2>&1 || echo "WARN manifest: $onnx"
  done
done
BENCH=(
  "acasxu_2023"
  "cifar100_2024"
  "tinyimagenet_2024"
  "cora_2024 dist_shift_2023 linearizenn_2024"
  "safenlp_2024 cersyve test"
  "nn4sys sat_relu"
  "relusplitter relusplitter_2026"
  "malbeware metaroom_2023 collins_rul_cnn_2022"
  "vit_2023 yolo_2023 vggnet16_2022"
  "cgan_2023 cgan2026 soundnessbench soundnessbench_2026 challenging_certified_training_2026"
  "lsnc_relu traffic_signs_recognition_2023 cctsdb_yolo_2023 collins_aerospace_benchmark ml4acopf_2024 tllverifybench_2023 adaptive_cruise_control_non_linear_2026 isomorphic_acasxu_2026 monotonic_acasxu_2026 smart_turn_multimodal_2026"
)
mkdir -p ~/sweep_logs
for i in "${!BENCH[@]}"; do
  G="${BENCH[$i]}"
  FOLDERS=$(echo "$G" | sed "s/[^ ]*/'&'/g" | tr ' ' ',')
  cat > ~/sweep_logs/grp$i.m <<EOF
cd('/home/ubuntu/nnv/code/nnv'); startup_nnv;
addpath('/home/ubuntu/nnv/code/nnv/examples/Submission/VNN_COMP2026');
cd('/home/ubuntu/nnv/code/nnv/examples/Submission/VNN_COMP2026');
run_all_benchmarks('/home/ubuntu/vnncomp2026_benchmarks/benchmarks', 300, {$FOLDERS}, 'all');
fprintf('GROUP-$i-DONE\n');
EOF
  nohup matlab -batch "run('/home/ubuntu/sweep_logs/grp$i.m')" > ~/sweep_logs/grp$i.log 2>&1 &
  sleep 3
done
echo "launched ${#BENCH[@]} workers; logs in ~/sweep_logs/; combine results_*.csv when all GROUP-*-DONE"
