#!/bin/bash
# Sweep reevaluation: confirm the merged fixes (#349 ReshapeLayer, #350
# dispatcher, #352 cgan) unlock the previously-erroring/timing-out benchmarks.
# Runs on FRESH origin/master (which now has the fixes) in a SINGLE copy (disk-
# careful after the earlier crisis). Compare verdicts vs the sweep's error rows.
set -u
TBX='clpmex:1.0:glnxa64\nglpkmex:1.0:glnxa64\nfourier:1.0:glnxa64\nsedumi:1.3:glnxa64\nmpt:3.2.1:all\nmptdoc:3.0.4:all\nlcp:1.0.3:glnxa64\nhysdel:2.0.6:glnxa64\ncddmex:1.0.1:glnxa64\n'
echo "disk before: $(df -h /|tail -1|awk '{print $4}') free"

[ -d ~/nnv-fix ] || cp -r ~/nnv ~/nnv-fix
cd ~/nnv-fix
git fetch -q origin master && git checkout -q -f origin/master
echo "fix-validation on $(git log --oneline -1)"
grep -q "isprop(layer, 'Vars')" code/nnv/engine/nn/layers/ReshapeLayer.m && echo "  ReshapeLayer fix present" || echo "  WARN: ReshapeLayer fix MISSING"
grep -q 'contains(category, "adaptive_cruise")' code/nnv/examples/Submission/VNN_COMP2026/run_vnncomp_instance.m && echo "  adaptive_cruise fix present" || echo "  WARN: adaptive_cruise fix MISSING"
grep -q 'contains(category, "smart_turn")' code/nnv/examples/Submission/VNN_COMP2026/run_vnncomp_instance.m && echo "  smart_turn fix present" || echo "  WARN: smart_turn fix MISSING"

RAB=code/nnv/examples/Submission/VNN_COMP2026/run_all_benchmarks.m
sed -i "s/datestr(now, 'yyyymmdd_HHMMSS')/datestr(now, 'yyyymmdd_HHMMSSFFF')/" "$RAB"
sed -i 's/cold_grace_s = 120;/cold_grace_s = 60;/' "$RAB"
printf "$TBX" > code/nnv/tbxmanager/tbxenabled.txt
echo "disk after copy: $(df -h /|tail -1|awk '{print $4}') free"

ROOT=$HOME/vnncomp2026_benchmarks/benchmarks
# benchmark:timeout -- the fixed ones. smart_turn abstains (instant unknown).
BENCH=( "relusplitter_2026:120" "relusplitter:120" "adaptive_cruise_control_non_linear_2026:60" "traffic_signs_recognition_2023:60" "smart_turn_multimodal_2026:30" "cgan_2023:60" "cgan2026:60" )
# regenerate manifests for the manifest-category fixed benchmarks (traffic/cgan)
for gz in "$ROOT"/traffic_signs_recognition_2023/*/onnx/*.onnx.gz "$ROOT"/cgan*/*/onnx/*.onnx.gz; do
  [ -f "$gz" ] || continue; gunzip -kf "$gz"; onnx="${gz%.gz}"
  [ -f "${onnx%.onnx}.nnv.mat" ] || python3 ~/nnv-fix/code/nnv/tools/onnx2nnv_python/onnx2nnv.py "$onnx" "${onnx%.onnx}.nnv.mat" >/dev/null 2>&1 || echo "WARN manifest: $onnx"
done
mkdir -p ~/fixval_logs
job=0
for bt in "${BENCH[@]}"; do
  bench="${bt%%:*}"; to="${bt##*:}"
  [ -d "$ROOT/$bench" ] || { echo "skip missing: $bench"; continue; }
  M=~/fixval_logs/fv${job}_${bench}.m
  cat > "$M" <<EOF
jraf = java.io.RandomAccessFile('/home/ubuntu/fixval_logs/startup.lock','rw'); jlk = jraf.getChannel().lock();
try
    cd('/home/ubuntu/nnv-fix/code/nnv'); startup_nnv;
catch ME
    jlk.release(); jraf.close(); rethrow(ME);
end
jlk.release(); jraf.close();
addpath('/home/ubuntu/nnv-fix/code/nnv/examples/Submission/VNN_COMP2026');
cd('/home/ubuntu/nnv-fix/code/nnv/examples/Submission/VNN_COMP2026');
run_all_benchmarks('$ROOT', $to, {'$bench'}, 'all');
fprintf('FV-${job}-${bench}-DONE\n');
EOF
  nohup matlab -batch "run('$M')" > ~/fixval_logs/fv${job}_${bench}.log 2>&1 &
  echo "launched fv$job: $bench @${to}s"
  job=$((job+1)); sleep 3
done
echo "launched $job fix-validation jobs; disk=$(df -h /|tail -1|awk '{print $4}') free"
