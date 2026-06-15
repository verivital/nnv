#!/bin/bash
# wait for install to finish in the licensed tmux session, then run the suite
for i in $(seq 1 120); do grep -q "INSTALL-DONE" /home/ubuntu/ml_session.log 2>/dev/null && break; sleep 15; done
if ! grep -q "INSTALL-DONE" /home/ubuntu/ml_session.log 2>/dev/null; then echo "INSTALL TIMEOUT" > /home/ubuntu/dryrun_status.txt; exit 1; fi
echo "install done; launching suite" > /home/ubuntu/dryrun_status.txt
tmux send-keys -t ml "run(\"/home/ubuntu/dryrun_suite.m\")" C-m
for i in $(seq 1 240); do grep -q "SUITE-DONE" /home/ubuntu/ml_session.log 2>/dev/null && break; sleep 15; done
grep -q "SUITE-DONE" /home/ubuntu/ml_session.log && echo "SUITE COMPLETE" > /home/ubuntu/dryrun_status.txt || echo "SUITE TIMEOUT" > /home/ubuntu/dryrun_status.txt
