#!/bin/bash
##############################################################################
# Authors: Ana Radulovic and Rok Blagus
#
# run_all.sh - runs the whole real-data analysis (complete subset and 300
# subsamples, all 22 methods) with N parallel workers (default 40).
#
# Usage (from the folder real_data, after slofit_subset.Rdata exists):
#   bash run_all.sh        # 40 workers
#   bash run_all.sh 30     # 30 workers
# Restart after an interruption: simply run it again (finished jobs are
# skipped). Do not start it while workers are still running, because the
# locks of unfinished jobs are removed at the start.
# Progress: Rscript 03_progress.R   |   Logs: logs/worker_<i>.log
##############################################################################
N=${1:-40}
mkdir -p logs results
rm -rf locks && mkdir locks
for i in $(seq 1 $N); do
  setsid nohup "$(R44 RHOME)/bin/Rscript" 02_worker.R $i > logs/worker_$i.log 2>&1 < /dev/null &
done
echo "Started $N workers. Progress: Rscript 03_progress.R"
