##############################################################################
# Authors: Ana Radulovic and Rok Blagus
# Part of the code for: "Estimating the number of changepoints in segmented
# regression: new sequential tests and comparative evaluation"
#
# 02_worker.R - one worker of the parallel analysis
#
# Usage: Rscript 02_worker.R <worker_id>   (started N times by run_all.sh)
# All workers go through the same list of jobs (00_setup.R: job_list(),
# longest jobs first). A worker skips jobs that are finished (result file
# exists) or taken by another worker (lock directory exists; creating a
# directory is atomic, so two workers can never take the same job), and runs
# the first free job. When no free jobs are left, the worker stops.
# Each job saves its result as soon as it is finished.
##############################################################################

args   <- commandArgs(trailingOnly = TRUE)
WORKER <- if (length(args) >= 1) as.integer(args[1]) else 1
source("00_setup.R")
jobs <- job_list()
cat(sprintf("Worker %d: %d jobs in total\n\n", WORKER, nrow(jobs)))

n_done <- 0
for (j in seq_len(nrow(jobs))) {
  job <- jobs[j, ]
  if (file.exists(result_file(job))) next                                   # finished
  if (!dir.create(file.path("locks", job$name), showWarnings = FALSE)) next  # taken
  cat(format(Sys.time(), "%Y-%m-%d %H:%M"), "start", job$name, "\n")
  res <- tryCatch(run_job(job), error = function(e) {
    cat("   ERROR in", job$name, ":", conditionMessage(e), "\n"); NULL })
  if (!is.null(res)) {
    n_done <- n_done + 1
    cat(format(Sys.time(), "%Y-%m-%d %H:%M"), "done ", job$name,
        if ("K" %in% names(res)) paste("| K =", paste(res$K, collapse = "/")) else "",
        "|", res$seconds[1], "s\n")
  }
}
cat(sprintf("\nWorker %d finished: %d jobs run, no free jobs left.\n", WORKER, n_done))
