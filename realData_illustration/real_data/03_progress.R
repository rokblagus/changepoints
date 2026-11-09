##############################################################################
# Authors: Ana Radulovic and Rok Blagus
#
# 03_progress.R - progress of the analysis: finished jobs by method and
# subsample size, and the jobs that are currently running
# Usage: Rscript 03_progress.R
##############################################################################
source("00_setup.R")
jobs <- job_list()
jobs$done    <- sapply(seq_len(nrow(jobs)), function(j) file.exists(result_file(jobs[j, ])))
jobs$running <- !jobs$done & dir.exists(file.path("locks", jobs$name))
jobs$size    <- ifelse(jobs$type == "complete", "complete", paste0("n", jobs$npy))
cat("\nFinished jobs (complete: out of 1 per method and bandwidth; subsamples: out of 100):\n")
print(with(jobs, tapply(done, list(what, size), sum)))
cat(sprintf("\nIn total: %d of %d jobs finished, %d running\n", sum(jobs$done), nrow(jobs), sum(jobs$running)))
if (any(jobs$running)) { cat("Running:\n"); print(jobs$name[jobs$running]) }
