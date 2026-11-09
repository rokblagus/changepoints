# Real-data illustration: BMI of Slovenian schoolchildren (SLOfit)

**Authors:** Ana Radulović and Rok Blagus

Code and anonymised data for the real-data illustration (Section 4 and
Supplement S2) of *Estimating the number of changepoints in segmented
regression: new sequential tests and comparative evaluation*.

## Data

`slofit_subset.Rdata` (object `slofit`; also `slofit_subset.csv`) is an
anonymised subset of the SLOfit database, a national surveillance system of
the physical fitness of Slovenian schoolchildren, with four variables:
`BMI` (with added noise), `SEX`, `AGE_YRS` and `YEAR` (calendar year of
measurement, 1989-2021). It was created from the complete SLOfit file by
`01_make_subset.R`; the complete file is not public.

**No repeated measurements by construction.** Children in SLOfit are measured
every year, but the anonymised data contain no child identifier. To ensure
that each child contributes at most one measurement (independent
observations, as assumed by all compared methods), each day of birth (day of
the year) is assigned to exactly one calendar year of measurement, cyclically
over the 33 years (day d to year 1989 + (d - 1) mod 33). In each year, only
children born on the days assigned to that year are eligible, so a child can
appear in at most one year. Within a year, at most one record per birth date,
sex and school is kept (possible double records of the same child), and records with unknown birth dates (recorded as 1 January)
are excluded. Repeated measurements could thus only occur if a birth date was
recorded inconsistently. Since the day of birth is unrelated to BMI, age, sex
and school, each year remains a random sample of children of all ages and
schools (and, unlike sampling a few schools per year, the design does not
induce clustering within schools). Then 1,000 records per year are sampled at
random.

**Anonymisation.** Strata (sex x age x year) with fewer than 5 children are
removed, and Gaussian noise N(0, 0.5^2) is added to each BMI value.

## Methods

As in the simulation study (22 methods): AIC and BIC; AIC-S and BIC-S
(`selgmented()` of the **segmented** package); two sequential
implementations of Davies' test (Davies-F with fixed and Davies-R with jointly
re-estimated changepoint locations); and the F-type and GoF tests in their
simulation-based (S, S2) and permutation-based (P, P-het, the latter with the
adjustment for heteroscedasticity) versions, each with and without the
Bonferroni adjustment (alpha = 0.05 / K_max). K_max = 5 and B = 1,000
permutations or simulations. Changepoint locations are estimated with
`segmented()` (bootstrap restarting with 50 samples, search restricted to the
middle 98% of the range of the predictor). For the heteroscedasticity-adjusted
tests, the bandwidth of the Nadaraya-Watson estimator is c/sqrt(n) with
c = (range of year)^2 / 2 = 512 (this gives 450/sqrt(n) for the predictor
range [0, 30] of the simulation study); it is multiplied by 0.5 and 2 in the
sensitivity analysis.

Table 1 is based on 300 stratified subsamples: 100 replicates for each of
200, 500 and 750 records per calendar year.

## Requirements

R >= 4.4.3 with the packages `segmented` (version 2.2-2), `data.table`,
`ggplot2`, and `foreign` (only for `01_make_subset.R`).

## Files

| File | Purpose |
|---|---|
| `00_setup.R` | common setup: data, constants, sequential procedure, subsampling design, list of all jobs |
| `01_make_subset.R` | creates the anonymised subset from the complete SLOfit file |
| `02_worker.R` | one worker of the parallel analysis (runs free jobs, longest first) |
| `run_all.sh` | starts N workers (default 40) that together run the whole analysis |
| `03_progress.R` | progress: finished and running jobs |
| `04_table1.R` | Table 1 and the full table for Supplement S2 |
| `05_figure4.R` | Figure 4 |
| `06_heteroscedasticity.R` | residual variability by year (justification of the heteroscedasticity-adjusted tests) |
| `07_coefficients.R` | coefficients with 95% CIs of the selected models; check of local optima |
| `functions/functions_ic_davies.R` | AIC, BIC, AIC-S, BIC-S, Davies-F, Davies-R |
| `functions/functions_tests.R` | F-type and GoF tests |
| `results/` | results: `complete_*.csv` for the complete subset, `sub/<method>/task*.csv` for the subsamples |

## How to run

```
Rscript 01_make_subset.R <path to 1989-2021_SLOfit_v2.sav>   # only with the complete file
bash run_all.sh 40              # the whole analysis with 40 workers
Rscript 03_progress.R           # progress
Rscript 04_table1.R             # when all jobs are finished
Rscript 05_figure4.R
Rscript 06_heteroscedasticity.R
Rscript 07_coefficients.R
```

The analysis consists of 2,713 jobs: all methods on the complete subset
(including the bandwidth sensitivity analysis) and every method on every
subsample. All workers go through the same list of jobs, longest first; a
worker takes the first job that is neither finished nor taken by another
worker (a lock directory is created atomically), so the load is balanced
automatically. Each job saves its own result file, so the analysis can be
interrupted and restarted at any time with `bash run_all.sh` (finished jobs
are skipped; do not restart while workers are still running).

**Computing time.** The permutation-based tests dominate: on the complete
subset up to about 12 hours per test; on the subsamples about 2,900
core-hours in total. With 40 cores, the whole analysis takes about 3 days.

## Reproducibility

Every subsample and every step of the sequential procedure has its own seed,
so the results do not depend on which worker or computer runs a job. With the
same versions of R and **segmented**, the results of the permutation-based
tests on an earlier subset were identical on two different computers.

Note on `selgmented()`: it re-creates the data from the `lm()` call by name,
which works only if the `data` argument is a simple name; the data are
therefore passed under a simple local name in `ic_rok()`.

## Correspondence with the paper

| Output (`results/`) | Paper |
|---|---|
| `table1_mean_sd.csv` | Table 1 (selected methods); Supplement S2 (all methods) |
| `table1_distribution.csv`, `s_vs_s2_agreement.csv` | Supplement S2 |
| `complete_*.csv` | Section 4; Supplement S2 (p-values, bandwidth sensitivity) |
| `figure4.pdf` | Figure 4 |
| `hetero_by_year.pdf`, `hetero_tests.csv` | Supplement S2 (heteroscedasticity) |
| `coefficients.csv`, `local_optima_check.csv` | Supplement S2 (coefficients) |
| `subset_steps.csv`, `subset_by_year.csv` | Supplement S2 (construction of the subset) |
