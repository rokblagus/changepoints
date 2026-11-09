##############################################################################
# Authors: Ana Radulovic and Rok Blagus
# Part of the code for: "Estimating the number of changepoints in segmented
# regression: new sequential tests and comparative evaluation"
#
# 01_make_subset.R - creates the anonymised subset used in the illustration
#
# Input : the complete SLOfit file (1989-2021_SLOfit_v2.sav), which is NOT
#         public; this script is run only by the data holders.
# Output: slofit_subset.Rdata (object 'slofit') and slofit_subset.csv, with
#         four variables only: BMI (with added noise), SEX, AGE_YRS, YEAR.
#
# Design (no repeated measurements of the same child by construction):
#  - each day of birth (day of the year, 1-366) is assigned to exactly one
#    calendar year of measurement, cyclically over the 33 years
#    (day d -> year 1989 + (d - 1) mod 33), and in each year only children
#    born on the days assigned to that year are eligible. Since a child has
#    only one birthday, he or she can appear in at most one year;
#  - within a year, at most one record per birth date, sex and school is
#    kept (possible double records of the same child in the same year);
#  - records with unknown birth dates (recorded as 1 January) are excluded.
# Repeated measurements could thus only occur if a birth date was recorded
# inconsistently. Since the day of birth is unrelated to BMI, age, sex and
# school, each year is still a random sample of children of all ages and
# schools, and the cyclic assignment includes children born in all months.
#
# Then 1,000 records per year are sampled at random, and the data are
# anonymised as before: strata (sex x age x year) with fewer than 5 children
# are removed, and Gaussian noise N(0, 0.5^2) is added to each BMI value.
##############################################################################

# path to the complete SLOfit file: Rscript 01_make_subset.R <path>
# (default: the file in the current folder)
args      <- commandArgs(trailingOnly = TRUE)
FULL_FILE <- if (length(args) >= 1) args[1] else "1989-2021_SLOfit_v2.sav"
dir.create("results", showWarnings = FALSE)
N_PER_YEAR  <- 1000
SEED_SAMPLE <- 2026
SEED_NOISE  <- 123

dd <- foreign::read.spss(FULL_FILE, to.data.frame = TRUE)
cat("Complete data:", nrow(dd), "records\n")

birth <- as.Date(dd$BIRTH / 86400, origin = "1582-10-14")   # SPSS date
d <- data.frame(YEAR = dd$LETO_MERITEV, SEX = dd$SEX, AGE_YRS = dd$AGE_YRS,
                BMI = dd$BMI, BIRTH = birth, SCHOOL = dd$SOLA)
rm(dd); invisible(gc())

steps <- data.frame(step = character(0), records = integer(0))
add_step <- function(label) steps[nrow(steps) + 1, ] <<- list(label, nrow(d))
add_step("complete data")

d <- d[complete.cases(d[, c("YEAR", "SEX", "AGE_YRS", "BMI", "BIRTH")]), ];                                add_step("non-missing BMI, sex, age, year, birth date")
d <- d[d$AGE_YRS >= 6 & d$AGE_YRS <= 18, ];                 add_step("age 6-18 years")
d <- d[format(d$BIRTH, "%m-%d") != "01-01", ];              add_step("birth date not 1 January")

years <- sort(unique(d$YEAR))
stopifnot(length(years) == 33, min(years) == 1989)
doy <- as.integer(format(d$BIRTH, "%j"))
d$ASSIGNED_YEAR <- min(years) + (doy - 1) %% length(years)
d <- d[d$YEAR == d$ASSIGNED_YEAR, ];                        add_step("born on a day assigned to the year of measurement")

set.seed(SEED_SAMPLE)
d <- d[sample(nrow(d)), ]                                   # random order
d <- d[!duplicated(d[, c("YEAR", "BIRTH", "SEX", "SCHOOL")]), ];      add_step("one record per birth date, sex and school within a year")
eligible <- table(d$YEAR)

d <- do.call(rbind, lapply(split(d, d$YEAR), function(z) z[seq_len(min(N_PER_YEAR, nrow(z))), ]))
add_step(sprintf("random sample of %d per year", N_PER_YEAR))

# anonymisation
strata <- interaction(d$SEX, d$AGE_YRS, d$YEAR, drop = TRUE)
small  <- names(which(table(strata) < 5))
cat(sprintf("Strata (sex x age x year): %d, of which %d with fewer than 5 children (%d records removed)\n",
            nlevels(strata), length(small), sum(strata %in% small)))
d <- d[!(strata %in% small), ];                             add_step("strata with at least 5 children")
d <- d[order(d$YEAR, d$SEX, d$AGE_YRS), ]

set.seed(SEED_NOISE)
slofit <- data.frame(BMI = d$BMI + rnorm(nrow(d), mean = 0, sd = 0.5),
                     SEX = d$SEX, AGE_YRS = d$AGE_YRS, YEAR = d$YEAR)
rownames(slofit) <- NULL

save(slofit, file = "slofit_subset.Rdata")
write.csv(slofit, "slofit_subset.csv", row.names = FALSE)
write.csv(steps, "results/subset_steps.csv", row.names = FALSE)
write.csv(data.frame(YEAR = names(eligible), eligible = as.integer(eligible),
                     sampled = as.integer(table(slofit$YEAR)[names(eligible)])),
          "results/subset_by_year.csv", row.names = FALSE)

cat("\n"); print(steps, row.names = FALSE)
cat("\nEligible records per year: min", min(eligible), "| max", max(eligible), "\n")
cat("Saved: slofit_subset.Rdata, slofit_subset.csv (", nrow(slofit), "records )\n")
