# =====================================================================
# Authors: Ana Radulovic and Rok Blagus
# Part of the code for: "Estimating the number of changepoints in segmented
# regression: new sequential tests and comparative evaluation"
#
# functions_ic_davies.R - AIC, BIC, AIC-S, BIC-S (ic_rok) and the two
# sequential implementations of Davies' test (Davies-F: davies_fixed,
# Davies-R: davies_reestimated), as in the simulation study, with the model
# formula Y ~ SEX + AGE_YRS + X.
#
# Note on selgmented(): it re-creates the data from the lm() call by name,
# which works only if the data argument is a simple name. Therefore, the data
# are passed to lm() under a simple local name (dat_ic) in ic_rok().
# =====================================================================

library(segmented)

ic_rok <- function(mydata, max_segments = 5) {

  # The data are passed to lm() under a simple local name.
  # selgmented() re-creates the data from the lm() call by name; with
  # data = mydata$data this fails ("only the first element is used as variable
  # name") and selgmented() then evaluates mydata$data in the global environment,
  # i.e. it fits the models with K >= 1 to whatever global 'mydata' exists.
  dat_ic <- mydata$data
  fit <- try(lm(Y ~ SEX + AGE_YRS + X, data = dat_ic), silent = TRUE)

  ###using selgmented (AIC-S, BIC-S) (as in the simulation study)
  selfit.aic <- try(selgmented(fit, seg.Z = ~X, alpha = 0.05, type = "aic",
                               control = seg.control(n.boot = 50, alpha = 0.01),
                               return.fit = TRUE, bonferroni = TRUE,
                               Kmax = max_segments, msg = TRUE), silent = TRUE)
  selfit.bic <- try(selgmented(fit, seg.Z = ~X, alpha = 0.05, type = "bic",
                               control = seg.control(n.boot = 50, alpha = 0.01),
                               return.fit = TRUE, bonferroni = TRUE,
                               Kmax = max_segments, msg = TRUE), silent = TRUE)
  if (class(selfit.aic)[1] != "try-error") sel.aic <- selfit.aic$selection.psi$npsi else sel.aic <- NA
  if (class(selfit.bic)[1] != "try-error") sel.bic <- selfit.bic$selection.psi$npsi else sel.bic <- NA

  ###"our" approaches (AIC, BIC) (as in the simulation study)
  bics <- aics <- rep(NA, max_segments + 1)
  npsi_actual  <- c(0, rep(NA, max_segments))          # added: diagnostics only
  psi_list     <- vector("list", max_segments + 1)     # added: diagnostics only
  bics[1] <- BIC(fit)
  aics[1] <- AIC(fit)
  for (ii in 1:max_segments) {
    segmented_fit <- try(segmented(fit, seg.Z = ~X, npsi = ii,
                                   control = seg.control(n.boot = 50, alpha = 0.01)),
                         silent = TRUE)
    if (class(segmented_fit)[1] != "try-error") {
      bics[ii + 1] <- BIC(segmented_fit)
      aics[ii + 1] <- AIC(segmented_fit)
      npsi_actual[ii + 1] <- if (is.null(segmented_fit$psi)) 0 else nrow(segmented_fit$psi)
      if (!is.null(segmented_fit$psi)) psi_list[[ii + 1]] <- sort(segmented_fit$psi[, "Est."])
    }
  }
  n.bic <- which.min(bics) - 1
  n.aic <- which.min(aics) - 1

  psi_sel <- function(sf) if (class(sf)[1] == "try-error" || is.null(sf$psi)) numeric(0) else sort(sf$psi[, "Est."])

  list(n.aic = n.aic, n.bic = n.bic, sel.aic = sel.aic, sel.bic = sel.bic,
       aics = aics, bics = bics, npsi_actual = npsi_actual, psi_list = psi_list,
       psi_selaic = psi_sel(selfit.aic), psi_selbic = psi_sel(selfit.bic))
}


davies_fixed <- function(mydata, max_segments, alpha.a = 0.05) {
  data.dav <- mydata$data
  mm <- 0; pd <- NA; chpnts <- numeric(0); pv <- numeric(0)
  fit.0 <- try(lm(Y ~ SEX + AGE_YRS + X, data = data.dav), silent = TRUE)
  if (inherits(fit.0, "try-error")) return(list(nm.davies = NA, chpnts = chpnts, pv = pv))
  fitlm <- fit.0
  while (mm < max_segments) {
    pd <- try(davies.test(fitlm, seg.Z = ~X), silent = TRUE)
    if (inherits(pd, "try-error")) { if (mm == 0) mm <- NA; break }
    pv <- c(pv, pd$p.value)
    if (pd$p.value > alpha.a) break
    fit1 <- try(segmented(fitlm, seg.Z = ~X, npsi = 1,
                          control = seg.control(n.boot = 50, alpha = 0.01)), silent = TRUE)
    if (inherits(fit1, "try-error") || is.null(fit1$psi) || nrow(fit1$psi) < 1) break
    psi.new <- fit1$psi[1, 2]
    GGi <- pmax(data.dav$X - psi.new, 0)
    if (is.null(data.dav$ZZ)) data.dav$ZZ <- GGi else data.dav$ZZ <- cbind(data.dav$ZZ, GGi)
    fit2 <- try(lm(Y ~ SEX + AGE_YRS + X + ZZ, data = data.dav), silent = TRUE)
    if (inherits(fit2, "try-error")) break
    mm <- mm + 1; chpnts <- c(chpnts, psi.new); fitlm <- fit2
  }
  list(nm.davies = mm, chpnts = chpnts, pv = pv)
}

davies_reestimated <- function(mydata, max_segments, alpha.a = 0.05) {
  data.dav <- mydata$data
  mm <- 0; pd <- NA; chpnts <- numeric(0); pv <- numeric(0)
  fit.0 <- try(lm(Y ~ SEX + AGE_YRS + X, data = data.dav), silent = TRUE)
  if (inherits(fit.0, "try-error")) return(list(nm.davies = NA, chpnts = chpnts, pv = pv))
  fitlm <- fit.0
  while (mm < max_segments) {
    pd <- try(davies.test(fitlm, seg.Z = ~X), silent = TRUE)
    if (inherits(pd, "try-error")) { cat("davies.test error:", as.character(pd)); if (mm == 0) mm <- NA; break }
    pv <- c(pv, pd$p.value)
    if (pd$p.value > alpha.a) break
    npsi.new <- mm + 1
    fit1 <- try(segmented(fit.0, seg.Z = ~X, npsi = npsi.new,
                          control = seg.control(n.boot = 50, alpha = 0.01)), silent = TRUE)
    if (inherits(fit1, "try-error") || is.null(fit1$psi)) break
    npsi.actual <- nrow(fit1$psi)
    if (is.null(npsi.actual) || npsi.actual < npsi.new) break
    mm <- npsi.actual; fitlm <- fit1; chpnts <- as.numeric(fit1$psi[, 2])
  }
  list(nm.davies = mm, chpnts = chpnts, pv = pv)
}


est.nochp <- function(fit) ifelse(is.null(nrow(fit$psi)), 0, nrow(fit$psi))
psi_of    <- function(f) if (is.null(f$psi)) numeric(0) else sort(f$psi[, "Est."])
