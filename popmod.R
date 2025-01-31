########
# popmod
########

# Brownsk populasjonsmodell


popmod <- function(pop.size, 
                   years = NULL, 
                   nboot = 10000, 
                   name = "", 
                   plot = TRUE, 
                   se = F, r = F, b1 = F, b2 = F, b3 = F, 
                   mean1 = F, mean2 = F, mean3 = F, 
                   bootout = TRUE) {
  # sd is demographic variance

#  catch <- NULL
  sd = 0.1
#  cov1 = NULL
#  cov2 <- NULL
#  cov3 <- NULL
#  nsim <- 0
  
  X <- ln(pop.size)
  
#¤#  if (is.vector(X)) X <- data.frame("lnN" = X)
  
  na.und <- function(x) ifelse(is.na(x), 0, x)
  
  log.gauss <- function(X, my, sig2) -0.5*log(2*pi*sig2)-0.5*((X-my)^2)/sig2
  
  ######### likelihood function

  lik <- function(par) {
    L <- 0
#    k  <- na.und(-exp(par[q[1]]))
    m  <- par[1]
    s  <- exp(par[2])
    w  <- log(exp(X))
#    b1 <- na.und(par[q[4]])
#    b2 <- na.und(par[q[5]])
#    b3 <- na.und(par[q[6]])
#    my <- w - 0.5 * sd * exp(-w) + m + k * exp(w)
#    s <- s + sd * exp(-w)
#    L  <- L - sum(log.gauss(X[-1], my[-length(my)], s[-length(s)]),
#                  na.rm=T)
    my <- w - 0.5 * sd * exp(-w) + m
    s <- s + sd * exp(-w)
    L  <- L - sum(log.gauss(X[-1], my[-length(my)], s[-length(s)]), na.rm=T)
    return(L)
  }
  
  ############################################
  
  estparam <- function(par) {
    test <- optim(par, lik, 
                 control = list(trace = F,
                                maxit = 12000,
                                beta = 0.1,
                                gamma = 1.2))
    return(test)
  }
  
  ########## set initial values for optimization and check for density dependence

  #if (!is.null(bootout))
  # {if (file.exists(bootout))
  #   {stop("ERROR: the specified file already exists.\n")}
  #  else
  #   {if(!file.create(bootout))
  #     {stop("ERROR: the specified file could not be created.\n")}
  #    else
  #     {? <- file.remove(bootout)}}}
  if (is.null(years)) stop(
    "timepoints of data must be specified, if no NA's, try years=1:length(X)")
  if (any(X > 30, na.rm = T)) {
    cat("\nOBS: Jeg har en viss mistanke om at du har glemt " %+%
         "å log-transformere bestandstallene.\n\n")
  }
  if (any(X < 0, na.rm = T)) stop ("N kan ikke være mindre enn 1!")
  if (length(na.omit(X)) < 6) stop ("For få år med data!")
#  if (is.null(catch)) catch <- rep(0, length(years))
#  if (is.null(cov1))   cov1 <- rep(0, length(years))
#  if (is.null(cov2))   cov2 <- rep(0, length(years))
#  if (is.null(cov3))   cov3 <- rep(0, length(years))
  if (years %!=% (min(years):max(years))) {
    cat("\nOBS: Årstallene er ikke sammenhengende. Er det meningen?\n\n")
    NAS <- cbind(years, X)
    NAS <- merge(years = min(years):max(years), NAS, all.x = T)
    years <- NAS[,1]
    X <- as.vector(NAS[, 2])
  }
  
  orig.years <- data.frame(years=min(years):max(years))
  
  bootres <- psim <- NULL
  predicted <- residual <- demcomp <- X
  pred.r <- lci <- uci <- X
  p.res.t <- p.res.N <- p.white <- 1
  N <- exp(X) #¤ tja, eller pop.size!
  names(N) <- years
  NAs <- is.na(X)

  X.orig <- X
  par <- c(pi / 100, -2)
#  q <- 0
#  q[2] <- 1
#  q[3] <- 2

  ######################### estimation of parameters

  estimat <- estparam(par)
  my.est  <- se.est <- b1.est <- b2.est <- b3.est <- 0
  my.est  <-     estimat$par[1]
  se.est  <- exp(estimat$par[2])
  lnL     <-    -estimat$value
  if (estimat$convergence) {
    cat("\n##########################################" %+%
        "\nOBS: Optimaliseringa konvergerte ikke!!!!!" %+%
        "\n##########################################\n")
  }
  
  ww <- log(exp(X))
  predicted <- c(NA, ww + my.est - 0.5 * sd * exp(-ww))
  residual <-  c(X, NA) - predicted
  residual <- residual[-1]
  predicted <- predicted[-length(predicted)]
  demcomp <- sd / (exp(X))
  df.res <- data.frame("years" = years,
                       "residual" = residual,
                       "demcomp" = demcomp, 
                       "predicted" = predicted)
  res <- merge(orig.years, df.res, all.x = T)
  residual <- res$residual
  demcomp <- res$demcomp
  predicted <- res$predicted
  names(predicted) <- names(demcomp) <- names(residual) <- 
    as.character(orig.years$years)
  
  #####################################################
  
  r2 <- NA
  pred.r <- my.est - 0.5 * sd * exp(-X)
  r2 <- summary(lm(diff(X) ~ pred.r[-length(X)]))$r.sq
  uci <- lci <- NA
  for (j in 1:(length(pred.r))) {
    lci[j] <- pred.r[j] + sqrt(se.est + sd * exp(-X[j])) * qnorm(0.025)
    uci[j] <- pred.r[j] + sqrt(se.est + sd * exp(-X[j])) * qnorm(0.975)
  }
  names(pred.r) <- names(lci) <- names(uci) <- as.character(orig.years$years)

  ########## simulate data based on estimated parameters
  # NY versjon av simdat
  # (skrevet av Hanno etter forbilde av Vidars logismodfit):
  
  simdat <- function() {
    sim  <- function() {
      erst <- min(which(!is.na(X.orig)))
      xsim <- rep(NA,   length(X.orig))
      xsim[erst] <- X.orig[erst]
      for (j in erst:(length(X.orig) - 1)) {
        xsim[j + 1] <- xsim[j] - 0.5 * sd * exp(-xsim[j]) + my.est -
          sqrt(se.est + sd * exp(-xsim[j])) * rnorm(1)
        if (xsim[j+1] < 0) xsim[j+1] <- 0
      }
      return(xsim)
    }
    xsim <- X.orig
    xsim <- sim()
    xsim[is.na(X.orig)] <- 12
    zahl <- max(sum(!is.na(X.orig)) / 12, any(X.orig == 0, na.rm = T))
    while (sum(xsim == 0) > zahl) {
      xsim <- sim()
      xsim[is.na(X.orig)] <- 12
    }
    xsim[is.na(X.orig)] <- NA
    return(xsim)
  }

  simdat2 <- function() {
    # IKKE tilpasset til logistisk modell!!!
    # og så gir det ikke mening per nå!
    sim <- function() {
      erst <- min(which(!is.na(X.orig)))
      lezt <- max(which(!is.na(X.orig)))
      xsim <- rep(NA, length(X.orig))
      xsim[erst] <- X.orig[lezt]
      for (j in (lezt - 1):erst) xsim[erst] <- xsim[erst] - my.est
      eps <- 5
      while (eps > 0.000001) {
        for (j in erst:(length(X.orig) - 1)) {
          xsim[j+1] <- xsim[j] - 0.5 * sd * exp(-xsim[j]) + my.est +
            sqrt(se.est + sd * exp(-xsim[j])) * rnorm(1)
          if (xsim[j+1] < 0) xsim[j+1] <- 0
        }
        if (xsim[length(xsim)] < X.orig[length(X.orig)]) {
          xsim[erst] <- xsim[erst] + eps
        } else {
          xsim[erst] <- xsim[erst] - eps
        }
      }
      return(xsim)
    }
    xsim <- X.orig
    xsim <- sim()
    xsim[is.na(X.orig)] <- 12
    zahl <- max(sum(!is.na(X.orig)) / 12, any(X.orig == 0, na.rm = T))
    while (sum(xsim == 0) > zahl) {
      xsim <- sim()
      xsim[is.na(X.orig)] <- 12
    }
    xsim[is.na(X.orig)] <- NA
    return(xsim)
  }
  
  ############################### parametric bootstrapping

  if (nboot > 0) {
    cat("bootstrapping\n")
    catval <- c(1, seq(0, nboot, nboot / 10))
    boottab <- matrix(NA, nboot, 2)
    for (j in 1:nboot) {
      if (any(j == catval)) cat("boot", j, "of", nboot, "\n")
      X <- simdat()
      estimat <- estparam(par)
      boottab[j, ] <- estimat$par
    }
    boottab[, 2] <- exp(boottab[, 2])
    bootres <- matrix(NA, 18, 2)
    rownames(bootres) <- c("min", "0.05%", "0.1%", "0.5%", "1%", "2.5%",
                           "5%", "25%", "median", "75%", "95%", "97.5%", "99%", 
                           "99.5%", "99.9%", "99.95%", "max", "mean")
    colnames(bootres) <- c("r", "sig2e")
    for (j in 1:ncol(bootres)) {
      bootres[1:17, j] <- quantile(boottab[, j], 
                                   c(0, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 
                                     0.25, 0.5, 0.75, 0.95, 0.975, 0.99, 0.995, 
                                     0.999, 0.9995, 1), na.rm = T)
      bootres[18, j] <- mean(boottab[, j], na.rm = T)
    }
  }
  if (!bootout) boottab <- NULL
  
  ####### 3 TESTS!!!

  lmfit <- lm(residual ~ orig.years$years, na.action = "na.exclude")
  p.res.t <- signif(anova(lmfit)[1,5],3)
    
  lmfit <- lm(residual ~ N, na.action = "na.exclude")
  p.res.N <- signif(anova(lmfit)[1, 5], 3)
    
  portman.Q <- function (X, K) {
      # portman.Q uses the cummulative ACF to test for whiteness  of a time series.
      # This is the Ljung-Box version of the the Portemanteau  test for
      # whiteness (Tong 1990). It may in particular be  usefull to test
      # for whiteness in the residuals from time series models.
      # 
      # A vector is returned consisting of the asymtpotic chi-square
      # value, the associated d.f. and asymptotic p.val for the test of
      # whiteness. p.val < 0.05 -> non-whiteness!
      # Tong, H. (1990) Non-linear time series : a dynamical  system approach. Clarendon Press, Oxford.
      # Author: Ottar N. Bjornstad onb1@psu.edu
    K <- min(K, length(na.omit(X)) - 1) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Q <- 0
    n <- length(X)
    p <- acf(X, plot = FALSE, lag.max = K, na.action=na.pass)$acf[2:(K + 1)]
    for (k in 1:K) Q <- Q + p[k]^2 / (n - k)
    Q <- n * (n + 2) * Q
    res <- list(chisq = signif(Q, 3), df = K, p.val = signif(1 - pchisq(Q, K), 3))
    unlist(res)
  }
  p.white <- portman.Q(residual, K = 10)[[3]]

  #############################
  
#  if (!is.null(names(N))) {
#    navn <- names(N)
#    names(k.est) <- names(my.est) <- names(se.est) <- names(r2) <- navn
#    names(b1.est) <- names(b2.est) <- names(b3.est) <- navn
#    names(p.res.t) <- names(p.res.N) <- names(p.white) <- navn
#  }
  
  n.par <- 2
  n.eff <- length(na.omit(predicted))
  aic   <- 2 * (n.par - lnL)
  aicc  <- aic + 2 * n.par * (n.par + 1) / (n.eff - n.par - 1)

  res <- list(
    "K" = NA,
    "s" = my.est,
    "sig2" = se.est,
    "sigd" = sd,
    "years" = orig.years$years,
    "N" = N,
    "residuals" = residual, "predicted" = predicted, "demcomp" = demcomp,
    "pred.r" = pred.r, "lci" = lci, "uci" = uci,
    "bootout" = boottab, "boot" = bootres,
    "model" = "Brownian", "name"=name,
    "p.res.t" = p.res.t, "p.res.N" = p.res.N, "p.white" = p.white,
    "lnL" = lnL, "r2" = r2, "aic" = aic, "aicc" = aicc)
  if (plot) {
    brownian.plot(res)
    cat("\nSee plot for diagnostics.\n")
  }
  cat("\nA Brownian (density-independent) model was fitted.\n")
  cat("\nThe demographic variance was assumed to be " %+% round(sd, 3) %+% ".\n")
  cat("\nParameter estimates:\n")
  partab <- matrix(NA, 2, 3)
  partab  [1,1]   <- res$s
  partab  [1,2:3] <- res$boot[c("2.5%", "97.5%"), "r"]
  partab  [2,1]   <- res$sig2
  partab  [2,2:3] <- res$boot[c("2.5%", "97.5%"), "sig2e"]
  rownames(partab) <- c("Popul. growth rate (r):", "Environmental variance:")
  colnames(partab) <- c("Mean", "Lower 95% CI", "Upper 95% CI")
  print(partab)
  return(res)
}

