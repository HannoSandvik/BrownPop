
########
BrownPop <- function(
########
                     pop.size,      # time series with population sizes
                     years = NULL,  # years of pop.size
                     name = ""      # name of dataset
                     ) {
  
  # Brownian population model
  # =========================
  # Code written by Hanno Sandvik, Norwegian Institute for Nature Research, 2025,
  # based on code by Vidar Grøtan, Norwegian University of Science and Technology 

    
  ##############################################################################
  # Preparations
  
  # Constants ------------------------------------------------------------------
  sd <- 0.1                               # demographic variance
  qntl <- c(0.05, 0.25, 0.5, 0.75, 0.95)  # quantiles to be used in the output
  nboot <- 100000                         # number of simulations
  seed <- NULL                            # seed for the random number generator

  # Functions ------------------------------------------------------------------

  # Natural logarithm should be abbreviated as "ln"!
  ln <- function(x) log(x)
  
  
  # Tests equality of arguments - insensitive to rounding error!
  "%=%" <- function(arg1, arg2) { 
    attributes(arg1) <- NULL
    attributes(arg2) <- NULL
    return(identical(all.equal(arg1, arg2), TRUE))
  }
  
  
  # Combines text variables into one string
  "%+%" <- function(string1, string2) paste0(string1, string2)
  
  
  # Log-likelihood function of normal distribution
  log.gauss <- function(X, my, sig2) -0.5*ln(2*pi*sig2)-0.5*((X-my)^2)/sig2
  

  lik <- function(par) { # likelihood function
    L <- 0
    m  <- par[1]
    s  <- exp(par[2])
    w  <- X
    my <- w - 0.5 * sd * exp(-w) + m
    s <- s + sd * exp(-w)
    L  <- L - sum(log.gauss(X[-1], my[-length(my)], s[-length(s)]), na.rm=T)
    return(L)
  }
  

  estparam <- function(par) { # parameter estimation
    test <- optim(par, lik, 
                  control = list(trace = F,
                                 maxit = 12000,
                                 beta = 0.1,
                                 gamma = 1.2))
    return(test)
  }


  simdat  <- function(i) { # simulation of population trajectories
    xsim <- matrix(NA, length(i), length(X.orig))
    xsim[, 1] <- X.orig[1]
    for (j in 2:length(X.orig)) {
      # A stochastic Brownian model is used to simulate population sizes
      xsim[, j] <- xsim[, j - 1] - 0.5 * sd * exp(-xsim[, j - 1]) + my -
        sqrt(se + sd * exp(-xsim[, j - 1])) * rnorm(length(i))
      xsim[xsim[, j] < 0, j] <- 0
    }
    xsim[, is.na(X.orig)] <- NA
    return(xsim)
  }
  
  
  portman.Q <- function (X, K = 10) { # Portemanteau test for whiteness
    # portman.Q uses the cumulative ACF to test for whiteness  of a time series.
    # This is the Ljung-Box version of the the Portemanteau  test for
    # whiteness (Tong 1990). It may in particular be  usefull to test
    # for whiteness in the residuals from time series models.
    # Tong, H. (1990) Non-linear time series: a dynamical  system approach. 
    # Clarendon Press, Oxford.
    # Author: Ottar N. Bjørnstad, onb1@psu.edu
    K <- min(K, length(na.omit(X)) - 1)
    Q <- 0
    n <- length(X)
    p <- acf(X, plot = FALSE, lag.max = K, na.action=na.pass)$acf[2:(K + 1)]
    for (k in 1:K) Q <- Q + p[k]^2 / (n - k)
    Q <- n * (n + 2) * Q
    res <- list(chisq = signif(Q, 3), df = K, p.val = signif(1 - pchisq(Q, K), 3))
    unlist(res)
  }
  
  
  # Handling of input ----------------------------------------------------------
  N <- pop.size
  while (is.na(N[1]))         N <- N[-1]
  while (is.na(N[length(N)])) N <- N[-length(N)]
  
  if (any(N < 1, na.rm = T)) stop ("Population sizes cannot be smaller than 1!")
  if (length(unique(na.omit(N[N > 1]))) < 9) stop ("Too few years of data!")
  if (is.null(years)) {
    years <- 1:length(N)
    cat("\nNB: Years were not provided and are simply numbered from 1 to" %+%
        length(N) %+% ".\n\n")
  }
  if (length(years) %=% 1) {
    years <- years + (1:length(N)) - 1
    cat("\nNB: Years are assumed to be " %+% min(years) %+% "-" %+% 
                                             max(years) %+% ".\n\n")
  }
  if (!(years %=% (min(years):max(years)))) {
    cat("\nNB: The years provided are not consecutive. Is that by purpose?\n\n")
    NAs <- cbind(years, N)
    NAs <- merge(data.frame(years = min(years):max(years)), NAs, all.x = T)
    years <- NAs[, 1]
    N <- as.vector(NAs[, 2])
  }

  # Definition of required variables -------------------------------------------
  X <- ln(N)
  bootres <- psim <- NULL
  predicted <- residual <- X
  p.res.t <- p.res.N <- p.white <- 1
  names(N) <- years
  X.orig <- X
  par <- c(pi / 100, -2)

  set.seed(seed) # initiates the random seed generator

  
  ##############################################################################
  # Parameter estimation

  estimat <- estparam(par)
  my  <- se <- 0
  my  <-     estimat$par[1]
  se  <- exp(estimat$par[2])
  lnL     <-    -estimat$value
  if (estimat$convergence) {
    cat("\nNB: The optimalisation did not converge!\n\n")
  }
  
  predicted <- c(NA, X + my - 0.5 * sd / N)
  residual <-  c(X, NA) - predicted
  residual <- residual[-1]
  predicted <- predicted[-length(predicted)]
  names(predicted) <- names(residual) <- years

  r2 <- NA
  pred.r <- my - 0.5 * sd / N
  r2 <- summary(lm(diff(X) ~ pred.r[-length(X)]))$r.sq

  
  ##############################################################################
  # Parametric bootstrapping

  # These simulations are needed to estimate the uncertainties of the parameters
  # (population growth rate and environmental variance)
  if (nboot > 0) {
    boottab <- matrix(NA, nboot, 2)
    # Simulate population trajectories based on the parameter estimates:
    xsim <- simdat(1:nboot)
    # If the simulated trajectories approach zero too soon, it become inmpossible
    # to estimate the parameters. To avoid this situation, new simulations are done:
    nuller   <- apply(xsim > 0, 1, sum, na.rm = T) < 6
    while(any(nuller)) {
      xsim[which(nuller), ] <- simdat(which(nuller))
      nuller <- apply(xsim > 0, 1, sum, na.rm = T) < 6
    }

    for (j in 1:nboot) {
      # Re-estimate the parameters for each of the simulated trajectories:
      X <- xsim[j, ]
      estimat <- estparam(par)
      boottab[j, ] <- estimat$par
      cat("Bootstrapping: " %+% floor(100 * j / nboot) %+% "% done.\r")
    }
    cat("\n")
    boottab[, 2] <- exp(boottab[, 2])
    bootres <- matrix(NA, length(qntl) + 1, 2)
    rownames(bootres) <- c((qntl * 100) %+% "-percentile", "mean")
    colnames(bootres) <- c("r", "sigma2e")
    for (j in 1:ncol(bootres)) {
      bootres[1:length(qntl), j] <- quantile(boottab[, j], qntl, na.rm = T)
      bootres[1+length(qntl), j] <- mean(boottab[, j], na.rm = T)
    }
  }
  X <- X.orig

  
  ##############################################################################
  # Model testing

  # Results of these tests are displayed only if a model assumption may be violated: 
  lmfit <- lm(residual ~ years, na.action = "na.exclude")
  p.res.t <- signif(anova(lmfit)[1,5],3)
  if (p.res.t < 0.05) {
    cat("\nNB: The residuals have a " %+% ifelse(p.res.t < 0.01, "marked ", "") %+% 
        ifelse(coef(lmfit)[[2]] > 0, "posi", "nega") %+% "tive trend (p = " %+%
        format(p.res.t, dig=2, sci=F, dec=",") %+% ").\n")
  }
  lmfit <- lm(residual ~ N, na.action = "na.exclude")
  p.res.N <- signif(anova(lmfit)[1, 5], 3)
  tetth <- coef(lmfit)[[2]] < 0
  if (p.res.N < 0.05) {
    cat("\nNB: The residuals are " %+% ifelse(p.res.t < 0.01, "strongly ", "") %+% 
        ifelse(tetth, "nega", "posi") %+% "tively correlated with population " %+%
        "size (p = " %+% format(p.res.N, dig=2, sci=F, dec=",") %+% ").\n")
    if (tetth) {
      cat("This indicates that a density-dependent model " %+% ifelse(
        p.res.t < 0.01, "is", "might be") %+% " a better choice!\n")
    }
  }
  p.white <- portman.Q(residual)[[3]]
  if (p.white < 0.05) {
    cat("\nNB: The residuals deviate " %+% ifelse(p.white<0.01,"strongly ","") %+% 
        "from the assumption of whiteness (p = " %+% 
        format(p.white, dig=2, sci=F, dec=".") %+% ").\n")
  }

  
  ##############################################################################
  # Output
  
  cat("\nBrownian population model:")
  cat("\n==========================\n\n")
  cat("* Population growth rate (r):\n\n")
  print(bootres[1:length(qntl), 1])
  cat("\n* Environmental variance:\n\n")
  print(bootres[1:length(qntl), 2])
  cat("\n")
  
  n.par <- 2
  n.eff <- length(na.omit(predicted))
  aic   <- 2 * (n.par - lnL)
  aicc  <- aic + 2 * n.par * (n.par + 1) / (n.eff - n.par - 1)

  res <- list(
    "K" = NA,                 # carrying capacity (NA for Brownian models)
    "s" = my,                 # estimate of intrinsic population growth rate
    "sig2e" = se,             # estimate of environmental variance
    "sig2d" = sd,             # assumed demographic variance
    "years" = years,          # years of observation
    "N" = N,                  # observed  population sizes
    "predicted" = predicted,  # predicted population sizes
    "residuals" = residual,   # residual  population sizes
    "bootout" = boottab,      # matrix with simulated values for each parameter
    "boot" = bootres,         # matrix of quantiles for parameters estimated
    "model" = "Brownian",     # type of population model
    "name" = name,            # name of dataset
    "p.res.t" = p.res.t,      # p-value of correlation between residuals and year
    "p.res.N" = p.res.N,      # p-value of correlation between residuals and N
    "p.white" = p.white,      # p-value of correlation of test for whiteness
    "lnL" = lnL,              # log-likelihood of the model
    "r2" = r2,                # variance explained by the model (R^2)
    "aic" = aic,              # Akaike's Information Criterion (AIC) of the model
    "aicc" = aicc             # AIC corrected for small sample size
  )

  # The output is a list (see above)
  invisible(res)
}





