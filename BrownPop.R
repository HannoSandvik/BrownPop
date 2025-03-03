##########
# BrownPop
##########

# Brownian population model



BrownPop <- function(pop.size,       # time series with population sizes
                     years = NULL,   # years of pop.size
                     nboot = 10000,  # number of simulations - to be fixated!
                     name = ""       # name of dataset
                     ) {
  
  ##############################################################################
  # Preparations
  
  # Constants ------------------------------------------------------------------
  # (These are former function arguments that have been fixated)
  sd <- 0.1         # demographic variance
  qntl <- c(0.05, 0.25, 0.5, 0.75, 0.95) # quantiles to be used in the output
  seed <- NULL      # seed for the random number generator

  
  # Functions ------------------------------------------------------------------

  # Natural logarithm should be abbreviated as "ln"!
  ln <- function(x) log(x)
  
  
  # Decadal logarithm should be abbreviated as "lg"!
  lg <- function(x) log10(x)
  
  
  # Tests equality of arguments - insensitive to rounding error!
  "%=%" <- function(arg1, arg2) { 
    attributes(arg1) <- NULL
    attributes(arg2) <- NULL
    return(identical(all.equal(arg1, arg2), TRUE))
  }
  
  
  # Tests inequality of arguments - insensitive to rounding error!
  "%!=%" <- function(arg1, arg2) !(arg1 %=% arg2)
  
  
  # Combines text variables into one string
  "%+%" <- function(string1, string2) paste0(string1, string2)
  
  
  # Pretty formatting of integers #¤
  heltall <- function(x, m = 6, n = 4) {
    if (is.numeric(x)) x <- round(x)
    x <-  as.character(x)
    while(nchar(x) < m)           x <- " " %+%  x
    while(nchar(x) < m + 1 + n)   x <-  x  %+% " "
    return(x)
  }
  
  
  # Pretty formatting of decimal numbers #¤
  destall <- function(x, m = 6, n = 4, z = 3, des = ",") {
    x <- round(as.numeric(x), z)
    z <- min(z, n)
    x <- unlist(strsplit(as.character(x), ".", fixed = T))
    en <- x[1]
    while(nchar(en) < m) en <- " " %+% en
    to <- x[2]
    if (is.na(to)) {
      if (z %=% 0) {
        to <- ""
      } else {
        to <- des %+% paste(rep("0", z), collapse = "")
      }
    } else {
      while(nchar(to) < z)   to <- to %+% "0"
      to <- des %+% to
    }
    while(nchar(to) < n + 1) to <- to %+% " "
    return(en %+% to)
  }
  
  
  # Log-likelihood function of normal distribution
  log.gauss <- function(X, my, sig2) -0.5*log(2*pi*sig2)-0.5*((X-my)^2)/sig2
  

  lik <- function(par) { # likelihood function
    L <- 0
    m  <- par[1]
    s  <- exp(par[2])
    w  <- log(exp(X))
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


  simdat  <- function(i) { #¤
    xsim <- matrix(NA, length(i), length(X.orig))
    xsim[, 1] <- X.orig[1]
    for (j in 2:length(X.orig)) {
      xsim[, j] <- xsim[, j - 1] - 0.5 * sd * exp(-xsim[, j - 1]) + my -
        sqrt(se + sd * exp(-xsim[, j - 1])) * rnorm(length(i))
      xsim[xsim[, j] < 0, j] <- 0
    }
    xsim[, is.na(X.orig)] <- NA
    return(xsim)
  }
  
  
  portman.Q <- function (X, K = 10) { #¤
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
  
  
  na.und <- function(x) ifelse(is.na(x), 1, ifelse(x < C, 1, x)) #¤ eller C?!
  
  
  #risk <- function(m) 100 - c((sum(m[, 10] > C)) * 100 / nrow(m),
  #                            (sum(m[, 20] > C)) * 100 / nrow(m), 
  #                            (sum(m[,100] > C)) * 100 / nrow(m))
  # utdøingssannsynlighet - brukes ikke per nå ¤
  

  # Handling of input ----------------------------------------------------------
  N <- pop.size #¤?!
  while (is.na(N[1]))         N <- N[-1]
  while (is.na(N[length(N)])) N <- N[-length(N)]
  
  if (any(N < 1, na.rm = T)) stop ("Bestandsstørrelser kan ikke være mindre enn 1!") #¤
  if (length(unique(na.omit(N[N > 1]))) < 9) stop ("For få år med data!") #¤
  if (is.null(years)) stop(
    "timepoints of data must be specified, if no NA's, try years=1:length(N)") #¤
  if (length(years) %=% 1) {
    years <- years + (1:length(N)) - 1
    cat("\nOBS: Årstallene er satt til " %+% min(years) %+% "-" %+% 
                                             max(years) %+% ".\n\n")
  }
  if (years %!=% (min(years):max(years))) {
    cat("\nOBS: Årstallene er ikke sammenhengende. Er det meningen?\n\n") #¤
    NAs <- cbind(years, N)
    NAs <- merge(data.frame(years = min(years):max(years)), NAs, all.x = T)
    years <- NAs[, 1]
    N <- as.vector(NAs[, 2])
  }

  # Definition of required variables -------------------------------------------
  X <- ln(N)
  last <- length(N)
  lastN <- N[last] #¤ + bruk den!
  ###x   <- ln(N) #¤ + trengs ikke??!
  bootres <- psim <- NULL
  predicted <- residual <- demcomp <- X
  pred.r <- lci <- uci <- X
  p.res.t <- p.res.N <- p.white <- 1
  names(N) <- years
  NAs <- is.na(X) #¤ er ikke i bruk!

  X.orig <- X
  par <- c(pi / 100, -2)

  set.seed(seed)

  ##############################################################################
  # Parameter estimation

  estimat <- estparam(par)
  my  <- se <- 0
  my  <-     estimat$par[1]
  se  <- exp(estimat$par[2])
  lnL     <-    -estimat$value
  if (estimat$convergence) { #¤
    cat("\n##########################################" %+%
        "\nOBS: Optimaliseringa konvergerte ikke!!!!!" %+%
        "\n##########################################\n")
  }
  
  ww <- log(exp(X)) #¤ = X?!
  predicted <- c(NA, ww + my - 0.5 * sd * exp(-ww))
  residual <-  c(X, NA) - predicted
  residual <- residual[-1]
  predicted <- predicted[-length(predicted)]
  demcomp <- sd / (exp(X)) #¤ N
#  res <- data.frame(years     = years,
#                    residual  = residual,
#                    demcomp   = demcomp, 
#                    predicted = predicted)
#  residual <- res$residual
#  demcomp <- res$demcomp
#  predicted <- res$predicted
  names(predicted) <- names(demcomp) <- names(residual) <- years

  # ----------------------------------------------------------------------------
  
  r2 <- NA
  pred.r <- my - 0.5 * sd * exp(-X)
  r2 <- summary(lm(diff(X) ~ pred.r[-length(X)]))$r.sq
  uci <- lci <- NA
  for (j in 1:(length(pred.r))) {
    lci[j] <- pred.r[j] + sqrt(se + sd * exp(-X[j])) * qnorm(0.025)
    uci[j] <- pred.r[j] + sqrt(se + sd * exp(-X[j])) * qnorm(0.975)
  }
  names(pred.r) <- names(lci) <- names(uci) <- years

  ##############################################################################
  # parametric bootstrapping

  if (nboot > 0) {
    cat("bootstrapping\n")
    catval <- c(1, seq(0, nboot, nboot / 10))
    boottab <- matrix(NA, nboot, 2)
    
    xsim <- simdat(1:nboot)
    # Om de simulerte forløpene går altfor fort mot null, blir det umulig å esti-
    # mere parameterne. For å forebygge krøll, erstattes de av nye simuleringer:
    nuller   <- apply(xsim > 0, 1, sum, na.rm = T) < 6
    while(any(nuller)) {
      xsim[which(nuller), ] <- simdat(which(nuller))
      nuller <- apply(xsim > 0, 1, sum, na.rm = T) < 6
    }

    for (j in 1:nboot) {
      if (any(j == catval)) cat("boot", j, "of", nboot, "\n")
      X <- xsim[j, ]
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
    } #¤ nei, bruk qntl!
  }
  X <- X.orig

  ##############################################################################
  # Model testing

  lmfit <- lm(residual ~ years, na.action = "na.exclude")
  p.res.t <- signif(anova(lmfit)[1,5],3)
  if (p.res.t < 0.05) {
    cat("\nOBS: Residualene har en " %+% ifelse(p.res.t < 0.01, "utprega ", "") %+% 
        ifelse(coef(lmfit)[[2]] > 0, "posi", "nega") %+% "tiv trend (p = " %+%
        format(p.res.t, dig=2, sci=F, dec=",") %+% ").\n")
  }
  lmfit <- lm(residual ~ N, na.action = "na.exclude")
  p.res.N <- signif(anova(lmfit)[1, 5], 3)
  tetth <- coef(lmfit)[[2]] < 0
  if (p.res.N < 0.05) {
    cat("\nOBS: Residualene er " %+% ifelse(p.res.t < 0.01, "sterkt ", "") %+% 
        ifelse(tetth, "nega", "posi") %+% "tivt korrelert med bestandsstø" %+%
        "rrelsen (p = " %+% format(p.res.N, dig=2, sci=F, dec=",") %+% ").\n")
    if (tetth) {
      cat("Det indikerer at en tetthetsavhengig modell " %+% ifelse(
        p.res.t < 0.01, "er et bedre valg!\n", "kan være et bedre valg.\n"))
    }
  }
  p.white <- portman.Q(residual)[[3]]
  if (p.white < 0.05) {
    cat("\nOBS: Residualene avviker " %+% ifelse(p.white < 0.01, "sterkt ", "") %+% 
        "fra en antagelse om hvit støy (p = " %+% 
        format(p.white, dig=2, sci=F, dec=",") %+% ").\n")
  }
  #¤ merk at her brukes komma som desimaltegn

  ##############################################################################
  # ¤¤ forhenværende avslutning av modellfunksjonen - må ryddes!
  
  n.par <- 2
  n.eff <- length(na.omit(predicted))
  aic   <- 2 * (n.par - lnL) #¤ enten vise eller droppe!
  aicc  <- aic + 2 * n.par * (n.par + 1) / (n.eff - n.par - 1)

  res <- list( #¤ den må gjenngås grundig!
    "K" = NA,
    "s" = my,
    "sig2e" = se,
    "sig2d" = sd,
    "years" = years,
    "N" = N,
    "residuals" = residual, "predicted" = predicted, "demcomp" = demcomp,
    "pred.r" = pred.r, "lci" = lci, "uci" = uci,
    "bootout" = boottab, "boot" = bootres,
    "model" = "Brownian", "name" = name,
    "p.res.t" = p.res.t, "p.res.N" = p.res.N, "p.white" = p.white,
    "lnL" = lnL, "r2" = r2, "aic" = aic, "aicc" = aicc
  )

  cat("\n                          Persentiler:")
  cat("\n                        ")
  for (i in qntl) {
    cat("  " %+% destall(i * 100, z = 1))
  }
  cat("\nPopulasjonsmodell...")
  cat("\nBestandsvekstrate (r):  ")
  for (i in c(6, 8, 9, 10, 12)) cat("  " %+% destall(    bootres[i, 1], 
                                    z = 2 - floor(lg(abs(bootres[9, 1])))))
  cat("\nMiljøvarians (sigma-e): ")
  for (i in c(6, 8, 9, 10, 12)) cat("  " %+% destall(    bootres[i, 2], 
                                    z = 2 - floor(lg(abs(bootres[9, 2])))))

  invisible(res)
}







