
########
BrownPVA <- function(
########  
                     pop.model,      # output from BrownPop
                     mortality = 0,  # proportion of population that dies in crash
                     C = 20          # quasi-extinction threshold - to be fixated!
                     ) {
  
  # Brownian population viability analysis
  # ======================================
  # Code written by Hanno Sandvik, Norwegian Institute for Nature Research, 2025
  
  
  ##############################################################################
  # Preparations
  
  # Constants ------------------------------------------------------------------
  # (including former function arguments that have been fixated)
  qntl <- c(0.05, 0.25, 0.5, 0.75, 0.95) # quantiles to be used in the output
  mort.year <- 0   # in which years of simulations does the collapse happen ¤ needed??!
  # C <- 20          # quasi-extinction threshold
  tmax <- 1000     # maximum number of years to simulate population trajectories
  tfig <- 100      # number of years of population trajectories shown in graph
  nlines <- 25     # number of random population trajectories shown in graph
  
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
  
  
  # Combines text variables into one string
  "%+%" <- function(string1, string2) paste0(string1, string2)
  
  
  # replaces population sizes that are smaller than C or missing by 1
  na.und <- function(x) ifelse(is.na(x), 1, ifelse(x < C, 1, x))
  
  
  # Handling of input ----------------------------------------------------------
  N         <- pop.model$N        # time series with population sizes
  years     <- pop.model$years    # years of N
  sd        <- pop.model$sig2d    # demographic variance
  boottab   <- pop.model$bootout  # bootstrapped values of r and env. variance
  name      <- pop.model$name     # name of dataset
  nboot     <- nrow(boottab)      # number of simulations
  mortality <- sort(mortality)    # sorts mortalities in ascending order
  nmort     <- length(mortality)  # number of mortalities to consider
  
  # Definition of required variables -------------------------------------------
  X <- ln(N)          # vector of real log-population sizes
  last <- length(N)   # (number of) last year with observation
  lastN <- N[last]    # last observed population size
  lnC <- ln(C)        # logarithm og quasi-extinction threshold
  ext <- rec <- array(NA, c(10, nmort, tmax))
  
  # Graph ----------------------------------------------------------------------
  graphics.off()
  title <- name %+% ", "
  if (!nchar(name)) title <- "Population viability analysis, "
  if (nmort > 1) {
    title <- title %+% "mass mortality " %+% (min(mortality) * 100) %+%
             "% to " %+% (max(mortality) * 100) %+% "%"
  } else {
    title <- title %+% ifelse(mortality > 0, 
                             (mortality * 100) %+% "%", "no")
    title <- title %+% " mass mortality"
  }
  plot(years, X, main = title,
       xlim = c(years[1], years[last] + tfig),
       ylim = c(lnC, max(1.2 * (X[last] - lnC) + lnC, X, na.rm=T)), yaxt = "n",
       type = "l", lwd = 2, xlab = "Year", ylab = "Population size")
  if (exp(max(1.2 * (X[last]-lnC) + lnC, X, na.rm=T)) >= 10^(ceiling(lg(C)) + 2)) {
    lab <- sort(unique(c(C, 10^(2:12))))
  } else {
    lab <- c(1, 2, 5) * 10^rep(0:6, each = 3)
  }
  for (i in 1:9) {
    axis(2, ln(i * 10^(0:12)), F, T, tcl = i / 20 - 0.55)
  }
  axis(2, ln(lab), lab, T, F)
  
  
  ##############################################################################
  # Simulation of future population trajectories
  
  x <- matrix(lastN, nmort, nboot)
  # (Starting point of simulations is the last observed population size)
  sample <- 1:min(nlines, nboot)
  selection <- array(0, c(nmort, nlines,       tmax))
  quantiles <- array(0, c(nmort, length(qntl), tmax))
  x[] <- round(x * rep(1 - mortality, nboot))
  # (The mass mortality event is applied to the last observed population size)
  boottab <- array(rep(boottab, each = nmort), c(nmort, dim(boottab)))
  for (t in 1:tmax) {
    alive <- x > C
    x <- ln(x)
    if (any(alive)) {
      # A stochastic Brownian model is used to forecast population sizes for the
      # following year:
      x <- round(exp(x - 0.5 * sd * exp(-x) + boottab[, , 1] + 
                     sqrt(boottab[, , 2] + sd * exp(-x)) * 
                     rep(rnorm(nboot), each = nmort)))
    } else {
      # If all simulated trajectories are extinct, things are easier:
      x <- rep(0, nboot)
    }
    x[is.na(x)] <- 0
    x[x <= C]   <- 0
    selection[, , t] <- x[, sample]
    quantiles[, , t] <- t(apply(x, 1, quantile, qntl))
    for (Z in 1:10) {
      # Dividing the simulations into 10 subsets allows estimates to be based
      # on means with 1 decimal.
      # Number of trajectories that are extinct or that have recovered:
      subset <- 1:(nboot / 10) + (Z - 1) * nboot / 10
      ext[Z, , t] <- apply(x[, subset, drop = FALSE] <      C, 1, sum)
      rec[Z, , t] <- apply(x[, subset, drop = FALSE] >= lastN, 1, sum)
    }
  }

  # Plot the results -----------------------------------------------------------
  y1 <- ln(round(lastN * (1 - mortality[nmort])))
  if (nmort %=% 1) {
    # Illustrate a sample of simulated trajectories (only shown if there is 
    # exactly one mortality):
    for (i in sample) {
      lines(years[last] + 0:tmax, c(y1, ln(na.und(selection[nmort, i, ]))), 
            col = grey(0.72))
    }
  }
  # The mass mortality will be shown as a vertical red line:
  if (max(mortality) > 0) {
    lines(rep(years[last], 2), c(X[last], y1), lwd = 2, col = "red")
  }
  # Illustrate the quantiles of the simulations (as black lines if there is
  # exactly one mortality, otherwise with blue lines for the lowest mortality
  # and with red lines for the highest mortality)
  for (i in 1:length(qntl)) {
    lines(years[last] + 0:tmax, c(y1, ln(na.und(quantiles[nmort, i, ]))), 
          col=ifelse(nmort > 1, "red", "black"))
  }
  if (nmort > 1) {
    for (i in 1:length(qntl)) {
      lines(years[last] + 0:tmax, c(ln(lastN), ln(na.und(quantiles[1, i, ]))), 
            col = "blue")
    }
  }

  
  ##############################################################################
  # Output
  
  cat("\nPopulation viability analysis:")
  cat("\n==============================\n\n")
  EXT <- REC <- array(0, c(10, nmort, length(qntl)),
                      list(1:10,
                           "Mortality." %+% (mortality * 100) %+% "%",
                           (qntl * 100) %+% "-percentile"))
  for (Z in 1:10) { # for each of the 10 subsets of the simulations
    extinct <- matrix(0, nmort, length(qntl), FALSE,
                      list("Mortality." %+% (mortality * 100) %+% "%",
                           (qntl * 100) %+% "-percentile"))
    recover <- matrix(0, nmort, length(qntl), FALSE,
                      list("Mortality." %+% (mortality * 100) %+% "%",
                           (qntl * 100) %+% "-percentile"))
    for (i in 1:length(qntl)) {
      for (j in 1:nmort) {
        if (any(ext[Z, j, ] / nboot * 10 >= max(qntl[i], 0.000001))) {
          # identifying the year in which the proportion of simulated trajectories 
          # that have gone extinct equals the quantile:
          extinct[j, i] <- min(which(ext[Z, j, ] / nboot * 10 >= 
                                     max(qntl[i], 0.000001)))
        } else {
          extinct[j, i] <- tmax
        }
        if (any(rec[Z, j, ] / nboot * 10 >= max(qntl[i], 0.000001))) {
          # identifying the year in which the proportion of simulated trajectories 
          # that have recovered equals the quantile:
          recover[j, i] <- min(which(rec[Z, j, ] / nboot * 10 >=
                                     max(qntl[i], 0.000001)))
        } else {
          recover[j, i] <- tmax
        }
      }
    }
    EXT[Z, , ] <- extinct
    REC[Z, , ] <- recover
  }
  extinct <- apply(EXT, 2:3, mean)  # mean across subsets
  recover <- apply(REC, 2:3, mean)  # mean across subsets
  extinct[apply(EXT, 2:3, max) == 1000] <- 1000  # to avoid spurious averages
  recover[apply(REC, 2:3, max) == 1000] <- 1000  # to avoid spurious averages

  cat("* Population lifetime (time to quasi-exinction in years):\n\n")
  print(extinct)
  cat("\n* Recovery time (time to population recovery in years):\n\n")
  print(recover)
  cat("\n")
  if (any(extinct == tmax) | any(recover == tmax)) {
    cat("(Note: Values given as " %+% tmax %+% " may be much larger.)\n\n")
  }
  Q <- list()
  for (i in 1:nmort) {
    Q[[i]] <- quantiles[i, , c(1, 1:tmax)]
    Q[[i]][, 1] <- round(lastN * (1 - mortality[i]))
    dimnames(Q[[i]]) <- list((100 * qntl) %+% "-percentile", years[last] + 0:tmax)
  }
  names(Q) <- "Mortality." %+% (mortality * 100) %+% "%"

  # The output is a list consisting of three objects:
  # 1. a matrix of the quantiles of population lifetime for each mortality 
  #    (identical to the table of population lifetimes printed on screen)
  # 2. a matrix of the quantiles of recovery time for each mortality 
  #    (identical to the table of recovery times printed on screen)
  # 3. a list of matrices, one matrix per mortality, providing the quantiles of the
  #    simulated population sizes for each year (corresponding to the red and blue
  #    lines in the graph; note that these matrices are rather big [1001 columns],
  #    so printing them in their entirety is normally inexpedient)
  invisible(list(PopulationLifetime = extinct,
                 RecoveryTime = recover,
                 Trajectories = Q))
}



