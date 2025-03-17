##########
# BrownPVA
##########

# Brownian population viability analysis



BrownPVA <- function(pop.model,       # output from BrownPop
                     crash.prop = 0,  # proportion of population that dies in crash
                     C = 20           # quasi-extinction threshold - to be fixated!
                     ) {
  
  ##############################################################################
  # Preparations
  
  # Constants ------------------------------------------------------------------
  # (including former function arguments that have been fixated)
  qntl <- c(0.05, 0.25, 0.5, 0.75, 0.95) # quantiles to be used in the output
  crash.year <- 0   # in which years of simulations does the collapse happen ¤ needed??!
  # C <- 20           # quasi-extinction threshold
  tmax <- 1000      # maximum number of years to simulate population trajectories
  tfig <- 100       # number of years of population trajectories shown in graph
  nlines <- 25      # number of random population trajectories shown in graph
  
  
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
  
  
  # Pretty formatting of integers
  intnum <- function(x, m = 6, n = 4) {
    if (is.numeric(x)) x <- round(x)
    x <-  as.character(x)
    while(nchar(x) < m)           x <- " " %+%  x
    while(nchar(x) < m + 1 + n)   x <-  x  %+% " "
    return(x)
  }
  
  
  # replaces population sizes that are smaller than C or missing by 1
  na.und <- function(x) ifelse(is.na(x), 1, ifelse(x < C, 1, x))
  
  
  # Handling of input ----------------------------------------------------------
  N       <- pop.model$N
  years   <- pop.model$years
  sd      <- pop.model$sig2d
  boottab <- pop.model$bootout
  name    <- pop.model$name
  nboot <- nrow(boottab)
  ncrash <- length(crash.prop)
  
  # Definition of required variables -------------------------------------------
  X <- ln(N)          # vector of real log-population sizes
  last <- length(N)   # (number of) last year with observation
  lastN <- N[last]    # last observed population size
  lnC <- ln(C)        # logarithm og quasi-extinction threshold
  ext <- rec <- matrix(NA, ncrash, tmax)  # extinction & recovery
  
  # Graph ----------------------------------------------------------------------
  graphics.off()
  title <- name %+% ", "
  if (!nchar(name)) title <- "Population viability analysis, "
  title <- title %+% ifelse(crash.prop[ncrash] > 0, 
                           (crash.prop[ncrash] * 100) %+% "%", "no")
  title <- title %+% " population crash"
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
  
  x <- matrix(lastN, ncrash, nboot)
  sample <- 1:min(nlines, nboot)
  selection <- array(0, c(ncrash, nlines,       tmax))
  quantiles <- array(0, c(ncrash, length(qntl), tmax))
  x[] <- round(x * rep(1 - crash.prop, nboot))
  boottab <- array(rep(boottab, each = ncrash), c(ncrash, dim(boottab)))
  for (t in 1:tmax) {
    alive <- x > C
    x <- ln(x)
    if (any(alive)) {
      x <- round(exp(x - 0.5 * sd * exp(-x) + boottab[, , 1] + 
                     sqrt(boottab[, , 2] + sd * exp(-x)) * 
                     rep(rnorm(nboot), each = ncrash)))
#    if (any(alive)) {
#      x[alive] <- round(exp(
#        x[alive] - 0.5 * sd * exp(-x[alive]) + boottab[, , 1][alive] +
#          sqrt(boottab[, , 2][alive] + sd * exp(-x[alive])) * rnorm(sum(alive))))
    } else {
      x <- rep(0, nboot)
    }
    x[is.na(x)] <- 0
    x[x <= C]   <- 0
    selection[, , t] <- x[, sample]
    quantiles[, , t] <- t(apply(x, 1, quantile, qntl))
    ext[, t] <- apply(x <      C, 1, sum)   # n of trajectories that are extinct
    rec[, t] <- apply(x >= lastN, 1, sum)   # n of trajectories that have recovered
  }
  
  y1 <- ln(round(lastN * (1 - crash.prop[ncrash])))
  for (i in sample) {
    lines(years[last] + 0:tmax, c(y1, ln(na.und(selection[ncrash, i, ]))), 
          col = grey(0.72))
  }
  if (crash.prop[ncrash] > 0) {
    lines(rep(years[last], 2), c(X[last], y1), lwd = 2, col = "red")
  }
  for (i in 1:length(qntl)) {
    lines(years[last] + 0:tmax, c(y1, ln(na.und(quantiles[ncrash, i, ]))))
  }
  

  ##############################################################################
  # Output
  
  output <- matrix("", length(qntl), 2)
  colnames(output) <- c("Pop. lifetime", "Time to recovery")
  rownames(output) <- (qntl * 100) %+% "%"
  for (i in 1:length(qntl)) {
    if (any          (ext[ncrash, ] / nboot >= max(qntl[i], 0.000001))) {
      lt <- min(which(ext[ncrash, ] / nboot >= max(qntl[i], 0.000001)))
    } else {
      lt <- "> " %+% format(tmax, sci = F)
    }
    if (any          (rec[ncrash, ] / nboot >= max(qntl[i], 0.000001))) {
      rt <- min(which(rec[ncrash, ] / nboot >= max(qntl[i], 0.000001)))
    } else {
      rt <- "> " %+% format(tmax, sci = F)
    }
    output[i, ] <- c(lt, rt)
  }
  cat("\nPopulation viability analysis...")
  cat("\n- Popul. lifetime  (a): ")
  for (i in 1:length(qntl)) cat("  " %+% intnum(output[i, 1]))
  cat("\n- Time to recovery (a): ")
  for (i in 1:length(qntl)) cat("  " %+% intnum(output[i, 2]))
  cat("\n")
  
  extinct <- matrix(0, ncrash, length(qntl), FALSE,
                    list("Crash." %+% (crash.prop * 100) %+% "%",
                         (qntl * 100) %+% "-percentile"))
  recover <- matrix(0, ncrash, length(qntl), FALSE,
                    list("Crash." %+% (crash.prop * 100) %+% "%",
                         (qntl * 100) %+% "-percentile"))
  for (i in 1:length(qntl)) {
    for (j in 1:ncrash) {
      if (any(ext[j, ] / nboot >= max(qntl[i], 0.000001))) {
        extinct[j, i] <- min(which(ext[j, ] / nboot >= max(qntl[i], 0.000001)))
      } else {
        extinct[j, i] <- tmax - 0.01
      }
      if (any(rec[j, ] / nboot >= max(qntl[i], 0.000001))) {
        recover[j, i] <- min(which(rec[j, ] / nboot >= max(qntl[i], 0.000001)))
      } else {
        recover[j, i] <- tmax - 0.01
      }
    }
  }
  
  invisible(list(PopulationLifetime = extinct, RecoveryTime = recover))
}






