##########
# BrownPVA
##########

# Brownian population viability analysis



BrownPVA <- function(pop.model,       # output from BrownPop
                     crash.prop = 0,  # fraction of population that dies
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
  fig <- 100        # number of years of population trajectories shown in graph
  fil <- NULL       # omit? ¤
  linjer <- 25      # number of random population trajectories shown in graph
  
  
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
  
  
  # Pretty formatting of integers
  heltall <- function(x, m = 6, n = 4) {
    if (is.numeric(x)) x <- round(x)
    x <-  as.character(x)
    while(nchar(x) < m)           x <- " " %+%  x
    while(nchar(x) < m + 1 + n)   x <-  x  %+% " "
    return(x)
  }
  
  
  # Pretty formatting of decimal numbers
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
  
  
  na.und <- function(x) ifelse(is.na(x), 1, ifelse(x < C, 1, x)) #¤ or C?!
  
  
  #risk <- function(m) 100 - c((sum(m[, 10] > C)) * 100 / nrow(m),
  #                            (sum(m[, 20] > C)) * 100 / nrow(m), 
  #                            (sum(m[,100] > C)) * 100 / nrow(m))
  # utdøingssannsynlighet - brukes ikke per nå ¤
  
  
  # Handling of input ----------------------------------------------------------
  N       <- pop.model$N
  years   <- pop.model$years
  sd      <- pop.model$sig2d
  se      <- pop.model$sig2e
  boottab <- pop.model$bootout
  name    <- pop.model$name
  nboot <- nrow(boottab)

  # Definition of required variables -------------------------------------------
  X <- ln(N)
  last <- length(N)
  lastN <- N[last] #¤ + bruk den!
  lnC <- ln(C)
  ext <- rec <- NA # extinction & recovery
  
  
  # sik <- list() #¤ fjernes om det ikke skal lagres noen resultater!
  
  # Graph ----------------------------------------------------------------------
  
  graphics.off()
  #  plot(years[last - 30:0], X[last - 30:0], 
  #       xlim=c(years[last] - 40, years[last] + fig),
  #       ylim=c(lnC, max(1.618 * (X[last] - lnC) + lnC, X[last - 30:0], na.rm=T)),
  #       type="l", lwd=2, xlab="Years", ylab="Ln population size", ...)
  #  plot(years, N, 
  #       xlim = c(years[1], years[last] + fig),
  #       ylim = c(0, 1.1 * max(1.1 * lastN, N, na.rm = T)), yaxs = "i",
  #       type = "l", lwd = 2, xlab = "År", ylab = "Bestandsstørrelse") # ...
  title <- name %+% ", "
  if (!nchar(name)) title <- "Population viability analysis, "
  title <- title %+% ifelse(crash.prop > 0, (crash.prop * 100) %+% "%", "no")
  title <- title %+% " population crash"
  plot(years, X, main = title,
       xlim = c(years[1], years[last] + fig),
       ylim = c(lnC, max(1.2 * (X[last] - lnC) + lnC, X, na.rm=T)), yaxt = "n",
       type = "l", lwd = 2, xlab = "Year", ylab = "Population size") # ...¤
  if (exp(max(1.2 * (X[last]-lnC) + lnC, X, na.rm=T)) >= 10^(ceiling(lg(C)) + 2)) {
    lab <- sort(unique(c(C, 10^(2:12))))
  } else {
    lab <- c(1, 2, 5) * 10^rep(0:6, each = 3)
  }
  for (i in 1:9) {
    axis(2, ln(i * 10^(0:12)), F, T, tcl = i / 20 - 0.55)
  }
  axis(2, ln(lab), lab, T, F)
  
  ### Density-independent model ¤ ----------------------------------------------
  x <- rep(lastN, nboot)
  stikkpr <- 1:min(linjer, nboot)
  utvalg  <- matrix(0, linjer, tmax)
  kvantiler <- matrix(0, length(qntl), tmax)
  #  est <- matrix(0, nrow(model$bootout), 5)
  #  est[, 1:ncol(model$bootout)] <- model$bootout #¤ dette går enklere! est <- bootout! (eller boottab!)
  x <- round(x * (1 - crash.prop))
  for (t in 1:tmax) {
    alive <- x > C
    x <- ln(x)
    if (any(alive)) {
      x[alive] <- round(exp(
        x[alive] - 0.5 * sd * exp(-x[alive]) + boottab[alive, 1] +
          sqrt(boottab[alive, 2] + sd * exp(-x[alive])) * rnorm(sum(alive))))
    } else {
      x <- rep(0, nboot)
    }
    x[x <= C] <- 0
    utvalg[,    t] <- x[stikkpr]
    kvantiler[, t] <- quantile(x, qntl)
    ext[t] <- sum(x <      C) # n of trajectories that  are  extinct  by time t
    rec[t] <- sum(x >= lastN) # n of trajectories that have recovered by time t
  }
  
  y1 <- ln(round(lastN * (1 - crash.prop)))
  for (i in stikkpr) {
    lines(years[last] + 0:tmax, c(y1, ln(na.und(utvalg[i, ]))), col = grey(0.72))
  }
  if (crash.prop > 0) {
    lines(rep(years[last], 2), c(X[last], y1), lwd = 2, col = "red")
  }
  for (i in 1:length(qntl)) {
    lines(years[last] + 0:tmax, c(y1, ln(na.und(kvantiler[i, ]))))
  }
  
  # Utmating -------------------------------------------------------------------
  utmating <- matrix("", length(qntl), 2)
  colnames(utmating) <- c("Levetid", "Restitueringstid")
  rownames(utmating) <- (qntl * 100) %+% "-persentil"
  for (i in 1:length(qntl)) {
    if (any(ext / nboot >= max(qntl[i], 0.000001))) {
      lt <- min(which(ext / nboot >= max(qntl[i], 0.000001)))
    } else {
      lt <- "> " %+% format(tmax, sci = F)
    }
    if (any(rec / nboot >= max(qntl[i], 0.000001))) {
      rt <- min(which(rec / nboot >= max(qntl[i], 0.000001)))
    } else {
      rt <- "> " %+% format(tmax, sci = F)
    }
    utmating[i, ] <- c(lt, rt)
  }
  #  print.table(utmating, right = T)
  cat("\nLevedyktighetsanalyse...")
  cat("\n- Tid til utdøing  (år):")
  for (i in 1:length(qntl)) cat("  " %+% heltall(utmating[i, 1]))
  cat("\n- Restitueringstid (år):")
  for (i in 1:length(qntl)) cat("  " %+% heltall(utmating[i, 2]))
  
  
  
  
  #    for (i in qntl) {
  #      if (any(ext / n >= max(i, 0.000001))) {
  #        lt <- min(which(ext / n >= max(i, 0.000001)))
  #      } else {
  #        lt <- ">1000"
  #      }
  #      cat("Levetidas " %+% (i * 100) %+% "-persentil: " %+% lt %+% " år\n")
  #    }
  #
  #    if (is.null(fil)) {
  #      j <- 1
  #      fil <- getwd() %+% "/pva" %+% j %+% ".R"
  #      while (file.exists(fil)) {
  #        j <- j + 1
  #        fil <- getwd() %+% "/pva" %+% j %+% ".R"
  #      }
  #    }
  #    for (q in c(0.05, 0.1, 0.2, 0.8, 0.9, 0.95)) {
  #      lines(years[last] - 1 + 1:fig, na.und(apply(m[,1:fig], 2, quantile, q)))
  #    }
  #    lines(years[last] - 1 + 1:fig, na.und(apply(m[, 1:fig], 2, quantile, 0.5)), 
  #          lwd = 2)
  #    save(sik, m, ext, file=fil)
  #    cat("\nThe PVA results have been saved as \"" %+% fil %+% "\".\n\n")
  #    extR <- matrix(NA, 13, 4, dimnames=list(
  #      c("extinction.risk", "expected.pop.lifetime",
  #        c(50, 20, 80, 10, 90, 5, 95, 2.5, 97.5, 0.5, 99.5) %+% ".pct.pop.lifetime"),
  #      c("10a", "20a", "100a", "exact")))
  #    extR[1, 1:3] <- risk(m)
  #    extR[2, 1:3] <- -c(10, 20, 100) / ln(1 - risk(m) / 100)
  #    if(all(ext < tmax)) extR[2,4] <- mean(ext)
  #    for (i in 3:13) {
  #      q <- c(0.5, 0.2, 0.8, 0.1, 0.9, 0.05, 0.95, 0.025, 0.975, 0.005, 0.995)[i - 2]
  #      extR[i, 1:3] <- -extR[2, 1:3] * ln(1 - q)
  #      if (sum(ext < tmax) > (n * q)) extR[i,4] <- quantile(ext, q)
  #    }
  #    cat("Extinction risk (percent within n years):\n")
  #    er <- as.vector(extR[1, 1:3])
  #    names(er) <- c("10", "20", "100")
  #    print(er)
  #    cat("\nPopulation lifetime (quantiles in years):\n")
  #    pl <- round(as.vector(extR[3:13, 4]))
  #    if (any(is.na(pl))) pl[is.na(pl)] <- ">" %+% tmax
  #    names(pl) <- c(50, 20, 80, 10, 90, 5, 95, 2.5, 97.5, 0.5, 99.5) %+% "%"
  #    print(pl, quote=F)
  #    if (crash.prop > 0) {
  #      return <- NA
  #      cat("\nReturn to the pre-crash population size (quantiles in years):\n")
  #      for (i in 1:length(sik)) {
  #        x  <- sik[[i]][-1]
  #        x1 <- sik[[i]][1]
  #        if (all(x < x1)) {
  #          return[i] <- tmax
  #        } else {
  #          return[i] <- which(x >= x1)[1]
  #        }
  #      }
  #      q <- round(quantile(return, c(0.5, 0.2, 0.8, 0.1, 0.9,
  #                                    0.05, 0.95, 0.025, 0.975, 0.005, 0.995)))
  #      if (any(q == tmax)) q[which(q == tmax)] <- ">" %+% tmax
  #      print(q, quote=F)
  #    }
  #    ext <- extR
  
  invisible(utmating)
}






