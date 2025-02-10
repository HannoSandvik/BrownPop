##########
# BrownPVA
##########

# Brownsk levedyktighetsanalyse
# (Brownian population viability analysis)



BrownPVA <- function(pop.size, 
                     years = NULL, 
                     nboot = 10000,  # fikseres til slutt!
                     name = "",      # ønskelig med tittel?
                     crash.prop = 0, # fikseres?
                     C = 20          # fikseres til slutt!
                     ) {
  
  ##############################################################################
  # Forberedelser
  
  # Konstanter -----------------------------------------------------------------
  # (Dette er tidligere funksjonsargumenter som har blitt fiksert)
  sd <- 0.1         # demografisk varians
  kvntl <- c(0.025, 0.25, 0.5, 0.75, 0.975) # kvantiler som skal vises i utmatinga
                    # eller c(0.025, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.975) ?!
  seed <- NULL      # frø til slumptallgeneratoren
  crash.year <- 1   # hvilket år av simuleringene skal bestanden kollapse
  C <- 20           # kvasiutdøingsterskel (¤¤¤!)
  tmax <- 1000      # maks. antall år å modellere bestandsforløp for
  fig <- 100        # antall år med modellerte bestandsforløp som vises i figuren
  fil <- NULL       # fjernes? ¤
  linjer <- 25      # antall tilfeldige bestandsforløp som skal vises i grafen
  
  
  # Funksjoner -----------------------------------------------------------------

  # Naturlig logaritme skal forkortes som "ln"!
  ln <- function(x) log(x)
  
  
  # Dekadisk logaritme skal forkortes som "lg"!
  lg <- function(x) log10(x)
  
  
  # Tester om argumentene er like - ikke følsom for avrundingsfeil!
  "%=%" <- function(arg1, arg2) { 
    attributes(arg1) <- NULL
    attributes(arg2) <- NULL
    return(identical(all.equal(arg1, arg2), TRUE))
  }
  
  
  # Tester om argumentene er ulike - ikke følsom for avrundingsfeil!
  "%!=%" <- function(arg1, arg2) !(arg1 %=% arg2)
  
  
  # Limer sammen tekstvariabler
  "%+%" <- function(string1, string2) paste0(string1, string2)
  
  
  # Pen utmating av heltall
  heltall <- function(x, m = 6, n = 4) {
    if (is.numeric(x)) x <- round(x)
    x <-  as.character(x)
    while(nchar(x) < m)           x <- " " %+%  x
    while(nchar(x) < m + 1 + n)   x <-  x  %+% " "
    return(x)
  }
  
  
  # Pen utmating av desimaltall
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
  
  
  # Normalfordelingas log-likelihood-funksjon
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
  

  estparam <- function(par) { # parameterestimering
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
  

  # Bearbeidelse av innmatinga -------------------------------------------------
  N <- pop.size #¤?!
  while (is.na(N[1]))         N <- N[-1]
  while (is.na(N[length(N)])) N <- N[-length(N)]
  
  if (any(N < 1, na.rm = T)) stop ("Bestandsstørrelser kan ikke være mindre enn 1!")
  if (length(unique(na.omit(N[N > 1]))) < 9) stop ("For få år med data!")
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

  # Definisjon av nødvendige variabler -----------------------------------------
  X <- ln(N)
  last <- length(N) #¤ sist?!
  lastN <- N[last] #¤ + bruk den!
  ###x   <- ln(N) #¤ + trengs ikke??!
  lnC <- ln(C)
  bootres <- psim <- NULL
  predicted <- residual <- demcomp <- X
  pred.r <- lci <- uci <- X
  p.res.t <- p.res.N <- p.white <- 1
  names(N) <- years
  NAs <- is.na(X) #¤ er ikke i bruk!
  ext <- rec <- NA # extinction & recovery

  X.orig <- X
  par <- c(pi / 100, -2)

  set.seed(seed)

  ##############################################################################
  # Parameterestimering

  estimat <- estparam(par)
  my  <- se <- 0
  my  <-     estimat$par[1]
  se  <- exp(estimat$par[2])
  lnL     <-    -estimat$value
  if (estimat$convergence) {
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
  res <- data.frame(years     = years,
                    residual  = residual,
                    demcomp   = demcomp, 
                    predicted = predicted)
  residual <- res$residual
  demcomp <- res$demcomp
  predicted <- res$predicted
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
    } #¤ nei, bruk kvntl!
  }
  X <- X.orig

  ##############################################################################
  # Modelltesting

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
  aic   <- 2 * (n.par - lnL)
  aicc  <- aic + 2 * n.par * (n.par + 1) / (n.eff - n.par - 1)

  res <- list(
    "K" = NA,
    "s" = my,
    "sig2" = se,
    "sigd" = sd,
    "years" = years,
    "N" = N,
    "residuals" = residual, "predicted" = predicted, "demcomp" = demcomp,
    "pred.r" = pred.r, "lci" = lci, "uci" = uci,
    "bootout" = boottab, "boot" = bootres,
    "model" = "Brownian", "name"=name,
    "p.res.t" = p.res.t, "p.res.N" = p.res.N, "p.white" = p.white,
    "lnL" = lnL, "r2" = r2, "aic" = aic, "aicc" = aicc
  )
  cat("\nParameter estimates:\n")
  partab <- matrix(NA, 2, 3)
  partab  [1,1]   <- res$s
  partab  [1,2:3] <- res$boot[c("2.5%", "97.5%"), "r"]
  partab  [2,1]   <- res$sig2
  partab  [2,2:3] <- res$boot[c("2.5%", "97.5%"), "sig2e"]
  rownames(partab) <- c("Popul. growth rate (r):", "Environmental variance:")
  colnames(partab) <- c("Mean", "Lower 95% CI", "Upper 95% CI")
  print(partab)
  
  cat("\n                          Persentiler:")
  cat("\n                        ")
  for (i in kvntl) {
    cat("  " %+% destall(i * 100, z = 1))
  }
  cat("\nPopulasjonsmodell...")
  cat("\nBestandsvekstrate (r):  ")
  for (i in c(6, 8, 9, 10, 12)) cat("  " %+% destall(    bootres[i, 1], 
                                    z = 2 - floor(lg(abs(bootres[9, 1])))))
  cat("\nMiljøvarians (sigma-e): ")
  for (i in c(6, 8, 9, 10, 12)) cat("  " %+% destall(    bootres[i, 2], 
                                    z = 2 - floor(lg(abs(bootres[9, 2])))))

  ##############################################################################
  # PVA #¤ @@@
  
  sik <- list() #¤ fjernes om det ikke skal lagres noen resultater!
  
  # Graf -----------------------------------------------------------------------
  
  graphics.off()
#  plot(years[last - 30:0], X[last - 30:0], 
#       xlim=c(years[last] - 40, years[last] + fig),
#       ylim=c(lnC, max(1.618 * (X[last] - lnC) + lnC, X[last - 30:0], na.rm=T)),
#       type="l", lwd=2, xlab="Years", ylab="Ln population size", ...)
#  plot(years, N, 
#       xlim = c(years[1], years[last] + fig),
#       ylim = c(0, 1.1 * max(1.1 * lastN, N, na.rm = T)), yaxs = "i",
#       type = "l", lwd = 2, xlab = "År", ylab = "Bestandsstørrelse") # ...
  plot(years, X, main = name,
       xlim = c(years[1], years[last] + fig),
       ylim = c(lnC, max(1.2 * (X[last] - lnC) + lnC, X, na.rm=T)), yaxt = "n",
       type = "l", lwd = 2, xlab = "År", ylab = "Bestandsstørrelse") # ...
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
  x <- rep(N[last], nboot)
  stikkpr <- 1:min(linjer, nboot)
  utvalg  <- matrix(0, linjer, tmax)
  kvantiler <- matrix(0, length(kvntl), tmax)
#  est <- matrix(0, nrow(model$bootout), 5)
#  est[, 1:ncol(model$bootout)] <- model$bootout #¤ dette går enklere! est <- bootout! (eller boottab!)
  for (t in 1:tmax) {
    alive <- x > C
    if (t %=% crash.year) x <- round(x * (1 - crash.prop))
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
    kvantiler[, t] <- quantile(x, kvntl)
    ext[t] <- sum(x <      C) # n of trajectories that  are  extinct  by time t
    rec[t] <- sum(x >= lastN) # n of trajectories that have recovered by time t
  }
  
  for (i in stikkpr) {
    lines(years[last] + 0:1, c(X[last], ln(na.und(utvalg[i, 1]))), col = grey(0.72))
    lines(years[last] + 1:tmax,         ln(na.und(utvalg[i,  ])),  col = grey(0.72))
  }
  for (i in 1:length(kvntl)) {
    lines(years[last] + 0:1, c(X[last], ln(na.und(kvantiler[i, 1]))),
          col = if (crash.prop > 0) "red" else "black")
    lines(years[last] + 1:tmax, ln(na.und(kvantiler[i, ])))
  }
  
  # Utmating -------------------------------------------------------------------
  utmating <- matrix("", length(kvntl), 2)
  colnames(utmating) <- c("Levetid", "Restitueringstid")
  rownames(utmating) <- (kvntl * 100) %+% "-persentil"
  for (i in 1:length(kvntl)) {
    if (any(ext / nboot >= max(kvntl[i], 0.000001))) {
      lt <- min(which(ext / nboot >= max(kvntl[i], 0.000001)))
    } else {
      lt <- "> " %+% format(tmax, sci = F)
    }
    if (any(rec / nboot >= max(kvntl[i], 0.000001))) {
      rt <- min(which(rec / nboot >= max(kvntl[i], 0.000001)))
    } else {
      rt <- "> " %+% format(tmax, sci = F)
    }
    utmating[i, ] <- c(lt, rt)
  }
#  print.table(utmating, right = T)
  cat("\nLevedyktighetsanalyse...")
  cat("\n- Tid til utdøing  (år):")
  for (i in 1:length(kvntl)) cat("  " %+% heltall(utmating[i, 1]))
  cat("\n- Restitueringstid (år):")
  for (i in 1:length(kvntl)) cat("  " %+% heltall(utmating[i, 2]))
  
  
  
    
  #    for (i in kvntl) {
  #      if (any(ext / n >= max(i, 0.000001))) {
  #        lt <- min(which(ext / n >= max(i, 0.000001)))
  #      } else {
  #        lt <- ">1000"
  #      }
  #      cat("Levetidas " %+% (i * 100) %+% "-persentil: " %+% lt %+% " år\n")
  #    }
  
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
  
  invisible(list(extinction = ext, recovery = rec))
}

