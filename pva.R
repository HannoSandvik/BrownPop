#####
# pva
#####



pva <- function(model, 
                crash.year=NULL, 
                crash.prop=0, 
                C=20, 
                seed=NULL, 
                ...) {
  # model: en fiks-ferdig populasjonsmodell (med bootout=TRUE!)
  # cov1-3: matriser med kovariatene (nrow=n, ncol=tmax) for fremtida

  na.und <- function(x) ifelse(is.na(x), C, ifelse(x < C, C, x))
  risk <- function(m) 100 - c((sum(m[, 10] > C)) * 100 / nrow(m),
                              (sum(m[, 20] > C)) * 100 / nrow(m), 
                              (sum(m[,100] > C)) * 100 / nrow(m)) #¤ brukes ikke per nå
  ok <- T
  N <- unlist(model$N)
  n <- nrow(model$bootout) #!¤ ikke tillat ulike n i popmodellen og pvaen!
  cov1 <- NULL #¤ skriv disse ut av skriptet ...
  cov2 <- NULL
  cov3 <- NULL
  tmax <- 1000
  fig <- 100
  fil <- NULL # ... og hit ¤
  linjer <- 25 # antall tilfeldige bestandsforløp som skal vises i grafen
  years <- model$years
  K <- NULL #¤
  sd <- model$sigd
  se <- model$sig2
  last <- max(which(!is.na(N)))
  x   <- ln(N)
  lnC <- ln(C)
#  if (is.null(tmax)) {
#    if (is.null(cov1)) {
#      tmax <- 1000
#    } else {
#      tmax <- ncol(cov1)
#    }
#  }
  if (is.null(crash.year)) crash.year <- Inf
  if (is.null(cov1)) cov1 <- matrix(0, n, tmax)
  if (is.null(cov2)) cov2 <- matrix(0, n, tmax)
  if (is.null(cov3)) cov3 <- matrix(0, n, tmax)
  if (nrow(cov1) < n) stop("ERROR1")
  if (nrow(cov2) < n) stop("ERROR2")
  if (nrow(cov3) < n) stop("ERROR3")
  if (nrow(model$bootout) < n) stop("ERROR4")
  sik <- list()
  m <- matrix(0, n, max(100, fig))
  ext <- NA
  if (ok) {
    graphics.off()
    plot(years[last - 30:0], x[last - 30:0], 
         xlim=c(years[last] - 40, years[last] + fig),
         ylim=c(lnC, max(1.618 * (x[last] - lnC) + lnC, x[last - 30:0], na.rm=T)),
         type="l", lwd=2, xlab="Years", ylab="Ln population size", ...)
    if (is.null(K)) {
      ### Density-independent model
      x <- forrige <- rep(N[last], n)
      kvntl <- c(0.025, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.975)
      kvantiler <- rep(0, length(kvntl))
      fkvant <- quantile(x, kvntl)
      stikkpr <- round((1:linjer - 0.5) * n / linjer)
      est <- matrix(0, nrow(model$bootout), 5)
      est[, 1:ncol(model$bootout)] <- model$bootout
      if (!is.null(seed)) set.seed(seed)
      for (t in 1:1000) {
        alive <- x > C
        if (t %=% crash.year) {
          x <- round(x * (1 - crash.prop))
        }
        x <- ln(x)
        forrige <- x
        #x[!alive] <- 0
        #t <- 1
        #x <- ln(N[last])
        #while(alive & (t <= tmax)) {
        if (any(alive)) {
          x[alive] <- round(exp(x[alive] - 0.5 * sd * exp(-x[alive]) + est[alive, 1] +
#                        est[i,3] * cov1[i,t] + 
#                        est[i,4] * cov2[i,t] + 
#                        est[i,5] * cov3[i,t] + 
                        sqrt(est[alive, 2] + sd * exp(-x[alive])) * rnorm(sum(alive))))
        } else {
          x <- rep(0, n)
        }
        x[x <= C] <- 0
        kvantiler <- quantile(x, kvntl)
        ext[t] <- sum(x < C) # number of trajectories that are extinct by time t
        #          if (x[t+1] <= C) {
#            x[t+1] <- 0
#            alive <- F
#          }
#          t <- t + 1
        #}
        #sik[[i]] <- x
        #m[i, 1:min(t, fig)] <- x[1:min(t, fig)]
        #ext[i] <- t - 1
        #if (!(i %% (ceiling(10000 / 6000) * 120))) {
        #  cat(i %+% " of " %+% n %+% " simulations done.\n  Extinction risk: ")
        #  cat(paste(round(risk(m[1:i,])), collapse="/") %+%
        #        "% (within 10/20/100 years.)\n")
        for (i in stikkpr) {
          lines(years[last] + t - 1:0, ln(c(na.und(exp(forrige[i])), na.und(x[i]))), 
                col=grey(0.72))
        }
        for (i in 1:length(kvntl)) {
          lines(years[last] + t - 1:0, ln(c(na.und(fkvant[i]), na.und(kvantiler[i]))))
        }
        fkvant <- kvantiler
        #}
      }
    } else {
      ### Density-dependent model
      cat("Sorry, the density-dependent version is not yet implemented!\n")
    }
    
    for (i in kvntl) {
      if (any(ext / n >= max(i, 0.000001))) {
        lt <- min(which(ext / n >= max(i, 0.000001)))
      } else {
        lt <- ">1000"
      }
      cat("Levetidas " %+% i %+% "-persentil: " %+% lt %+% " år\n")
    }
    
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
  }
  invisible(ext)
}

