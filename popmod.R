########
# popmod
########



popmod <- function(pop.size, 
                   years = NULL, 
                   catch = NULL, 
                   sd = 0.1, 
                   cov1 = NULL, cov2 = NULL, cov3 = NULL, 
                   nboot = 10000, 
                   nsim = 0, 
                   name = "", 
                   plot = TRUE, 
                   se = F, r = F, b1 = F, b2 = F, b3 = F, 
                   mean1 = F, mean2 = F, mean3 = F, 
                   bootout = TRUE, 
                   K = F) {
  # sd is demographic variance
  # if K:  logistic model
  # if !K: Brownian model
  
  X <- ln(pop.size)
  
  if (is.vector(X)) X <- data.frame("lnN" = X)
  
  na.und <- function(x) ifelse(is.na(x), 0, x)
  
  log.gauss <- function(X, my, sig2) -0.5*log(2*pi*sig2)-0.5*((X-my)^2)/sig2
  
  ######### likelihood function

  lik <- function(par) {
    L <- 0
    for (i in 1:npop) {
      k  <- na.und(-exp(par[q[1,i]]))
      m  <- ifelse(K, exp(par[q[2,i]]), par[q[2,i]])
      s  <- exp(par[q[3,i]])
      w  <- log(exp(X[[i]])-catch[[i]])
      b1 <- na.und(par[q[4,i]])
      b2 <- na.und(par[q[5,i]])
      b3 <- na.und(par[q[6,i]])
      my <- w - 0.5 * sd * exp(-w) + m + k * exp(w) + 
        b1*cov1[[i]] + b2*cov2[[i]] + b3*cov3[[i]]
      s <- s + sd * exp(-w)
      L  <- L - sum(log.gauss(X[[i]][-1], my[-length(my)], s[-length(s)]),
                na.rm=T)
    }
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
  if (is.null(sd))  stop ("demographic variance (sd) not specified, try sd=0")
  if (any(X > 30, na.rm = T)) {
    cat("\nOBS: Jeg har en viss mistanke om at du har glemt " %+%
         "å log-transformere bestandstallene.\n\n")
  }
  if (any(X < 0, na.rm = T)) stop ("N kan ikke være mindre enn 1!")
  for (i in 1:ncol(X)) {
    if (length(na.omit(X[, i])) < 6) {
      if (ncol(X) %=% 1) {
        stop ("For få år med data!")
      } else {
        stop ("For få år med data i kolonne " %+% "i" %+% "!")
      }
    }
  }
  if (is.null(catch)) catch <- rep(0, length(years))
  if (is.null(cov1))   cov1 <- rep(0, length(years))
  if (is.null(cov2))   cov2 <- rep(0, length(years))
  if (is.null(cov3))   cov3 <- rep(0, length(years))
  if (years %!=% (min(years):max(years))) {
    cat("\nOBS: Årstallene er ikke sammenhengende. Er det meningen?\n\n")
    NAS <- cbind(years, cov1, cov2, cov3, X)
    NAS <- merge(years = min(years):max(years), NAS, all.x = T)
    years <- NAS[,1]
    cov1 <- NAS[,2]
    cov2 <- NAS[,3]
    cov3 <- NAS[,4]
    X <- as.data.frame(NAS[, 5:ncol(NAS)])
  }
  
  ncov <- sum(c(any(cov1 != 0), any(cov2 != 0), any(cov3 != 0)))
  npop <- ncol(X)
  nk   <- K * npop
  nr   <- ifelse(r,  npop, 1)
  nse  <- ifelse(se, npop, 1)
  nb1  <- ifelse(b1,  npop, 1) * (ncov > 0.5)
  nb2  <- ifelse(b2,  npop, 1) * (ncov > 1.5)
  nb3  <- ifelse(b3,  npop, 1) * (ncov > 2.5)
  
  orig.years <- data.frame(years=min(years):max(years))
  
  bootres <- psim <- NULL
  Catch <- catch
  Year <- years
  Cov1 <- as.list(as.data.frame(cov1))
  Cov2 <- as.list(as.data.frame(cov2))
  Cov3 <- as.list(as.data.frame(cov3))
  if ((length(Cov1) > 1) & !b1 & mean1) Cov1 <- list(apply(cov1, 1, mean))
  if ((length(Cov2) > 1) & !b2 & mean2) Cov2 <- list(apply(cov2, 1, mean))
  if ((length(Cov3) > 1) & !b3 & mean3) Cov3 <- list(apply(cov3, 1, mean))
  X <- N <- catch <- years <- cov1 <- cov2 <- cov3 <- as.list(X)
  predicted <- residual <- demcomp <- as.list(X)
  pred.r <- lci <- uci <- as.list(X)
  p.res.t <- p.res.N <- p.white <- vector("numeric", npop)
  for (i in 1:length(X)) {
    N[[i]] <- exp(X[[i]])
    names(N[[i]]) <- Year
    years[[i]] <- Year
    cov1[[i]] <- Cov1[[min(i, length(Cov1))]]
    cov2[[i]] <- Cov2[[min(i, length(Cov2))]]
    cov3[[i]] <- Cov3[[min(i, length(Cov3))]]
    catch[[i]] <- Catch
    names(catch[[i]]) <- Year
    NAS <- is.na(X[[i]] + cov1[[i]] + cov2[[i]] + cov3[[i]])
    cov1[[i]][NAS] <- NA
    cov2[[i]][NAS] <- NA
    cov3[[i]][NAS] <- NA
    if(any(na.omit(cov1[[i]])!=0)) cov1[[i]] <- cov1[[i]] - mean(cov1[[i]], na.rm=T)
    if(any(na.omit(cov2[[i]])!=0)) cov2[[i]] <- cov2[[i]] - mean(cov2[[i]], na.rm=T)
    if(any(na.omit(cov3[[i]])!=0)) cov3[[i]] <- cov3[[i]] - mean(cov3[[i]], na.rm=T)
  }
  
  X.orig <- X
  par <- rep(-2, nk + nr + nse + nb1 + nb2 + nb3)
  if (K) {
    est.r <- est.k <- rep(0, npop)
    for (i in 1:npop) {
      fit <- lm(diff(X[[i]]) ~ exp(X[[i]][-length(X[[i]])]))
      est.r[i] <- ln(max(1e-12, coef(fit)[1]))
      est.k[i] <- est.r[i] - ln(max(10, min(N[[i]], na.rm = T) / 10, 
                                    -coef(fit)[1] / coef(fit)[2]))
    }
    par[1:nk] <- est.k
    if (nr > 1) {
      par[nk + 1:nr] <- est.r
    } else {
      par[nk + 1] <- mean(est.r)
    }
  } else {
    par[nk + 1:nr] <- pi / 100
  }
  if (ncov > 0) par[nk + nr + nse + 1:nb1] <- 0
  if (ncov > 1) par[nk + nr + nse + nb1 + 1:nb2] <- 0
  if (ncov > 2) par[nk + nr + nse + nb1 + nb2 + 1:nb3] <- 0
  q <- matrix(length(par) + 12, 6, npop)
  if (K) q[1,] <- 1:nk
  q[2,] <- nk + 1:nr
  q[3,] <- nk + nr + 1:nse
  if (ncov > 0) q[4,] <- nk + nr + nse + 1:nb1
  if (ncov > 1) q[5,] <- nk + nr + nse + nb1 + 1:nb2
  if (ncov > 2) q[6,] <- nk + nr + nse + nb1 + nb2 + 1:nb3
  
  ######################### estimation of parameters

  estimat <- estparam(par)
  k.est <- my.est <- se.est <- b1.est <- b2.est <- b3.est <- rep(0, npop)
  my.est[1:npop] <- ifelse(K, exp(estimat$par[q[2,]]), estimat$par[q[2,]])
  if (K) {
    k.est[1:npop]  <- my.est / exp(estimat$par[q[1,]])
  } else {
    k.est[1:npop] <- 120
  }
  se.est[1:npop] <-    exp(estimat$par[q[3,]])
  b1.est[1:npop] <- na.und(estimat$par[q[4,]])
  b2.est[1:npop] <- na.und(estimat$par[q[5,]])
  b3.est[1:npop] <- na.und(estimat$par[q[6,]])
  lnL <- -estimat$value
  if (estimat$convergence)
  {cat("\n############################################\n" %+%
         "OBS: Optimaliseringen konvergerte ikke!!!!!" %+%
         "\n############################################\n")}
  
  for (i in 1:npop) {
    ww <- log(exp(X[[i]]) - catch[[i]])
    predicted[[i]] <- c(NA, ww + my.est[i] -
                        K * exp(ww) * my.est[i] / k.est[i] + 
                        b1.est[i] * cov1[[i]] + 
                        b2.est[i] * cov2[[i]] +
                        b3.est[i] * cov3[[i]] -
                        0.5 * sd * exp(-ww))
  residual[[i]] <-  c(X[[i]], NA) - predicted[[i]]
  residual[[i]] <- residual[[i]][-1]
  predicted[[i]] <- predicted[[i]][-length(predicted[[i]])]
  demcomp[[i]] <- sd / (exp(X[[i]]) - catch[[i]])
  df.res <- data.frame("years" = years[[i]],
                       "residual" = residual[[i]],
                       "demcomp" = demcomp[[i]], 
                       "predicted" = predicted[[i]])
  res <- merge(orig.years, df.res, all.x = T)
  residual[[i]] <- res$residual
  demcomp[[i]] <- res$demcomp
  predicted[[i]] <- res$predicted
  names(predicted[[i]]) <- names(demcomp[[i]]) <-
    names(residual[[i]]) <- as.character(orig.years$years)}
  
  #####################################################
  
  r2 <- NA
  for (i in 1:npop) {
    pred.r[[i]] <- my.est[i] - 0.5 * sd * exp(-X[[i]]) -
      K * exp(X[[i]]) * my.est[i] / k.est[i] +
      b1.est[i] * na.und(cov1[[i]]) +
      b2.est[i] * na.und(cov2[[i]]) +
      b3.est[i] * na.und(cov3[[i]])
    r2[i] <- summary(lm(diff(X[[i]]) ~ pred.r[[i]][-length(X[[i]])]))$r.sq
    uci[[i]] <- lci[[i]] <- NA
    for (j in 1:(length(pred.r[[i]]))) {
      lci[[i]][j] <- pred.r[[i]][j] + 
        sqrt(se.est[i] + sd * exp(-X[[i]][j])) * qnorm(0.025)
      uci[[i]][j] <- pred.r[[i]][j] +
        sqrt(se.est[i] + sd * exp(-X[[i]][j])) * qnorm(0.975)
    }
    names(pred.r[[i]]) <- names(lci[[i]]) <- names(uci[[i]]) <-
      as.character(orig.years$years)
  }
  r2.tot <- summary(lm(diff(
    as.vector(unlist(X)))[-(length(X[[1]]) * 1:npop)] ~
    as.vector(unlist(pred.r))[-(length(X[[1]]) * 1:npop)]))$r.sq
  
  ########## simulate data based on estimated parameters
  # NY versjon av simdat
  # (skrevet av Hanno etter forbilde av Vidars logismodfit):
  
  simdat <- function(b1=T, b2=T, b3=T) {
    if (b1) b1 <- b1.est else b1 <- rep(0, length(b1.est))
    if (b2) b2 <- b2.est else b2 <- rep(0, length(b2.est))
    if (b3) b3 <- b3.est else b3 <- rep(0, length(b3.est))
    sim <- function(i) {
      erst <- min(which(!is.na(X.orig[[i]])))
      xsim <- rep(NA, length(X.orig[[i]]))
      xsim[erst] <- X.orig[[i]][erst]
      for (j in erst:(length(X.orig[[i]])-1)) {
        xsim[j + 1] <- xsim[j] - 0.5 * sd * exp(-xsim[j]) + my.est[i] -
          K * exp(xsim[j]) * my.est[i] / k.est[i] +
          b1[i] * na.und(cov1[[i]][j]) + 
          b2[i] * na.und(cov2[[i]][j]) + 
          b3[i] * na.und(cov3[[i]][j]) +
          sqrt(se.est[i] + sd * exp(-xsim[j])) * rnorm(1)
        if (xsim[j+1] < 0) xsim[j+1] <- 0
      }
      return(xsim)
    }
    xsim <- X.orig
    for (i in 1:length(X.orig)) {
      xsim[[i]] <- sim(i)
      xsim[[i]][is.na(X.orig[[i]])] <- 12
      zahl <- max(sum(!is.na(X.orig[[i]])) / 12, any(X.orig[[i]] == 0, na.rm = T))
      while (sum(xsim[[i]] == 0) > zahl) {
        xsim[[i]] <- sim(i)
        xsim[[i]][is.na(X.orig[[i]])] <- 12
      }
      xsim[[i]][is.na(X.orig[[i]])] <- NA
    }
    return(xsim)
  }

  simdat2 <- function(b1=T, b2=T, b3=T) {
    # IKKE tilpasset til logistisk modell!!!
    if (b1) b1 <- b1.est else b1 <- rep(0, length(b1.est))
    if (b2) b2 <- b2.est else b2 <- rep(0, length(b2.est))
    if (b3) b3 <- b3.est else b3 <- rep(0, length(b3.est))
    sim <- function(i) {
      erst <- min(which(!is.na(X.orig[[i]])))
      lezt <- max(which(!is.na(X.orig[[i]])))
      xsim <- rep(NA, length(X.orig[[i]]))
      xsim[erst] <- X.orig[[i]][lezt]
      for (j in (lezt - 1):erst) {
        xsim[erst] <- xsim[erst] - my.est[i] -
          b1[i] * na.und(cov1[[i]][j]) -
          b2[i] * na.und(cov2[[i]][j]) -
          b3[i] * na.und(cov3[[i]][j])
      }
      eps <- 5
      while (eps > 0.000001) {
        for (j in erst:(length(X.orig[[i]]) - 1)) {
          xsim[j+1] <- xsim[j] - 0.5 * sd * exp(-xsim[j]) + my.est[i] +
            b1[i] * na.und(cov1[[i]][j]) + 
            b2[i] * na.und(cov2[[i]][j]) + 
            b3[i] * na.und(cov3[[i]][j]) +
            sqrt(se.est[i] + sd * exp(-xsim[j])) * rnorm(1)
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
    for (i in 1:length(X.orig)) {
      xsim[[i]] <- sim(i)
      xsim[[i]][is.na(X.orig[[i]])] <- 12
      zahl <- max(sum(!is.na(X.orig[[i]])) / 12, any(X.orig[[i]] == 0, na.rm = T))
      while (sum(xsim[[i]] == 0) > zahl) {
        xsim[[i]] <- sim(i)
        xsim[[i]][is.na(X.orig[[i]])] <- 12
      }
      xsim[[i]][is.na(X.orig[[i]])] <- NA
    }
    return(xsim)
  }
  
  ############################### parametric bootstrapping

  if (nboot>0) {
    cat("bootstrapping\n")
    catval <- c(1, seq(0, nboot, nboot / 10))
    boottab <- matrix(NA, nboot, nk + nr + nse + nb1 + nb2 + nb3)
    for (j in 1:nboot) {
      if (any(j==catval)) cat("boot",j, "of", nboot, "\n")
      X <- simdat()
      estimat <- estparam(par)
      boottab[j,] <- estimat$par
    }
    if (K) {
      boottab[,nk + 1:nr] <- exp(boottab[,nk + 1:nr])
      boottab[,1:nk] <- boottab[,nk + 1:nr] / exp(boottab[,1:nk])}
      boottab[,nk + nr + 1:nse] <- exp(boottab[,nk + nr + 1:nse])
      bootres <- matrix(NA, 18, nk + nr + nse + nb1 + nb2 + nb3)
      rownames(bootres) <- c("min", "0.05%", "0.1%", "0.5%", "1%", "2.5%",
                             "5%", "25%", "mean", "75%", "95%", "97.5%", "99%", 
                             "99.5%", "99.9%", "99.95%", "max", "median")
      colnames(bootres) <- 1:ncol(bootres)
      if (K) colnames(bootres)[1:nk] <- "k." %+% 1:nk
      colnames(bootres)[nk + 1:nr] <- "r." %+% 1:nr
      colnames(bootres)[nk + nr + 1:nse] <- "sig2." %+% 1:nse
      if (nb1) colnames(bootres)[nk + nr + nse +   1:nb1] <- "b1." %+% 1:nb1
      if (nb2) colnames(bootres)[nk + nr + nse+nb1+1:nb2] <- "b2." %+% 1:nb2
      if (nb3) colnames(bootres)[nk+nr+nse+nb1+nb2+1:nb3] <- "b3." %+% 1:nb3
      for (j in 1:ncol(bootres)) {
        bootres[1:17,j] <- quantile(boottab[,j], 
                                    c(0, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 
                                      0.25, 0.5, 0.75, 0.95, 0.975, 0.99, 0.995, 
                                      0.999, 0.9995, 1), na.rm = T)
      bootres[18,j] <- median(boottab[,j], na.rm=T)
    }
  }
  if (!bootout) boottab <- NULL
  
  ##################### simulation test of covariates

  if (nsim > 0) {
    cat("simulering\n")
    catval <- c(1,seq(0,nsim,nsim/10))
    simres <- matrix(NA,nsim,3*npop)
    #dyn.load(rot.dir() %+% "uni\\simdat_brownian.dll","sim", "cdecl") 
    #par <- c(0.04, -2, rep(0, ncov))
    #q <- matrix(rep(1:(2 + ncov), npop), 2 + ncov, npop)
    
    if (ncov > 0) {
      for (j in 1:nsim) {
        if (any(j==catval)) {
          cat("sim",j, "of", nsim, "covariate 1 of ",ncov,  "\n")
        }
        X <- simdat(b1=F)
        estimat <- estparam(par)
        simres[j,1:npop] <- na.und(estimat$par[q[3,]])
      }
    }
    if (ncov > 1) {
      for (j in 1:nsim) {
        if (any(j==catval)) {
          cat("sim",j, "of", nsim, "covariate 2 of ",ncov,  "\n")
        }
        X <- simdat(b2=F)
        estimat <- estparam(par)
        simres[j,npop+1:npop] <- na.und(estimat$par[q[4,]])
      }
    }
    if (ncov > 2) {
      for (j in 1:nsim) {
        if (any(j==catval)) {
          cat("sim",j, "of", nsim, "covariate 3 of ",ncov,  "\n")
        }
        X <- simdat(b3=F)
        estimat <- estparam(par)
        simres[j,2*npop+1:npop] <- na.und(estimat$par[q[5,]])
      }
    }
    
    psim <- matrix(NA, 3*npop, 3,
                   dimnames=list(rep(c("b1","b2","b3"), each=npop) %+% "." %+% 1:npop,
                                 c("p.less", "p.greater", "p.diff.0")))
    for (i in 1:(3*npop)) {
      psim[i,1] <- sum(simres[,i]<b1.est[1 + (i - 1) %% npop])/nsim
      psim[i,2] <- sum(simres[,i]>b1.est[1 + (i - 1) %% npop])/nsim
      psim[i,3] <- sum(abs(simres[,i])>abs(b1.est[1 + (i - 1) %% npop]))/nsim
    }
    
    #save(simres, file="d:\\r31" %+% i %+% ".rda")
    
    #   dyn.unload(rot.dir() %+% "uni\\simdat_brownian.dll")
  }

  ####### 3 TESTS!!!

  for (i in 1:npop) {
    lmfit <- lm(residual[[i]] ~ orig.years$years, na.action="na.exclude")
    p.res.t[i] <- signif(anova(lmfit)[1,5],3)
    
    lmfit <- lm(residual[[i]] ~ N[[i]], na.action="na.exclude")
    p.res.N[i] <- signif(anova(lmfit)[1,5],3)
    
    portman.Q <- function (X, K) {
      # portman.Q uses the cummulative ACF to test for whiteness  of a time series.
      # This is the Ljung-Box version of the the Portemanteau  test for
      # whiteness (Tong 1990). It may in particular be  usefull to test
      # for whiteness in the residuals from time  series models.
      # 
      # A vector is returned consisting of the asymtpotic  chi-square
      # value, the associated d.f. and asymptotic  p.val for the test of
      # whiteness. p.val<0.05 -> non-whiteness!
      # Tong, H. (1990) Non-linear time series : a dynamical  system approach. Clarendon Press, Oxford.
      # Author: Ottar N. Bjornstad onb1@psu.edu
      K <- min(K, length(na.omit(X))-1) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Q <- 0
      n <- length(X)
      p <- acf(X, plot = FALSE, lag.max = K,na.action=na.pass)$acf[2:(K + 1)]
      for (k in 1:K) Q <- Q + p[k]^2/(n - k)
      Q <- n * (n + 2) * Q
      res <- list(chisq = signif(Q, 3), df = K, p.val = signif(1 - pchisq(Q, K), 3))
      unlist(res)
    }
    p.white[i] <- portman.Q(residual[[i]],K=10)[[3]]
  }
  
  #############################
  
  if (!is.null(names(N))) {
    navn <- names(N)
    names(k.est) <- names(my.est) <- names(se.est) <- names(r2) <- navn
    names(b1.est) <- names(b2.est) <- names(b3.est) <- navn
    names(p.res.t) <- names(p.res.N) <- names(p.white) <- navn
  }
  
  n.par <- nk + nr + nse + nb1 + nb2 + nb3
  n.eff <- 0
  for (i in 1:npop) {
    n.eff <- n.eff + length(na.omit(predicted[[i]]))
  }
  aic <- 2 * (n.par - lnL)
  aicc<- aic + 2 * n.par * (n.par + 1) / (n.eff - n.par - 1)

  if (!K) k.est <- rep(NA, npop)
  res <- list(
    "K"=k.est,
    "s"=my.est,
    "sig2"=se.est,
    "sigd"=sd,
    "b1"=b1.est, "b2"=b2.est, "b3"=b3.est,
    "years"=orig.years$years,
    "N"=N, "catch"=Catch,
    "ncov"=ncov, "cov1"=cov1, "cov2"=cov2, "cov3"=cov3,
    "residuals"=residual, "predicted"=predicted, "demcomp"=demcomp,
    "pred.r"=pred.r, "lci"=lci, "uci"=uci,
    "bootout"=boottab, "boot"=bootres, "p.cov"=psim,
    "model"=ifelse(K, "logistic", "Brownian"), "name"=name,
    "p.res.t"=p.res.t, "p.res.N"=p.res.N, "p.white"=p.white,
    "lnL"=lnL, "r2"=r2, "r2.tot"=r2.tot, "aic"=aic, "aicc"=aicc)
  if (plot & (npop %=% 1)) {
    if (K) {
      logistic.plot(res)
    } else {
      brownian.plot(res)
    }
    cat("\nSee plot for diagnostics.\n")
  }
  cat("\nA " %+% ifelse(K, "logistic (density-dependent)",
                        "Brownian (density-independent)") %+% " model was fitted.\n")
  cat("\nThe demographic variance was assumed to be " %+%
      rund(sd, 3) %+% ".\n")
  cat("\nParameter estimates:\n")
  partab <- matrix(NA, 6, 3)
  if (!is.na(res$K)) {
    partab[1,1]   <- res$K
    partab[1,2:3] <- res$boot[c("2.5%", "97.5%"), "K.1"]
  }
  partab  [2,1]   <- res$s
  partab  [2,2:3] <- res$boot[c("2.5%", "97.5%"), "r.1"]
  partab  [3,1]   <- res$sig2
  partab  [3,2:3] <- res$boot[c("2.5%", "97.5%"), "sig2.1"]
  if (ncov > 0) {
    partab[4,1]   <- res$b1
    partab[4,2:3] <- res$boot[c("2.5%", "97.5%"), "b1.1"]
  }
  if (ncov > 1) {
    partab[5,1]   <- res$b2
    partab[5,2:3] <- res$boot[c("2.5%", "97.5%"), "b2.1"]
  }
  if (ncov > 2) {
    partab[6,1]   <- res$b3
    partab[6,2:3] <- res$boot[c("2.5%", "97.5%"), "b3.1"]
  }
  rownames(partab) <- c("Carrying capacity  (K):", "Popul. growth rate (r):",
                        "Environmental variance:", "Slope of  covariate  1:",
                        "Slope of  covariate  2:", "Slope of  covariate  3:")
  colnames(partab) <- c("Mean", "Lower 95% CI", "Upper 95% CI")
  if (ncov < 3) partab <- partab[-6,]
  if (ncov < 2) partab <- partab[-5,]
  if (ncov < 1) partab <- partab[-4,]
  if (is.na(res$K)) partab <- partab[-1,]
  print(partab)
  return(res)
}

