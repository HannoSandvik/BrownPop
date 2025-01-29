#####
# pva
#####

# ikke rydda ennå!


pva <- function(model, crash.year=NULL, crash.prop=0, C=20, cov1=NULL, cov2=NULL, cov3=NULL, n=10000, tmax=NULL, fig=100, fil=NULL, seed=NULL, ...)
{#.....
  # model: en fiks-ferdig populasjonsmodell (med bootout=TRUE!)
  # cov1-3: matriser med kovariatene (nrow=n, ncol=tmax) for fremtida
  na.und <- function(x) ifelse(is.na(x), C, ifelse(x < C, C, x))
  risk <- function(m) 100 - c((sum(m[,10] > C)) * 100 / nrow(m),
                              (sum(m[,20] > C)) * 100 / nrow(m), (sum(m[,100] > C)) * 100 / nrow(m))
  ok <- T
  N <- unlist(model$N)
  years <- model$years
  K <- NULL #???????
  sd <- model$sigd
  se <- model$sig2
  last <- max(which(!is.na(N)))
  x <- ln(N)
  C <- ln(C)
  if (is.null(tmax))
  {if (is.null(cov1))
  {tmax <- 1000}
    else
    {tmax <- ncol(cov1)}}
  if (is.null(crash.year))
  {crash.year <- Inf}
  if (is.null(cov1))
  {cov1 <- matrix(0, n, tmax)}
  if (is.null(cov2))
  {cov2 <- matrix(0, n, tmax)}
  if (is.null(cov3))
  {cov3 <- matrix(0, n, tmax)}
  if (nrow(cov1) < n)
  {stop("ERROR1")}
  if (nrow(cov2) < n)
  {stop("ERROR2")}
  if (nrow(cov3) < n)
  {stop("ERROR3")}
  if (nrow(model$bootout) < n)
  {stop("ERROR4")}
  sik <- list()
  m <- matrix(0, n, max(100,fig))
  ext <- NA
  if (ok)
  {graphics.off()
    plot(years[last - 30:0], x[last - 30:0], xlim=c(years[last] - 40, years[last] + fig),
         ylim=c(C, max(1.618 * (x[last] - C) + C, x[last - 30:0], na.rm=T)),
         type="l", lwd=2, xlab="Years", ylab="Ln population size", ...)
    if (is.null(K))
    {### Density-independent model
      est <- matrix(0, nrow(model$bootout), 5)
      est[, 1:ncol(model$bootout)] <- model$bootout
      for (i in 1:n)
      {alive <- T
      t <- 1
      if (!is.null(seed))
      {set.seed(seed + i)}
      x <- ln(N[last])
      while(alive & (t <= tmax))
      {x[t+1] <- exp(x[t] - 0.5 * sd * exp(-x[t]) + est[i,1] +
                       est[i,3] * cov1[i,t] + est[i,4] * cov2[i,t] + est[i,5] * cov3[i,t] +
                       sqrt(est[i,2] + sd * exp(-x[t])) * rnorm(1))
      if (t %=% crash.year)
      {x[t+1] <- x[t+1] * (1 - crash.prop)}
      x[t+1] <- ln(round(x[t+1]))
      if (x[t+1] <= C)
      {x[t+1] <- 0
      alive <- F}
      t <- t + 1}
      sik[[i]] <- x
      m[i, 1:min(t,fig)] <- x[1:min(t,fig)]
      ext[i] <- t - 1
      if (!(i %% (ceiling(10000/6000)*120)))
      {cat(i %+% " of " %+% n %+% " simulations done.\n  Extinction risk: ")
        cat(paste(round(risk(m[1:i,])), collapse="/") %+% "% (within 10/20/100 years.)\n")
        lines(years[last] - 1 + 1:fig, na.und(x[1:fig]), col=grey(0.72))}}}
    
    else
    {### Density-dependent model
      cat("Sorry, the density-dependent version is not yet implemented!\n")}
    
    if (is.null(fil))
    {j <- 1
    fil <- getwd() %+% "/pva" %+% j %+% ".R"
    while (file.exists(fil))
    {j <- j + 1
    fil <- getwd() %+% "/pva" %+% j %+% ".R"}}
    for (q in c(0.05, 0.1, 0.2, 0.8, 0.9, 0.95))
    {lines(years[last] - 1 + 1:fig, na.und(apply(m[,1:fig], 2, quantile, q)))}
    lines(years[last] - 1 + 1:fig, na.und(apply(m[,1:fig], 2, quantile, 0.5)), lwd=2)
    save(sik, m, ext, file=fil)
    cat("\nThe PVA results have been saved as \"" %+% fil %+% "\".\n\n")
    extR <- matrix(NA, 13, 4, dimnames=list(c("extinction.risk", "expected.pop.lifetime",
                                              c(50, 20, 80, 10, 90, 5, 95, 2.5, 97.5, 0.5, 99.5) %+% ".pct.pop.lifetime"),
                                            c("10a", "20a", "100a", "exact")))
    extR[1,1:3] <- risk(m)
    extR[2,1:3] <- -c(10,20,100) / ln(1 - risk(m)/100)
    if(all(ext < tmax))
    {extR[2,4] <- mean(ext)}
    for (i in 3:13)
    {q <- c(0.5, 0.2, 0.8, 0.1, 0.9, 0.05, 0.95, 0.025, 0.975, 0.005, 0.995)[i - 2]
    extR[i,1:3] <- -extR[2,1:3] * ln(1 - q)
    if (sum(ext < tmax) > (n * q))
    {extR[i,4] <- quantile(ext, q)}}
    cat("Extinction risk (percent within n years):\n")
    er <- as.vector(extR[1,1:3])
    names(er) <- c("10", "20", "100")
    print(er)
    cat("\nPopulation lifetime (quantiles in years):\n")
    pl <- round(as.vector(extR[3:13,4]))
    if (any(is.na(pl)))
    {pl[is.na(pl)] <- ">" %+% tmax}
    names(pl) <- c(50, 20, 80, 10, 90, 5, 95, 2.5, 97.5, 0.5, 99.5) %+% "%"
    print(pl, quote=F)
    if (crash.prop > 0)
    {return <- NA
    cat("\nReturn to the pre-crash population size (quantiles in years):\n")
    for (i in 1:length(sik))
    {x  <- sik[[i]][-1]
    x1 <- sik[[i]][1]
    if (all(x < x1))
    {return[i] <- tmax}
    else
    {return[i] <- which(x >= x1)[1]}}
    q <- round(quantile(return, c(0.5, 0.2, 0.8, 0.1, 0.9,
                                  0.05, 0.95, 0.025, 0.975, 0.005, 0.995)))
    if (any(q == tmax))
    {q[which(q == tmax)] <- ">" %+% tmax}
    print(q, quote=F)}
    ext <- extR}
  invisible(ext)}

