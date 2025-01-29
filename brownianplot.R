###############
# brownian.plot
###############

# ikke rydda ennå!


brownian.plot <- function (fit) 
{
  if (tolower(fit$model) != "brownian") 
    stop("The fit is not based on the brownian model!")
  for (i in 1:length(fit)) {
    if (is.list(fit[[i]])) {
      fit[[i]] <- fit[[i]][[1]]
    }
  }
  lok <- paste(fit$name)
  year <- as.numeric(names(fit$N))
  obs.growth <- log(c(fit$N[-1], NA)) - log((fit$N - fit$catch))
  Nmax <- max(c(fit$K, fit$N), na.rm = T)
  if (Nmax > 2 * max(fit$N, na.rm = T)) 
    Nmax <- 1.2 * max(fit$N, na.rm = T)
  Nmin <- min(fit$N, na.rm = T)
  Nrange <- seq(trunc(Nmin * 0.9), trunc(Nmax * 1.1), by = 0.1)
  Nrange <- Nrange[Nrange > 0.9]
  pred <- NA
  pred <- -0.5 * fit$sigd * exp(-log(Nrange)) + fit$s
  is.R <- !is.null(version$language)
  if (is.R) {
    portman.Q <- function(X, K) {
      Q <- 0
      n <- length(X)
      p <- acf(X, plot = FALSE, lag.max = K, na.action = na.pass)$acf[2:(K + 
                                                                           1)]
      for (k in 1:K) Q <- Q + p[k]^2/(n - k)
      Q <- n * (n + 2) * Q
      res <- list(chisq = round(Q, 4), df = K, p.val = round(1 - 
                                                               pchisq(Q, K), 4))
      unlist(res)
    }
    par(mfrow = c(3, 2), mai = c(0.7, 0.7, 0.3, 0.2))
    plot(Nrange[is.finite(pred) & (pred < max(3 * abs(obs.growth), 
                                              na.rm = T))], pred[is.finite(pred) & (pred < max(3 * 
                                                                                                 abs(obs.growth), na.rm = T))], type = "n", ylab = "predicted/observed r", 
         xlab = "N", ylim = range(c(obs.growth, pred[is.finite(pred) & 
                                                       (pred < max(3 * abs(obs.growth), na.rm = T))]), 
                                  na.rm = T), sub = paste("s: ", round(fit$s, 3), 
                                                          "  (sig_e)^2: ", round(fit$sig2, 3)))
    lines(Nrange[is.finite(pred)], pred[is.finite(pred)], 
          col = "red", lwd = 2)
    points(fit$N, obs.growth)
    sel <- which((abs(scale(fit$residuals)) > 1.5))
    if (length(sel) > 0) 
      text(fit$N[sel], obs.growth[sel], labels = year[sel], 
           pos = 4, cex = 0.84)
    if (!is.null(lok)) 
      title(lok)
    plot(year, log(fit$N), type = "n", ylim = range(log(fit$N), 
                                                    na.rm = T) - c(0.1, -0.1), ylab = "observed/predicted ln N")
    points(year, log(fit$N))
    lines(year, log(fit$N))
    points(year, fit$predicted, col = "red", pch = 16)
    plot(year, obs.growth, type = "n", ylim = range(c(obs.growth, 
                                                      diff(fit$predicted)), na.rm = T), xlim = range(year), 
         ylab = "observed/predicted r")
    points(year, obs.growth)
    lines(year, obs.growth)
    points(year, -0.5 * fit$sigd * exp(-log(fit$N - fit$catch)) + 
             fit$s, col = "red", pch = 16)
    abline(h = 0, lty = 2)
    lmfit <- lm(fit$residuals ~ year, na.action = "na.exclude")
    prob <- anova(lmfit)[1, 5]
    out <- prob
    plot(year, fit$residuals, ylab = "residuals", sub = paste("p-val: ", 
                                                              round(prob, 3)))
    lines(year, fit$residuals)
    abline(h = 0, lty = 3, col = "blue")
    if (prob < 0.05) 
      abline(lmfit, lwd = 2, col = "red")
    else abline(lmfit, lwd = 1, lty = 2)
    lmfit <- lm(fit$residuals ~ fit$N, na.action = "na.exclude")
    prob <- anova(lmfit)[1, 5]
    out[2] <- prob
    port <- portman.Q(fit$residuals, K = 10)
    out[3] <- port[3]
    acf(fit$residuals, na.action = na.pass, lag.max = 10, 
        main = "", sub = paste("whiteness p-val: ", round(port[3], 
                                                          3)))
    plot(fit$N, fit$residuals, ylab = "residuals", xlab = "N", 
         sub = paste("p-val: ", round(prob, 3)))
    if (prob < 0.05) 
      abline(lmfit, lwd = 2, col = "red")
    else abline(lmfit, lwd = 1, lty = 2)
  }
}



