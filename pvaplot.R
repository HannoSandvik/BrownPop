
# OBS: Det er (minst!) to versjon av funksjonen i omløp.
# Foreløpig kaller jeg dem helt vilkårlig "pva.plot" og "pva.plot2"

# ikke rydda ennå!


pva.plot <- function(model, pva, years, traj=100, tradj=0, log=T, C=20, CI, xlab, ylab, lwd, col, lty, ...)
{N  <- model$N[[1]]
yr <- as.numeric(names(N))
names(N) <- NULL
load(pva)
if (missing(CI))   CI <- c(0.6, 0.8, 0.9)
if (missing(years)) years <- min(yr):(max(yr) + 100)
if (min(years) > max(yr))
{N  <- numeric(0)}
else
{N  <- N[(yr %A% years) - yr[1] + 1]}
yrsi <- years %A% (max(yr):max(years))
Y    <- length(yrsi)
traj <- 1:min(traj, length(sik)) + tradj
aar  <- 1:length(yr)
m    <- matrix(0, length(sik),  Y)
ci   <- matrix(0, length(CI)*2, Y)
for (i in 1:length(sik)) m[i, 1:min(Y, length(sik[[i]]))] <- sik[[i]][1:min(Y, length(sik[[i]]))]
rm(sik)
if (missing(xlab)) xlab <- "År"
if (log)
{N <- ln(N)
m[m <= ln(C)] <- ln(C)
if (missing(ylab)) ylab <- "Logaritme av bestandsstørrelsen"
ylim <- c(ln(C), max(N, m[traj,], na.rm=T))}
else
{m <- exp(m)
m[m < 0.12] <- 0
if (missing(ylab)) ylab <- "Bestandsstørrelse"
ylim <- c(0, max(N, m[traj,], na.rm=T))}
for (i in 1:length(CI))
{ci[i*2-1,] <- apply(m, 2, quantile, 0.5 - CI[i] / 2)
ci[i*2,]   <- apply(m, 2, quantile, 0.5 + CI[i] / 2)}
med  <- apply(m, 2, median)
m    <- m[traj,]
if (missing(lwd)) lwd <- c(2,1,2,1)
if (missing(col)) col <- c("black", grey(0.72), "black", "black")
if (missing(lty)) lty <- rep(1,4)
plot(yr %A% years, N, xlim=range(years), ylim=ylim, xlab=xlab, ylab=ylab,
     type="l", lwd=lwd[1], col=col[1], lty=lty[1], ...)
for (i in traj) lines(yrsi,  m[i,], lwd=lwd[2], col=col[2], lty=lty[2])
lines(yrsi, med, lwd=lwd[3], col=col[3], lty=lty[3])
for (i in 1:nrow(ci)) lines(yrsi, ci[i,], lwd=lwd[4], col=col[4], lty=lty[4])}



pva.plot2 <- function(model, pva, years, traj=100, log=T, C=20, CI, xlab, ylab, lwd, col, lty, ...)
{N  <- model$N[[1]]
yr <- as.numeric(names(N))
names(N) <- NULL
load(pva)
if (missing(CI))   CI <- c(0.6, 0.8, 0.9)
if (missing(years)) years <- min(yr):(max(yr) + 100)
yrsi <- years %A% (max(yr):max(years))
Y    <- length(yrsi)
traj <- min(traj, length(sik))
m    <- matrix(0, length(sik),  Y)
ci   <- matrix(0, length(CI)*2, Y)
for (i in 1:length(sik)) m[i, 1:min(Y, length(sik[[i]]))] <- sik[[i]][1:min(Y, length(sik[[i]]))]
rm(sik)
if (missing(xlab)) xlab <- "?r"
if (log)
{N <- ln(N)
m[m <= ln(C)] <- ln(C)
if (missing(ylab)) ylab <- "Logaritme av bestandsst?rrelsen"
ylim <- c(ln(C), max(N, m[1:traj,], na.rm=T))}
else
{m <- exp(m)
m[m < 0.12] <- 0
if (missing(ylab)) ylab <- "Bestandsst?rrelse"
ylim <- c(0, max(N, m[1:traj,], na.rm=T))}
for (i in 1:length(CI))
{ci[i*2-1,] <- apply(m, 2, quantile, 0.5 - CI[i] / 2)
ci[i*2,]   <- apply(m, 2, quantile, 0.5 + CI[i] / 2)}
med  <- apply(m, 2, median)
m    <- m[1:traj,]
if (missing(lwd)) lwd <- c(2,1,2,1)
if (missing(col)) col <- c("black", grey(0.72), "black", "black")
if (missing(lty)) lty <- rep(1,4)
plot(yr %A% years, N[(yr %A% years) - yr[1] + 1], xlim=range(years), ylim=ylim, xlab=xlab, ylab=ylab,
     type="l", lwd=lwd[1], col=col[1], lty=lty[1], ...)
for (i in 1:traj) lines(yrsi,  m[i,], lwd=lwd[2], col=col[2], lty=lty[2])
lines(yrsi, med, lwd=lwd[3], col=col[3], lty=lty[3])
for (i in 1:nrow(ci)) lines(yrsi, ci[i,], lwd=lwd[4], col=col[4], lty=lty[4])}



