## ----include = FALSE----------------------------------------------------------
knitr::knit_hooks$set(pngquant = knitr::hook_pngquant)
knitr::opts_chunk$set(
  dev = "png",
  dev.args = list(type = if (Sys.info()["sysname"] == "Darwin") "quartz" else "cairo-png"),
  fig.width = 640 / 72,
  fig.height = 480 / 72,
  dpi = 72,
  fig.retina = 1,
  collapse = TRUE,
  comment = "#>",
  pngquant = "--speed=1 --quality=50-60"
)

## ----setup--------------------------------------------------------------------
library(smoothemplik)

quiet <- function(expr) {  # Suppressing all output
  sink(if (.Platform$OS.type == "windows") "NUL" else "/dev/null")
  ret <- expr
  sink()
  return(ret)
}

## -----------------------------------------------------------------------------
x  <- c(-3, -2, 2, 3, 4)
ct <- c(10, 4:1)
grid.full <- seq(-3.5, 5.5, 0.05)
grid.keep <- grid.full < -2.25 | grid.full > 3.25
selr0 <- sapply(grid.full, function(m) -2*EL0(x, ct, mu = m, chull.fail = "none")$logelr)
selr1 <- sapply(grid.full, function(m) -2*EL0(x, ct, mu = m, chull.fail = "taylor")$logelr)
selr2 <- sapply(grid.full, function(m) -2*EL0(x, ct, mu = m, chull.fail = "wald")$logelr)
plot(grid.full, selr0, ylim = c(0, 120), xlim = c(-3.5, 4.5), bty = "n",
     main = "-2 * weighted log-EL", ylab = "")
points(grid.full[grid.keep], selr1[grid.keep], col = 4, pch = 0)
wm <- weighted.mean(x, ct)
wv <- weighted.mean((x - wm)^2, ct) / sum(ct)
points(grid.full[grid.keep], selr2[grid.keep], col = 2, pch = 2)
lines(grid.full, (grid.full - wm)^2 / wv, col = 2, pch = 2)
rug(x, lwd = 2)
abline(v = c(-2.5, 3.5), lty = 3)

## -----------------------------------------------------------------------------
pickMinBW <- function(X, d = 2 / sqrt(length(X))) {
  n <- length(X)
  nd <- ceiling(d*n)
  b <- numeric(n)
  gaps <- pmin(c(NA, diff(X)), c(diff(X), NA))
  gaps[1] <- X[2] - X[1]
  gaps[n] <- X[n] - X[n-1]
  for (i in 1:n) {
    bw.start <- gaps[i] * 2
    obs.neigh <- function(bw) sum(abs(X[i] - X) < bw) - 1
    b[i] <- uniroot(function(bw) obs.neigh(bw) - nd, c(bw.start/2, bw.start/1.5), extendInt = "upX")$root
  }
  b
}

## -----------------------------------------------------------------------------
n <- 50
set.seed(1)
X <- sort(rchisq(n, df = 3))
Y <- 1 + X + (rchisq(n, df = 3) - 3) * (1 + X)
mod0 <- lm(Y ~ X)
vhat0 <- kernelSmooth(X, mod0$residuals^2, bw = max(diff(X))*1.2, kernel = "epanechnikov")
mod1 <- lm(Y ~ X, weights = 1 / vhat0)
cbind(OLS = mod0$coefficients, WOLS = mod1$coefficients)

## -----------------------------------------------------------------------------
bw0 <- bw.CV(X, kernel = "epanechnikov")
bw0 <- max(bw0, max(diff(X))*1.1)
wF <- kernelWeights(X, bw = bw0, kernel = "epanechnikov")

# Assuming Gaussian CDF, which is not true
Xs <- scale(X)
XP <- pnorm(Xs)
wP <- kernelWeights(XP, bw = bw.CV(XP), kernel = "epanechnikov")
rowMeans(wP > 0) - 1/n

bw.adapt <- pickMinBW(X, d = 0.15)
plot(X, bw.adapt, bty = "n", main = "Adaptive bandwidth ensuring 15% non-zero weights")
wA <- kernelWeights(X, bw = bw.adapt, kernel = "epanechnikov")
rowMeans(wA > 0) - 1/n

rX <- rank(X)
wNN <- kernelWeights(rX/max(rX), bw = 0.09, kernel = "epanechnikov")
rowMeans(wNN > 0) - 1/n

## -----------------------------------------------------------------------------
g  <- function(theta, ...) Y - theta[1] - theta[2]*X

wF <- wF / rowSums(wF)
wP <- wP / rowSums(wP)
wNN <- wNN / rowSums(wNN)
wA <- wA / rowSums(wA)

g1 <- function(theta) smoothEmplik(rho = g, theta = theta, data = NULL, sel.weights = wF, minus = TRUE)
g2 <- function(theta) smoothEmplik(rho = g, theta = theta, data = NULL, sel.weights = wP, minus = TRUE)
g3 <- function(theta) smoothEmplik(rho = g, theta = theta, data = NULL, sel.weights = wNN, minus = TRUE)
g4 <- function(theta) smoothEmplik(rho = g, theta = theta, data = NULL, sel.weights = wA, minus = TRUE)
  
th1  <- optim(mod1$coefficients, g1, method = "BFGS", control = list(ndeps = rep(1e-5, 2)))
th2  <- optim(mod1$coefficients, g2, method = "BFGS", control = list(ndeps = rep(1e-5, 2)))
th3  <- optim(mod1$coefficients, g3, method = "BFGS", control = list(ndeps = rep(1e-5, 2)))
th4  <- optim(mod1$coefficients, g4, method = "BFGS", control = list(ndeps = rep(1e-5, 2)))

cbind(WOLS = mod1$coefficients, Fixed = th1$par, PIT = th2$par, NNeighb = th3$par, Adapt = th4$par)

