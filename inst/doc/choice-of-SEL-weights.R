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

## -----------------------------------------------------------------------------
x  <- c(-3, -2, 2, 3, 4)
ct <- c(10, 4:1)
grid.full <- c(seq(-3.5, 5.5, length.out = 201))
grid.keep <- grid.full <= -2.5 | grid.full >= 3.5
selr0 <- sapply(grid.full, function(m) -2*EL0(x, m, ct)$logelr)
selr1 <- -2*ExEL1(x, mu = grid.full, ct = ct, exel.control = list(xlim = c(-2.5, 3.5)))
selr2 <- -2*ExEL2(x, mu = grid.full, ct = ct, exel.control = list(xlim = c(-2.5, 3.5)))
plot(grid.full, selr0, ylim = c(0, 120), xlim = c(-3.5, 4.5), bty = "n",
     main = "-2 * weighted log-EL", ylab = "", type = "l")
points(grid.full[grid.keep], selr1[grid.keep], col = 4, pch = 0, cex = 0.75)
points(grid.full[grid.keep], selr2[grid.keep], col = 2, pch = 2, cex = 0.75)
rug(x, lwd = 2)
abline(v = c(-2.5, 3.5), lty = 3)
legend("top", c("Taylor", "Wald"), title = "Extrapolation type", col = c(4, 2),
       pch = c(0, 2), bty = "n")

## -----------------------------------------------------------------------------
pickMinBW <- function(X, d = 2 / sqrt(NROW(X)), tol = 1e-12) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  k <- ceiling(d * n)
  k <- max(1, min(k, n - 1))

  # has.ms <- requireNamespace("matrixStats", quietly = TRUE)
  b <- sapply(seq_len(n), function(i) {
    # L_inf distances from X_i
    A  <- abs(sweep(X, 2, X[i, ], "-", check.margin = FALSE))
    # Di <- if (has.ms) matrixStats::rowMaxs(A) else do.call(pmax, as.data.frame(A))
    Di <- do.call(pmax, as.data.frame(A))
    Di[i] <- Inf  # Exclude self
    if (k < n - 1) {
      return(sort(Di, partial = k + 1)[k + 1])
    } else {
      return(sort(Di, partial = n - 1)[n - 1])
    }
  })
  b * (1 - tol)
}

## -----------------------------------------------------------------------------
n <- 50
set.seed(1)
X <- sort(rchisq(n, df = 3))
Y <- 1 + X + (rchisq(n, df = 3) - 3) * (1 + X)
mod0 <- lm(Y ~ X)
vhat0 <- kernelSmooth(X, mod0$residuals^2, bw = max(diff(X))*1.2, kernel = "epanechnikov")
mod1 <- lm(Y ~ X, weights = 1 / vhat0)
plot(X, Y, bty = "n")
abline(c(1, 1), lty = 2)
abline(mod0, col = 1); abline(mod1, col = 2)
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

bw.adapt <- pickMinBW(X, d = 0.16)
plot(X, bw.adapt, bty = "n", main = "Adaptive bandwidth ensuring 15% non-zero weights",
     ylim = c(0, max(bw.adapt)))
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

g1 <- function(theta) smoothEmplik(rho = g, theta = theta, data = NULL, sel.weights = wF, minus = TRUE, type = "EuL")
g2 <- function(theta) smoothEmplik(rho = g, theta = theta, data = NULL, sel.weights = wP, minus = TRUE, type = "EuL")
g3 <- function(theta) smoothEmplik(rho = g, theta = theta, data = NULL, sel.weights = wNN, minus = TRUE, type = "EuL")
g4 <- function(theta) smoothEmplik(rho = g, theta = theta, data = NULL, sel.weights = wA, minus = TRUE, type = "EuL")
  
th1  <- optim(mod1$coefficients, g1, method = "BFGS", control = list(ndeps = rep(1e-5, 2), reltol = 1e-5))
th2  <- optim(mod1$coefficients, g2, method = "BFGS", control = list(ndeps = rep(1e-5, 2), reltol = 1e-5))
th3  <- optim(mod1$coefficients, g3, method = "BFGS", control = list(ndeps = rep(1e-5, 2), reltol = 1e-5))
th4  <- optim(mod1$coefficients, g4, method = "BFGS", control = list(ndeps = rep(1e-5, 2), reltol = 1e-5))

cbind(WOLS = mod1$coefficients, Fixed = th1$par, PIT = th2$par, NNeighb = th3$par, Adapt = th4$par)

## -----------------------------------------------------------------------------
n <- 200
set.seed(1)
X <- matrix(rchisq(n*3, df = 3), ncol = 3)
Y <- 1 + rowSums(X) + (rchisq(n, df = 3) - 3) * (1 + rowSums(X))
mod0 <- lm(Y ~ X)
vhat0 <- kernelSmooth(X, mod0$residuals^2, PIT = TRUE, bw = 0.2,
                      kernel = "epanechnikov", no.dedup = TRUE)
mod1 <- lm(Y ~ X, weights = 1 / vhat0)
cbind(OLS = mod0$coefficients, WOLS = mod1$coefficients)

## -----------------------------------------------------------------------------
bw0 <- 1
atleast2 <- FALSE
while(!atleast2) {
 wF <- kernelWeights(X, bw = bw0, kernel = "epanechnikov")
 bw0 <- bw0 * 1.1
 atleast2 <- min(rowSums(wF > 0)) > 1
}
rowMeans(wF > 0) - 1/n
min(rowSums(wF > 0))

Xs <- scale(X)
XP <- apply(Xs, 2, pnorm)
wP <- kernelWeights(XP, bw = bw.rot(XP)*2, kernel = "epanechnikov")
rowMeans(wP > 0) - 1/n
min(rowSums(wP > 0))

bw.adapt <- pickMinBW(X, d = 0.16)
plot(rowMeans(X), bw.adapt, bty = "n", main = "Adaptive bandwidth ensuring 15% non-zero weights",
     ylim = c(0, max(bw.adapt)))
wA <- kernelWeights(X, bw = bw.adapt, kernel = "epanechnikov")
rowMeans(wA > 0) - 1/n
min(rowSums(wA > 0))

## -----------------------------------------------------------------------------
g  <- function(theta, ...) Y - drop(cbind(1, X) %*% theta)

wF <- wF / rowSums(wF)
wP <- wP / rowSums(wP)
wA <- wA / rowSums(wA)

g1 <- function(theta) smoothEmplik(rho = g, theta = theta, data = NULL, sel.weights = wF,
                                   minus = TRUE, type = "EuL")
g2 <- function(theta) smoothEmplik(rho = g, theta = theta, data = NULL, sel.weights = wP,
                                   minus = TRUE, type = "EuL")
g4 <- function(theta) smoothEmplik(rho = g, theta = theta, data = NULL, sel.weights = wA,
                                   minus = TRUE, type = "EuL")

## -----------------------------------------------------------------------------
# th1  <- optim(mod1$coefficients, g1, method = "BFGS", control = list(ndeps = rep(1e-5, 4)))
# th2  <- optim(mod1$coefficients, g2, method = "BFGS", control = list(ndeps = rep(1e-5, 4)))
# th4  <- optim(mod1$coefficients, g4, method = "BFGS", control = list(ndeps = rep(1e-5, 4)))

# round(cbind(True = rep(1, 4), WOLS = mod1$coefficients, Fixed = th1$par, PIT = th2$par, Adapt = th4$par), 3)
#             True  WOLS Fixed    PIT Adapt
# (Intercept)    1 0.166 2.755  3.668 0.110
# X1             1 0.059 0.132  0.241 0.790
# X2             1 0.991 1.036  0.694 0.929
# X3             1 0.752 0.625 -0.161 0.863

