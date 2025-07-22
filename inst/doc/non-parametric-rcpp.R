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

clock <- function(expr, n = 20) {
  # Capturing the output into the warm-up iteration to retuce the timing
  # overhead from spinning up the processor from any sleep or idle state
  # (idea from the microbenchmark package)
  ret <- expr
  tic0 <- proc.time()
  replicate(n, suppressWarnings(suppressMessages(expr)))
  tm <- proc.time() - tic0
  et <- as.numeric(tm["elapsed"]) * 1000
  times <- c(median = et, mean = mean(et), sd = sd(et))
  ftimes <- ifelse(times > 1, sprintf("%1.1f", times),
                   ifelse(times == 1, "~0.5-1.5", "<0.5"))
  ftimes[3] <- sprintf("%1.1f", times[3])
  cat("Median [mean +- SD] time per evaluation: ", ftimes[1], " [",
      ftimes[2], "+-",  "]", " ms\n", sep = "")
  attr(ret, "time") <- times
  return(ret)
}

quiet <- function(expr) {  # Suppressing all output
  sink(if (.Platform$OS.type == "windows") "NUL" else "/dev/null")
  ret <- expr
  sink()
  return(ret)
}

## -----------------------------------------------------------------------------
# x = numeric vector, xout = numeric vector, bw = scalar, kfun = function
kernelWeights1R <- function(x, xout, bw, kfun) kfun(outer(xout, x, "-") / bw)

## -----------------------------------------------------------------------------
kernelWeightsR <- function(x, xout = NULL, bw = 1, kfun = stats::dnorm) {
  if (is.null(dim(x))) x <- as.matrix(x)
  if (is.null(xout)) xout <- x
  if (is.null(dim(xout))) xout <- as.matrix(xout)
  d <- ncol(x)
  if (d > 1 && length(bw) == 1) bw <- rep(bw, d)
  pk <- kernelWeights1R(x = x[, 1], xout = xout[, 1], bw = bw[1], kfun = kfun)
  if (d > 1) { # Accumulating the product kernel
    for (i in 2:d)
      pk <- pk * kernelWeights1R(x = x[, i], xout = xout[, i], bw = bw[i], kfun = kfun)
  }
  return(pk)
}

## -----------------------------------------------------------------------------
kernelDensityR <- function(x, xout = NULL, bw = 1, kfun = stats::dnorm) {
  x1d <- is.null(dim(x))
  n <- if (x1d) length(x) else nrow(x)
  if (isTRUE(ncol(x) > 1) && length(bw) == 1) bw <- rep(bw, ncol(x))
  pk <- kernelWeightsR(x = x, xout = xout, bw = bw, kfun = kfun)
  return(rowSums(pk) / (n * prod(bw)))
}

## -----------------------------------------------------------------------------
kernelSmoothR <- function(x, y, xout = NULL, bw = 1,
                          kfun = stats::dnorm, LOO = FALSE) {
  pk <- kernelWeightsR(x = x, xout = xout, bw = bw, kfun = kfun)
  if (LOO) diag(pk) <- 0
  return(rowSums(sweep(pk, 2, y, "*")) / rowSums(pk))
}

## -----------------------------------------------------------------------------
set.seed(1)
X <- sort(rchisq(300, 3))
xg <- seq(0, max(X), length.out = 100)
Y <- sin(X) + rnorm(300)
b <- 0.5

## ----message=FALSE------------------------------------------------------------
wCPP <- clock(kernelWeights(X, xout = xg, bw = b))
wR   <- clock(kernelWeightsR(X, xout = xg, bw = b))

fCPP <- clock(kernelDensity(X, xout = xg, bw = b))
fR   <- clock(kernelDensityR(X, xout = xg, bw = b))

mCPP <- clock(kernelSmooth(X, Y, xout = xg, bw = b))
mR   <- clock(kernelSmoothR(X, Y, xout = xg, bw = b))

## -----------------------------------------------------------------------------
all.equal(wR, wCPP, tolerance = 1e-16)
all.equal(fR, fCPP, tolerance = 1e-16)
all.equal(mR, mCPP, tolerance = 1e-16)

## -----------------------------------------------------------------------------
oldpar <- par(mfrow = c(2, 1), mar = c(4, 4, 2, 3))
plot(X, Y, bty = "n", main = "Non-parametric regression (black) and density (blue) estimate", las = 1)
lines(xg, mCPP, lwd = 2)
par(new = TRUE)
plot(xg, fCPP, type = "l", lwd = 2, yaxt = "n", bty = "n", xaxt = "n", col = 4, xlab = "", ylab = "", las = 1)
axis(4, col = 4, col.axis = 4, las = 1)
plot(xg, fR - fCPP, bty = "n", ylab = "", xlab = "X",
     main = "Discrepancy between the R and C++ implementation")

## ----message=FALSE------------------------------------------------------------
par(mfrow = c(1, 1))
fCPPunif <- kernelDensity(X, xout = xg, bw = b * sqrt(3), kernel = "uniform")
fCPPtrng <- kernelDensity(X, xout = xg, bw = b * sqrt(3), kernel = "triangular")
fCPPepan <- kernelDensity(X, xout = xg, bw = b * sqrt(5), kernel = "epanechnikov")
plot(xg, fCPP, ylim = range(0, fCPP, fCPPunif, fCPPtrng, fCPPepan), xlab = "X",
     type = "l", bty = "n", main = "Various kernel shapes", ylab = "Density")
rug(X)
lines(xg, fCPPunif, col = 2)
lines(xg, fCPPtrng, col = 3)
lines(xg, fCPPepan, col = 4)
legend("topright", c("Gaussian", "Uniform", "Triangular", "Epanechnikov"), lwd = 1, col = 1:4)

## ----message=FALSE------------------------------------------------------------
ns <- c(500, 1000, 2000)
timings <- sapply(ns, function(n) {
  set.seed(1)
  X <- rchisq(n, 3)
  b <- bw.rot(X)
  tR <- quiet(clock(fR <- kernelDensityR(X, bw = b), n = 3))
  tCPP <- quiet(clock(fCPP <- kernelDensity(X, bw = b), n = 3))
  c(R = attr(tR, "time"), CPP = attr(tCPP, "time"))
})
colnames(timings) <- paste0("n=", ns)
print(timings, 2)

## -----------------------------------------------------------------------------
set.seed(1)
n <- 100
ng <- 30
X <- matrix(rchisq(n*2, 3), ncol = 2)
xg0 <- seq(0, 13, length.out = ng)
xg <- expand.grid(X1 = xg0, X2 = xg0)
Y <- sin(X[, 1]) + sin(X[, 2]) + rnorm(n)
b <- c(0.7, 0.8)

wCPP <- clock(kernelWeights(X, xout = xg, bw = b), n = 10)
wR <- clock(kernelWeightsR(X, xout = xg, bw = b), n = 10)
all.equal(wR, wCPP, tolerance = 1e-16)

## -----------------------------------------------------------------------------
max.dw <- apply(wR - wCPP, 1, function(x) max(abs(x), na.rm = TRUE))
filled.contour(xg0, xg0, log10(matrix(max.dw, ng, ng)), xlab = "X1", ylab = "X2",
               main = "Log10(discrepancy) between R and C++ kernel weights", asp = 1)
points(X[, 1], X[, 2], pch = 16)

## ----message=FALSE------------------------------------------------------------
fCPP <- clock(kernelDensity(X, xout = xg, bw = b), n = 5)
fR <- clock(kernelDensityR(X, xout = xg, bw = b), n = 5)
all.equal(fR, fCPP, tolerance = 1e-16)

## ----message=FALSE------------------------------------------------------------
fCPPunif <- kernelDensity(X, xout = xg, bw = b * sqrt(3), kernel = "uniform")
fCPPtrng <- kernelDensity(X, xout = xg, bw = b * sqrt(3), kernel = "triangular")
fCPPepan <- kernelDensity(X, xout = xg, bw = b * sqrt(5), kernel = "epanechnikov")
par(mfrow = c(2, 2), mar = c(0.5, 0.5, 2, 0.5))
persp(xg0, xg0, matrix(fCPP, nrow = ng), theta = 120, phi = 20, main = "Gaussian", xlab = "X1", ylab = "X2", zlab = "Density")
persp(xg0, xg0, matrix(fCPPunif, nrow = ng), theta = 120, phi = 20, main = "Uniform", xlab = "X1", ylab = "X2", zlab = "Density")
persp(xg0, xg0, matrix(fCPPtrng, nrow = ng), theta = 120, phi = 20, main = "Triangular", xlab = "X1", ylab = "X2", zlab = "Density")
persp(xg0, xg0, matrix(fCPPepan, nrow = ng), theta = 120, phi = 20, main = "Epanechnikov", xlab = "X1", ylab = "X2", zlab = "Density")

## ----message=FALSE------------------------------------------------------------
ns <- c(125, 500)
dims <- c(2, 6)
nd <- expand.grid(n = ns, dim = dims)

timings <- t(sapply(seq_len(nrow(nd)), function(i) {
  set.seed(1)
  X <- matrix(rchisq(nd$n[i] * nd$dim[i], 3))
  b <- bw.rot(X)
  aR <- system.time(kernelDensityR(X, bw = b))["elapsed"]
  aCPP <- system.time(kernelDensity(X, bw = b))["elapsed"]
  c(R = aR, CPP = aCPP)
}))
timings <- cbind(nd, timings)
print(timings, 2)

## -----------------------------------------------------------------------------
bwgrid <- seq(0.2, 3, .1)
bws0 <- suppressWarnings(LSCV(X, Y, bw = bwgrid, degree = 0, same = TRUE))
bws1 <- suppressWarnings(LSCV(X, Y, bw = bwgrid, degree = 1, same = TRUE))
bws2 <- suppressWarnings(LSCV(X, Y, bw = bwgrid, degree = 2, same = TRUE))
bw.cv <- cbind(bws0, bws1, bws2)
bw.opt <- bwgrid[apply(bw.cv, 2, which.min)]
par(mar = c(4, 4, 0.5, 0.5))
matplot(bwgrid, bw.cv, lty = 1, type = "l", bty = "n", lwd = 2, col = 1:3,
        xlab = "Bandwidth", ylab = "Out-of-sample prediction error", log = "y")
legend("topright", paste("Degree", 0:2), lwd = 2, lty = 1, col = 1:3, bty = "n")

## ----message=FALSE------------------------------------------------------------
xg0 <- seq(0, 8, length.out = ng)
xg <- expand.grid(X1 = xg0, X2 = xg0)
yhat <- sapply(1:3, function(i) kernelSmooth(x = X,  y = Y, xout = xg, bw = bw.opt[i], degree = i-1))
par(mfrow = c(2, 2), mar = c(0.2, 0.2, 2, 0.2))
for (i in 1:3) {
  p <- persp(xg0, xg0, matrix(yhat[, i], ng, ng), ticktype = "detailed",
             main = paste0("Degree ", i-1), theta = 30, phi = 25,
             xlab = "X1", ylab = "X2", zlab = "Y", zlim = range(yhat, Y))
  points(trans3d(X[, 1], X[, 2], Y, p), col = 2)
}

## ----message=FALSE------------------------------------------------------------
ys <- sapply(1:3, function(i) kernelSmooth(x = X,  y = Y, bw = bw.opt[i], degree = i-1))
colnames(ys) <- paste("Degree", 0:2)
print(cor(cbind(Y, ys))[1, 2:4], 3)

## -----------------------------------------------------------------------------
set.seed(1)
n <- 1000
X <- sort(rnorm(n))
w1 <- kernelWeights(X, bw = 0.1)
w2 <- kernelWeights(X, bw = 0.1, sparse = TRUE)
print(c(object.size(w1), object.size(w2)))
print(c(class(w1)[1], class(w2)[1]))

## -----------------------------------------------------------------------------
print(c(object.size(kernelWeights(X, bw = 0.5)),
        object.size(kernelWeights(X, bw = 0.5, sparse = TRUE))))

## -----------------------------------------------------------------------------
print(c(object.size(kernelWeights(X, kernel = "epanechnikov", bw = 0.6)),
        object.size(kernelWeights(X, kernel = "epanechnikov", bw = 0.6, sparse = TRUE))))

## -----------------------------------------------------------------------------
set.seed(1)
nz <- c(1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1)
m <- matrix(runif(16) * nz, nrow = 4, byrow = TRUE)
print(m)
print(apply(m, 1, sparseVectorToList, renormalise = FALSE))

## -----------------------------------------------------------------------------
kw <- kernelWeights(X, kernel = "epanechnikov", bw = 0.6)
print(c(object.size(kw),
        object.size(kernelWeights(X, kernel = "epanechnikov", bw = 0.6, sparse = TRUE)),
        object.size(apply(kw, 1, sparseVectorToList, renormalise = FALSE))))

## ----include=FALSE------------------------------------------------------------
par(oldpar)

