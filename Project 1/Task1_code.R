# ==============================================================================
# Question 1: KDE Comparison
# ==============================================================================

# Data Generation
n <- 500
sim.data <- rnorm(n, mean = 0, sd = 1)

# Kernel Functions
kern.gauss <- function(u) { (1 / sqrt(2 * pi)) * exp(-0.5 * u^2) }
kern.rect <- function(u) { ifelse(abs(u) <= 1, 0.5, 0) }
kern.tri <- function(u) { ifelse(abs(u) <= 1, 1 - abs(u), 0) }
kern.epan <- function(u) { ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0) }
kern.biweight <- function(u) { ifelse(abs(u) <= 1, (15/16) * (1 - u^2)^2, 0) }

# KDE Function
kde <- function(x, data, h, fn) {
  n <- length(data)
  est <- numeric(length(x))
  
  for (i in 1:length(x)) {
    u <- (x[i] - data) / h
    est[i] <- sum(fn(u)) / (n * h)
  }
  return(est)
}

# Bandwith using Silverman Equation
iqr.val <- IQR(sim.data)
sigma.hat <- min(sd(sim.data), iqr.val / 1.34)
h <- 1.06 * sigma.hat * n^(-1/5)

# True density of data
true.dens.data <- dnorm(sim.data)

est.gauss.data <- kde(sim.data, sim.data, h, kern.gauss)
est.rect.data  <- kde(sim.data, sim.data, h, kern.rect)
est.tri.data   <- kde(sim.data, sim.data, h, kern.tri)
est.epan.data  <- kde(sim.data, sim.data, h, kern.epan)
est.biwt.data  <- kde(sim.data, sim.data, h, kern.biweight)

# Plotting 
x.grid <- seq(-4, 4, length.out = 200)
est.gauss.grid <- kde(x.grid, sim.data, h, kern.gauss)
est.biwt.grid  <- kde(x.grid, sim.data, h, kern.biweight)
est.rect.grid  <- kde(x.grid, sim.data, h, kern.rect)
est.tri.grid  <- kde(x.grid, sim.data, h, kern.tri)
est.epan.grid <- kde(x.grid, sim.data, h, kern.epan)
true.dens.grid <- dnorm(x.grid) # True density for plotting

hist(sim.data, freq = FALSE, breaks = 30, col = "gray90", border = "gray",
     main = "KDE Comparison", xlab = "Value", ylim = c(0, 0.5))

# Add True Density
lines(x.grid, true.dens.grid, col = "darkgrey", lwd = 3, lty = 2)
lines(x.grid, est.gauss.grid, col = "blue", lwd = 2)
lines(x.grid, est.biwt.grid, col = "green", lwd = 2)
lines(x.grid, est.rect.grid, col = "red", lwd = 2, lty = 3)
lines(x.grid, est.epan.grid, col = "yellow", lwd = 2, lty = 3)
lines(x.grid, est.tri.grid, col = "black", lwd = 2, lty = 3)
legend("topright", 
       legend = c("True Density", "Gaussian", "Biweight", "Rectangular", "Epanechnikov", "Triweight"),
       lwd = c(3, 2, 2, 2, 2, 2), 
       lty = c(2, 1, 1, 3, 3, 3), 
       col = c("darkgrey", "blue", "green", "red", "yellow", "black"),
       cex = 0.5)

# Histogram Estimation for MSE and ISE
hist.obj <- hist(sim.data, breaks = 30, plot = FALSE)
bin.indices <- findInterval(sim.data, hist.obj$breaks, all.inside = TRUE)
est.hist.data <- hist.obj$density[bin.indices]

grid.bin.indices <- findInterval(x.grid, hist.obj$breaks)
est.hist.grid <- numeric(length(x.grid))
n_bins <- length(hist.obj$density)
valid.idx <- grid.bin.indices > 0 & grid.bin.indices <= n_bins
est.hist.grid[valid.idx] <- hist.obj$density[grid.bin.indices[valid.idx]]

# MSE Calculations & Comparison
mse.hist  <- mean((est.hist.data - true.dens.data)^2) # Histogram MSE
mse.gauss <- mean((est.gauss.data - true.dens.data)^2)
mse.rect  <- mean((est.rect.data - true.dens.data)^2)
mse.tri   <- mean((est.tri.data - true.dens.data)^2)
mse.epan  <- mean((est.epan.data - true.dens.data)^2)
mse.biwt  <- mean((est.biwt.data - true.dens.data)^2)

# Integrated Squared Error (ISE) Calculation 
dx <- x.grid[2] - x.grid[1] 

# Helper function for ISE
calc_ise <- function(est, true, dx) {
  return(sum((est - true)^2) * dx)
}

ise.hist  <- calc_ise(est.hist.grid, true.dens.grid, dx)
ise.gauss <- calc_ise(est.gauss.grid, true.dens.grid, dx)
ise.rect  <- calc_ise(est.rect.grid, true.dens.grid, dx)
ise.tri   <- calc_ise(est.tri.grid, true.dens.grid, dx)
ise.epan  <- calc_ise(est.epan.grid, true.dens.grid, dx)
ise.biwt  <- calc_ise(est.biwt.grid, true.dens.grid, dx)

# Comparison table
results <- data.frame(
  Estimator = c("Histogram", "Rectangular", "Triangular", "Epanechnikov", "Biweight", "Gaussian"),
  MSE = c(mse.hist, mse.rect, mse.tri, mse.epan, mse.biwt, mse.gauss),
  ISE = c(ise.hist, ise.rect, ise.tri, ise.epan, ise.biwt, ise.gauss)
)
results <- results[order(results$MSE), ] # Sort by best performance (lowest MSE)

cat("\n--- Performance Comparison (MSE) ---\n")
print(results, row.names = FALSE)

cat("\nNow we will perform a bootstrap procedure to compare Histogram vs KDE estimators.\n")

# Bootstrap Procedure to Compare Histogram vs KDE
B = 1000 # Number of Bootstrap Samples
bootstrap.kernel <- function(sim.data, kernel, B, m = 200, grid.min = -4, grid.max = 4) {
  
  n <- length(sim.data)
  
  x.grid <- seq(grid.min, grid.max, length.out = m)
  f.true <- dnorm(x.grid)
  
  ISE.hist <- numeric(B)
  ISE.kde  <- numeric(B)
  
  for (b in 1:B) {
    
    # Bootstrap resample
    x.boot <- sample(sim.data, size = n, replace = TRUE)
    
    # Histogram (Freedman–Diaconis)
    h.fd <- 2 * IQR(x.boot) / (n^(1/3))
    breaks <- seq(grid.min, grid.max, by = h.fd)
    
    hist.est <- hist(x.boot, breaks = breaks, plot = FALSE)
    
    f.hist <- approx(hist.est$mids, hist.est$density,
                     xout = x.grid,
                     rule = 1, yleft = 0, yright = 0)$y
    
    # Kernel Density Estimator
    f.kde <- kde(x.grid, x.boot, h, kernel)
    
    # Integrated Squared Error
    ISE.hist[b] <- mean((f.hist - f.true)^2)
    ISE.kde[b]  <- mean((f.kde  - f.true)^2)
  }
  
  list(
    kernel = kernel,
    kde.better.count = sum(ISE.kde < ISE.hist),
    ISE.hist = ISE.hist,
    ISE.kde  = ISE.kde
  )
}

res.gaussian <- bootstrap.kernel(sim.data, kern.gauss, B)
res.epanechnikov <- bootstrap.kernel(sim.data, kern.epan, B)
res.uniform <- bootstrap.kernel(sim.data, kern.rect, B)
res.triangular <- bootstrap.kernel(sim.data, kern.tri, B)
res.biweight <- bootstrap.kernel(sim.data, kern.biweight, B)

cat("Gaussian KDE better than Histogram:",
    res.gaussian$kde.better.count, "times out of", B, "\n")

cat("Epanechnikov KDE better than Histogram:",
    res.epanechnikov$kde.better.count, "times out of", B, "\n")

cat("Rectangular KDE better than Histogram:",
    res.uniform$kde.better.count, "times out of", B, "\n")

cat("Triangular KDE better than Histogram:",
    res.triangular$kde.better.count, "times out of", B, "\n")

cat("Biweight KDE better than Histogram:",
    res.biweight$kde.better.count, "times out of", B, "\n")

# Boxplots of ISE Comparisons
par(mfrow = c(2, 3))
boxplot(res.gaussian$ISE.hist, res.gaussian$ISE.kde,
        names = c("Histogram", "Gaussian KDE"),
        col = c("gray", "lightblue"),
        main = "Bootstrap ISE Comparison")

boxplot(res.epanechnikov$ISE.hist, res.epanechnikov$ISE.kde,
        names = c("Histogram", "Epanechnikov KDE"),
        col = c("gray", "lightblue"),
        main = "Bootstrap ISE Comparison")

boxplot(res.uniform$ISE.hist, res.uniform$ISE.kde,
        names = c("Histogram", "Rectangular KDE"),
        col = c("gray", "lightblue"),
        main = "Bootstrap ISE Comparison")

boxplot(res.triangular$ISE.hist, res.triangular$ISE.kde,
        names = c("Histogram", "Triangular KDE"),
        col = c("gray", "lightblue"),
        main = "Bootstrap ISE Comparison")    

boxplot(res.biweight$ISE.hist, res.biweight$ISE.kde,
        names = c("Histogram", "Biweight KDE"),       
        col = c("gray", "lightblue"),
        main = "Bootstrap ISE Comparison")

# ==============================================================================
# Question 2: Nadaraya-Watson vs Local Linear Estimator
# ==============================================================================

# Data Generation
n <- 100
x <- sort(runif(n, 0, 1))
true.m <- exp(x)

# Errors
eps.norm <- rnorm(n, 0, 1)
u.unif <- runif(n, -0.5, 0.5)
eps.lap <- -sign(u.unif) * log(1 - 2 * abs(u.unif))

y.norm <- true.m + eps.norm
y.lap  <- true.m + eps.lap

# Estimators
# Nadaraya-Watson
nadaraya.watson <- function(x, X, Y, h) {
  res <- numeric(length(x))
  for(i in 1:length(x)) {
    x0 <- x[i]
    w <- kern.gauss((X - x0) / h)
    res[i] <- sum(w * Y) / sum(w)
  }
  return(res)
}

# LLE Optim (Argmax)
loc.lin.opt <- function(x, X, Y, h) {
  res <- numeric(length(x))
  # Objective: Minimize weighted Sum of Squared Errors
  obj.fn <- function(beta, x0, X, Y, h) {
    w <- kern.gauss((X - x0) / h)
    preds <- beta[1] + beta[2] * (X - x0)
    return(sum(w * (Y - preds)^2))
  }
  for(i in 1:length(x)) {
    x0 <- x[i]
    opt.res <- optim(par = c(mean(Y), 0), fn = obj.fn, x0 = x0, X = X, Y = Y, h = h)
    res[i] <- opt.res$par[1]
  }
  return(res)
}

# LLE Closed Form (Actual Parameters)
loc.lin.closed <- function(x, X, Y, h) {
  res <- numeric(length(x))
  for(i in 1:length(x)) {
    x0 <- x[i]
    w <- kern.gauss((X - x0) / h)
    W.mat <- diag(w)
    X.mat <- cbind(1, X - x0)
    
    # Solve (X'WX)beta = X'WY
    XtW <- t(X.mat) %*% W.mat
    lhs <- XtW %*% X.mat
    rhs <- XtW %*% Y
    
    beta.hat <- tryCatch({ solve(lhs, rhs) }, error = function(e) c(NA, NA))
    res[i] <- beta.hat[1]
  }
  return(res)
}

# Execution
h <- 0.1
est.nw.norm  <- nadaraya.watson(x, x, y.norm, h)
est.lle.norm <- loc.lin.closed(x, x, y.norm, h)
est.nw.lap   <- nadaraya.watson(x, x, y.lap, h)
est.lle.lap  <- loc.lin.closed(x, x, y.lap, h)

# MSE
mse.nw.norm  <- mean((est.nw.norm - true.m)^2)
mse.lle.norm <- mean((est.lle.norm - true.m)^2)
mse.nw.lap   <- mean((est.nw.lap - true.m)^2)
mse.lle.lap  <- mean((est.lle.lap - true.m)^2)

cat("\n--- Regression MSE Results ---\n")
cat("NW (Normal):  ", mse.nw.norm, "\n")
cat("LLE (Normal): ", mse.lle.norm, "\n")
cat("NW (Laplace): ", mse.nw.lap, "\n")
cat("LLE (Laplace):", mse.lle.lap, "\n")

# Plotting
par(mfrow = c(1, 2))
plot(x, y.norm, pch = 20, col = "gray", main = "Normal Errors")
lines(x, true.m, lty = 2); lines(x, est.nw.norm, col = "green"); lines(x, est.lle.norm, col = "orange")
legend("topleft", legend = c("True", "NW", "LLE"), col = c("black", "green", "orange"), 
       lty = c(2, 1, 1), cex = 0.6)

plot(x, y.lap, pch = 20, col = "gray", main = "Laplace Errors")
lines(x, true.m, lty = 2); lines(x, est.nw.lap, col = "green"); lines(x, est.lle.lap, col = "orange")
legend("topleft", legend = c("True", "NW", "LLE"), col = c("black", "green", "orange"),
       lty = c(2, 1, 1), cex = 0.6)

bootstrap.compare <- function(X, Y, true.m, h, B = 1000) {
  
  n <- length(X)
  
  mse.nw  <- numeric(B)
  mse.ll  <- numeric(B)
  
  for (b in 1:B) {
    idx <- sample(1:n, size = n, replace = TRUE)
    
    Xb <- X[idx]
    Yb <- Y[idx]
    true.mb <- true.m[idx]
    
    # Sort (important for smoothness)
    ord <- order(Xb)
    Xb <- Xb[ord]
    Yb <- Yb[ord]
    true.mb <- true.mb[ord]
    
    est.nw <- nadaraya.watson(Xb, Xb, Yb, h)
    est.ll <- loc.lin.closed(Xb, Xb, Yb, h)
    
    mse.nw[b] <- mean((est.nw - true.mb)^2, na.rm = TRUE)
    mse.ll[b] <- mean((est.ll - true.mb)^2, na.rm = TRUE)
  }
  
  list(
    mse.nw = mse.nw,
    mse.ll = mse.ll,
    prob.ll.better = mean(mse.ll < mse.nw),
    prob.nw.better = mean(mse.nw < mse.ll)
  )
}
boot.norm <- bootstrap.compare(x, y.norm, true.m, h, B = 1000)
boot.lap <- bootstrap.compare(x, y.lap, true.m, h, B = 1000)


cat("\n--- Bootstrap Results (10000 replications) ---\n")

cat("\nNormal errors:\n")
cat("P(LL better than NW): ", boot.norm$prob.ll.better, "\n")
cat("P(NW better than LL): ", boot.norm$prob.nw.better, "\n")

cat("\nLaplace errors:\n")
cat("P(LL better than NW): ", boot.lap$prob.ll.better, "\n")
cat("P(NW better than LL): ", boot.lap$prob.nw.better, "\n")

par(mfrow = c(1, 2))

# Normal errors
boxplot(
  boot.norm$mse.nw,
  boot.norm$mse.ll,
  names = c("NW", "LL"),
  main = "MSE Comparison (Normal Errors)",
  ylab = "MSE",
  col = c("lightblue", "lightgreen")
)

# Laplace errors
boxplot(
  boot.lap$mse.nw,
  boot.lap$mse.ll,
  names = c("NW", "LL"),
  main = "MSE Comparison (Laplace Errors)",
  ylab = "MSE",
  col = c("lightblue", "lightgreen")
)

