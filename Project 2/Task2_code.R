
# ==========================================
#       HW 2 Question 1
# ==========================================
# Draw the regression quantile lines for various choices of quantile indexes.


#loading required libraries
if(!require(quantreg)) install.packages("quantreg")
library(quantreg)
library(ggplot2)
library(dplyr)
library(plotly)

#load data set
walmart <- read.csv("Walmart_Store_sales.csv")
str(walmart)
summary(walmart)

#checking for outliers
boxplot(walmart$Weekly_Sales, main = "Boxplot of Weekly Sales", ylab = "Weekly Sales", col = "lightblue")

#OLS(baseline)
ols_model <- lm(Weekly_Sales ~ Temperature + Fuel_Price + CPI + Unemployment, data = walmart)
summary(ols_model)

#quantile regression models
taus <- c(0.05, 0.25, 0.5, 0.75, 0.95)

quantile_regression_models <- lapply(taus, function(tau) {
  rq(Weekly_Sales ~ Temperature + Fuel_Price + CPI + Unemployment, tau = tau, data = walmart)
})

#comparing coefficients across different regression quantiles
coef_table <- sapply(quantile_regression_models, coef)
colnames(coef_table) <- paste0("tau=", taus)
round(coef_table, 2)

#calculating mse
compute_mse <- function(y_true, y_pred) {
  mean((y_true - y_pred)^2)
}

#mse for ols
ols_pred <- predict(ols_model, newdata = walmart)
mse_ols <- compute_mse(walmart$Weekly_Sales, ols_pred)
mse_ols

#mse for quantile regressions
mse_quantiles <- sapply(seq_along(taus), function(i) {
  
  qr_pred <- predict(quantile_regression_models[[i]],
                     newdata = walmart)
  
  compute_mse(walmart$Weekly_Sales, qr_pred)
})

#comparison of mse
mse_comparison <- data.frame(
  Model = c("OLS", rep("Quantile Regression", length(taus))),
  Quantile = c(NA, taus),
  MSE = c(mse_ols, mse_quantiles)
)
mse_comparison

#2d plotting function
plot_quantile_2d <- function(var_name) {
  
  # Build formulas
  f_qr  <- as.formula(paste("Weekly_Sales ~", var_name))
  f_ols <- as.formula(paste("Weekly_Sales ~", var_name))
  
  # Fit models
  qr_model  <- rq(f_qr, tau = taus, data = walmart)
  ols_model <- lm(f_ols, data = walmart)
  
  # Extract coefficients (2 × length(taus))
  coef_mat <- coef(qr_model)
  
  # Plot
  ggplot(walmart, aes_string(x = var_name, y = "Weekly_Sales")) +
    geom_point(alpha = 5, color = "black") +
    geom_abline(
      intercept = coef(ols_model)[1],
      slope = coef(ols_model)[2],
      color = "red",
      linewidth = 1.2
    ) +
    lapply(seq_along(taus), function(i) {
      geom_abline(
        intercept = coef_mat[1, i],
        slope = coef_mat[2, i],
        linetype = "dashed",
        color = scales::hue_pal()(length(taus))[i],
        linewidth = 0.8
      )
    }) +
    labs(
      title = paste("Quantile Regression Lines for", var_name),
      subtitle = "Red: OLS | Dashed: Quantile Regression",
      x = var_name,
      y = "Weekly Sales"
    )
}

#2d plots
par(mfrow = c(2,2))
plot_quantile_2d("Temperature")
plot_quantile_2d("Fuel_Price")
plot_quantile_2d("CPI")
plot_quantile_2d("Unemployment")

#3d visualization function
plot_quantile_3d <- function(var1, var2) {
  
  taus_3d <- taus[c(2, 3, 4)]   # 0.25, 0.5, 0.75
  surface_colors <- c("blue", "black", "red")
  
  f_qr <- as.formula(paste("Weekly_Sales ~", var1, "+", var2))
  f_ols <- f_qr
  
  qr_3d  <- rq(f_qr, tau = taus_3d, data = walmart)
  ols_3d <- lm(f_ols, data = walmart)
  
  x_seq <- seq(min(walmart[[var1]]),
               max(walmart[[var1]]),
               length.out = 30)
  
  y_seq <- seq(min(walmart[[var2]]),
               max(walmart[[var2]]),
               length.out = 30)
  
  grid_3d <- expand.grid(x_seq, y_seq)
  colnames(grid_3d) <- c(var1, var2)
  
  # Predict quantile surfaces
  pred_surfaces <- lapply(seq_along(taus_3d), function(i) {
    matrix(
      predict(qr_3d, newdata = grid_3d, tau = taus_3d[i]),
      nrow = length(x_seq),
      ncol = length(y_seq)
    )
  })
  
  # Small vertical shift for visibility (purely visual)
  eps <- sd(walmart$Weekly_Sales) * 0.05
  pred_surfaces[[1]] <- pred_surfaces[[1]] - eps
  pred_surfaces[[3]] <- pred_surfaces[[3]] + eps
  
  # OLS surface
  ols_surface <- matrix(
    predict(ols_3d, newdata = grid_3d),
    nrow = length(x_seq),
    ncol = length(y_seq)
  )
  
  fig <- plot_ly()
  
  # observed data
  fig <- fig %>%
    add_markers(
      x = walmart[[var1]],
      y = walmart[[var2]],
      z = walmart$Weekly_Sales,
      marker = list(size = 2, opacity = 0.25),
      name = "Observed data"
    )
  
  # quantile surfaces
  for (i in seq_along(taus_3d)) {
    fig <- fig %>%
      add_surface(
        x = x_seq,
        y = y_seq,
        z = pred_surfaces[[i]],
        opacity = 0.7,
        showscale = FALSE,
        colorscale = list(
          c(0, surface_colors[i]),
          c(1, surface_colors[i])
        ),
        name = paste0("Quantile τ = ", taus_3d[i])
      )
  }
  
  # OLS plane
  fig <- fig %>%
    add_surface(
      x = x_seq,
      y = y_seq,
      z = ols_surface,
      opacity = 0.4,
      colorscale = list(c(0, "gray"), c(1, "gray")),
      showscale = FALSE,
      name = "OLS mean surface"
    )
  
  fig
}

#3d plots
plot_quantile_3d("Temperature", "Unemployment")
plot_quantile_3d("Temperature", "Fuel_Price")
plot_quantile_3d("CPI", "Unemployment")
plot_quantile_3d("CPI", "Fuel_Price")
plot_quantile_3d("Temperature", "CPI")
plot_quantile_3d("Fuel_Price", "Unemployment")


# ==========================================
#       HW 2 Question 2
# ==========================================
# Draw the regression quantile curve for various choices of quantile indexes.


# Dataset loading
df <- read.csv("Updated Quality of Life Data.csv")
ind <- sample(1:dim(df)[1], 3e3)

# Consider avg sleep hours as input variable (taking only 3000 rows for faster visulization)
x <- df$avg_sleep_hours_per_day[ind]

# Age at Death as Output column
y <- df$age_at_death[ind]
n <- length(x) 

data <- data.frame(x, y)

# gaussian kernel 
gaussian_kernel <- function(u) {
  dnorm(u)
}

# epanechnikov kernel
epan_kernel <- function(u) {
  0.75 * (1 - u^2) * (abs(u) <= 1)
}

# Uniform Kernel
uniform_kernel <- function(u) {
  0.5 * (abs(u) <= 1)
}

# regression Quantile
regression_quantile <- function(x, y, x0, h, alpha, kernel_function) 
{
  quantile_est <- numeric(length(x0))
  
  for (i in 1:length(x0)) {
    weights <- kernel_function((x - x0[i]) / h) 
    
    if (sum(weights) == 0) {
      quantile_est[i] <- NA
    } else {
      quantile_est[i] <- weighted.quantile(y, weights, alpha)
    }
  }
  return(quantile_est)
}

weighted.quantile <- function(x, w, alpha) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  
  cw <- cumsum(w) / sum(w)
  return(x[which(cw >= alpha)[1]])
}


h_values <- c(0.5, 1, 1.5, 2) # Bandwidth

# Choice of kernel
kernels <- list(
  Gaussian = gaussian_kernel,
  Epanechnikov = epan_kernel,
  Uniform = uniform_kernel
)

plot_regression_quantiles <- function(x, y, 
                                      tau_values, 
                                      h_values, 
                                      kernels, 
                                      x_grid = NULL) {
  
  # Create grid once
  if(is.null(x_grid)) {
    x_grid <- seq(min(x), max(x), length.out = 200)
  }
  
  par(mfrow = c(length(h_values), length(kernels)),
      mar = c(2, 2, 2, 1), 
      oma = c(0, 0, 2, 0)) 
  
  line_colors <- rainbow(length(tau_values), start = 0, end = 0.8) 
  
  # Loop of Bandwidth (h) 
  for(h in h_values) {
    
    # Loop of Kernels
    for(kname in names(kernels)) {
      
      plot(x, y, 
           pch = 16, 
           col = rgb(0.5, 0.5, 0.5, 0.2), 
           main = paste("Kernel:", kname, "| h =", h),
           xlab = "", ylab = "")
      
      # Loop over Quantile indices
      for(i in seq_along(tau_values)) {
        alpha <- tau_values[i]
        
        # Calculate the estimate
        q_est <- regression_quantile(
          x, y, x_grid, 
          h = h, 
          alpha = alpha, 
          kernel_function = kernels[[kname]]
        )
        
        lines(x_grid, q_est, 
              col = line_colors[i], 
              lwd = 2)
      }
      
      if(h == h_values[1] && kname == names(kernels)[1]) {
        legend("topleft", 
               legend = paste("Tau =", tau_values),
               col = line_colors, 
               lwd = 2, 
               cex = 0.7, 
               bg = "white")
      }
    }
  }
  
  par(mfrow = c(1,1))
}

tau_values <- c(0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9)

plot_regression_quantiles(
  x = x, 
  y = y, 
  tau_values = tau_values, 
  h_values = h_values, 
  kernels = kernels
)

mse_results <- data.frame()

for(alpha in tau_values) {
  
  for(h in h_values) {
    
    for(kname in names(kernels)) {
      
      q_hat <- regression_quantile(
        x, y, x,
        h = h,
        alpha = alpha,
        kernel_fun = kernels[[kname]]
      )
      
      mse_val <- mean((y - q_hat)^2, na.rm = TRUE)
      
      mse_results <- rbind(mse_results,
                           data.frame(
                             Tau = alpha,
                             Bandwidth = h,
                             Kernel = kname,
                             MSE = mse_val
                           ))
    }
  }
}

mse_results <- mse_results[order(mse_results$MSE),]

print(mse_results)


# ==========================================
#       HW 2 Question 3
# ==========================================
# Draw the data depth and multivariate quantile contour lines


# Load required libraries
if(!require(ddalpha)) install.packages("ddalpha")
if(!require(mvtnorm)) install.packages("mvtnorm")
library(ddalpha)
library(mvtnorm)

# Function to generate plots for Symmetric and Asymmetric data
plot_depth_comparison <- function(data, type_label) {
  
  # Grid Setup
  grid_size <- 50
  x_range <- seq(min(data[,1]) - 0.5, max(data[,1]) + 0.5, length.out = grid_size)
  y_range <- seq(min(data[,2]) - 0.5, max(data[,2]) + 0.5, length.out = grid_size)
  grid <- expand.grid(x = x_range, y = y_range)
  
  # Define Depth Methods
  methods <- list(
    "Mahalanobis" = depth.Mahalanobis,
    "Halfspace" = depth.halfspace,
    "Simplicial" = depth.simplicial,
    "Simplicial Volume" = depth.simplicialVolume,
    "Spatial" = depth.spatial
  )
  
  # Set up 2x3 layout (5 plots + 1 legend)
  par(mfrow = c(2, 3), mar = c(3, 3, 3, 1), oma = c(1, 1, 2, 1))
  
  for (name in names(methods)) {
    # Calculate Depths
    sample_depths <- as.numeric(methods[[name]](data, data))
    grid_depths <- as.numeric(methods[[name]](grid, data))
    z <- matrix(grid_depths, nrow = grid_size)
    
    # Thresholds for Quantile Contours
    thresholds <- quantile(sample_depths, probs = c(0.25, 0.5, 0.75))
    
    # Plotting
    plot(data, pch = 16, col = rgb(0,0,0,0.5), main = name, xlab="", ylab="", cex=0.5)
    
    # Depth Curves
    contour(x_range, y_range, z, add = TRUE, col = "black", nlevels = 8, lwd = 0.5)
    
    # Quantile Contours
    contour(x_range, y_range, z, levels = thresholds, 
            add = TRUE, col = c("#1b9e77", "#d95f02", "#7570b3"), lwd = 1.5, labels = "")
  }
  
  # Legend
  plot.new()
  legend("center", 
         legend = c("75% Quantile", "50% Quantile", "25% Quantile", "General Depth Curves"),
         col = c("#1b9e77", "#d95f02", "#7570b3", "black"), 
         lwd = c(2, 2, 2, 0.5), 
         lty = 1, 
         cex = 1.2, 
         bty = "n", 
         title = "Quantile Contours")
  
  mtext(paste("Comparison of 5 Depth Functions -", type_label), outer = TRUE, cex = 1.2, font = 2)
}


# Generate the data and run
set.seed(123)
n <- 500

# Generate Symmetric Data
mu <- c(0,0)
sigma <- matrix(c(1, 0.6, 0.6, 1), 2, 2)
data_sym <- rmvnorm(n, mean = mu, sigma = sigma)

# Run for Symmetric
plot_depth_comparison(data_sym, "Symmetric Distribution")


# Generate Asymmetric Data (Mixture of two Gaussians)
# This creates a "tail" or "skew" that is mathematically stable
group1 <- matrix(rnorm(n*0.7, mean=0, sd=1), ncol=2)
group2 <- matrix(rnorm(n*0.3, mean=3, sd=1.5), ncol=2) # The "skew" group
data_asym <- rbind(group1, group2)

# Run for Asymmetric
plot_depth_comparison(data_asym, "Asymmetric Distribution")



