# QUESTION 1

# ==========================================
# Install and Load Required Library
# ==========================================
if(!require(ddalpha)) install.packages("ddalpha")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(gridExtra)) install.packages("gridExtra")
library(ddalpha)
library(ggplot2)
library(gridExtra)
library(MASS)

# Clean Data & Calculate parameters to simulate a matching Normal distribution
data(Boston)
n <- dim(Boston)[1]
ind <- sample(1:n,100)
boston_clean <- as.matrix(Boston[ind, c("crim","rm","medv")])
mu <- colMeans(boston_clean)
sigma <- cov(boston_clean)

# Simulate multivariate normal data with the same properties - mean and covariance 
set.seed(123)
sim_normal <- mvrnorm(n = nrow(boston_clean), mu = mu, Sigma = sigma)

# To make a DD-plot, we must evaluate a combined dataset (Z) 
# against both the original data (F) and the theoretical data (G)
Z <- rbind(boston_clean, sim_normal)

# Helper function for plotting
create_dd_plot <- function(depth_data, depth_normal, title) {
  df <- data.frame(Data_Depth = depth_data, Normal_Depth = depth_normal)
  ggplot(df, aes(x = Data_Depth, y = Normal_Depth)) +
    geom_point(alpha = 0.6, color = "darkblue") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
    ggtitle(title) +
    theme_minimal() +
    labs(x = "Depth w.r.t Iris Data", y = "Depth w.r.t Normal Data")
}

plots <- list()

# Calculate Depths
# Mahalanobis
d_data_mah <- depth.Mahalanobis(Z, boston_clean)
d_norm_mah <- depth.Mahalanobis(Z, sim_normal)
plots[[1]] <- create_dd_plot(d_data_mah, d_norm_mah, "Mahalanobis Depth")

# Spatial
d_data_spa <- depth.spatial(Z, boston_clean)
d_norm_spa <- depth.spatial(Z, sim_normal)
plots[[2]] <- create_dd_plot(d_data_spa, d_norm_spa, "Spatial Depth")

# Halfspace (Tukey)
d_data_hs <- depth.halfspace(Z, boston_clean, exact = FALSE)
d_norm_hs <- depth.halfspace(Z, sim_normal, exact = FALSE)
plots[[3]] <- create_dd_plot(d_data_hs, d_norm_hs, "Halfspace (Tukey) Depth")

# Simplicial
d_data_sim <- depth.simplicial(Z, boston_clean, exact = FALSE)
d_norm_sim <- depth.simplicial(Z, sim_normal, exact = FALSE)
plots[[4]] <- create_dd_plot(d_data_sim, d_norm_sim, "Simplicial Depth")

# Simplicial Volume
d_data_sv <- depth.simplicialVolume(Z, boston_clean, exact = FALSE)
d_norm_sv <- depth.simplicialVolume(Z, sim_normal, exact = FALSE)
plots[[5]] <- create_dd_plot(d_data_sv, d_norm_sv, "Simplicial Volume Depth")

# Display all plots
grid.arrange(grobs = plots, ncol = 3)

# If the plot points align closely to a straight diagonal line, 
# it suggests the data follows a multivariate normal distribution.

# ==========================================
#KS Test and CVM test
# ==========================================

# ==========================================
# Inbuilt Function: KS Statistic
# ==========================================
get_ks_stat <- function(depths_sample, depths_ref) {
  # Empirical CDF of the reference (theoretical) and the sample
  F_ref <- ecdf(depths_ref)
  F_samp <- ecdf(depths_sample)
  
  # Evaluate over all unique points in both sets
  grid <- sort(unique(c(depths_sample, depths_ref)))
  
  # KS Stat is the maximum absolute difference
  ks_stat <- max(abs(F_samp(grid) - F_ref(grid)))
  return(ks_stat)
}

# ==========================================
# Inbuilt Function: Cramér-von Mises Statistic
# ==========================================
get_cvm_stat <- function(depths_sample, depths_ref) {
  n <- length(depths_sample)
  F_ref <- ecdf(depths_ref)
  
  # Sort the sample depths
  sorted_samp <- sort(depths_sample)
  
  # Evaluate theoretical CDF at sorted sample points
  U <- F_ref(sorted_samp)
  
  # Computational formula for CvM statistic
  cvm_stat <- sum((U - (2 * (1:n) - 1) / (2 * n))^2) + 1 / (12 * n)
  return(cvm_stat)
}

# ==========================================
# Main Function: Depth-Based MVN Test
# ==========================================
depth_mvn_test <- function(data, method, B = 100, alpha = 0.05) {
  # data: The multivariate dataset to test
  # method: Depth method (e.g., "Mahalanobis", "spatial")
  # B: Number of bootstrap iterations (recommend 200-500 for stable results)
  # alpha: Significance level
  # ...: Extra arguments for depth function (like exact=FALSE, k=50000)
  
  data <- as.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  
  # Step 1: Estimate parameters from the original data
  mu_hat <- colMeans(data)
  Sigma_hat <- cov(data)
  
  # Step 2: Generate a large reference MVN to represent the "Theoretical" distribution
  set.seed(123)
  n_ref <- max(2000, 10 * n) # Ensure reference is large enough
  Z_ref <- mvrnorm(n_ref, mu_hat, Sigma_hat)
  
  # Step 3: Calculate depths for the original sample
  d_sample <- depth.(data, Z_ref, notion = method, exact = FALSE)
  d_ref    <- depth.(Z_ref, Z_ref, notion = method,exact = FALSE)
  
  # Step 4: Calculate Observed Test Statistics
  ks_obs  <- get_ks_stat(d_sample, d_ref)
  cvm_obs <- get_cvm_stat(d_sample, d_ref)
  
  # Step 5: Parametric Bootstrap for Critical Values
  ks_boots  <- numeric(B)
  cvm_boots <- numeric(B)
  
  cat(sprintf("Running %d bootstrap iterations... This may take a moment depending on the depth method.\n", B))
  
  for (b in 1:B) {
    if (b %% 25 == 0) cat("Iteration", b, "\n")
    
    # 5a. Simulate a new dataset from the estimated MVN parameters
    X_boot <- mvrnorm(n, mu_hat, Sigma_hat)
    
    # 5b. Re-estimate parameters for the bootstrap sample (Lilliefors effect)
    mu_boot <- colMeans(X_boot)
    Sigma_boot <- cov(X_boot)
    
    # 5c. Generate a new reference for this bootstrap iteration
    Z_ref_boot <- mvrnorm(n_ref, mu_boot, Sigma_boot)
    
    # 5d. Calculate depths
    d_samp_boot <- depth.(X_boot, Z_ref_boot, notion = method, exact = FALSE)
    d_ref_boot  <- depth.(Z_ref_boot, Z_ref_boot, notion = method, exact = FALSE)
    
    # 5e. Calculate bootstrap statistics
    ks_boots[b]  <- get_ks_stat(d_samp_boot, d_ref_boot)
    cvm_boots[b] <- get_cvm_stat(d_samp_boot, d_ref_boot)
  }
  
  # Step 6: Calculate Critical Values and P-values
  ks_crit  <- quantile(ks_boots, 1 - alpha)
  cvm_crit <- quantile(cvm_boots, 1 - alpha)
  
  ks_p  <- mean(ks_boots >= ks_obs)
  cvm_p <- mean(cvm_boots >= cvm_obs)
  
  # Return cleanly formatted results
  list(
    Observed_KS = ks_obs,
    Critical_Value_KS = ks_crit,
    P_Value_KS = ks_p,
    Observed_CvM = cvm_obs,
    Critical_Value_CvM = cvm_crit,
    P_Value_CvM = cvm_p
  )
}

# ==========================================
# Execution
# ==========================================
# Test 1: Mahalanobis

res_mah <- depth_mvn_test(boston_clean, method = "Mahalanobis", B = 100)
print(res_mah)
# Test 2: Spatial

res_spa <- depth_mvn_test(boston_clean, method = "spatial", B = 100)
print(res_spa)
# Test 3: Halfspace

res_hal <- depth_mvn_test(boston_clean, method = "halfspace", B = 100)
print(res_hal)

# ==========================================
# ==========================================


# QUESTION 2

# ==========================================
# Data Preparation
# ==========================================
# ==========================================
# Setup and Library Loading
# ==========================================

library(ddalpha)
library(ggplot2)

data(iris) # Load dataset

# ==========================================
# Exploratory Data Analysis: Symmetry
# ==========================================
cat("--- Checking Elliptical Symmetry ---\n")

elliptically_symmetric <- function(data_subset, species_name) {
  n <- nrow(data_subset)
  p <- ncol(data_subset)
  
  # Calculate depth to find the Multivariate Median
  depth_vals <- depth.spatial(as.matrix(data_subset), as.matrix(data_subset))
  center_point <- colMeans(data_subset) 
  
  # Calculation of Mahalanobis Distances (MD^2)
  S <- cov(data_subset)
  md2 <- mahalanobis(data_subset, center = center_point, cov = S)
  
  # Prepare Q-Q Plot Data
  expected_q <- qchisq(ppoints(n), df = p)
  observed_q <- sort(md2)
  df_qq <- data.frame(Theoretical = expected_q, Sample = observed_q)
  
  # Generate Plot
  p_plot <- ggplot(df_qq, aes(x = Theoretical, y = Sample)) +
    geom_point(color = "blue") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(title = paste("Elliptical Symmetry Diagnostic:", species_name),
         subtitle = "Linear alignment indicates Elliptical Symmetry",
         x = "Theoretical Chi-Square Quantiles",
         y = "Sample Mahalanobis Distances") +
    theme_minimal()
  
  return(p_plot)
}

# Create a list of data frames for each species
species_list <- list(
  setosa = iris[iris$Species == "setosa", 1:4],
  versicolor = iris[iris$Species == "versicolor", 1:4],
  virginica = iris[iris$Species == "virginica", 1:4]
)

# Print symmetry plots for each species
for(name in names(species_list)) {
  print(elliptically_symmetric(species_list[[name]], name))
}

# ==========================================
# Exploratory Data Analysis: Outliers
# ==========================================
cat("\n--- Detecting Multivariate Outliers ---\n")

X_all <- iris[, 1:4] 
mu <- colMeans(X_all)
S  <- cov(X_all)

# Calculate Mahalanobis Distance (MD^2)
md_values <- mahalanobis(X_all, center = mu, cov = S)

# Threshold setup (Chi-square distribution with p degrees of freedom)
p <- ncol(X_all)
cutoff <- qchisq(0.975, df = p)

outlier_indices <- which(md_values > cutoff)
cat("Identified Outlier Indices:", outlier_indices, "\n")

# Visualization: Distance Plot
plot(md_values, 
     main = "Multivariate Outlier Detection (Mahalanobis)",
     xlab = "Observation Index", 
     ylab = "Mahalanobis Distance (MD^2)",
     pch = 19, col = ifelse(md_values > cutoff, "red", "black"))

# Add the Chi-square threshold line
abline(h = cutoff, col = "red", lwd = 2, lty = 2)
text(x = 20, y = cutoff + 2, labels = "97.5% Threshold", col = "red")


# ==========================================
# Data Splitting for Classification
# ==========================================
# Split: First 25 of each class for training, remaining 25 for testing
train_idx <- c(1:25, 51:75, 101:125)
test_idx  <- c(26:50, 76:100, 126:150)

X_train <- as.matrix(iris[train_idx, 1:4])
y_train <- iris[train_idx, 5]

X_test  <- as.matrix(iris[test_idx, 1:4])
y_test  <- iris[test_idx, 5]

classes <- unique(y_train)

# ==========================================
# Depth-Based Classification (5 Methods)
# ==========================================
methods <- c("halfspace", "Mahalanobis", "simplicial", "spatial", "simplicialVolume")
method_names <- c("Half-space (Tukey)", "Mahalanobis", "Simplicial", "Spatial (L1)", "Oja (Simplicial Volume)")

misclass_rates <- numeric(length(methods))
n_test <- nrow(X_test)
n_classes <- length(classes)

cat("\nCalculating depths and classifying... (Simplicial and Oja may take a few seconds)\n")

for (m in seq_along(methods)) {
  method <- methods[m]
  depth_matrix <- matrix(0, nrow = n_test, ncol = n_classes)
  colnames(depth_matrix) <- as.character(classes)
  
  # Calculate depth of each test point w.r.t each training class
  for (i in 1:n_classes) {
    c <- classes[i]
    X_train_c <- X_train[y_train == c, ]
    
    # Unified wrapper for all depth functions
    depth_matrix[, i] <- depth.(X_test, X_train_c, notion = method)
  }
  
  # Assign to the class where the test point is "deepest"
  pred_indices <- max.col(depth_matrix, ties.method = "first")
  y_pred <- classes[pred_indices]
  
  # Record the misclassification rate
  misclass_rates[m] <- mean(y_pred != y_test)
}

# Print results table
results_df <- data.frame(
  Depth_Function = method_names,
  Misclassification_Rate = round(misclass_rates, 4),
  Accuracy = round(1 - misclass_rates, 4)
)

print("\n--- Classification Results ---")
print(results_df)


# ==========================================
# Visualizations
# ==========================================
cat("\nGenerating analytical plots...\n")

# Subset to Petal.Length and Petal.Width for 2D visual mapping
X_train_2D <- as.matrix(iris[train_idx, 3:4])
X_test_2D  <- as.matrix(iris[test_idx, 3:4])
X_train_vers_2D <- X_train_2D[y_train == "versicolor", ]
X_train_virg_2D <- X_train_2D[y_train == "virginica", ]

# --- Plot A: Compare Depth Contours ---
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))

for (m in seq_along(methods)) {
  depth.contours(X_train_vers_2D, depth = methods[m], 
                 main = method_names[m],
                 xlab = "Petal Length", ylab = "Petal Width",
                 col = "blue",
                 exact = TRUE) # Prevents "static" noise on Oja/Simplicial
  points(X_train_vers_2D, pch = 16, col = "blue")
}

# --- Plot B: DD-Plot and Boxplots ---
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

# 1. DD-Plot (Spatial Depth)
depth_wrt_vers <- depth.spatial(X_test_2D, X_train_vers_2D)
depth_wrt_virg <- depth.spatial(X_test_2D, X_train_virg_2D)

plot(depth_wrt_vers, depth_wrt_virg, 
     col = as.numeric(as.factor(y_test)), 
     pch = 16,
     xlab = "Depth w.r.t Versicolor",
     ylab = "Depth w.r.t Virginica",
     main = "DD-Plot (Spatial Depth)")
abline(a = 0, b = 1, lty = 2, col = "red", lwd = 2)
legend("topright", legend = levels(as.factor(y_test)), 
       col = 1:3, pch = 16, title = "True Class", cex = 0.8)

# 2. Depth Boxplots (Testing Spatial depth behavior for true Virginica points)
X_test_virg <- as.matrix(iris[126:150, 1:4]) 

depths_in_setosa <- depth.spatial(X_test_virg, as.matrix(iris[1:25, 1:4]))
depths_in_vers   <- depth.spatial(X_test_virg, as.matrix(iris[51:75, 1:4]))
depths_in_virg   <- depth.spatial(X_test_virg, as.matrix(iris[101:125, 1:4]))

boxplot(depths_in_setosa, depths_in_vers, depths_in_virg,
        names = c("Setosa", "Versicolor", "Virginica"),
        col = c("#FF9999", "#99FF99", "#9999FF"),
        main = "Depth of True Virginica Points",
        ylab = "Spatial Depth Score")

# Reset plotting layout to default
par(mfrow = c(1, 1))

# ==========================================
# ==========================================

