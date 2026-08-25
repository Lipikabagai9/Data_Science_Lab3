# ======================================
# ======================================
# install.packages(c("fda","ddalpha", "quantmod"))

# ======================================
# ======================================
# QUESTION 1
# ======================================
# ======================================
library(MASS) 

# Load the built-in dataset
data(iris)

# Extract features and scale them
df <- scale(iris[, 1:4])

# Extract the true labels for the confusion matrix
true_labels <- iris$Species

# Set k to 3 classes
k <- 3

# K-means Clustering
set.seed(123)
km_res <- kmeans(df, centers = k, nstart = 25)

# Data Depth Clustering
depth_clustering <- function(data, k, max_iter = 100) {
  n <- nrow(data)
  p <- ncol(data)
  set.seed(123)
  centers <- data[sample(1:n, k), ]
  clusters <- rep(0, n)
  cov_matrices <- list()
  for(j in 1:k) cov_matrices[[j]] <- cov(data)
  
  for (iter in 1:max_iter) {
    old_clusters <- clusters
    depths <- matrix(0, nrow = n, ncol = k)
    
    for (j in 1:k) {
      md <- mahalanobis(data, center = centers[j,], cov = cov_matrices[[j]])
      depths[, j] <- 1 / (1 + md)
    }
    
    clusters <- apply(depths, 1, which.max)
    if (all(clusters == old_clusters)) break
    
    for (j in 1:k) {
      cluster_data <- data[clusters == j, , drop = FALSE]
      if (nrow(cluster_data) > p) {
        centers[j, ] <- colMeans(cluster_data)
        cov_matrices[[j]] <- cov(cluster_data) + diag(1e-5, p)
      }
    }
  }
  return(list(cluster = clusters))
}

depth_res <- depth_clustering(df, k = k)


# Optimal Cluster Alignment (Solving Label Switching)
# This function maps predicted clusters to the true labels 
# by finding the permutation that maximizes accuracy.

align_clusters <- function(true_labels, pred_clusters) {
  true_levels <- unique(true_labels)
  k <- length(true_levels)
  true_numeric <- match(true_labels, true_levels)
  
  # A simple recursive function to generate all permutations in base R
  get_perms <- function(v) {
    if (length(v) == 1) return(matrix(v))
    X <- NULL
    for (i in 1:length(v)) {
      X <- rbind(X, cbind(v[i], get_perms(v[-i])))
    }
    return(X)
  }
  
  all_perms <- get_perms(1:k)
  best_acc <- -1
  best_mapped <- NULL
  
  # Evaluate every permutation to find the highest overlap
  for (i in 1:nrow(all_perms)) {
    mapped_numeric <- all_perms[i, ][pred_clusters]
    acc <- sum(mapped_numeric == true_numeric)
    
    if (acc > best_acc) {
      best_acc <- acc
      best_mapped <- mapped_numeric
    }
  }
  
  # Return the predicted clusters properly renamed as the true label strings
  return(true_levels[best_mapped])
}

# Apply the optimal alignment to both sets of results
aligned_km <- align_clusters(true_labels, km_res$cluster)
aligned_depth <- align_clusters(true_labels, depth_res$cluster)

# Clean 3-Panel Visual Comparison via PCA
pca <- princomp(df)
pca_data <- pca$scores[, 1:2]

# Function to calculate and draw 95% confidence ellipses via eigen decomposition
draw_ellipse <- function(data, cluster_id, col) {
  c_data <- data[cluster_id, , drop = FALSE]
  if(nrow(c_data) > 2) {
    center <- colMeans(c_data)
    cov_mat <- cov(c_data)
    eig <- eigen(cov_mat)
    
    # 95% confidence interval for 2D is roughly chisq with 2df = 5.991
    theta <- seq(0, 2*pi, length.out = 100)
    circle <- cbind(cos(theta), sin(theta))
    ellipse <- t(center + sqrt(5.991) * t(circle %*% diag(sqrt(eig$values)) %*% t(eig$vectors)))
    lines(ellipse, col = col, lwd = 2, lty = 2)
  }
}

# Set up a 1x3 plotting grid
par(mfrow = c(1, 3), mar = c(5, 4, 4, 1))

# Plot 1: The Ground Truth
plot(pca_data, col = as.numeric(true_labels), pch = 19,
     main = "Ground Truth (True Labels)", xlab = "PC 1", ylab = "PC 2", 
     bty = "n", cex = 1.2)
for(i in 1:k) draw_ellipse(pca_data, as.numeric(true_labels) == i, i)


# Plot 2: k-means Assignment
plot(pca_data, col = aligned_km, pch = 19,
     main = "k-means Clustering", xlab = "PC 1", ylab = "", 
     bty = "n", cex = 1.2)
for(i in 1:k) draw_ellipse(pca_data, as.numeric(aligned_km) == i, i)


# Plot 3: Mahalanobis Depth Assignment
plot(pca_data, col = aligned_depth, pch = 19,
     main = "Depth-Based Clustering", xlab = "PC 1", ylab = "", 
     bty = "n", cex = 1.2)
for(i in 1:k) draw_ellipse(pca_data, as.numeric(aligned_depth) == i, i)


# True Misclassification Rates & Final Tables
cat("\n--- True Misclassification Rates ---\n")

# Calculate and print for k-means
km_misclass <- sum(aligned_km != true_labels) / length(true_labels)
cat(sprintf("k-means Misclassification Rate: %.2f%%\n", km_misclass * 100))

cat("\nAligned k-means Confusion Matrix:\n")
print(table(True_Class = true_labels, Predicted = aligned_km))

# Calculate and print for Mahalanobis Depth
depth_misclass <- sum(aligned_depth != true_labels) / length(true_labels)
cat(sprintf("\nDepth-Based Misclassification Rate: %.2f%%\n", depth_misclass * 100))

cat("\nAligned Depth Confusion Matrix:\n")
print(table(True_Class = true_labels, Predicted = aligned_depth))


#  Create a data frame matching states to their assigned clusters
comparison_table <- data.frame(
  K_Means_Cluster = aligned_km,
  Depth_Cluster = aligned_depth,
  True_Cluster = true_labels
)

# Add columns to highlight differences in assignment
comparison_table$Kmeans_Match <- ifelse(
  comparison_table$K_Means_Cluster == comparison_table$True_Cluster, 
  "Yes", 
  "No"
)
comparison_table$Depth_Match <- ifelse(
  comparison_table$Depth_Cluster == comparison_table$True_Cluster, 
  "Yes", 
  "No"
)

# View the full table in R's viewer
View(comparison_table)
print(head(comparison_table, 15))









# ======================================
# ======================================
# QUESTION 2
# ======================================
# ======================================








library(fda)
library(ddalpha)

# ==============================================================================
# 0. FROM-SCRATCH STATISTICAL FUNCTIONS
# ==============================================================================

# Kolmogorov-Smirnov (KS) Permutation Test
ks_test_scratch <- function(x, y, n_perm = 1e4) {
  n <- length(x); m <- length(y)
  pooled <- c(x, y)
  
  calc_ks <- function(g1, g2) {
    grid <- sort(c(g1, g2))
    e1 <- sapply(grid, function(t) sum(g1 <= t) / length(g1))
    e2 <- sapply(grid, function(t) sum(g2 <= t) / length(g2))
    return(max(abs(e1 - e2)))
  }
  
  obs_stat <- calc_ks(x, y)
  perm_stats <- replicate(n_perm, {
    shuf <- sample(pooled)
    calc_ks(shuf[1:n], shuf[(n+1):(n+m)])
  })
  
  p_val <- mean(perm_stats >= obs_stat)
  return(list(statistic = obs_stat, p.value = p_val))
}

# Cramér-von Mises (CvM) Permutation Test
cvm_test_scratch <- function(x, y, n_perm = 1e4) {
  n <- length(x); m <- length(y)
  N <- n + m
  pooled <- c(x, y)
  
  calc_cvm <- function(g1, g2) {
    grid <- sort(c(g1, g2))
    e1 <- sapply(grid, function(t) sum(g1 <= t) / length(g1))
    e2 <- sapply(grid, function(t) sum(g2 <= t) / length(g2))
    return((length(g1) * length(g2) / (N^2)) * sum((e1 - e2)^2))
  }
  
  obs_stat <- calc_cvm(x, y)
  perm_stats <- replicate(n_perm, {
    shuf <- sample(pooled)
    calc_cvm(shuf[1:n], shuf[(n+1):(n+m)])
  })
  
  p_val <- mean(perm_stats >= obs_stat)
  return(list(statistic = obs_stat, p.value = p_val))
}

# PCA from scratch
pca_scratch <- function(data_matrix, n_comp = 3) {
  if(nrow(data_matrix) > ncol(data_matrix)) data_matrix <- t(data_matrix)
  X_c <- sweep(data_matrix, 2, colMeans(data_matrix))
  cov_mat <- cov(X_c)
  eigen_decomp <- eigen(cov_mat)
  scores <- X_c %*% eigen_decomp$vectors[, 1:n_comp]
  return(scores)
}

# Load Canadian Weather Data
data("CanadianWeather", package = "fda")
temp_data <- CanadianWeather$dailyAv[,, "Temperature.C"]
precip_data <- CanadianWeather$dailyAv[,, "Precipitation.mm"]


# ==============================================================================
# 1. FIRST SCENARIO: SAME DATA (NULL HYPOTHESIS TRUE)
# ==============================================================================
cat("========================================================\n")
cat("SCENARIO 1: SAME DATA (Randomly split Temperature data)\n")
cat("========================================================\n")

set.seed(123)
n_stations <- ncol(temp_data)
idx_g1 <- sample(1:n_stations, size = 17)
idx_g2 <- setdiff(1:n_stations, idx_g1)

# Apply PCA to the full temperature dataset
scores_same <- pca_scratch(t(temp_data), n_comp = 3)

# Split the scores into the two random groups
scores_same_g1 <- scores_same[idx_g1, ]
scores_same_g2 <- scores_same[idx_g2, ]

cat("\n--- (a) PCA on raw data -> CvM & KS on PC1 Scores ---\n")
ks_same_pca <- ks_test_scratch(scores_same_g1[,1], scores_same_g2[,1])
cvm_same_pca <- cvm_test_scratch(scores_same_g1[,1], scores_same_g2[,1])

cat(sprintf("KS Test  | Stat: %.4f, P-value: %.3f\n", ks_same_pca$statistic, ks_same_pca$p.value))
cat(sprintf("CvM Test | Stat: %.4f, P-value: %.3f\n", cvm_same_pca$statistic, cvm_same_pca$p.value))


cat("\n--- (b) Depth of PCA -> CvM & KS on Depth values ---\n")
# Calculate spatial depth of all pooled scores
depths_same <- depth.spatial(x = scores_same, data = scores_same)

# Split the 1D depth values into the two random groups
depths_same_g1 <- depths_same[idx_g1]
depths_same_g2 <- depths_same[idx_g2]

ks_same_depth <- ks_test_scratch(depths_same_g1, depths_same_g2)
cvm_same_depth <- cvm_test_scratch(depths_same_g1, depths_same_g2)

cat(sprintf("KS Test  | Stat: %.4f, P-value: %.3f\n", ks_same_depth$statistic, ks_same_depth$p.value))
cat(sprintf("CvM Test | Stat: %.4f, P-value: %.3f\n", cvm_same_depth$statistic, cvm_same_depth$p.value))

cat("-> High p-values. We fail to reject the null.They come from the same distribution\n\n")


# ==============================================================================
# 2. SECOND SCENARIO: DIFFERENT DATA (NULL HYPOTHESIS FALSE)
# ==============================================================================
cat("\n========================================================\n")
cat("SCENARIO 2: DIFFERENT DATA (Temperature vs Precipitation)\n")
cat("========================================================\n")

# Scale both datasets to isolate "shape" differences rather than pure scale/magnitude
temp_scaled <- scale(as.vector(temp_data))
precip_scaled <- scale(as.vector(precip_data))
dim(temp_scaled) <- dim(temp_data)
dim(precip_scaled) <- dim(precip_data)

# Combine into one dataset (Rows 1:35 = Temp, Rows 36:70 = Precip)
mixed_data <- cbind(temp_scaled, precip_scaled) 

# Apply PCA to the mixed dataset
scores_diff <- pca_scratch(t(mixed_data), n_comp = 3)

# Split the scores into the two logical groups
scores_diff_temp <- scores_diff[1:35, ]
scores_diff_precip <- scores_diff[36:70, ]

cat("\n--- (a) PCA on raw data -> CvM & KS on PC1 Scores ---\n")
ks_diff_pca <- ks_test_scratch(scores_diff_temp[,1], scores_diff_precip[,1])
cvm_diff_pca <- cvm_test_scratch(scores_diff_temp[,1], scores_diff_precip[,1])

cat(sprintf("KS Test  | Stat: %.4f, P-value: %.3f\n", ks_diff_pca$statistic, ks_diff_pca$p.value))
cat(sprintf("CvM Test | Stat: %.4f, P-value: %.3f\n", cvm_diff_pca$statistic, cvm_diff_pca$p.value))

cat("\n--- (b) Depth of PCA -> CvM & KS on Depth values ---\n")
# Calculate spatial depth of all pooled scores
depths_diff <- depth.spatial(x = scores_diff, data = scores_diff)

# Split the 1D depth values into the two logical groups
depths_diff_temp <- depths_diff[1:35]
depths_diff_precip <- depths_diff[36:70]

ks_diff_depth <- ks_test_scratch(depths_diff_temp, depths_diff_precip)
cvm_diff_depth <- cvm_test_scratch(depths_diff_temp, depths_diff_precip)

cat(sprintf("KS Test  | Stat: %.4f, P-value: %.3f\n", ks_diff_depth$statistic, ks_diff_depth$p.value))
cat(sprintf("CvM Test | Stat: %.4f, P-value: %.3f\n", cvm_diff_depth$statistic, cvm_diff_depth$p.value))

cat("-> P-values < 0.05 for the CVM test and close to 0.05 for KS Test. We reject the null. They come from different distributions\n")









# ======================================
# ======================================
# QUESTION 3
# ======================================
# ======================================









# ---- 1. Load Libraries ----
library(quantmod)
library(ddalpha)

# ---- 2. Load Stock Data ----
getSymbols(c("RELIANCE.NS", "TCS.NS"), from = "2024-01-01")

rel <- Cl(RELIANCE.NS)
tcs <- Cl(TCS.NS)

# Log returns (Brownian-like)
X_ret <- na.omit(diff(log(rel)))
Y_ret <- na.omit(diff(log(tcs)))

# ---- 3. Create Functional Data ----
create_curves <- function(series, window_size = 30) {
  n <- length(series)
  curves <- list()
  
  for(i in 1:(n - window_size)){
    curves[[i]] <- as.numeric(series[i:(i + window_size - 1)])
  }
  
  do.call(rbind, curves)
}

X_curves <- create_curves(X_ret)
Y_curves <- create_curves(Y_ret)

# ---- 4. Visualization ----
par(mfrow=c(1,2))
matplot(t(X_curves[1:10, ]), type="l", lty=1,
        main="Reliance Functional Curves",
        xlab="Time", ylab="Returns")

matplot(t(Y_curves[1:10, ]), type="l", lty=1,
        main="TCS Functional Curves",
        xlab="Time", ylab="Returns")

# ---- 5. Joint Data ----
joint_curves <- cbind(X_curves, Y_curves)

# ---- 6. Depth Computation ----

# Half-space depth (ddalpha)
Dx_hs <- depth.halfspace(x = X_curves, data = X_curves)
Dy_hs <- depth.halfspace(x = Y_curves, data = Y_curves)
D_joint_hs <- depth.halfspace(x = joint_curves, data = joint_curves)

# Spatial depth (DepthProc)
Dx_sp <- as.numeric(depth.spatial(x = X_curves, data = X_curves))
Dy_sp <- as.numeric(depth.spatial(x = Y_curves, data = Y_curves))
D_joint_sp <- as.numeric(depth.spatial(x = joint_curves, data = joint_curves))

# ==============================================================================
# EMPIRICAL CDF INDEPENDENCE TEST (KS and CvM) FROM SCRATCH
# ==============================================================================

independence_test <- function(x, y, n_perm = 1e3) {
  
  # Remove NAs
  valid <- complete.cases(x, y)
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  n <- length(x)
  
  # --- Helper Function for Fast Computation ---
  calc_stats_fast <- function(x_vec, y_vec) {
    
    # 1. Create logical grids
    # X_mat[k, i] is TRUE if x_vec[k] <= x_vec[i]
    X_mat <- outer(x_vec, x_vec, "<=")
    Y_mat <- outer(y_vec, y_vec, "<=")
    
    # 2. Marginal CDFs (averaging down the columns)
    Fx <- colMeans(X_mat)
    Fy <- colMeans(Y_mat)
    
    # 3. Joint CDF at the specific data points (for CvM)
    # Checks where BOTH X_k <= X_i AND Y_k <= Y_i
    Fxy_pts <- colMeans(X_mat & Y_mat)
    
    # CvM Statistic: Mean of squared differences at observed points
    cvm_stat <- mean((Fxy_pts - Fx * Fy)^2)
    
    # 4. Joint CDF across the ENTIRE grid (for KS)
    # The matrix dot product instantly counts how many times 
    # (X_k <= X_i) AND (Y_k <= Y_j) occur across all k points.
    Fxy_grid <- crossprod(X_mat, Y_mat) / n
    
    # Independent Grid (Fx * Fy)
    Fx_Fy_grid <- outer(Fx, Fy, "*")
    
    # KS Statistic: Maximum absolute difference anywhere on the grid
    ks_stat <- max(abs(Fxy_grid - Fx_Fy_grid))
    
    return(c(cvm = cvm_stat, ks = ks_stat))
  }
  
  # --- Calculate Observed Statistics ---
  obs_stats <- calc_stats_fast(x, y)
  
  # --- Permutation Test Loop ---
  perm_cvm <- numeric(n_perm)
  perm_ks <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    # Shuffle Y to simulate the null hypothesis of independence
    y_shuffled <- sample(y)
    
    perm_res <- calc_stats_fast(x, y_shuffled)
    perm_cvm[i] <- perm_res["cvm"]
    perm_ks[i]  <- perm_res["ks"]
  }
  
  # Calculate p-values (+1 ensures we include the observed state)
  p_val_cvm <- sum(perm_cvm >= obs_stats["cvm"]) / n_perm
  p_val_ks  <- sum(perm_ks >= obs_stats["ks"]) / n_perm
  
  return(list(
    CvM = list(statistic = obs_stats["cvm"], p_value = p_val_cvm),
    KS  = list(statistic = obs_stats["ks"],  p_value = p_val_ks)
  ))
}

# ---- 8. Apply Tests ----
res_hs   <- independence_test(Dx_hs, Dy_hs)
res_sp   <- independence_test(Dx_sp, Dy_sp)

# ---- 9. Results Table ----
results <- data.frame(
  Method = c("Half-space", "Spatial"),
  Statistic = c(res_hs$KS$statistic, res_sp$KS$statistic),
  P_value = c(res_hs$KS$p_value, res_sp$KS$p_value)
)

print(results)

# ---- 10. Visualization: Depth Relationship ----
par(mfrow=c(1,3))

plot(Dx_hs * Dy_hs, D_joint_hs, main="Half-space",
     xlab="D(X)D(Y)", ylab="D(X,Y)")
abline(lm(D_joint_hs ~ I(Dx_hs * Dy_hs)))

plot(Dx_sp * Dy_sp, D_joint_sp, main="Spatial",
     xlab="D(X)D(Y)", ylab="D(X,Y)")
abline(lm(D_joint_sp ~ I(Dx_sp * Dy_sp)))

# ================================
# FINAL INTERPRETATION
# ================================
# If p-value < 0.05:
# Reject H0: D(X,Y) = D(X)*D(Y)
# ⇒ Dependence exists
# ================================

