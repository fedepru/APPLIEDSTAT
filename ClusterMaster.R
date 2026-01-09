# ==============================================================================
# MASTER CLUSTERING SCRIPT (Exam Optimized)
# ==============================================================================

# ----------------------------- 1. LIBRARIES -----------------------------------
library(cluster)    # clustering algorithms
library(dbscan)     # DBSCAN, kNNdistplot
library(MASS)       # isoMDS (optional)
library(car)        # plotting helpers
library(mvtnorm)    # multivariate normal tools
library(rgl)        # 3D plots (if needed)

# ----------------------------- 2. DATA LOADING --------------------------------
FILENAME <- "dataset.txt"   # e.g., "folk_music.txt", "sharks.txt"
IS_DIST_MATRIX <- TRUE      # TRUE for distance matrix (Folk/Metro), 
                            #FALSE for raw data (Sharks/Kapok)

#two main settings: distance matrix vs raw data

if (IS_DIST_MATRIX) {
  # Load square distance matrix
  # row.names=1 ensures row labels match column labels
  dm <- read.table(FILENAME, header=TRUE, row.names=1) 
  D  <- as.dist(dm)
  cat("Loaded DISTANCE MATRIX. Dimensions:", dim(dm), "\n")
} else {
  # Load raw multivariate data
  X_raw <- read.table(FILENAME, header=TRUE)
  D     <- dist(X_raw, method = "euclidean")
  cat("Loaded RAW DATA. Dimensions:", dim(X_raw), "\n")
}

# ----------------------------- 3. THEORY RECAP --------------------------------
# (A) HIERARCHICAL CLUSTERING: 
#     - Agglomerative: starts with n clusters, merges closest pair iteratively.
#     - Result: Dendrogram. Requires cutting to get specific k clusters.
#     - Linkages determine how distance between *clusters* is defined:
#       1. Single Linkage (Nearest Neighbor): Min dist between points. 
#          PRO: Can find non-spherical shapes. CON: "Chaining" effect (long thin clusters).
#       2. Complete Linkage (Farthest Neighbor): Max dist between points.
#          PRO: Avoids chaining. CON: Forces compact, spherical clusters; sensitive to outliers.
#       3. Average Linkage: Avg dist between all pairs.
#          PRO: Compromise between Single/Complete. Robust.
#       4. Ward Linkage: Minimizes total within-cluster variance (error sum of squares).
#          PRO: Compact, spherical, often equal-sized clusters. BEST for separating Gaussian-like blobs.
#          NOTE: Ward requires Squared Euclidean distance (usually handled internally by hclust^2).
#          also: might not be adequate in case of distance matrix: if the distance matrix is not Euclidean,
#           Ward linkage may produce misleading results.
# (B) K-MEANS:
#     - Partitional: Moves k centroids to minimize Within-Cluster Sum of Squares (WCSS).
#     - PRO: Fast, efficient.
#     - CON: Requires specifying k; Assumes spherical clusters; Sensitive to outliers;
#            REQUIRES COORDINATES (cannot run directly on a distance matrix without MDS).

# (C) DBSCAN (Density-Based Spatial Clustering of Applications with Noise):
#     - Groups points in high-density regions.
#     - Params: eps (radius), minPts (min points to form core).
#     - PRO: Finds arbitrary shapes; Detects outliers (noise = cluster 0); No need to specify k.
#     - CON: Sensitive to parameter tuning; Struggles with varying densities.

# (D) MULTIDIMENSIONAL SCALING (MDS):
#     - Visualization technique. Projects high-dim distances into 2D (or 3D) while preserving pairwise distances.
#     - Metric MDS (cmdscale): Tries to preserve actual distance values (Principal Coordinate Analysis).
#     - "Exaggerate vs Compress": Check if 2D dists < or > original dists.


#  SOLVING THE EXAM QUESTIONS (Step-by-Step)



# Q1: "Which methods are appropriate?" (Scenario: Distance Matrix) --------


if (IS_DIST_MATRIX) {
  cat("\n[Q1 Analysis]\n")
  cat("Since input is a distance matrix:\n")
  cat("- Hierarchical Clustering: YES (works directly on distances).\n")
  cat("- DBSCAN: YES (works directly on distances).\n")
  cat("- K-Means: NO (requires raw coordinates to compute means. Only possible if we run MDS first).\n")
}




# Q2: Hierarchical Clustering & Linkage Selection (Cophenetic Coeff) --------

# Calculate all main linkages
hc_single   <- hclust(D, method = "single")
hc_complete <- hclust(D, method = "complete")
hc_average  <- hclust(D, method = "average")
hc_ward     <- hclust(D, method = "ward.D2")

# Compute Cophenetic Correlation Coefficient (CPCC)
# Measure of how faithfully the dendrogram preserves the original pairwise distances.
# Higher CPCC (~0.7-0.9) = Better fit.
coph_single   <- cor(D, cophenetic(hc_single))
coph_complete <- cor(D, cophenetic(hc_complete))
coph_average  <- cor(D, cophenetic(hc_average))
coph_ward     <- cor(D, cophenetic(hc_ward))

cat("\n[Q2 - Linkage Comparison (Cophenetic Correlation)]\n")
res_link <- data.frame(Linkage = c("Single", "Complete", "Average", "Ward"),
                       CPCC = c(coph_single, coph_complete, coph_average, coph_ward))
print(res_link)
best_link <- res_link$Linkage[which.max(res_link$CPCC)]
cat("-> Best linkage based on CPCC:", best_link, "\n")
cat("-> Linkages to exclude: Usually Single (due to chaining) or those with very low CPCC.\n")

# --- PLOT DENDROGRAM (Choose the best one or the one asked in text) ---
# Example: Using Average as default or best
hc_final <- hc_average 
plot(hc_final, main = paste("Dendrogram -", best_link), sub = "", xlab = "")

# --- CUT DENDROGRAM ---
# Visually inspect the plot to choose k (e.g., horizontal cut across longest branches)
k_choice <- 3  # <--- MODIFY THIS BASED ON DENDROGRAM INSPECTION
rect.hclust(hc_final, k = k_choice, border = "red")
clusters_hc <- cutree(hc_final, k = k_choice)
# --- Centroids and Sizes of Clusters ---
cat("\n[Hierarchical Clusters]\n")
print(table(clusters_hc)) # Size of clusters

if (!IS_DIST_MATRIX) {
  # If Raw Data: Calculate Centroids
  cat("\n[Cluster Centroids]\n")
  print(aggregate(X_raw, by=list(cluster=clusters_hc), FUN=mean))
}


# Q3: DBSCAN Tuning & Execution -------------------------------------------


# Parameter setup
k_minPts <- 4 # As often requested in exams (e.g., 4 or 5)
cat("\n[Q3 - DBSCAN]\n")
cat("Using minPts =", k_minPts, "\n")

# A) Tune epsilon (eps)
# Look for the "knee" or "elbow" in the k-NN distance plot
kNNdistplot(D, k = k_minPts)
abline(h = 0.25, col = "red", lty = 2) # <--- ADJUST 'h' visually to the knee
cat("Inspect plot. The y-value of the 'elbow' is your optimal eps.\n")

# B) Run DBSCAN
eps_choice <- 0.25 # <--- INSERT VALUE FROM PLOT HERE
dbs <- dbscan(D, eps = eps_choice, minPts = k_minPts)

cat("DBSCAN Results (0 = Noise/Outliers):\n")
print(table(dbs$cluster))

# Satisfaction check:
# - Too much noise (0)? Increase eps.
# - Everything is one cluster? Decrease eps.
# - "Satisfied" if clusters correspond to dense regions without excessive noise.


# Q4: Visualization (MDS) & Distortion ------------------------------------

cat("\n[Q4 - Visualization (MDS)]\n")

# Perform Classical MDS (Principal Coordinate Analysis)
# Transforms distances D into 2D coordinates
mds_fit <- cmdscale(D, k = 2) 
x_mds <- mds_fit[,1]
y_mds <- mds_fit[,2]

# Plot with Hierarchical Clusters
par(mfrow=c(1,2))
plot(x_mds, y_mds, col = clusters_hc + 1, pch = 19, 
     main = "MDS Map (Hierarchical)", xlab = "Dim 1", ylab = "Dim 2")
# Plot with DBSCAN Clusters (0 is usually black/white)
plot(x_mds, y_mds, col = dbs$cluster + 1, pch = 19, 
     main = "MDS Map (DBSCAN)", xlab = "Dim 1", ylab = "Dim 2")
par(mfrow=c(1,1))

# Check for Distortion (Shepard-like logic)
# Compare Euclidean distance in 2D map vs Original Distances
dist_2d <- dist(mds_fit)
# Visualizzazione dell'errore (Stress)
plot(D, dist_2d, main = "Shepard Plot: Map Dist vs True Dist", 
     xlab = "Original Distance", ylab = "MDS 2D Distance")
abline(0, 1, col = "red")

cat("Interpretation of Map Distortion:\n")
cat("- Points ABOVE line: Map exaggerates distance (points look farther than reality).\n")
cat("- Points BELOW line: Map compresses distance (points look closer than reality).\n")


# Q5: Bonferroni Intervals (ONLY IF RAW DATA) -----------------------------

if (!IS_DIST_MATRIX) {
  cat("\n[Q5 - Bonferroni Intervals for Mean & Variance]\n")
  
  # Parameters
  alpha_global <- 0.10             # Global significance level (usually 10% or 5%)
  p_vars <- ncol(X_raw)            # Number of variables
  k_cl   <- k_choice               # Number of clusters chosen
  n_families <- 2 * p_vars * k_cl  # Mean CIs + Var CIs for every variable in every cluster
  
  alpha_bonf <- alpha_global / n_families
  cat("Bonferroni Correction: alpha_adj =", alpha_bonf, "\n")
  
  # Loop through clusters
  for (g in 1:k_cl) {
    cat(paste("\n--- Cluster", g, "---\n"))
    subset_data <- X_raw[clusters_hc == g, ]
    n_g <- nrow(subset_data)
    
    # Loop through variables (columns)
    for (v in names(subset_data)) {
      x_vec <- subset_data[[v]]
      
      # 1. Normality Assumption Check
      # shapiro.test(x_vec) # Optional check
      
      # Statistics
      x_bar <- mean(x_vec)
      s2    <- var(x_vec)
      
      # Mean CI (t-distribution)
      t_crit <- qt(1 - alpha_bonf/2, df = n_g - 1)
      ci_mean <- c(x_bar - t_crit * sqrt(s2/n_g), 
                   x_bar + t_crit * sqrt(s2/n_g))
      
      # Variance CI (Chi-square distribution)
      chi_low  <- qchisq(alpha_bonf/2, df = n_g - 1)
      chi_high <- qchisq(1 - alpha_bonf/2, df = n_g - 1)
      ci_var   <- c((n_g - 1) * s2 / chi_high, 
                    (n_g - 1) * s2 / chi_low)
      
      cat(sprintf("Var: %s | Mean CI: [%.2f, %.2f] | Var CI: [%.2f, %.2f]\n", 
                  v, ci_mean[1], ci_mean[2], ci_var[1], ci_var[2]))
    }
  }
  cat("\nAssumptions: Data within clusters is i.i.d. Gaussian.\n")
} else {
  cat("\n[Q5 - Bonferroni Intervals skipped (Input is Distance Matrix)]\n")
}