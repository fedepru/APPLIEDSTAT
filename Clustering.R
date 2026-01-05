# Required libraries
library(mvtnorm)
library(MVN)
library(rgl)
library(car)
library(dbscan)
library(cluster)
library(fields)

# -------------------------------------------------------------------------
# 1. Load Data
# -------------------------------------------------------------------------
# Replace the filename depending on the dataset you want to analyze
# Examples:
#   data <- read.table("metro.txt", h=TRUE)     # city distances
#   data <- read.table("sharks.txt", h=TRUE)    # shark measurements
data <- read.table("yourfile.txt", h=TRUE)

matplot(data, main="Raw Data Plot")

# -------------------------------------------------------------------------
# 2. Hierarchical Clustering (Euclidean + Average linkage)
# -------------------------------------------------------------------------
data.dist <- dist(data, method="euclidean")
data.hclust_avg <- hclust(data.dist, method="average")

plot(data.hclust_avg, main="Hierarchical Clustering Dendrogram (Average Linkage)")

# Choose number of clusters (set k accordingly after inspecting the dendrogram)
k <- 4
clusters <- cutree(data.hclust_avg, k=k)
clusters

# -------------------------------------------------------------------------
# 3. Centroids and Sizes of Clusters
# -------------------------------------------------------------------------
centroids <- lapply(1:k, function(i) colMeans(data[clusters==i,]))
sizes <- sapply(1:k, function(i) sum(clusters==i))

centroids
sizes

# -------------------------------------------------------------------------
# 4. Bonferroni Intervals for Means and Variances
#    (if the dataset is measurement-based, e.g. shark body size)
# -------------------------------------------------------------------------
alpha <- 0.1
p <- ncol(data)   # number of variables
k_total <- 2 * p * k   # 2 statistics (mean+var) × variables × clusters
alpha.bf <- alpha / k_total

bonferroni_results <- list()

for (i in 1:k) {
  ni <- sizes[i]
  cov_i <- cov(data[clusters==i,])
  mean_i <- colMeans(data[clusters==i,])
  
  # Means
  ICmean <- cbind(
    inf = mean_i - sqrt(diag(cov_i)/ni) * qt(1 - alpha.bf/2, ni-1),
    center = mean_i,
    sup = mean_i + sqrt(diag(cov_i)/ni) * qt(1 - alpha.bf/2, ni-1)
  )
  
  # Variances
  ICvar <- cbind(
    inf = diag(cov_i)*(ni-1) / qchisq(1 - alpha.bf/2, ni-1),
    center = diag(cov_i),
    sup = diag(cov_i)*(ni-1) / qchisq(alpha.bf/2, ni-1)
  )
  
  bonferroni_results[[paste0("Cluster_", i)]] <- list(
    Means = ICmean,
    Variances = ICvar
  )
}

bonferroni_results

# -------------------------------------------------------------------------
# 5. Interpretation
# -------------------------------------------------------------------------
# - Number of clusters chosen
# - Centroids (approximate "typical" points per cluster)
# - Cluster sizes
# - Bonferroni confidence intervals highlight uncertainty in means and variances
# - In application (e.g. sharks), cluster count might correspond to # of seismic events
# - In cities, clusters may represent regional groupings
