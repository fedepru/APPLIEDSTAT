## ===================== CLUSTERING TEMPLATE (minimal) =====================
## Supports:
## - RAW features data frame  -> Euclidean dist + Average HClust + DBSCAN (+ Bonferroni CIs)
## - DISTANCE matrix (square, symmetric with row/col names) -> HClust + DBSCAN + MDS plot
## =========================================================================

library(cluster)
library(dbscan)

## -------------------- 1) LOAD DATA --------------------------------------
## Set ONE of the two below:

# (A) RAW features (e.g., sharks)
# dat <- read.table("sharks.txt", header = TRUE)
# is_dist <- FALSE

# (B) DISTANCE matrix (e.g., metro)
# m <- as.matrix(read.table("metro.txt", header = TRUE, row.names = 1, check.names = FALSE))
# dat <- as.dist(m)
# is_dist <- TRUE

## -------------------- 2) PREP -------------------------------------------
if (!exists("dat")) stop("Load 'dat' as data.frame (raw) or as.dist (distance).")
if (!exists("is_dist")) stop("Set is_dist = FALSE (raw) or TRUE (distance).")

if (!is_dist) {
  X <- as.data.frame(dat)
  D <- dist(X, method = "euclidean")
} else {
  D <- dat
}

## -------------------- 3) HIERARCHICAL (Average linkage) -----------------
hc <- hclust(D, method = "average")
plot(hc, main = "Hierarchical Clustering (Average linkage)")

## choose k by dendrogram (edit as needed)
k <- 4
cl_hc <- cutree(hc, k = k)
table(cl_hc)

## centroids & sizes (only for RAW features)
if (!is_dist) {
  sizes <- table(cl_hc)
  cents <- aggregate(X, list(cluster = cl_hc), mean)
  print(sizes); print(cents)
}

## -------------------- 4) DBSCAN -----------------------------------------
## For distance input, dbscan() accepts a 'dist' object directly.
## For raw features, it will compute internal distances (metric = "euclidean").

minPts <- 5

if (!is_dist) {
  ## quick eps heuristic (optional): 5-NN distance “knee”
  # kNNdistplot(X, k = minPts); abline(h = ..., col = 2, lty = 2)
  eps <- 0.25   # <- set after inspecting kNNdistplot (or grid-search)
  dbs <- dbscan(X, eps = eps, minPts = minPts)
} else {
  ## eps grid on distances (rounded by 0.05 as in your style)
  eps_grid <- seq(0.05, 0.50, by = 0.05)
  scan <- sapply(eps_grid, function(e) {
    fit <- dbscan(D, eps = e, minPts = minPts)
    c(n_clusters = length(setdiff(unique(fit$cluster), 0)),
      noise = sum(fit$cluster == 0))
  })
  colnames(scan) <- sprintf("eps=%.2f", eps_grid); print(t(scan))
  eps <- 0.25  # <- pick from scan
  dbs <- dbscan(D, eps = eps, minPts = minPts)
}

table(dbs$cluster)   # 0 = noise

## -------------------- 5) 2D VISUALIZATION -------------------------------
if (!is_dist) {
  ## raw: use PCA for a quick 2D map
  p <- prcomp(scale(X, center = TRUE, scale = TRUE))
  S <- p$x[, 1:2]
  cols <- cl_hc
  plot(S, col = cols, pch = 19, asp = 1,
       main = "2D map (PCA) colored by HClust")
  ## switch to DBSCAN colors:
  # cols <- dbs$cluster + 1L
  # plot(S, col = cols, pch = 19, asp = 1, main = "2D map (PCA) colored by DBSCAN")
} else {
  ## distances: use classical MDS
  S <- cmdscale(D, k = 2)
  cols <- cl_hc
  plot(S, col = cols, pch = 19, asp = 1,
       main = "2D map (MDS) colored by HClust", xlab = "MDS1", ylab = "MDS2")
  ## switch to DBSCAN colors:
  # cols <- dbs$cluster + 1L
  # plot(S, col = cols, pch = 19, asp = 1, main = "2D map (MDS) colored by DBSCAN")
  
  ## distortion check (under/over-estimation)
  plot(D, dist(S), xlab = "Original distances", ylab = "Map distances",
       main = "Distance preservation (ideal: y = x)")
  abline(0, 1, lty = 2)
}

## -------------------- 6) BONFERRONI CIs (optional; RAW features only) ---
## Global 90% for means & variances of a chosen variable per cluster.
## Example for 'circumference' column.
if (!is_dist) {
  alpha <- 0.10
  var_name <- grep("circ", names(X), value = TRUE, ignore.case = TRUE)
  if (length(var_name) == 0) var_name <- names(X)[2]  # fallback
  
  K  <- length(unique(cl_hc))
  m_tot <- 2 * K  # (mean + variance) per cluster
  a_bf <- alpha / m_tot
  
  CI_mean <- CI_var <- vector("list", K)
  names(CI_mean) <- names(CI_var) <- paste0("cl_", 1:K)
  
  for (g in 1:K) {
    x <- X[cl_hc == g, var_name]
    n <- length(x); xb <- mean(x); s2 <- var(x); s <- sqrt(s2)
    # mean
    tcrit <- qt(1 - a_bf/2, df = n - 1)
    CI_mean[[g]] <- c(inf = xb - tcrit * s / sqrt(n),
                      center = xb,
                      sup = xb + tcrit * s / sqrt(n))
    # variance
    qL <- qchisq(1 - a_bf/2, df = n - 1)
    qU <- qchisq(a_bf/2,     df = n - 1)
    CI_var[[g]] <- c(inf = (n - 1) * s2 / qL,
                     center = s2,
                     sup = (n - 1) * s2 / qU)
  }
  CI_mean; CI_var
}
## =================== END TEMPLATE =======================================
