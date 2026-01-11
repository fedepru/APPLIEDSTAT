# ========================= MASTER FDA TEMPLATE (R + fda) =========================
# Supports:
#   A) Functional sample on a shared grid (e.g., load by hour, temperatures by hour)
#      -> smoothing (B-spline or Fourier) + FPCA + grouping/classification
#   B) Single noisy curve y(x) (e.g., road elevation)
#      -> smoothing (B-spline) + derivatives + steepest slope + grade flags
#
# Each section contains concise "Interpretation" comments near the relevant code.

suppressPackageStartupMessages({
  library(fda)
  library(ggplot2)
  library(dplyr)
  library(MASS)        # LDA for quick 1-D discriminants
  library(effectsize)  # eta^2 for ANOVA interpretations
})

# ------------------------------- 0) MODE & INPUTS --------------------------------
MODE <- "A"  # "A" = functional sample; "B" = single NOISY curve

# MODE A: functional sample ------------------------------------------------------------------
A_DATA_PATH   <- "data.txt"                     # filename

A_GRID_LABELS <- paste0("h", a:b)              # grid columns in data file (e.g., h0, h1, ..., h24)
A_ARGVALS     <- a:b                           # numeric grid in DOMAIN (e.g. [0,24])

A_GROUP_VAR   <- "factor"                      #factor variable that indicates groups,
                                                #if not present set to NULL
# Choose basis for MODE A:
BASIS_TYPE  <- "fourier"    # "bspline" or "fourier"

# If B-splines: (order=LFD+2 rule of thumb)
BS_ORDER_A  <- 3                # order of B-spline (Quadratic ->3,Cubic -> 4)
BS_BREAKS_A <- A_ARGVALS        # knots at each point in domain (e.g., 0,1,2,...,24)
                                # (when using bspline)

# If Fourier (often for daily periodic data):
N_BASIS_FOURIER <- 13           #specified in exercise

# Penalty (applies to A):
LFD_ORDER_A <- 1                # order of derivative to penalize: usually specified in text
                                # if not specified:
                                # good rule: order - 2 (e.g., for order=4, use LFD=2)
                                # Fourier: usually LFD=1 (penalize roughness)
LAMBDA_A    <- 1                # lambda: specified by text

# FPCA settings (A)
NHARM       <- 10               # try up to this many harmonics (upper bound for fdpcs)
CUM_PVE_TGT <- 0.95             # target cumulative variance explained (e.g., 95%)


# MODE B: single NOISY CURVE ----------------------------------------------

B_DATA_PATH <- "data.txt"       # filename
BS_ORDER_B  <- 4                # order of B-spline (Quadratic ->3,Cubic -> 4)
LFD_ORDER_B <- 2                # order of derivative to penalize: usually specified in text
                                # if not specified: use same rule: LFD = order - 2
LAMBDA_B    <- 1                # fixed λ provided by text
GRADE_FLAGS <- c(0.08, 0.10)    # optional grade thresholds (8%, 10%)

# ------------------------------- 1) LOAD DATA ------------------------------------
if (MODE == "A") {
  datA <- read.table(A_DATA_PATH, header = TRUE, check.names = FALSE)
  stopifnot(all(A_GRID_LABELS %in% names(datA)))
  Ymat <- as.matrix(datA[, A_GRID_LABELS])  # N x T
  Y    <- t(Ymat)                           # T x N (fda expects rows = argvals)
  arg  <- A_ARGVALS
  group <- if (!is.null(A_GROUP_VAR) && A_GROUP_VAR %in% names(datA)) {
    factor(datA[[A_GROUP_VAR]])
  } else NULL
  cat(sprintf("[A] Loaded %d curves on a %d-point grid.\n", ncol(Y), nrow(Y)))
}

if (MODE == "B") {
  datB <- read.table(B_DATA_PATH, header = TRUE)
  x_vec <- datB$distance   # km
  y_vec <- datB$altitude   # m
  stopifnot(is.numeric(x_vec), is.numeric(y_vec))
  cat(sprintf("[B] Loaded single curve with %d points over [%.2f, %.2f].\n",
              length(x_vec), min(x_vec), max(x_vec)))
}

# ------------------------------- 2) BASIS SETUP ----------------------------------
#choose if using b-splines or fourier basis
mk_basis_bspline <- function(rng, order, breaks) {
  create.bspline.basis(rangeval = rng, norder = order, breaks = breaks)
}
mk_basis_fourier <- function(rng, nbasis) {
  create.fourier.basis(rangeval = rng, nbasis = nbasis)
}

# ------------------------------- 3) SMOOTHING ------------------------------------
if (MODE == "A") {
  rng <- range(arg)
  if (BASIS_TYPE == "bspline") {
    basis <- mk_basis_bspline(rng, BS_ORDER_A, BS_BREAKS_A)
  } else {
    basis <- mk_basis_fourier(rng, N_BASIS_FOURIER)
  }
  nbasis <- basis$nbasis
  fdpar  <- fdPar(basis, Lfdobj = LFD_ORDER_A, lambda = LAMBDA_A)
  fit    <- smooth.basis(arg, Y, fdpar)
  
  # --- Report (A): number of basis functions, mean GCV, mean effective df ---
  cat(sprintf("[A:SMOOTH] nbasis = %d | mean GCV = %.4f | mean eff. df ≈ %.2f\n",
              nbasis, mean(fit$gcv), mean(fit$df)))
  
  # --- Plot smoothed curves ---
  # Interpretation: Higher λ -> smoother (lower eff.df). For daily load, Fourier K=13 captures daily cycles well.
  plot(fit$fd, main = "Smoothed curves", xlab = "x", ylab = "y(x)") #number of splines = nbasis
  # mean eff. df = approx dim of the space in which the smoothed curves live
  
} else if (MODE == "B") {
  rngB   <- range(x_vec)
  basisB <- mk_basis_bspline(rngB, BS_ORDER_B, breaks = x_vec)    # knot at every observation
  nbasisB <- basisB$nbasis
  fdparB <- fdPar(basisB, Lfdobj = LFD_ORDER_B, lambda = LAMBDA_B)
  fitB   <- smooth.basis(x_vec, y_vec, fdparB)
  
  # --- Report (B) ---
  cat(sprintf("[B:SMOOTH] nbasis = %d | GCV = %.4f | eff. df ≈ %.2f\n",
              nbasisB, fitB$gcv, fitB$df))
  
  # Evaluate curve and derivatives; plots with interpretations
  yhat   <- eval.fd(x_vec, fitB$fd, Lfd = 0)
  yprime <- eval.fd(x_vec, fitB$fd, Lfd = 1)
  
  par(mfrow = c(1,2))
  plot(x_vec, yhat, type = "l", lwd = 2, col = 2,
       xlab = "x", ylab = "y(x)", main = "Fitted curve")
  points(x_vec, y_vec, pch = 16, cex = 0.5)
  plot(x_vec, yprime, type = "l", lwd = 2, col = 4,
       xlab = "x", ylab = "y'(x)", main = "First derivative")
  par(mfrow = c(1,1))
}

# ------------------------------- 4) FPCA (MODE A) --------------------------------
if (MODE == "A") {
  nh <- min(NHARM, basis$nbasis) #upper bound on number of harmonics
  pca <- pca.fd(fit$fd, nharm = nh, centerfns = TRUE)
  
  # Proportion of variance explained by PC2, and PCs for 95% cum. PVE
  pve <- pca$varprop
  cum_pve <- cumsum(pve)
  #pcs needed to reach target cumulative pve 
  k95 <- which(cum_pve >= CUM_PVE_TGT)[1]
  cat(sprintf("[A:FPCA] PC2 PVE = %.2f%% | PCs to reach %.0f%% = %d (cum = %.2f%%)\n",
              100*pve[2], 100*CUM_PVE_TGT, k95, 100*cum_pve[k95]))
  
  # --- Scree / PVE plot ---
  # Interpretation: Choose smallest k with cum PVE ≥ target; elbow = diminishing returns.
  plot(100*pve, type = "b", xlab = "PC", ylab = "PVE (%)", main = "FPCA variance explained")
  abline(h = 95, lty = 2, col = "grey60")
  
  # --- Interpretation plots for PC1 and PC2 (mean ± c*phi_j) ---
  # Interpretation: Where mean± shifts are large,
  #that PC modulates that time segment (e.g., morning/evening peaks).
  par(mfrow = c(1,2))
  plot(pca, harm = 1, nx = 101, pointplot = TRUE)
  plot(pca, harm = 2, nx = 101, pointplot = TRUE)
  par(mfrow = c(1,1))
  
  # Convenience: scores dataframe
  scores <- as.data.frame(pca$scores)
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  if (!is.null(group)) scores$group <- group
}

# plot degli scores dei primi due PC divisi per gruppo: classificazione visiva
scores <- pca$scores
plot(scores[, 1], scores[, 2], col = as.factor(group),
     pch = 19, xlab = "PC1", ylab = "PC2", main = "Scores of First Two Principal Components")
legend("topright", legend = as.factor(levels(as.factor(group))), col = as.factor(levels(as.factor(group))), pch = 19)

#HOW TO INTERPRET: IF THE COLOURS ARE NEATLY SEPARATED,
#IT MEANS THAT THE FIRST TWO PCS ARE USEFUL TO DISTINGUISH THE GROUPS
#diagonal line separates groups: we need bothm PC1 and PC2 to separate them
#horizontal/vertical line: only the respective PC is enough to separate them

plot(scores[, 1] ~ as.factor(group)) 
plot(scores[, 2] ~ as.factor(group)) 



# -------- 5) Group separation & quick classification (MODE A, if groups exist) ---
if (MODE == "A" && !is.null(group)) {
  # --- PC1 vs PC2 scatter ---
  # Interpretation: Visible clustering suggests these PCs carry factor information.
  print(
    ggplot(scores, aes(x = PC1, y = PC2, color = group)) +
      geom_point(alpha = 0.8) +
      theme_minimal() + labs(title = "FPCA scores by group", x = "PC1", y = "PC2")
  )
  
  # --- PC1 and PC2 by group: ANOVA + eta^2 ---
  # Interpretation: Significant ANOVA with non-trivial η² => component is useful for separation.
  for (j in 1:2) {
    g <- ggplot(scores, aes(group, .data[[paste0("PC", j)]], fill = group)) +
      geom_boxplot(outlier.alpha = 0.3) +
      theme_minimal() + guides(fill = "none") +
      labs(title = paste0("PC", j, " across groups"))
    print(g)
    aovj <- aov(as.formula(paste0("PC", j, " ~ group")), data = scores)
    print(summary(aovj))
    suppressWarnings(print(eta_squared(aovj)))
  }
  
  # --- “Is PC1+PC2 enough to distinguish groups?” (2D check) ---
  # Interpretation: If scatter is well-separated and η² sizable for PC1 and/or PC2, 2D scores are informative.
  
  # --- 1D representation *better than PC1* (requested in LOAD exercise) ---
  # Fisher/LDA 1-D discriminant using top K PCs; compare to PC1 alone.
  # Interpretation: If LDA-1D outperforms PC1 (accuracy/AUC), it’s a strictly better 1-D summary for discrimination.
  K <- min(10, ncol(pca$scores))     # use up to 10 PCs (or adjust)
  dat_lda <- data.frame(y = scores$group, pca$scores[, 1:K, drop=FALSE])
  lda_fit <- lda(y ~ ., data = dat_lda, CV = TRUE)  # LOOCV predictions
  # Calculate the discriminant scores using predict()
  lda_pred <- predict(lda_fit, dat_lda)
  lda_proj <- as.numeric(lda_pred$x)
  acc_lda <- mean(lda_fit$class == dat_lda$y)
  
  # Baseline: PC1-only thresholding via LDA on PC1 (still 1-D)
  lda_pc1  <- lda(y ~ PC1, data = scores, CV = TRUE)
  acc_pc1  <- mean(lda_pc1$class == scores$group)
  
  cat(sprintf("[A:LDA-1D] LOOCV acc (PC1-only) = %.1f%% | LOOCV acc (LDA on top-%d PCs) = %.1f%%\n",
              100*acc_pc1, K, 100*acc_lda))
  
  # Plot the 1-D discriminant vs. group (violin/box) for interpretability
  tmp <- data.frame(group = scores$group, LDA1 = lda_proj)
  print(
    ggplot(tmp, aes(group, LDA1, fill = group)) +
      geom_violin(trim = FALSE, alpha = 0.3) +
      geom_boxplot(width = 0.15, outlier.alpha = 0.3) +
      theme_minimal() + guides(fill = "none") +
      labs(title = "1-D Fisher/LDA projection (better than PC1?)", y = "LDA score")
  )
}

# ------------------- 6) Steepest slope & grade flags (MODE B) --------------------
if (MODE == "B") {
  yprime <- eval.fd(x_vec, fitB$fd, Lfd = 1)
  
  # Steepest ascent
  idx_max <- which.max(yprime)
  x_max   <- x_vec[idx_max]
  slope_m_per_km <- yprime[idx_max]
  grade_pct      <- slope_m_per_km / 10  # m/km -> %
  
  cat(sprintf("[B:STEEPEST] max slope ≈ %.1f m/km (≈ %.1f%%) at x ≈ %.2f km\n",
              slope_m_per_km, grade_pct, x_max))
  
  # Flag contiguous sections ≥ thresholds (optional interpretation for cyclists/operations)
  runs_where <- function(idx) { if (!length(idx)) return(list()); split(idx, cumsum(c(1, diff(idx) != 1))) }
  grade_vec <- yprime / 10  # decimal
  for (thr in GRADE_FLAGS) {
    idx <- which(grade_vec >= thr)
    segs <- runs_where(idx)
    if (length(segs) == 0) {
      cat(sprintf("[B:FLAGS] No sections ≥ %.0f%% grade.\n", 100*thr))
    } else {
      cat(sprintf("[B:FLAGS] Sections ≥ %.0f%% grade:\n", 100*thr))
      for (sg in segs) {
        cat(sprintf("  from %.2f km to %.2f km | median ≈ %.1f%% | max ≈ %.1f%%\n",
                    x_vec[min(sg)], x_vec[max(sg)],
                    100*median(grade_vec[sg]), 100*max(grade_vec[sg])) )
      }
    }
  }
}

# ------------------------------- 7) λ tuning (Optional) --------------------------
# Use ONLY when the exercise asks to pick λ by GCV. Keep separate from fixed-λ answers.
#adjust the log10_min and log10_max to cover the interval given by the exercise ex: 0.1 ,1000 => -1,3
lambda_grid_search <- function(arg, y, basis_obj, Lfd, log10_min = -3, log10_max = 3, by = 0.5) {
  lam_grid <- 10^seq(log10_min, log10_max, by = by)
  gcv_grid <- numeric(length(lam_grid))
  for (i in seq_along(lam_grid)) {
    fdpar_tmp <- fdPar(fdobj = basis_obj, Lfdobj = Lfd, lambda = lam_grid[i])
    gcv_grid[i] <- smooth.basis(arg, y, fdpar_tmp)$gcv
  }
  list(lambda = lam_grid, gcv = gcv_grid)
}
# Example (B only):  # res <- lambda_grid_search(x_vec, y_vec, basisB, LFD_ORDER_B, -1, 3, 0.5)
# =================================================================================
