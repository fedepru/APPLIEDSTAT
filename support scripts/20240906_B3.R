library(fda)
library(ggplot2)
library(dplyr)
library(MASS)        # LDA for quick 1-D discriminants
library(effectsize)

### Apply penalized smoothing to the altitude data using a basis of cubic B-splines, with knots placed
### at each horizontal distance point, penalizing the second-order derivative and using a smoothing
### parameter of lambda = 1. Report the number of splines used and the GCV error

B_DATA_PATH <- "stelvio.txt"    # columns: distance (km), altitude (m)
BS_ORDER_B  <- 4                # cubic B-splines
LFD_ORDER_B <- 2                # penalize curvature
LAMBDA_B    <- 1                # fixed λ (unless asked to tune)
GRADE_FLAGS <- c(0.08, 0.10)

datB <- read.table(B_DATA_PATH, header = TRUE)
x_vec <- datB$distance   # km
y_vec <- datB$altitude   # m
stopifnot(is.numeric(x_vec), is.numeric(y_vec))
cat(sprintf("[B] Loaded single curve with %d points over [%.2f, %.2f].\n",
            length(x_vec), min(x_vec), max(x_vec)))
mk_basis_bspline <- function(rng, order, breaks) {
  create.bspline.basis(rangeval = rng, norder = order, breaks = breaks)
}

rngB   <- range(x_vec)
basisB <- mk_basis_bspline(rngB, BS_ORDER_B, breaks = x_vec)    # knot at every observation
nbasisB <- basisB$nbasis
fdparB <- fdPar(basisB, Lfdobj = LFD_ORDER_B, lambda = LAMBDA_B)
fitB   <- smooth.basis(x_vec, y_vec, fdparB)

cat(sprintf("[B:SMOOTH] nbasis = %d | GCV = %.4f | eff. df ≈ %.2f\n",
            nbasisB, fitB$gcv, fitB$df))
### nbasis = 128 | GCV = 1490.3615

### Estimate the approximate dimension of the space in which the fitted curve lives.
### Provide a plot of the fitted curve along with its first derivative

### approx dim of the space where fitted curve lives: eff. df ≈ 14.32

# Evaluate curve and derivatives; plots with interpretations
yhat   <- eval.fd(x_vec, fitB$fd, Lfd = 0)
yprime <- eval.fd(x_vec, fitB$fd, Lfd = 1)

par(mfrow = c(1,2))
plot(x_vec, yhat, type = "l", lwd = 2, col = 2,
     xlab = "distance (km)", ylab = "altitude (m)", main = "Fitted curve")
points(x_vec, y_vec, pch = 16, cex = 0.5)
plot(x_vec, yprime, type = "l", lwd = 2, col = 4,
     xlab = "distance (km)", ylab = "slope (m/km)", main = "First derivative")
par(mfrow = c(1,1))

### Determine the value of lambda that minimizes the GCV error, and report the corresponding GCV error.
## Use a grid search with log10(lambda) values ranging from [-1,3] ....
### Refit the smoothed curve using this optimal value

lambda_grid_search <- function(arg, y, basis_obj, Lfd, log10_min = -1, log10_max = 3, by = 0.5) {
  lam_grid <- 10^seq(log10_min, log10_max, by = by)
  gcv_grid <- numeric(length(lam_grid))
  for (i in seq_along(lam_grid)) {
    fdpar_tmp <- fdPar(fdobj = basis_obj, Lfdobj = Lfd, lambda = lam_grid[i])
    gcv_grid[i] <- smooth.basis(arg, y, fdpar_tmp)$gcv
  }
  list(lambda = lam_grid, gcv = gcv_grid)
}
res <- lambda_grid_search(x_vec, y_vec, basisB, LFD_ORDER_B, -1, 3, 0.5)
opt_idx <- which.min(res$gcv)
opt_lambda <- res$lambda[opt_idx]
### GCV opt = 1460.233    lamda opt = 3.162278

fdparB_opt <- fdPar(basisB, Lfdobj = LFD_ORDER_B, lambda = opt_lambda)
fitB_opt   <- smooth.basis(x_vec, y_vec, fdparB_opt)

cat(sprintf("[B:SMOOTH] nbasis = %d | GCV = %.4f | eff. df ≈ %.2f\n",
            nbasisB, fitB_opt$gcv, fitB_opt$df))
### opt: nbasis = 128 | GCV = 1460.2334 | eff. df ≈ 10.99

### Calculate the slope1 at the steepest point of the ascent

yprime <- eval.fd(x_vec, fitB_opt$fd, Lfd = 1)

# Steepest ascent
idx_max <- which.max(yprime)
x_max   <- x_vec[idx_max]
slope_m_per_km <- yprime[idx_max]
grade_pct      <- slope_m_per_km / 10  # m/km -> %

cat(sprintf("[B:STEEPEST] max slope ≈ %.1f m/km (≈ %.1f%%) at x ≈ %.2f km\n",
            slope_m_per_km, grade_pct, x_max))
### [B:STEEPEST] max slope ≈ 107.4 m/km (≈ 10.7%) at x ≈ 11.20 km