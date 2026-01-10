library(fda)
library(ggplot2)
library(dplyr)
library(MASS)        # LDA for quick 1-D discriminants
library(effectsize)

### Apply penalized smoothing to the temperature data using a basis of quadratic B-splines, with knots
### placed at each observed time point (i.e., every hour). Penalize the first-order derivative and use
### a smoothing parameter of lambda = 1. Report the number of splines used and the mean generalized
### cross-validation (GCV) error. Provide a plot of the smoothed temperature curves. What is the
### approximate dimension of the space in which the smoothed curves live?

MODE <- "A"  # "A" = functional sample; "B" = single curve

# A) Functional sample
# Load (hourly GW over 2023; groups = daytype)  :contentReference[oaicite:3]{index=3}
A_DATA_PATH   <- "temperature.txt"              
A_GRID_LABELS <- paste0("h", 0:24)              # grid columns
A_ARGVALS     <- 0:24                           # numeric grid in [0,24]
A_GROUP_VAR   <- "season"                      
# Choose basis for MODE A:
BASIS_TYPE  <- "bspline"
# If B-splines (often for temps): order = degree+1 (Quadratic -> 3; Cubic -> 4)
BS_ORDER_A  <- 3
BS_BREAKS_A <- A_ARGVALS        # knots at each hour (when using bspline)
# Penalty (applies to A): first-derivative penalty is common for daily profiles
LFD_ORDER_A <- 1
LAMBDA_A    <- 1                # ω (fixed unless the exercise asks to tune)

datA <- read.table(A_DATA_PATH, header = TRUE, check.names = FALSE)
stopifnot(all(A_GRID_LABELS %in% names(datA)))
Ymat <- as.matrix(datA[, A_GRID_LABELS])  # N x T
Y    <- t(Ymat)                           # T x N (fda expects rows = argvals)
arg  <- A_ARGVALS
group <- if (!is.null(A_GROUP_VAR) && A_GROUP_VAR %in% names(datA)) {
  factor(datA[[A_GROUP_VAR]])
} else NULL
cat(sprintf("[A] Loaded %d curves on a %d-point grid.\n", ncol(Y), nrow(Y)))

mk_basis_bspline <- function(rng, order, breaks) {
  create.bspline.basis(rangeval = rng, norder = order, breaks = breaks)
}

rng <- range(arg)
basis <- mk_basis_bspline(rng, BS_ORDER_A, BS_BREAKS_A)
nbasis <- basis$nbasis
fdpar  <- fdPar(basis, Lfdobj = LFD_ORDER_A, lambda = LAMBDA_A)
fit    <- smooth.basis(arg, Y, fdpar)

# --- Report (A): number of basis fns, mean GCV, mean effective df ---
cat(sprintf("[A:SMOOTH] nbasis = %d | mean GCV = %.4f | mean eff. df ≈ %.2f\n",
            nbasis, mean(fit$gcv), mean(fit$df)))

# Plot smoothed curves
# Interpretation: Higher λ -> smoother (lower eff.df). For daily load, Fourier K=13 captures daily cycles well.
plot(fit$fd, main = "Smoothed curves", xlab = "hour", ylab = "y(t)")
# mean eff. df = approx dim of the space in which the smoothed curves live
### mean eff. df ≈ 9.83 => 10

### Conduct a FPCA on the smoothed functions. What proportion of the total variance is explained by the
### second principal component? From a dimensionality reduction perspective, how many principal
### components would you retain? Justify your choice.

# FPCA settings (A)
NHARM       <- 10               # try up to this many harmonics
CUM_PVE_TGT <- 0.95             # target cumulative variance explained (e.g., 95%)

nh <- min(NHARM, basis$nbasis)
pca <- pca.fd(fit$fd, nharm = nh, centerfns = TRUE)

# Proportion of variance explained by PC2, and PCs for 95% cum. PVE
pve <- pca$varprop
cum_pve <- cumsum(pve)
### 2nd PC explains the 0.0487272849  of the total variance
### I would retain 2 PC, they are able to explain 0.9820061  of the total variance

### Provide a plot showing the effect of the second principal component
par(mfrow = c(1,2))
plot(pca, harm = 1, nx = 101, pointplot = TRUE)
plot(pca, harm = 2, nx = 101, pointplot = TRUE)
par(mfrow = c(1,1))
### plot(pca, harm = 2, nx = 101, pointplot = TRUE)

### Is the representation given by the first principal component satisfying for distinguishing seasons?
### Support your answer with a plot.

# plot degli scores dei primi due PC divisi per gruppo
scores <- pca$scores
plot(scores[, 1], scores[, 2], col = as.factor(group),
     pch = 19, xlab = "PC1", ylab = "PC2", main = "Scores of First Two Principal Components")
legend("topright", legend = c("Autumn", "Spring", "Summer", "Winter"), col = as.factor(levels(as.factor(season))), pch = 19)

plot(scores[, 1] ~ as.factor(season)) # by PC1 we see difference among Summer, Winter and Autumn wiht Spring
plot(scores[, 2] ~ as.factor(season)) # by PC2 it's difficlut to ifnd difference among seasons

### Could we successfully classify seasons based on the representation given by the second PC?

## No, the 2nd PC is not able to distinguishing season