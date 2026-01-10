library(fda)
library(ggplot2)
library(dplyr)
library(MASS)        # LDA for quick 1-D discriminants
library(effectsize)

A_DATA_PATH   <- "load.txt"                     # or "temperature.txt" (:contentReference[oaicite:4]{index=4})
A_GRID_LABELS <- paste0("h", 0:24)              # grid columns
A_ARGVALS     <- 0:24                           # numeric grid in [0,24]
A_GROUP_VAR   <- "daytype"                      # e.g., "daytype" (load) or "season" (temperature)
# Choose basis for MODE A:
BASIS_TYPE  <- "fourier"
# If Fourier (often for daily periodic data like load):
N_BASIS_FOURIER <- 13           # <-- per LOAD exercise
# Penalty (applies to A): first-derivative penalty is common for daily profiles
LFD_ORDER_A <- 1
LAMBDA_A    <- 1 

### Smooth the data using a Fourier basis with 13 basis functions.
### Provide a plot of the smoothed observations

datA <- read.table(A_DATA_PATH, header = TRUE, check.names = FALSE)
stopifnot(all(A_GRID_LABELS %in% names(datA)))
Ymat <- as.matrix(datA[, A_GRID_LABELS])  # N x T
Y    <- t(Ymat)                           # T x N (fda expects rows = argvals)
arg  <- A_ARGVALS
group <- if (!is.null(A_GROUP_VAR) && A_GROUP_VAR %in% names(datA)) {
  factor(datA[[A_GROUP_VAR]])
} else NULL
cat(sprintf("[A] Loaded %d curves on a %d-point grid.\n", ncol(Y), nrow(Y)))
mk_basis_fourier <- function(rng, nbasis) {
  create.fourier.basis(rangeval = rng, nbasis = nbasis)
}
rng <- range(arg)
if (BASIS_TYPE == "bspline") {
  basis <- mk_basis_bspline(rng, BS_ORDER_A, BS_BREAKS_A)
} else {
  basis <- mk_basis_fourier(rng, N_BASIS_FOURIER)
}
nbasis <- basis$nbasis

matplot(Ymat, type='l',main='Kwat consumption', xlab='Hours', ylab='365 days')
plot(basis)
data_W.fd <- Data2fd(y=Ymat, argvals=A_ARGVALS, basisobj=basis)

plot.fd(data_W.fd)
data_W.fd$coefs[1:3, 1] # 116.344769 -12.747130  -2.863551
data_W.fd$coefs[, 1]
# Same procedure:
Xsp <- smooth.basis(argvals=time, y=data, fdParobj=basis) # 
Xsp0 = eval.fd(time, Xsp$fd)
plot(Xsp) # same as previously
coef(Xsp$fd)[1:3, 1]








fdpar  <- fdPar(basis, Lfdobj = LFD_ORDER_A, lambda = LAMBDA_A)
fit    <- smooth.basis(arg, Y, fdpar)

# --- Report (A): number of basis fns, mean GCV, mean effective df ---
cat(sprintf("[A:SMOOTH] nbasis = %d | mean GCV = %.4f | mean eff. df ≈ %.2f\n",
            nbasis, mean(fit$gcv), mean(fit$df)))

# --- Plot smoothed curves ---
# Interpretation: Higher λ -> smoother (lower eff.df). For daily load, Fourier K=13 captures daily cycles well.
plot(fit$fd, main = "Smoothed curves", xlab = "hour", ylab = "y(t)") #number of splines = nbasis
# mean eff. df = approx dim of the space in which the smoothed curves live