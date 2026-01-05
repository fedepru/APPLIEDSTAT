#### 0) Packages---------------------------------------------------------------

# Core
library(car)        # Anova(), linearHypothesis(), boxM()
library(MVN)        # mvn() diagnostics for MVN
library(MASS)       # lda(), qda()
library(heplots)
# Optional (for robustness / alternatives)
# install.packages(c("vegan","MANOVA.RM","rrcov"))
# library(vegan)     # adonis2(): permutation MANOVA
# library(MANOVA.RM) # RM(): robust rank-based MANOVA for factorial designs
# library(rrcov)     # Wilks.test(): robust one-way MANOVA variants

#### 1) Data & variables---------------------------------------------------------------

# Replace with your file and columns
data <- read.table("your_data.txt", header = TRUE)

# Responses (multivariate)
Ynames <- c("Y1", "Y2")            # <-- replace with your response columns
Y <- as.matrix(data[ , Ynames])

# Two binary (or K-level) factors
F1 <- factor(data$Factor1)         # e.g., air_purifier (0/1)
F2 <- factor(data$Factor2)         # e.g., auto_watering (0/1)

# Interaction (often useful for diagnostics/classification)
F12 <- interaction(F1, F2, drop = TRUE)

# Design check (balanced/unbalanced)
table(F1, F2)

#### 2) Fit the 2-way MANOVA (with interaction)---------------------------------------------------------------
# Full model: includes main effects and interaction
fit <- manova(Y ~ F1 * F2)
summary.manova(fit)

#extra
{# Report multiple test statistics (Wilks, Pillai, HL, Roy)
  man_wilks  <- summary(fit, test = "Wilks")
  man_pillai <- summary(fit, test = "Pillai")
  man_hl     <- summary(fit, test = "Hotelling-Lawley")
  man_roy    <- summary(fit, test = "Roy")
  
  man_wilks
  man_pillai
  man_hl
  man_roy
}

# Exam write-up guide:
# - State H0 for each effect (e.g., no F1 effect on the mean vector; similarly for F2 and F1:F2).
# - Report the chosen statistic (commonly Wilks’ Lambda) with df and p-value.
# - Conclude at α = 0.05 (or as given).

#### 3) Assumptions---------------------------------------------------------------

# (A) Multivariate Normality (MVN) in each cell of F1:F2
Ps <- sapply(levels(F12), function(lv) {
  mvn(Y[F12==lv, ])$multivariateNormality$`p value`
})

# (B) Homogeneity of covariance matrices (Box’s M)
## One covariance matrix per F12 level
S_list <- lapply(levels(F12), function(lv) cov(Y[F12 == lv, , drop = FALSE]))
names(S_list) <- paste0("S_", levels(F12))  # e.g., S_0.0, S_1.0, ...

S_list

summary(boxM(Y, F12))

# (C) Balance check matters:
# - Severe imbalance + heteroscedasticity => prefer Pillai’s trace and/or robust/permutation approaches.

#### 4) Should I modify the model?---------------------------------------------------------------
# Use this commented checklist to justify changes in the exam.

# 4.1 If MVN is *approximately* OK (no extreme outliers, QQ plots ~linear)
#     AND Box’s M is NOT significant:
#     -> Assumptions reasonable. Keep the standard MANOVA.
#
# 4.2 If MVN is borderline (some non-normality) OR Box’s M is significant but
#     sample sizes are reasonably similar across cells:
#     -> Still report all tests, but put more weight on Pillai’s trace
#        (most robust to assumption violations and imbalance).
#
# 4.3 If MVN clearly fails (heavy tails/outliers) and transformations help:
#     -> Consider transforming responses (e.g., log or Box–Cox per-variable),
#        refit MANOVA, and re-check assumptions.
#        Example:
#        Ytr <- cbind(log1p(data$Y1), log1p(data$Y2))
#        fit_tr <- manova(Ytr ~ F1 * F2)
#        summary(fit_tr, test="Pillai")
#
# 4.4 If heteroscedasticity is strong (Box’s M highly significant) + imbalance:
#     -> Report MANOVA but emphasize Pillai’s trace.
#     -> Consider *nonparametric / permutation MANOVA* as a sensitivity check:
#        (requires distance-based approach)
#        library(vegan)
#        df4 <- data.frame(Y, F1, F2)
#        adonis2_res <- adonis2(Y ~ F1 * F2, data = df4, method = "euclidean", permutations = 4999)
#        adonis2_res
#        (Explain it tests differences in multivariate centroids via permutations—robust to MVN.)
#
# 4.5 If design is factorial and you want a robust rank-based inference:
#     -> Use MANOVA.RM::RM() (rank-based, heteroscedasticity-robust):
#        library(MANOVA.RM)
#        dfRM <- data.frame(Y1 = Y[,1], Y2 = Y[,2], F1, F2)
#        rmfit <- RM(cbind(Y1, Y2) ~ F1 * F2, data = dfRM, iter = 10000)
#        summary(rmfit)
#        (State this as a robustness check, complementary to classical MANOVA.)
#
# 4.6 If interaction (F1:F2) is significant:
#     -> Interpret simple effects (effect of F1 within each level of F2, and vice versa).
#        You can subset and run one-way MANOVAs, or use linearHypothesis on fitted lm():
#        fit_lm <- lm(Y ~ F1*F2)
#        # Example contrast (replace with your coding)
#        # linearHypothesis(fit_lm, c("F11 - F12 = 0", "F11:Y1 - F12:Y1 = 0", ...)) # advanced
#     -> For follow-ups, proceed to univariate ANOVAs + adjusted multiple comparisons.

#### 5) Univariate follow-ups---------------------------------------------------------------
####  (which variables drive the effect?)
# After a significant multivariate test, examine univariate ANOVAs
summary.aov(fit)   # Lists per-response ANOVAs for each effect

# If multiple univariate tests are reported, mention multiplicity control
# (e.g., Bonferroni or Holm) when interpreting them jointly.

#### 6) Branch A: Discriminant Classification---------------------------------------------------------------
####  (if the exam asks for a classifier)

# Choose the target(s). Typical exam variants:
#   - Predict F1 from Y
#   - Predict F2 from Y
#   - Predict the joint class F12 from Y (especially if interaction is important)

# LDA (works best with MVN + equal covariances; still commonly reported)
lda_F1  <- lda(Y, F1, CV = TRUE)
lda_F2  <- lda(Y, F2, CV = TRUE)
lda_F12 <- lda(Y, F12, CV = TRUE)

err_F1  <- mean(lda_F1$class  != F1)
err_F2  <- mean(lda_F2$class  != F2)
err_F12 <- mean(lda_F12$class != F12)

err_F1; err_F2; err_F12

# If Box’s M suggested unequal covariances, also report QDA as a sensitivity check:
qda_F1  <- qda(Y, F1, CV = TRUE)
qda_F2  <- qda(Y, F2, CV = TRUE)
qda_F12 <- qda(Y, F12, CV = TRUE)

errq_F1  <- mean(qda_F1$class  != F1)
errq_F2  <- mean(qda_F2$class  != F2)
errq_F12 <- mean(qda_F12$class != F12)

errq_F1; errq_F2; errq_F12

# Exam write-up:
# - Report LOOCV error rates (actual error estimates).
# - Comment on LDA vs QDA in light of covariance homogeneity.
# - If interaction was significant, justify using F12 classifier.

#### 7) Branch B: Bonferroni CIs for mean differences---------------------------------------------------------------
#### (d) Bonferroni-adjusted CIs for factor effects on each response
####     Generic variables: Y (matrix/data.frame), F1, F2, F12

alpha <- 0.05

# Fit the additive MANOVA model (since interaction is not retained for main-effect CIs)
fit_add <- manova(as.matrix(Y) ~ F1 + F2)

# Residual SS matrix (pooled across groups) and df
W <- summary.manova(fit_add)$SS$Residuals
N <- nrow(Y); L1 <- nlevels(F1); L2 <- nlevels(F2)
df_resid <- N - L1 - L2 + 1

# Number of responses
p <- ncol(Y)

# Family size for Bonferroni:
#   - all pairwise level-contrasts for F1 and F2, across ALL responses
#   - for binary factors, choose(L,2)=1, so m = p*(1 + 1) = 2p
m <- p * (choose(L1, 2) + choose(L2, 2))
qT <- qt(1 - alpha/(2*m), df = df_resid)

# Helper: build CIs for all pairwise contrasts of a factor (marginal means)
build_CIs_for_factor <- function(Y, F, factor_name, W, df_resid, qT) {
  levs <- levels(F)
  # counts and marginal means by factor level
  n_by_lev <- sapply(levs, function(lv) sum(F == lv))
  mean_by_lev <- lapply(levs, function(lv) colMeans(Y[F == lv, , drop = FALSE]))
  names(mean_by_lev) <- levs
  
  # variance (per response) part common to all pairs: diag(W)/df_resid
  sigma2 <- diag(W) / df_resid
  
  # all pairs
  pairs <- t(combn(levs, 2))
  out <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    a <- pairs[i, 1]; b <- pairs[i, 2]
    mdiff <- mean_by_lev[[a]] - mean_by_lev[[b]]                 # length p
    se    <- sqrt(sigma2 * (1/n_by_lev[a] + 1/n_by_lev[b]))      # length p
    ci_lo <- mdiff - qT * se
    ci_hi <- mdiff + qT * se
    out[[i]] <- data.frame(
      factor     = factor_name,
      contrast   = paste0(a, " - ", b),
      response   = colnames(as.matrix(Y)),
      estimate   = as.numeric(mdiff),
      lower      = as.numeric(ci_lo),
      upper      = as.numeric(ci_hi)
    )
  }
  do.call(rbind, out)
}

# Compute CIs for both factors and bind
CI_F1 <- build_CIs_for_factor(Y, F1, "F1", W, df_resid, qT)
CI_F2 <- build_CIs_for_factor(Y, F2, "F2", W, df_resid, qT)
Bonf_CIs <- rbind(CI_F1, CI_F2)

# Add meta info
Bonf_CIs$alpha_global <- alpha
Bonf_CIs$family_size  <- m
Bonf_CIs$per_CI_level <- 1 - alpha/m

Bonf_CIs