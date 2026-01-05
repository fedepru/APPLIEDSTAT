library(MASS)     # lda
library(MVN)      # mvn()
library(heplots)  # boxM()

## ======================= PARAMS (edit only this) =========================
file_path   <- "health.txt"                      # e.g. "food.txt"
features    <- c("hr","activity","sleep","spo2") # e.g. c("taste","odor","smooth","bright")
label_name  <- "risk"                            # e.g. "result"
pos_label   <- "yes"                             # e.g. "low"
c_fp        <- 500                               # cost: predict pos when true neg
c_fn        <- 100000                            # cost: predict neg when true pos
pi_pos_true <- 0.001                             # true prevalence from text
cohort_n    <- 1000                              # size of next cohort (number of sample on which fitted LDA will be used)
## ========================================================================

# load & split -------------------------------------------------------------
dat <- read.table(file_path, h = TRUE)
X   <- dat[, features]
y   <- factor(dat[[label_name]])
stopifnot(all(pos_label %in% y))
neg_label <- setdiff(levels(y), pos_label)

# (a) checks ---------------------------------------------------------------
mvn(X[y == pos_label, , drop = FALSE], )
mvn(X[y == neg_label, , drop = FALSE], )


boxM(X, y) #if p-value > 0.05 => same covariance structure => not QDA => LDA

S_pos <- cov(X[y == pos_label, , drop = FALSE])
S_neg <- cov(X[y == neg_label, , drop = FALSE])
max(abs(S_pos/S_neg)); min(abs(S_pos/S_neg))  # rule-of-thumb check

# (b) priors & cost-adjusted priors ---------------------------------------
p_pos <- pi_pos_true; p_neg <- 1 - p_pos
p_adj_pos <- p_pos * c_fn / (p_pos * c_fn + p_neg * c_fp)
p_adj_neg <- p_neg * c_fp / (p_pos * c_fn + p_neg * c_fp)

# (c) LDA; APER & LOOCV AER -----------------------------------------------
lda_cv <- lda(X, y, prior = c(p_adj_neg, p_adj_pos), CV = TRUE)
lda_m  <- lda(X, y, prior = c(p_adj_neg, p_adj_pos))

# confusion matrices
tab_APP <- table(class.true = y, class.assigned = predict(lda_m)$class)
tab_CV  <- table(class.true = y, class.assignedCV = lda_cv$class)
tab_APP; tab_CV

# helper for rates wrt TRUE priors
err_wrt_truepriors <- function(tab, pos = pos_label, neg = neg_label, ppos = p_pos, pneg = p_neg){
  n_neg <- sum(tab[neg, ])
  n_pos <- sum(tab[pos, ])
  fp <- if (pos %in% colnames(tab)) tab[neg, pos] else 0
  fn <- if (neg %in% colnames(tab)) tab[pos, neg] else 0
  fp_rate <- if (n_neg > 0) fp / n_neg else 0
  fn_rate <- if (n_pos > 0) fn / n_pos else 0
  aer <- pneg * fp_rate + ppos * fn_rate
  c(aer = aer, fp_rate = fp_rate, fn_rate = fn_rate)
}

APER <- err_wrt_truepriors(tab_APP)  # apparent error (wrt true priors)
AER  <- err_wrt_truepriors(tab_CV)   # LOOCV error (wrt true priors)
APER; AER

# (d) budget for interventions/tests in next cohort ------------------------
# predicted-positive rate from CV:
ptd <- mean(lda_cv$class == pos_label)  # P(pred = positive)
budget <- cohort_n * ptd * c_fp
budget

# (e) expected cost: new vs test/treat-all --------------------------------
# Expected misclassification cost per item using CV conditional error rates:
EMC_per <- p_neg * AER["fp_rate"] * c_fp + p_pos * AER["fn_rate"] * c_fn
cost_new_per   <- EMC_per + ptd * c_fp
cost_new_total <- cohort_n * cost_new_per

cost_all_total <- cohort_n * c_fp
savings <- cost_all_total - cost_new_total

cost_new_total; cost_all_total; savings
