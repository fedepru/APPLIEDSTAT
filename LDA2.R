## PACKAGES ----------------------------------------------------------------
library(MVN)        # mvn()
library(car)        # Box's M, ellipse if you want plots later
library(heplots)    # boxM()
library(MASS)       # lda()

## DATA --------------------------------------------------------------------
# Edit the path if needed
health <- read.table("health.txt", header = TRUE)
str(health)
table(health$risk)  # counts by class

# Features and labels
X <- health[, c("hr", "activity", "sleep", "spo2")]
y <- health$risk

## (a) METHOD & ASSUMPTIONS ------------------------------------------------
# Multivariate normality within classes
mvn_yes <- mvn(X[y == "yes", , drop = FALSE], mvnTest = "mardia")
mvn_no  <- mvn(X[y == "no",  , drop = FALSE], mvnTest = "mardia")
print(mvn_yes$multivariateNormality); print(mvn_no$multivariateNormality)

# Equality of covariance matrices (LDA vs QDA)
boxM_out <- boxM(X, y)  #if p-value > 0.05 => same covariance structure => not QDA => LDA
print(boxM_out)

# Optional quick check: element-wise covariance ratios (rule-of-thumb)
S_yes <- cov(X[y == "yes", , drop = FALSE])
S_no  <- cov(X[y == "no",  , drop = FALSE])
ratios <- abs(S_yes / S_no)
cat("Max|S_yes/S_no|:", max(ratios), "   Min|S_yes/S_no|:", min(ratios), "\n")

## (b) COSTS, PRIORS, COST-ADJUSTED PRIORS --------------------------------
c_no_yes <- 100000   # cost of FN: classify "no" when true is "yes"
c_yes_no <- 500      # cost of FP: classify "yes" when true is "no"

p_yes <- 0.001       # prevalence from text (NOT sample prior)
p_no  <- 1 - p_yes

# Cost-adjusted priors (minimize expected misclassification cost)
p_adj_yes <- p_yes * c_no_yes / (p_yes * c_no_yes + p_no * c_yes_no)
p_adj_no  <- p_no  * c_yes_no / (p_yes * c_no_yes + p_no * c_yes_no)
priors_adj <- c("no" = p_adj_no, "yes" = p_adj_yes)

cat("Natural priors:      no =", p_no,  " yes =", p_yes,  "\n")
cat("Cost-adjusted priors no =", round(p_adj_no, 3), " yes =", round(p_adj_yes, 3), "\n")

## (c) CLASSIFIER (LDA), APER & LOOCV AER ----------------------------------
# Fit LDA with cost-adjusted priors (decision rule adjusted for costs)
fit_lda   <- lda(x = X, grouping = y, prior = priors_adj)
pred_appr <- predict(fit_lda)$class             # apparent (training predictions)

# LOOCV predictions for estimated actual error rate
fit_cv <- lda(x = X, grouping = y, prior = priors_adj, CV = TRUE)
pred_cv <- fit_cv$class

# Confusion matrices
cm_APP <- table(truth = y, pred = pred_appr)
cm_CV  <- table(truth = y, pred = pred_cv)
cm_APP; cm_CV

# Error-rate helpers (use TRUE priors p_no, p_yes for evaluation)
err_rate <- function(cm, p_no, p_yes) {
  n_no  <- sum(cm["no", ])
  n_yes <- sum(cm["yes", ])
  fp <- if ("yes" %in% colnames(cm)) cm["no", "yes"] else 0
  fn <- if ("no"  %in% colnames(cm)) cm["yes","no"]  else 0
  # class-conditional error rates
  fp_rate <- if (n_no  > 0) fp / n_no  else 0   # P(pred=yes | true=no)
  fn_rate <- if (n_yes > 0) fn / n_yes else 0   # P(pred=no  | true=yes)
  # overall error wrt true priors
  aer <- p_no * fp_rate + p_yes * fn_rate
  list(aer = aer, fp_rate = fp_rate, fn_rate = fn_rate)
}

APP  <- err_rate(cm_APP, p_no, p_yes)  # apparent error
CV   <- err_rate(cm_CV,  p_no, p_yes)  # LOOCV error

cat(sprintf("APER (wrt true priors): %.6f\n", APP$aer))
cat(sprintf("AER  (LOOCV, true priors): %.6f\n", CV$aer))
cat(sprintf("FP_rate (CV): %.4f   FN_rate (CV): %.4f\n", CV$fp_rate, CV$fn_rate))

# Comment (what to observe):
# - If APER << AER(CV), model is optimistic on training set.
# - With strong cost adjustment and few positives, boundary skews; LOOCV is crucial.

## (d) BUDGET FOR 1000 PATIENTS -------------------------------------------
# Predicted-positive rate from LOOCV (fraction flagged "yes")
r_hat <- mean(pred_cv == "yes")
budget_interventions <- 1000 * r_hat * 500   # $500 per flagged patient
cat(sprintf("Budget for interventions (1000 pts): $%.2f\n", budget_interventions))

## (e) EXPECTED COST: NEW vs TREAT-ALL ------------------------------------
# Expected misclassification cost per patient (from CV conditional error rates)
# EMC_per = p_no * FP_rate * c_yes_no + p_yes * FN_rate * c_no_yes
EMC_per <- p_no * CV$fp_rate * c_yes_no + p_yes * CV$fn_rate * c_no_yes

# BUT the system also spends $500 on every predicted "yes":
# Expected intervention spend per patient = r_hat * 500
cost_per_new <- EMC_per + r_hat * 500
cost_1000_new <- 1000 * cost_per_new

# Old strategy: intervene everyone (no FN loss, but 100% * $500)
cost_1000_all <- 1000 * 500

cat(sprintf("Expected total cost (new, 1000 pts): $%.2f\n", cost_1000_new))
cat(sprintf("Cost (treat-all, 1000 pts): $%.2f\n", cost_1000_all))
cat(sprintf("Savings with new approach: $%.2f\n",
            cost_1000_all - cost_1000_new))
