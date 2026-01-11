# ==============================================================================
# MASTER CLASSIFICATION SCRIPT (LDA, QDA, KNN, FISHER)
# ==============================================================================

# ----------------------------- 1. LIBRARIES -----------------------------------
library(MASS)       # lda, qda
library(class)      # knn
library(mvtnorm)    # Multivariate normal distributions
library(car)        # Anova, linearHypothesis
library(heplots)    # boxM test (Homoscedasticity)
library(MVN)        # Multivariate Normality Tests

# ----------------------------- 2. THEORY RECAP --------------------------------
# A) BAYES CLASSIFIER:
#    - Minimizes Total Probability of Misclassification (if costs equal) or Expected Cost.
#    - Allocates x to group g that maximizes posterior P(g|x).
#    - Rule: Assign to group k if: f_k(x) * pi_k > f_j(x) * pi_j  (for all j != k)
#    - "Maximum Likelihood Classifier": Special case where priors pi_k are equal.

# B) LDA (Linear Discriminant Analysis):
#    - Assumptions: 
#      1. Normality (Populations are Gaussian).
#      2. Homoscedasticity (Sigma_1 = Sigma_2 = ... = Sigma_g = Sigma).
#    - Boundary: Linear. Good for small samples if assumption holds.

# C) QDA (Quadratic Discriminant Analysis):
#    - Assumptions: 
#      1. Normality.
#      2. Heteroscedasticity (Sigma_1 != Sigma_2). Allows different covariances.
#    - Boundary: Quadratic (Curved). Requires more parameters (n_g > p).

# D) FISHER DA:
#    - Non-parametric (does not assume Normality).
#    - Finds linear combination 'a' maximizing between-group var / within-group var.
#    - For g=2 groups, Fisher's direction is proportional to LDA direction.

# E) KNN (K-Nearest Neighbors):
#    - Non-parametric. No assumptions on distribution.
#    - Assigns class based on majority vote of 'k' nearest neighbors.
#    - Sensitive to scale (scale data first!).

# ==============================================================================
#  3. DATA & INPUTS ( !!! MODIFY THIS SECTION !!! )
# ==============================================================================

FILE_NAME <- "health.txt"      # e.g., "doping.txt", "health.txt"
DATA      <- read.table(FILE_NAME, header=TRUE)

# Define Variables
VAR_RESP  <- "risk"            # The grouping variable (Label)
VAR_PREDS <- c("hr", "activity", "sleep", "spo2") # Predictor columns

# Define Classes (Check your data levels!)
# usually: Positive = "yes"/"doping", Negative = "no"/"clean"
LABEL_POS <- "yes"             # The "Bad" outcome (Doping, Risk, Fail)
LABEL_NEG <- "no"              # The "Good" outcome (Clean, No Risk, Pass)

# EXAM PARAMETERS (Set NULL if not given)
# Note: "Cost of FP" is usually the cost of the test/intervention itself.
#       "Cost of FN" is the loss if we miss the bad case.
COST_FP    <- 500              # Cost if Truth=Neg, Pred=Pos (False Alarm)
COST_FN    <- 100000           # Cost if Truth=Pos, Pred=Neg (Missed Detection)
PRIOR_TRUE <- c(no=0.999, yes=0.001) # True population priors (sum to 1)
# Make sure names match data levels!
N_COHORT   <- 1000             # Number of new people to test (for Budget Q)

# ----------------------------- 4. PRE-PROCESSING ------------------------------
y <- as.factor(DATA[[VAR_RESP]])
X <- DATA[, VAR_PREDS]
g <- length(levels(y))
p <- ncol(X)

cat("\n[Data Loaded] Response:", VAR_RESP, "Levels:", levels(y), "\n")

# ----------------------------- 5. ASSUMPTIONS ---------------------------------
cat("\n================ ASSUMPTIONS CHECK ================\n")

# A) Normality (Henze-Zirkler Test)
cat("--- 1. Multivariate Normality (H0: Normal) ---\n")
# We check normality WITHIN each group
for(lev in levels(y)) {
  cat("Group", lev, ":\n")
  # If sample size is small (< 4), this might error. Skip if needed.
  tryCatch({
    print(mvn(data = X[y == lev, ], mvnTest = "hz")$multivariateNormality)
  }, error=function(e){cat(" Sample too small for MVN test\n")})
}

# B) Homoscedasticity (Box's M Test)
cat("\n--- 2. Homogeneity of Covariances (Box's M) ---\n")
# H0: Sigma_1 = Sigma_2 = ...
box_res <- boxM(X, y)
print(box_res)

choice_method <- if(box_res$p.value > 0.05) "LDA" else "QDA"
cat(sprintf("-> p-value = %.4f. Recommendation: %s\n", box_res$p.value, choice_method))

# ==============================================================================
#  6. CLASSIFIER BUILD (Cost-Adjusted)
# ==============================================================================
cat("\n================ MODEL BUILDING ================\n")

# A) Calculate Adjusted Priors (Minimizing Expected Cost)
# If no costs given, use PRIOR_TRUE. If PRIOR_TRUE not given, use Sample Priors.
if(!is.null(PRIOR_TRUE) && !is.null(COST_FN)) {
  pi_pos <- PRIOR_TRUE[LABEL_POS]
  pi_neg <- PRIOR_TRUE[LABEL_NEG]
  
  # Formula: pi_adj_k proportional to pi_k * Cost(Misclassifying k)
  # Cost(Misclassifying Pos) = COST_FN (We said it was Neg, but it was Pos)
  # Cost(Misclassifying Neg) = COST_FP (We said it was Pos, but it was Neg)
  
  term_pos <- pi_pos * COST_FN
  term_neg <- pi_neg * COST_FP
  
  pi_adj_pos <- term_pos / (term_pos + term_neg)
  pi_adj_neg <- term_neg / (term_pos + term_neg)
  
  priors_fit <- setNames(c(pi_adj_neg, pi_adj_pos), c(LABEL_NEG, LABEL_POS))
  
  cat("[Strategy] Minimizing Expected Cost.\n")
  cat("True Priors:    ", PRIOR_TRUE, "\n")
  cat("Adjusted Priors:", priors_fit, "\n")
  
} else {
  priors_fit <- PRIOR_TRUE # Or rep(1/g, g) for Max Likelihood
  if(is.null(priors_fit)) priors_fit <- table(y)/length(y) # Sample priors
  cat("[Strategy] Standard Bayes / ML (No Cost Adjustment).\n")
}

# B) Fit Model
if (choice_method == "LDA") {
  fit <- lda(X, grouping = y, prior = priors_fit)
  # For Fisher: fit$scaling gives the coefficients of the linear discriminant
} else {
  fit <- qda(X, grouping = y, prior = priors_fit)
}
print(fit)

# ==============================================================================
#  7. EVALUATION (AER & CV)
# ==============================================================================
cat("\n================ PERFORMANCE ================\n")

# A) APER (Apparent Error Rate - Training Data)
pred_train <- predict(fit)$class
conf_matrix <- table(Truth = y, Predicted = pred_train)
cat("--- Confusion Matrix (Apparent) ---\n")
print(conf_matrix)

# B) LOOCV (Leave-One-Out Cross Validation) - Actual Error Rate Est.
# Note: lda/qda have built-in CV=TRUE argument
if (choice_method == "LDA") {
  fit_cv <- lda(X, grouping = y, prior = priors_fit, CV = TRUE)
} else {
  fit_cv <- qda(X, grouping = y, prior = priors_fit, CV = TRUE)
}
class_cv <- fit_cv$class
conf_matrix_cv <- table(Truth = y, Predicted = class_cv)

cat("\n--- Confusion Matrix (LOOCV / Actual Estimate) ---\n")
print(conf_matrix_cv)

# Calculate rates manually for the exam (Sensitivity, Specificity)
# Assuming 2x2 matrix
TP <- conf_matrix_cv[LABEL_POS, LABEL_POS]
TN <- conf_matrix_cv[LABEL_NEG, LABEL_NEG]
FP <- conf_matrix_cv[LABEL_NEG, LABEL_POS] # Truth=Neg, Pred=Pos
FN <- conf_matrix_cv[LABEL_POS, LABEL_NEG] # Truth=Pos, Pred=Neg

# Rates conditioned on TRUTH (The rows)
fn_rate <- FN / (TP + FN)  # P(Pred Neg | True Pos)
fp_rate <- FP / (TN + FP)  # P(Pred Pos | True Neg)
sens    <- TP / (TP + FN)  # Recall
spec    <- TN / (TN + FP)  # Specificity

cat(sprintf("\nFP Rate (alpha): %.4f (False Alarm Rate)", fp_rate))
cat(sprintf("\nFN Rate (beta):  %.4f (Missed Detection Rate)\n", fn_rate))

# AER (Actual Error Rate)
# IMPORTANT: If Population Priors are given, weight errors by those!
# AER = pi_pos * FN_rate + pi_neg * FP_rate
if(!is.null(PRIOR_TRUE)) {
  aer_global <- PRIOR_TRUE[LABEL_POS] * fn_rate + PRIOR_TRUE[LABEL_NEG] * fp_rate
  cat(sprintf("AER (weighted by True Priors): %.6f\n", aer_global))
} else {
  aer_global <- mean(class_cv != y)
  cat(sprintf("AER (unweighted sample mean): %.6f\n", aer_global))
}

# ==============================================================================
#  8. BUDGET & SAVINGS (Exam Specific)
# ==============================================================================
cat("\n================ ECONOMICS / BUDGET ================\n")

if(!is.null(N_COHORT) && !is.null(COST_FN)) {
  
  # 1. New Strategy Cost
  # Cost = (Cost_Test * N_tested) + (Cost_FN * N_missed_cases)
  # Here "Test" implies the intervention/2nd-stage-test triggered by a Positive Prediction.
  # If the "1st stage" (this classifier) is free/cheap, we only pay COST_FP for Predicted Positives.
  
  # Probability of Predicting Positive (P_hat_Pos)
  # P(Pred Pos) = P(Pred Pos | True Pos) * P(True Pos) + P(Pred Pos | True Neg) * P(True Neg)
  #             = Sensitivity * pi_pos + FP_rate * pi_neg
  
  p_pred_pos <- (sens * PRIOR_TRUE[LABEL_POS]) + (fp_rate * PRIOR_TRUE[LABEL_NEG])
  
  # Expected Number of Interventions (Tests)
  n_interventions <- N_COHORT * p_pred_pos
  cost_interventions <- n_interventions * COST_FP
  
  # Expected Number of Missed Cases (False Negatives)
  # P(Miss) = FN_rate * pi_pos
  p_miss <- fn_rate * PRIOR_TRUE[LABEL_POS]
  n_missed <- N_COHORT * p_miss
  cost_missed <- n_missed * COST_FN
  
  total_cost_new <- cost_interventions + cost_missed
  
  cat(sprintf("Cohort Size: %d\n", N_COHORT))
  cat(sprintf("Exp. Interventions (Pred Pos): %.1f -> Cost: $%.2f\n", n_interventions, cost_interventions))
  cat(sprintf("Exp. Missed Cases  (FN):       %.4f -> Cost: $%.2f\n", n_missed, cost_missed))
  cat(sprintf("TOTAL COST (New Strategy):     $%.2f\n", total_cost_new))
  
  # 2. Old Strategy (Treat Everyone / Test Everyone)
  # If we test everyone, we pay COST_FP for N_COHORT.
  # We assume the test is perfect, so 0 Missed Cases cost.
  total_cost_old <- N_COHORT * COST_FP
  
  cat(sprintf("TOTAL COST (Old Strategy - Test All): $%.2f\n", total_cost_old))
  cat(sprintf("SAVINGS: $%.2f\n", total_cost_old - total_cost_new))
  
} else {
  cat("Skipping budget calculation (Inputs missing).\n")
}

# ==============================================================================
#  9. METRICS EVALUATION (Bonus / If asked)
# ==============================================================================
# Calculate F1, Precision, etc. based on LOOCV
precision <- TP / (TP + FP)
recall    <- sens
f1_score  <- 2 * (precision * recall) / (precision + recall)

cat("\n--- Additional Metrics ---\n")
cat(sprintf("Precision (PPV): %.4f\n", precision))
cat(sprintf("Recall (Sens):   %.4f\n", recall))
cat(sprintf("Specificity:     %.4f\n", spec))
cat(sprintf("F1-Measure:      %.4f\n", f1_score))

# ==============================================================================
#  OPTIONAL: KNN SECTION (If asked)
# ==============================================================================
# cat("\n--- KNN (k=5) ---\n")
# X_scaled <- scale(X) # KNN requires scaling
# knn_pred <- knn.cv(train = X_scaled, cl = y, k = 5)
# print(table(Truth=y, Pred=knn_pred))