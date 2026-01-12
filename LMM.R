library(nlmeU) ## --> for the dataset
library(nlme)
library(insight)
library(glmnet)

# setup: data loading -----------------------------------------------------


rm(list=ls()); graphics.off()
#data= read.table(".txt", header = TRUE)
#WARNING: AVOID AS MUCH AS POSSIBLE TO CONVERT VARIABLES TO FACTOR TYPE
#data$"factorvar"= as.factor(data$"factorvar")

#Modifying dataset: linear combination between variables => new variable
# a and b are the 2 coeff of the new variable, we can then use it directly in the model

data_new = data.frame(data, newvar = (a*data$"var1" + b*training$"var2"))


# classic LM and diagnostics ----------------------------------------------
model = lm("response" ~ "Covariates", data= )
#be careful with handling factor type variables
#if they are already in text format there is no need for as.factor()

#DIAGNOSTICS PLOTS
# plot to confirm that residuals are normal
qqnorm(model$residuals) 
qqline(resid(model), col='red', lwd=2)
plot(model, which = 1) # plot to confirm that residuals are normal
plot(model, which = 2) # plot showing homoscedasticity,
#BE CAREFUL WITH SHAPIRO TEST: if n is large even small deviations from normality
#can lead to rejection of H0 of normality
shapiro.test(model$residuals) # So we reject H0 of residuals normality

colori = rainbow(4)
boxplot(model$factor ~ storage$factor, col=colori,
        xlab='Time.f', ylab='Residuals') 


# test of whether a var has influence -------------------------------------
#ex:
# p-value of test that Social_connections has negative effect:
#in the model summary
alpha= 0.01 #adjust as needed
confint(lm_model, "Social_connections", level=1-alpha)
#if 0 is in the interval 
# => no statistical evidence to affirm that it has influence at level 1-alpha
#confint upper bound <0 => negative effect
#confint lower bound >0 => positive effect


# LASSO REGRESSION --------------------------------------------------------
# 1. Define Response Variable
y <- data$"Response"

# 2. Define Predictor Matrix
# model.matrix automatically converts categorical vars (like Gender) to Dummy variables (0/1).
# [,-1] removes the Intercept column automatically created by model.matrix, 
# because glmnet adds its own intercept.
x <- model.matrix("Response" ~ "Covariates", data = data)[, -1]


# PART (B): LASSO REGRESSION (FIXED LAMBDA = find it in exam)
# "Perform variable selection... setting lambda Report significant coefficients."
# alpha = 1 -> Lasso (L1 penalty)
# alpha = 0 -> Ridge (L2 penalty)
lasso_fixed <- glmnet(x, y, alpha = 1, lambda = 0.3)

cat("\n--- (b) Lasso Coefficients (lambda = 0.3) ---\n")
beta_fixed <- coef(lasso_fixed)
print(beta_fixed)

# Extract only non-zero (significant) coefficients
# Lasso performs selection by shrinking irrelevant coefficients exactly to zero.
sig_coeffs <- beta_fixed[beta_fixed[,1] != 0, ]
cat("\nSignificant (Non-zero) Variables:\n")
print(names(sig_coeffs))

# PART (C): CROSS-VALIDATION FOR LAMBDA
# "Optimize lambda... interval [0.01, 10]. Random state fixed."

# 1. Set Seed as requested
set.seed(20231108)

# 2. Define the lambda grid interval [0.01, 10]
lambda_grid <- seq(0.01, 10, length = 1000)

# 3. Perform Cross-Validation (default is 10-fold)
cv_lasso <- cv.glmnet(x, y, alpha = 1, lambda = lambda_grid)

# 4. Identify Optimal Lambda (lambda.min)
best_lambda <- cv_lasso$lambda.min
cat(sprintf("\n--- (c) Optimization Results ---\n"))
cat(sprintf("Optimal Lambda: %.4f\n", best_lambda))

# 5. Report Associated Mean Cross Validation MSE
# cv_lasso$cvm contains the Mean Squared Errors for every lambda in the grid.
# We find the index of the best lambda and extract the corresponding MSE.
min_mse_index <- which(cv_lasso$lambda == best_lambda)
min_cv_mse    <- cv_lasso$cvm[min_mse_index]
cat(sprintf("Minimum CV MSE: %.4f\n", min_cv_mse))

# 6. Report Coefficients for the Optimal Lambda
cat("\nCoefficients at Optimal Lambda:\n")
best_coefs <- coef(cv_lasso, s = "lambda.min")
print(best_coefs)

# number of parameters estimated in the model -----------------------------
#use summary and count the coefficients, intercept counts as one
#if random intercept is used we add 2 more: 
# Random intercept - sigma_b
# Residual - eps   (RESIDUAL STD ERROR)


# expected response based on group (RI) ----------------------------------------
#random intercept needed: we want to find the group with  highest/lowest response
ranef(lme_model) #then choose highest or lowest as asked



# random intercept & PVRE --------------------------------------------------------

lme_model <- lme("response Var" ~ "Covariates", random = ~1|as.factor("factor"), data = )

summary(lme_model)

var_b <- get_variance_random(lme_model)
var_eps <- get_variance_residual(lme_model)
var_b #this is not sigma_b but sigma_b^2 !!!
sqrt(var_b) #sigma_b
var_eps

PVRE <- var_b/(var_b+var_eps) 
PVRE 


# This creates a "Caterpillar Plot" (Dot plot with 95% Confidence Intervals)

re_plot <- plot(ranef(lme_model, level = 1)) 
print(re_plot)
#The Zero Line (Vertical): Represents the Fixed Intercept beta_0
#The Dots: Represent the estimated Random Intercept bi for each group
#They show the deviation from the average
#Right of 0: This group starts with a higher value than the average
#Left of 0: This group starts with a lower value than the average
#The Horizontal Lines (Intervals): These are the 95% confidence intervals
#Crossing 0: The group is not significantly different from the population average.
#Not Crossing 0: The group is significantly different from the average
#=>(statistically significant random effect)

# random intercept + slope ------------------------------------------------

lme_model_slope <- lme("response Var" ~ "Covariates",
                       random = ~1+"covariate that interacts with factor"|as.factor("factor"), data = )

summary(lme_model_slope)
#To answer questions about the interaction between factor and the covariate chosen
#to have random slope
random_effects <- ranef(lme_model_slope)
view(random_effects)

# fixed slope
fixef(lme_model_slope)

# factor-specific random effects
ranef_slopes <- ranef(lme_model_slope)

# factor-specific total slopes
diet_slopes <- fixef(lme_model_slope) + ranef_slopes
diet_slopes


# varPower: structure for residuals ---------------------------------------
# Fit an extended model M2, allowing heteroscedastic residuals: 
# εij = N (0,ϑ2ij) with

# ϑij= sigma*|variableij|^delta
# for individual j → {1, . . . , ni}, in factor= i
# Estimate delta for model M2.
#WARNING: IF THERE ARE NO RANDOM EFFECTS IN THE MODEL RUN GLS INSTEAD OF LME
#VARIABLE HAS TO BE NUMERIC IN ORDER FOR IT TO WORK

lme_model2 <- lme("response Var" ~ "Covariates",
                 weights = varPower(form = ~"variable"), data = )
#this has just the residual structure asked (no random intercept)

summary(lme_model2)
# Variance function:
# Structure: Power of variance covariate
# Formula: ~variable
# Parameter estimates:
#  power 
# 0.5026533 = delta (example)

# VarIdent: structure for residuals ---------------------------------------

# (with varIdent heteroscedscity):
#if there is no specified structure we assume varIdent: 
#different levels => different variances


lme_VarIdent <- lme("Response" ~  "Covariates",
                  weights = varIdent(form =~1|as.factor("Factor var")),
                  data= )
summary(lme_VarIdent)

# Compound symmetry:structure for residuals -------------------------------------------------------
lme_cs <- lme("Response" ~ "Covariates",
              correlation = corCompSymm(form = ~ 1 |as.factor("Factor var")),
              data = )

summary(lme_cs)
intervals(lme_cs, which = "var-cov",level=0.99)
# With the estimates of rho, sigma and delta we can estimate the var-cov matrix
print(intervals)
# Extract ρ (on natural scale)
rho <- coef(fit_cs$modelStruct$corStruct, unconstrained = FALSE)
rho

# Variance components as usual
VarCorr(fit_cs)


# correlation structure for Residuals: AR1 --------------------------------
lme_ar1 <- lme("Response" ~ "Covariates",
               correlation = corAR1(form = ~ 1 | as.factor("Factor var")),
               data = )
# Extract intervals for correlation and variance parameters
intervals(lme_ar1, which = "var-cov")
# With the estimates of rho, sigma and delta we can estimate the var-cov matrix
print(intervals)

# Specifically for Phi (rho):
cat("95% CI for rho:", round(intervals$corStruct["Phi",], 4), "\n")
#INTERPRETATION CHECK IF 0 IS CONTAINED IN THE INTERVAL:
# if yes => no statistical evidence of that correlation structure at level alpha

# comparison between models -----------------------------------------------
# Should M2 be preferred over M1? Support your answer with a test.
anova(model1, model2) 
#example:
# Model df      AIC      BIC    logLik   Test           L.Ratio p-value
# lme_model      1  6 6642.329 6675.923 -3315.165                        
# lme_model2     2  7 6573.207 6612.399 -3279.604 1 vs 2 71.12235  <.0001
# Report the value of the test statistic: 71.122 - L.Ratio 

#how to choose:
#Low p-value: the more complex model is better
#High p-value: the less complex model is better

#AIC/BIC: lower = better.
#L.Ratio: likelihood ratio statistic.
#p-value: significance of improvement.

# predicting response for new data (classic or intervals) ----------------------------------------

# d) Estimate (using M2) the response given new data
test_data = data.frame("variables (no response var)")
predict(lm_model, test_data, level = FALSE) 
# level = FALSE - for new group, TRUE - for existing group

#INTERVALS
alpha=0.05 #adjust as needed
# predict(lm_model, test_data, interval="confidence", level=1-alpha)  # conf. interval
predict(lm_model, test_data, interval="prediction", level=1-alpha) # prediction interval



var_b2 <- get_variance_random(lme_model2)
var_eps2 <- get_variance_residual(lme_model2)
sqrt(var_eps2)
