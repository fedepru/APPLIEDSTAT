library(nlmeU) ## --> for the dataset
library(nlme)
library(insight)


# setup: data loading -----------------------------------------------------


rm(list=ls()); graphics.off()
#data= read.table(".txt", header = TRUE)

#data$"factorvar"= as.factor(data$"factorvar") needed only if its numeric

#Modifying dataset: linear combination between variables => new variable
# a and b are the 2 coeff of the new variable, we can then use it directly in the model

data_new = data.frame(data, newvar = (a*data$"var1" + b*training$"var2"))


# classic LM and diagnostics ----------------------------------------------
model = lm("response" ~ "Covariates", data= )
#be careful with handling factor type variables
#if they are already in text format there is no need for as.factor()



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

# diet-specific random effects
ranef_slopes <- ranef(lme_model_slope)

# diet-specific total slopes
diet_slopes <- fixef(lme_model_slope) + ranef_slopes
diet_slopes


# varPower: structure for residuals ---------------------------------------
# Fit an extended model M2, allowing heteroscedastic residuals: 
# εij = N (0,ϑ2ij) with

# ϑij= sigma*|variableij|^delta
# for individual j → {1, . . . , ni}, in factor= i
# Estimate delta for model M2.

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


# Extract ρ (on natural scale)
rho <- coef(fit_cs$modelStruct$corStruct, unconstrained = FALSE)
rho

# Variance components as usual
VarCorr(fit_cs)


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
