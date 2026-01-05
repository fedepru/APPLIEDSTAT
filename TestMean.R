#________________ Applied Statistics – Template 1 ________________
####      TEST MEAN, MVN, Bonferroni CI's   ####


# setup: data loading, housekeeping and libraries -------------------------

rm(list=ls()); graphics.off()
#data= read.table(".txt", header = TRUE)

## Libraries
suppressPackageStartupMessages({
  library(MVN)
  library(car)
  library(mvtnorm)
})

# helpers ---------------------------------------------------------------
safe_solve <- function(S) {
  # try a stable inverse; fallback to Moore–Penrose if needed
  out <- try(solve(S), silent = TRUE)
  if (inherits(out, "try-error")) {
    suppressPackageStartupMessages(require(MASS))
    MASS::ginv(S)
  } else out
}

choose_large_n <- function(n, p, k = 30) {
  # Rule used in lectures: large-n if n > k * p^2
  n > k * p^2
}

# difference builder (only if needed) ------------------------------------
build_differences <- FALSE  # <- set TRUE if the task is "paired differences"

if (build_differences) {
  D <- data.frame(
    t12 = data$t1 - data$t2,
    t13 = data$t1 - data$t3
    # add other contrasts as needed
  )
  data_diff <- D
}

#framework: check if n>30*P^2: large n (no need for Gauss) 
# or small n (we need gauss assumption -> mvn test)
#then we choose the appropriate test for the mean:
#we can either compare the mean to a fixed mu0
#or take the difference between paired obs (different people measure the same thing)
#and compare it against the null vect (are they equal or not)


# MVN tests ---------------------------------------------------------------
result <- mvn(data) # By default the Henze-Zirkler's test is used (preferred)
result$multivariateNormality 

#I have small n and no Gauss => box cox transformation (never happened)





# tests on the Mean -------------------------------------------------------


# Hotelling's T2-test ----------------------------------------------------
## Classic with n<30*p^2 ----------------------------------------------------------------
#ASSUMPTIONS: multivariate normality check with mvn
#the following assumption are often satisfied by the nature of the problem itself
#Multiple labs/people measuring the same object
#for two-sample(paired) tests, equal covariance matrices across groups.
#Scaling: variables should be on comparable scales (or scale them) because 
#the test uses the sample covariance matrix.

#difference test: if we need to compare we need to build the appropriate df
D <- data.frame(
  t12 = data$t1 - data$t2,
  t13 = data$t1 - data$t3
) #this is just an example, adjust as needed/skip if it's not the case

#if it's a difference test substitute data with data_diff or run the following:
data = D #skip if not a diff test

n <- dim(data)[1]
p <- dim(data)[2]
mu0   <- c(0, 0) #adjust mu0 as needed
alpha = 0.05  #adjust alpha as needed

M   <- sapply(data, mean)
S    <- cov(data)
invS <- solve(S)
T2  <- n * (M - mu0) %*%  invS %*% (M - mu0)
cfr.fisher <- ((n - 1) * p / (n - p)) * qf(1 - alpha, p, n - p)
#we use fisher quantile for small n (with gaussianity satisfied)
T2 < cfr.fisher
# Compute the p-value
P <- 1 - pf(T2 * (n - p) / (p * (n - 1)), p, n - p)
P


## Asymptotic Test (n>30*p^2) --------------------------------------------------------
#No need for normality assumption: it's actually the same test
#we compensate the missing normality assumption with large n
#we compare the same T2 statistic against a chisq


#difference test: if we need to compare we need to build the appropriate df
D <- data.frame(
  t12 = data$t1 - data$t2,
  t13 = data$t1 - data$t3
) #this is just an example, adjust as needed/skip if it's not the case

#if it's a difference test substitute data with data_diff or run the following:
data = D #skip if not a diff test

n <- dim(data)[1]
p <- dim(data)[2]
mu0   <- c(0, 0) #adjust mu0 as needed
alpha = 0.05  #adjust alpha as needed

M   <- sapply(data, mean)
S    <- cov(data)
invS <- solve(S)
T2A  <- n * (M - mu0) %*%  invS  %*% (M - mu0) #A: Asymptotic
cfr.chisq <- qchisq(1 - alpha, p)
T2A < cfr.chisq
# Compute the p-value
PA <- 1 - pchisq(T2A, p)
PA



## Test for repeated measures -----------------------------------------------------------
# Example of context: the same measurement on the same unit repeated through time
#small n: go to mvn section to check normality
#setup: usual stuff, we suppose data have been loaded into data

n <- dim(data)[1]
p <- dim(data)[2]
mu0   <- c(0, 0) #adjust mu0 as needed
alpha = 0.05  #adjust alpha as needed

M   <- sapply(data, mean)
S    <- cov(data)
invS <- solve(S)

# we build one of the possible contrast matrices to answer
# the question: it has to be adjusted, we need the same reference point
#in this case the first measurement then make the difference with the others:
#as an ex. this should give us:t2-t1,t3-t1,t4-t1
C <- matrix(c(-1, 1, 0, 0,
              -1, 0, 1, 0,
              -1, 0, 0, 1), 3, 4, byrow=T)
C #check that is properly formulated

# Test: H0: C %*% mu == c(0, 0, 0) vs H1: C %*% mu != c(0, 0, 0)
alpha   <- 0.05
mu0 <- c(0, 0, 0) #adjust as needed

Md <- C %*% M # Sample mean of the "contrasted" observations
Sd <- C %*% S %*% t(C) # Sample covariance of the contrasted observations
Sdinv <- solve(Sd)

# Hotelling T2 statistics
T2 <- n * t(Md - mu0) %*% Sdinv %*% (Md - mu0)

# (p-1)*(n-1)/(n-(p-1)) times the 1-alpha Fisher quantile with p-1 and n-p+1 df
cfr.fisher <- ((p - 1) * (n - 1) / (n - (p - 1))) * qf(1 - alpha, (p - 1), n - (p - 1)) 

T2 < cfr.fisher # Testing if we are in the rejection region
T2
cfr.fisher
#Calculate the p-value:
P <- 1 - pf(T2 * (n - (p - 1)) / ((p- 1) * (n - 1)), (p - 1), n - (p - 1))
P


### asymptotic variant (of rep measures)------------------------------------------------------
cfr.chisq <- qchisq(1 - alpha, p)
T2 < cfr.chisq
# Compute the p-value
PA <- 1 - pchisq(T2, p)
PA






## Test for the mean of two independent Gaussian populations -----------------------------------

t1 <- #data for population 1
t2 <- #data for population 2
  
n1 <- dim(t1)[1] # n1 = 3
n2 <- dim(t2)[1] # n2 = 4
p  <- dim(t1)[2] # p = 2 must be equal to dim(t2)[2]

# Test for multivariate normality: if sample size is enough
mvn(t1)
mvn(t2)
# we compute the sample mean, covariance matrices and the matrix Spooled

t1.mean <- sapply(t1, mean)
t2.mean <- sapply(t2, mean)
t1.cov  <-  cov(t1)
t2.cov  <-  cov(t2)
Sp      <- ((n1 - 1) * t1.cov + (n2 - 1) * t2.cov) / (n1 + n2 - 2)

# We compare the matrices -> here, using rule of thumb:
# we don't reject equality of covariance matrices if s1_ii and s2_ii differ from
# less than a factor ~4 (see J-W p.291)
list(S1 = t1.cov, S2 = t2.cov, Spooled = Sp)

# Test H0: mu1 == mu2  vs  H1: mu1 != mu2
# i.e.,
# Test H0: mu1 - mu2 == c(0, 0)  vs  H1: mu1 - mu2 != c(0, 0)

alpha   <- 0.05 #adjust if needed
delta.0 <- c(0, 0) #adjust if needed
Spinv   <- solve(Sp)

T2 <- n1 * n2 / (n1 + n2) * (t1.mean - t2.mean - delta.0) %*% Spinv %*% (t1.mean - t2.mean - delta.0)

cfr.fisher <- (p * (n1 + n2 - 2) / (n1 + n2 - 1 - p)) * qf(1 - alpha, p, n1 + n2 - 1 - p)
T2 < cfr.fisher 
#calculate p-value
P <- 1 - pf(T2 / (p * (n1 + n2 - 2) / (n1 + n2 - 1 - p)), p, n1 + n2 - 1 - p)
P










# Confidence Region Plots -------------------------------------------------
#setup of parameters:adjust them as needed,
#if the test is already done they should already be in the environment
#just assign to the proper letter/change the letters in the code
#mu0=c(0,0) 
#alpha= 0.05
#M=  #mean of difference/other
#S=  #covariance matrix, either pooled(2 populations) or other use the same one from the test
#n= #number of obs
#p= #dim of R.Vector
#cfr.fisher= #if we have small n and gauss
#cfr.chisq= #if we have large n and use asymptotic test

#plot ellipse Fisher
ellipse(center= M, shape= S/n, radius= sqrt(cfr.fisher), col = 'red', lty = 1, center.pch = 4,
          center.cex=1.5, lwd=2)
#plot ellipse Chisq (Asymptotic)
ellipse(center= M,shape= S/n,radius= sqrt(cfr.chisq), col = 'orange', lty = 1, center.pch = 4,
        center.cex=1.5, lwd=2)

#Plot sample Mean and mu0
points(M[1], M[2], pch = 16, col ='orange')
points(mu0[1], mu0[2], pch = 16, col ='blue')  

#Region of rejection:same structure: just switch center= M to center= mu0
 

# Bonferroni CI's ---------------------------------------------------------
#setup of parameters:adjust them as needed,
#if the test is already done they should already be in the environment
#just assign to the proper letter/change the letters in the code
alpha=0.05 #adjust if needed: global confidence level
#mu0=c(0,0) 
#M=  #mean of difference/other
#S=  #covariance matrix, either pooled(2 populations) or other use the same one from the test
#n= #number of obs
#p= #dim of R.Vector
#We need to first set the number of CI's (usually 2=p)
k=2 #adjust if needed
cfr.t <- qt(1-alpha/(2*k), n-1) # Student quantile
#here we choose as directions (1,0) and (0,1) (most common request)
#but it can be adjusted:
a1= c(1,0)
a2= c(0,1)
#if necessary (k>2) build other vectors (they have to be std ||a||=1)

IC.BF1 <- c(a1%*%M - cfr.t*sqrt((a1%*%S%*%a1)/n),
            a1%*%M,
            a1%*%M + cfr.t*sqrt((a1%*%S%*%a1)/n))

IC.BF2 <- c(a2%*%M - cfr.t*sqrt((a2%*%S%*%a2)/n),
            a2%*%M,
            a2%*%M + cfr.t*sqrt((a2%*%S%*%a2)/n))

Bf <- rbind(IC.BF1, IC.BF2) #add other ones if needed
dimnames(Bf)[[2]] <- c('inf','center','sup')
Bf

#to set up variance Bonferroni intervals: standard directions, also variance
k=4

ICmean <- cbind(inf=M - sqrt(diag(S)/n) * qt(1 - alpha/(2*k), n-1),
                center= M,
                sup= M + sqrt(diag(S)/n) * qt(1 - alpha/(2*k), n-1))

ICvar <- cbind(inf=diag(S)*(n-1) / qchisq(1 - alpha/(2*k), n-1),
               center=diag(S),
               sup=diag(S)*(n-1) / qchisq(alpha/(2*k), n-1))




# extra -------------------------------------------------------------------








## One-sample t-test: H0: mu = mu0
# t.test(x, mu = mu0)

## Two-sample t-test (independent)
# t.test(x ~ g, var.equal = FALSE)     # Welch (default)
# var.test(x ~ g)                      # test equal variances (optional)

## Paired t-test
# t.test(before, after, paired = TRUE)

## Nonparametric alternatives
# wilcox.test(x ~ g)                   # Mann–Whitney (independent)
# wilcox.test(before, after, paired=TRUE)  # Wilcoxon signed-rank

## Proportion test (one or two)
# prop.test(x = c(x1, x2), n = c(n1, n2))  # chi-square approx
# binom.test(x1, n1)                       # exact CI for 1 prop

## k-sample ANOVA
# fit <- aov(y ~ f, data=dat)
# summary(fit)
# TukeyHSD(fit)                            # multiple comparisons

## Two-way ANOVA (with interaction)
# fit2 <- aov(y ~ A * B, data=dat)
# summary(fit2)
# interaction.plot(dat$A, dat$B, dat$y)

## Assumptions diagnostics
# shapiro.test(residuals(fit))             # normality (small n)
# plot(fit)                                # residuals vs fitted, QQ
# bartlett.test(y ~ f, data=dat)           # equal variances (normal)
# fligner.test(y ~ f, data=dat)            # robust variance test




