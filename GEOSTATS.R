library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat)        ## Geostatistics
library(nlme)


# load data & define coords ------------------------------------------------

data <- read.table("stats.txt", h=TRUE) #change file name

## Define the sample coordinates (transform Table to SpatialPointsDataFrame)
coordinates(data) <- c('x', 'y')  # if needed change coords

# plot of fitted variogram, estimates of range and sill -------------------

#If I have a factor I should be careful evaluating the model with/without intercept:
#the betas0 can vary and I have to report the correct ones.
v1 <- variogram("Response" ~  "Covariates", data)
v1
plot(v1)


## how to initialize vgm ---------------------------------------------------
#ideally the initialization call for vgm is contained in the exercise text
#vgm(psill=  ,"Model", range= , nugget=  )
#if the model is asked "without nugget" just set nugget param to zero

#possible choices for Model: "Sph" Spherical, "Exp" Exponential, "Gau" Gaussian
# "Mat", "Lin"

#In the unfortunate case the call for vgm isn't included start with all blank param
#and try to estimate them after first iteration: nugget is a special case:
#it will tell you with/without: with nugget => set nugget to zero
#  without nugget  don't set nugget, omit it in the vgm call ex: vgm(3,"Sph", 2000)

# Fit the variogram model --------------------------------------------------
v1.fit <- fit.variogram(v1, vgm("fill with intialisation data"), fit.method = 7) # 7 - by default (Weighted least squares) 
v1.fit # here we see psill and range (estimated)
#then repeat the call with estimated param in vgm
#ex:
#model    psill    range
#1   Nug 1.451283    0.000
#2   Sph 1.550279 2001.374  psill=1.45+1.55=3 range=2000
#=> vgm(3,"Sph", 2000, 1)  
#nugget is more difficult to estimate, 
#go kinda lower than the y of the first point on the variogram: 
#look at the plot: we look for a fit in the data


plot(v1, v1.fit)
v1.fit
# Range is in v1.fit
# sill = SUM OF psill from Nugget and Spherical/Exp part


# Estimate the parameters using GLS ---------------------------------------------------------------------


g.tr <- gstat(formula = "Response" ~ "Covariates", data = data, model = v1.fit)
summary(g.tr) 
#How to write the model correctly: look carefully at the pedices of the coefficients
#PAY A LOT OF ATTENTION TO FACTOR VARIABLES: IN ALMOST EVERY CASE
#WE WILL REMOVE INTERCEPT WHEN DEALING WITH FACTORS

#ex1: y= beta0i + beta1*covariate  i factor with levels => formula= Response ~ -1 + covariate + factor
#ex2: y= beta0ij + beta1j*covariate 2 factors i and j => formula = Reasponse ~ -1 + factor1*factor2 + covariate*factor2
#THE EXAMPLE 2 IS CRUCIAL, WHEN WE SEE A beta0,ij WE KNOW WE HAVE TO REMOVE INTERCEPT and fit interaction between factors


#To estimate parameters of the model use a series of points (all in x,y = 0,0)
#in this example we have no intercept, one factor (2 levels) and one covariate
#for these BLUE should be TRUE: y= beta0i + beta1*moisture i=0,1 (false/true)
#we make prediction for multiple points and solve the system:


s0.new=data.frame(x=0.0, y=0.0, moisture=0.0, season=FALSE) # UTM coordinates
coordinates(s0.new)=c('x','y')

predict(g.tr, s0.new, BLUE = TRUE) # var1.pred =  value to consider for system to solve

s1.new=data.frame(x=0.0, y=0.0, moisture=0.0, season=TRUE) # UTM coordinates
coordinates(s1.new)=c('x','y')

predict(g.tr, s1.new, BLUE = TRUE) # var1.pred = 

s2.new=data.frame(x=0.0, y=0.0, moisture=1.0, season=FALSE) # UTM coordinates
coordinates(s2.new)=c('x','y')
#you can't put both beta0 to zero (one level will always be true)

predict(g.tr, s2.new, BLUE = TRUE) # var1.pred = 1.12994 
# here we predict as b_0_0 + b_1 = 1.12994 =>
# b_1 = 1.12994 - 1.117702 = 0.012238




# estimate response for new data ------------------------------------------
#update with data in the same order and same var names as the original df
sc.new=data.frame(x=0.0, y=0.0, "Other vars") # UTM coordinates (adjust if needed)
coordinates(sc.new)=c('x','y')

#BLUE = FALSE if we're making a "classic" prediction: we have both data and location
#otherwise we're just making inference on average (trends)
#be careful of how the question is phrased

#good rule: given coords => BLUE = FALSE

predict(g.tr, sc.new, BLUE = ) # var1.pred = 1.213959

#example: there might be transformation needed for the data
1.213959 * 31 # Question was about total growth during July (so all days in July)
# there for I had to multiply by 31 (because our data reports DAILY growth)





# new model/model update --------------------------------------------------

#which model should be preferred?
#longer range: better
#lower sill: better
#smoother plot => better

#we need for the plot to run: v2,v3 and v2.fit, v3.fit, just assign them
#v2=variogram(...)
#v2.fit= fit.variogram(...)


plot.new()
plot(v2$dist, v2$gamma, pch = 19, main = "Empirical and Fitted Variograms",
     xlab = "Distance", ylab = "Semivariance", ylim = c(0, max(v2$gamma)*1.2))

d.seq <- seq(0, max(v2$dist), length.out = 100)
gamma1 <- variogramLine(v2.fit, dist_vector = d.seq)$gamma
lines(d.seq, gamma1, col = "blue", lwd = 2)

gamma2 <- variogramLine(v3.fit, dist_vector = d.seq)$gamma
lines(d.seq, gamma2, col = "red", lwd = 2)

#legend("bottomright", legend = c("Common urban factor", "Individual intercept"),
      # col = c("blue", "red"), lwd = 2, bty = "n")

# I would prefer model 3, because it little bit smoother, because has higher range than model 2
# (while their sill is equal)
# Model 1 worse, because have much higher sill and not smooth enough



#example: we need to repeat the steps done before:
#fit the variogram, estimate parameters
#predict many points with BLUE=TRUE and solve the system

#ex: interpretation question

# How would you describe the effect of the canopy variable (during early summer)?
# It has positive effect, вue to coefficent > 0, so canopy helps trees growth 
# also Coefficents with canopy (TRUE) has higher value, so it has better effect for
# trees growth, compared to open-airs

