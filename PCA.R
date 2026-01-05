library(MVN)        # for mvn()
library(car)        # ellipse()
library(tidyverse)

# setup: load and inspect data --------------------------------------------
data <- read.table("stats.txt", h = TRUE)  # modify file name if needed

n <- nrow(data)
p <- ncol(data)

#before performing PCA we need to separate eventual factors
raw = data
data = data[,1:8] #in this ex: the factor is in the 9th column

# homogeneity of variances check: boxplot
boxplot(data)
# standardize
data.sd <- scale(data, center = TRUE, scale = TRUE)
data.sd <- data.frame(data.sd)
# new boxplot
boxplot(data.sd)

# perform PCA -------------------------------------------------------------
# PCA on standardized variables
pca <- princomp(data.sd, cor = FALSE)

# Explained variance & how many PCs for ≥80%
eig     <- pca$sdev^2
prop    <- eig / sum(eig)
cumprop <- cumsum(prop)
out_var <- tibble(PC = paste0("PC", seq_along(prop)),
                  PropVar = prop,
                  CumVar  = cumprop)
print(out_var, n = Inf)

# Number of PCs to reach ≥ 80%
K80 <- which(cumprop >= 0.80)[1]
cat("PCs needed for ≥80% variance:", K80, "\n")

# plot of loadings --------------------------------------------------------
# princomp loadings are named Comp.1, Comp.2, ...
loadings <- as_tibble(as.matrix(pca$loadings), rownames = "variable") |>
  pivot_longer(cols = starts_with("Comp"), names_to = "PC", values_to = "loading") |>
  filter(PC %in% c("Comp.1","Comp.2","Comp.3")) #modify if we need more/less

ggplot(loadings, aes(x = reorder(variable, loading), y = loading)) +
  geom_col() +
  facet_wrap(~ PC, ncol = 1, scales = "free_y") +
  coord_flip() +
  labs(x = NULL, y = "Loading", title = "Loadings for PC1–PC3") +
  theme_minimal(base_size = 12)




# scatterplot along PC's --------------------------------------------------
# color by a factor column
cols <- as.factor(raw$`Factor Var`) #modify this for the grouping variable (a factor)


plot(pca$scores[,1], pca$scores[,2],
     col = cols, pch = 19, xlab = "PC1", ylab = "PC2",
     main = "scatterplot of data along PC's")
legend("topright", legend = levels(cols), col = 1:length(levels(cols)), pch = 19)

# CR for mean of vector containing PC's -----------------------------------
# filter for a specific level of the grouping
#here PC1 & PC2 selected, adjust indexes for other PC's
df <- data.frame(
  pc1 = pca$scores[,1][raw$`Factor Var` == "Desired Level"],
  pc2 = pca$scores[,2][raw$`Factor Var` == "Desired Level"]
)


# check Normality
mvn(df)


# construct CR (Hotelling for the MEAN in 2D)
M <- colMeans(df)
S <- cov(df)
n <- nrow(df)
p <- 2 
alpha <- 0.05
c_thr <- ((n - 1) * p / (n - p)) * qf(1 - alpha, p, n - p)

# principal directions & semi-axes (for the mean ellipse)
ev    <- eigen(S)                # eigen of sample covariance
V     <- ev$vectors              # directions (columns)
lam   <- ev$values
lengths <- sqrt(lam * c_thr / n) # semi-axes lengths

# Ellipsoidal confidence region with confidence level 95%
plot(df, asp = 1, pch = 1)
points(M[1], M[2], pch = 4, lwd = 2, cex = 1.3)
# plotting ellipse for the MEAN: shape = S/n, radius = sqrt(c_thr)
car::ellipse(center = M, shape = S / n, radius = sqrt(c_thr), lwd = 2, add = TRUE)

arrows(M[1], M[2],
       M[1] + lengths[1]*V[1,1], M[2] + lengths[1]*V[2,1],
       col = "red", lwd = 2)
arrows(M[1], M[2],
       M[1] + lengths[2]*V[1,2], M[2] + lengths[2]*V[2,2],
       col = "darkgreen", lwd = 2)

abline(a = M[2] - (V[2,1]/V[1,1]) * M[1], b = V[2,1]/V[1,1],
       lty = 2, col = 'darkred', lwd = 2)
abline(a = M[2] - (V[2,2]/V[1,2]) * M[1], b = V[2,2]/V[1,2],
       lty = 2, col = 'red', lwd = 2)

# biplots -----------------------------------------------------------------
# Base biplot uses PC1–PC2 by default; for PC2–PC3 set choices = c(2,3)
biplot(pca, cex = 0.7)
biplot(pca, choices = c(2, 3), cex = 0.7)




#plot a new point on the biplot
# IMPORTANT: supply the new athlete with the SAME variable names and order as 'dat'
#this is an example, mirror the dataset given in the exam
new_data <- tibble(
  sprint_speed   = NA_real_,  # <-- fill with Table 1 value
  endurance      = NA_real_,
  vertical_jump  = NA_real_,
  agility        = NA_real_,
  strength       = NA_real_,
  reaction_time  = NA_real_,
  accuracy       = NA_real_,
  flexibility    = NA_real_,
  throwing_power = NA_real_
)

# Coordinates (scores) of the new point in the PCA space
new_scores <- predict(pca, newdata = new_data)
new_scores_1to3 <- new_scores[, 1:3, drop = FALSE] #to plot in PC1,PC2,PC3 adjust if needed
print(new_scores_1to3)

# Base R biplot
biplot(pca, choices = c(2, 3), cex = 0.7)

# Add new data
points(new_scores[ ,2], new_scores[ ,3], col = "blue", pch = 19, cex = 1.4)
#in this case: 2 and 3 are the coord in PC2 and PC3

#Recall what the biplot shows: example with athletes
#Points = athletes, located by their PC scores (linear combinations of original skills).
#Arrows = variables (skills), direction shows correlation with the PCs, length shows strength of association.
#To interpret:
#If a point lies in the direction of an arrow, that athlete tends to score high on that variable.
#If a point lies opposite an arrow, they tend to score low.
#If close to the origin → average performance, not strongly distinguished by that PC pair.

