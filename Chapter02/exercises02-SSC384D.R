


#### In Galton's footsteps

heights = read.csv("heights.csv")
lmM = lm(SHGT ~ MHGT, data=heights)
lmF = lm(SHGT ~ FHGT, data=heights)

# Get the means and standard deviations
colMeans(heights)

sd(heights)
# R will do this, but will object!  The following command
# accomplishes the same thing without R yelling at you
sapply(heights,sd)

# The attach() command allows R to see the variable names in heights
# without adding the "data = heights" option to all other commands
attach(heights)

par(mfrow=c(1,2))
plot(MHGT, SHGT, pch=19)
abline(lmM$coefficients)
plot(FHGT, SHGT, pch=19)
abline(lmF$coefficients)

par(mfrow=c(1,2))
plot(MHGT, lmM$residuals, pch=19, main="MHGT")
abline(lmM$coefficients)
plot(FHGT, lmF$residuals, pch=19, main="FHGT")
abline(lmF$coefficients)




## Orchard sprays and honeybees

data(OrchardSprays)

# Calculate the group means.
# You must load the mosaic package to use the following command
library(mosaic)
mean(decrease~treatment, data=OrchardSprays)

# The above command shows you the groups, their means
# and their sample sizes
# Notice the means themselves are in the "S" column.
# Use the $-sign operator to extract them:

GroupMeans = mean(decrease~treatment, data=OrchardSprays)$S

# Make a dot plot/strip chart showing the group-wise data
stripchart(decrease~treatment,
	data=OrchardSprays, vertical=TRUE,
	xlab="Concentration of Lime Sulphur (A = highest; H = none)",
	ylab="Decrease in volume of solution")

# Add the group means in a different color and label
points(GroupMeans, pch=19, col='blue', cex=1.5)


# Fit a model and look at the differences between means
lm2 = lm(decrease~treatment, data=OrchardSprays)
coef(lm2)

# Compare these coefficients to the group means:
# How are they related?
# Hint: try adding the "intercept" to each coefficient
# and compare the result to the corresponding group mean!
rbind(GroupMeans, coef(lm2))

