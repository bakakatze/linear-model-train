
library(faraway)

data(gavote)

head(gavote)

# calculating undercount (ballot issued - total votes receiced)
# this can be due to equipment failure, voters not showing up, or invalid votes

gavote$undercount = (gavote$ballots-gavote$votes)/gavote$ballots
summary(gavote$undercount)

hist(gavote$undercount, main="Undercount", xlab ="Percent")

gavote$pergore = gavote$gore/gavote$votes

# plot the proportion of AA with the prop who voted for Gore
plot_pergore = plot(pergore ~ perAA, gavote, 
                    xlab= "prop of African American", ylab = "prop for Gore")


# checking for correlations
nix = c(3,10,11,12)
cor(gavote[,nix])


#### fitting for a linear model (LSE) ####
lmod = lm(undercount ~ pergore+perAA, gavote)
coef(lmod)

# residual sum of squared (RSS) = deviance in linear models
deviance(lmod)

# standard deviation of residual errors
# sqrt(RSS/df)
sqrt(deviance(lmod) / df.residual(lmod))

lmodsum = summary(lmod)
lmodsum$sigma # a compact command to call for the SD of the residuals

# we can compute the R squared (coefficient of determination OR
# percentage of variance explained)
lmodsum$r.squared
# 5% fit... ugh...
# you can also think of R^2 as the correlation between predicted values and
# the response squared

# R^2 cannot be used as a criterion for choosing models
# so that's where adjusted R^2 comes into picture
# adj R^2 = 1 - [ (1-R^2) (n-1) / (n-k-1) ]
# k = the numbers of predictors


# summary lmod
summary(lmod)


# contrast matrix
contr.treatment(5)


# center the percentage who support Gore and percentage of African American
# to the mean
# the mean centering allow us to interpret the effect of urban/rural
# but beware of interaction term

gavote$cpergore = gavote$pergore - mean(gavote$pergore)
gavote$cperAA = gavote$perAA - mean(gavote$perAA)

lmodi = lm(undercount ~ cperAA + cpergore*rural + equip, gavote)
summary(lmodi)


# we then can compute F statistics by comparing the two models
anova(lmod,lmodi)

# we can also compute test the significance of the model
# by dropping one variable at a time
drop1(lmodi, test = "F")

# we can also compute 95% CI
confint(lmodi)

# a collection of useful diagnostic plots
plot(lmodi)

# we can look for extreme outliers using cook's distance
gavote[cooks.distance(lmodi) > 0.1,]

# another technique to detect high leverage
halfnorm(influence(lmodi)$hat)

gavote[influence(lmodi)$hat >0.3,]
# note that only two counties use paper ballot, so they have high leverage
# high leverage does not mean the are influential in terms of cook's distance

# we can also plot for partial residual plots where it removes
# all X except the one we want to examine, in this case: cperAA

termplot(lmodi, partial = TRUE, terms=1)
# this plot shows a snapshot of the marginal relationship between cperAA
# and the response. In this case: there is a linear relationship
# no need to do transformation

#### robust regression ####

# least squares method works well when the errors are normal, but is poor
# if the errors are long-tailed

