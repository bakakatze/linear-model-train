## Tutorial for Linear Model and its extension by Faraway (2006)

library(dplyr)
library(tidyr)
library(faraway)
library(MASS)
library(splines)
library(ggplot2)
library(survival)

#
##### CH1: Introduction & data overview #####
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


## plots for undercounts
boxplot(gavote$undercount) # 2 extreme outliers

gavote %>%
  filter(undercount > 0.12)


## summary of counties' sizes
summary(gavote$votes) # varies greatly


#
#### fitting a linear model (Least Squared Estimate) ####
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

#### Robust Regression ####

# least squares method works well when the errors are normal, but is poor
# if the errors are long-tailed.
# Alternative: downweights the effect of larger errors (Huber Method (default))

require(MASS)

rlmodi = rlm(undercount ~ cperAA + cpergore*rural + equip, gavote)
summary(rlmodi)

# not equipOS-PC is about half the size of the "lmodi" model

#### Weighted Least Squares ####
# because the number of population across different counties varies greatly and
# we expect smaller counties to have higher variance, we can assign weights based on ballots size (pop size)

wlmodi = lm(undercount ~ cperAA + cpergore*rural + equip, data = gavote, weights = ballots)
summary(wlmodi)
# the residual standard error becomes too large
# may not be a good idea to use WLS

#### Transformation ####

# different types of transformation (i.e. box-cox)
# but fancy transformation makes interpretation difficult
# However, transforming predictors variable is less problematic and we'll see here:

# transform the proportion of african american using orthogonal polynomial transformation
plmodi = lm(undercount ~ poly(cperAA,4) + cpergore*rural + equip, gavote)
summary(plmodi)

# with standard polynomials, elimination of one term would cause a change in the values
# of the remaining coefficients.
# However, this is not the case with orthogonal polynomials. The lower order term did not change
# even if we remove the higher order term.

termplot(plmodi, partial = TRUE, terms = 1)

termplot(lmodi, partial = TRUE, terms = 1)

# these show that quadratic polynomial is not so different from constant fit
# explaining the lack of significance

# next we can try piecewise polynomials (splines)
# this is more stable compared with orthogonal polynomials
require(splines)

# use cubic B-splines:
blmodi = lm(undercount ~ cperAA + bs(cpergore,4) + rural + equip, gavote)

termplot(blmodi, partial = TRUE, terms = 2)
# again this is just an example using 4 degrees of freedom, in fact it's not that different
# from a constant fit.


#### Variables selection: AIC ####

# let's build a big linear model:
# (all main effects + all two-way between qualitative + all two-way between a quali and a quanti)
biglm = lm(undercount ~ (equip+econ+rural+atlanta)^2 +
             (equip+econ+rural+atlanta) * (perAA+pergore), gavote)

# Then, using the "step" command. This will implement a stepwise search strategy
# to minimise the AIC

smallm = step(biglm, trace = F)

# automated approach does not always give the best answer
# we know that the proportion for Gore is correlated with the undercount
# but it is eliminated from the model.
# If the idea is to get an explanation, one can use a more manual approach
# that takes into account background information


#### Variables selection: testing based approach ####
# the F-test approach is considered more inferior than using the criterion-based method
# but it's applicability across a wide class of models is nice.

# let's compare the model from AIC-method above with the model with one fewer term:
drop1(smallm, test = "F")

# we can remove rural:perAA
finalm = lm(undercount ~ equip + econ + perAA +
              equip:econ + equip:perAA, gavote)
summary(finalm)
# eek... the interactions make it difficult to interpret

#### Interpreting Interactions ####

# it's often helpful to construct predictions for all levels of variables involved
# so, here we generate all combinations of equip and con for a median proportion of perAA:
pdf = data.frame(econ = rep(levels(gavote$econ), 5),
                 equip = rep(levels(gavote$equip), rep(3,5)), perAA = 0.233)

# we now compute predicted undercount for all 15 combi
pp = predict(finalm, new=pdf)

xtabs(round(pp,3) ~ econ + equip, pdf)
# poorer counties had higher undercount
# paper prediction in rich counties is nonsense, because there is no such counties


# we can do the same approach for perAA:equip
pdf = data.frame(econ=rep("middle",15),
           equip=rep(levels(gavote$equip),rep(3,5)),
           perAA=rep(c(.11,0.23,0.35),5))

pp = predict(finalm, new=pdf)

propAA <- gl(3,1,15,labels=c("low","medium", "high"))
xtabs(round(pp,3) ~ propAA + equip, pdf)
# geez... this is so messy


#### CH1 Exercise ####


data(swiss)

head(swiss)

lfert = lm(Fertility ~ Agriculture + Examination + Education + Catholic + Infant.Mortality,
           data = swiss)
summary(lfert)
plot(lfert)


drop1(lfert, test = "F") # the examination variable look horrible

termplot(lfert, partial = TRUE, terms=4) # lol... catholics...

# drop Examination
lfert1 = lm(Fertility ~ Agriculture + Education + Infant.Mortality + Catholic, data = swiss)
summary(lfert1)

rfert = rlm(Fertility ~ Agriculture + Education + Infant.Mortality + Catholic, data = swiss)
summary(rfert)

#
##### CH2: Binomial Dat ####

data(orings)
head(orings)

plot(damage/6 ~ temp, orings)

# we want to predict the probability of failure in a given O-ring
# in relation to launch temperature and predict the outcome probability when
# T = 32F

# a naive way is to use linear model (LOL):
lmod = lm(damage/6 ~ temp, orings)
abline(lmod)

# not a good idea


#### binomial logistic regression (applying link function) ####

orings$not_damaged = 6-orings$damage

logitmod = glm(cbind(damage, not_damaged) ~ temp, family = binomial, orings)
summary(logitmod)

plot(damage/6 ~ temp, orings, xlim=c(25,85), ylim = c(0,1))
x = seq(25,85,1)
lines(x, ilogit(11.6630-0.2162*x))

#
#### compare this with Probit link #####

probitmod = glm(cbind(damage,not_damaged) ~ temp, family = binomial(link=probit), orings)
summary(probitmod)

lines(x, pnorm(5.5915 - 0.1058*x), lty = 2)

# both are quire similar
# we can predict the response at 31 F:
ilogit(11.6630 - 0.2162*31)
pnorm(5.5915-0.1058*31)

# the probability of failure at 31 F is too high


#### but we need to test for deviance first #####

# recall
summary(logitmod)

pchisq(deviance(logitmod), df.residual(logitmod), lower = FALSE)
# >0.05 the more saturated model is better

pchisq(38.898, 22, lower = FALSE)
# the null model not so much

## we can also compare the two model:
## this is a better approach in testing which model is better over sigle deviance test
pchisq(38.9 - 16.9, 1, lower = FALSE)
# which is very small = thus we can conclude that launch temperature is statistically significant

# deviance-based test is preferred over z-test
# with sparse data, the SE can be overestimated and the z-test becomes too small
# This is known as Hauck-Donner effect


## the confidence interval can be constructed using normal approximations:
confint(logitmod) 
# WARNING!!: need to load MASS package first, otherwise the default method for
# ordinary linear models will be used (which is not quite right)

# you can also do the same thing for probit model
# they both are strikingly similar model (no need to go deeper at this stage)

#### odds ratio ####

data(babyfood)
xtabs(disease/(disease+nondisease) ~ sex + food, babyfood)
# this xtab show the difference in the prob of respiratory conditions
# between: breasfed, bottle-fed, and breast with supplements

## let's fit the model
mdl = glm(cbind(disease, nondisease) ~ sex + food, family = binomial, babyfood)
summary(mdl)
# as we can see from the residual deviance (small), no evidence of interaction

drop1(mdl, test = "Chi")

# let's compute the OR and the confidence interval
exp(coef(mdl))
exp(confint(mdl))

#
#### link functions ####
# usuall the choice is made based on physical knowledge or simple convenience
# let's contrusct the three link functions:

data(bliss)
bliss # insect dying at different concentration of insectiside

modl = glm(cbind(dead, alive) ~ conc, family = binomial, data = bliss)
modp = glm(cbind(dead, alive) ~ conc, family = binomial(link = probit), data = bliss)
modc = glm(cbind(dead, alive) ~ conc, family = binomial(link = cloglog), data = bliss)

fitted(modl) # this is the inverse logit of the linear predictors, just so that we got the probability back

# let's compare the logit, probit and complementary log-log fits:
cbind(fitted(modl), fitted(modp), fitted(modc))

# not much different, but let's see at a wider range
x = seq(-2,8,0.2)

pl = ilogit(modl$coef[1] + modl$coef[2]*x)
pp = pnorm(modp$coef[1] + modp$coef[2]*x)
pc = 1-exp(-exp(modc$coef[1] + modc$coef[2]*x))

plot(x, pl, type = "l", ylab = "Probability", xlab = "Dose")
lines(x,pp,lty=2)
lines(x,pc,lty=4)

# let's look at the relative different
# black lines = lower tail ratio
# red dashed lines = upper tail ratio
matplot(x,cbind(pp/pl,(1-pp)/(1-pl)),type="l",xlab="Dose",ylab="Ratio")

matplot(x, cbind(pc/pl, (1-pc)/(1-pl)), type = "l", xlab = "Dose", ylab = "Ratio")

# the problem is that we don't have a tool that can reliably predict
# at the lower and upper tail of the probability, regardless of the link function
# For example, it is easy to predict asbestosis where the mine workers are exposed to high concentration
# however, if we want to predict asbetosis in low-level exposure
# it is impossible to accurately predict the risk.

# the default choice is the logit link, because:
# simpler mathematics, easier to interpret and easier analysis of retrospectively sampled data


#### Estimation Problem ####
data(hormone) # urinary adrosterone & etiocholanolone in 26 healthy males

plot(estrogen ~ androgen, data = hormone, pch = as.character(orientation))
# there seem to be a convergence fail in the middle :O


# anyway, let's try to predict orientation based on two urinary hormones
modl = glm(orientation ~ estrogen + androgen, hormone, family = binomial)

summary(modl)
# the residual deviance is very small = extremely good fit
# but none of the predictors are significant due to high SE

#### goodness of fit ####

# using Pearson's chi-squared statistics:
modl = glm(cbind(dead,alive) ~ conc, family=binomial, data = bliss)

sum(residuals(modl, type = "pearson")^2)
deviance(modl)
# as shown here, there is little difference between pearson and deviance statistics
# although, we need to be careful because the model is fit to minimise the deviance and not the Pearson's X2

# the proportion of variane explained R^2 is a popular measuer of fit for normal linear models.
# we can apply the same concept to binomial regression by using the proportion of deviance explained
# however, a better statsitics is due to Naglekerke (1991):

(1-exp((modl$dev-modl$null)/150))/(1-exp(-modl$null/150))
# n=150, as there are 5 covariate class with 30 observation each


#### prediction and effective doses ####

modl = glm(cbind(dead,alive) ~ conc, family=binomial, data = bliss)
  
# let's predict the probability of death at dose = 2.5 + it's CI
pred.modl = predict(modl, newdata = data.frame(conc=2.5), se=T)

ilogit(c(pred.modl$fit - 1.96*pred.modl$se.fit, pred.modl$fit + 1.96*pred.modl$se.fit))
# confidence interval

## calculate LD50
ld50 = -modl$coef[1] / modl$coef[2]

# compute the SE using delta method
dr = c(-1/modl$coef[2],modl$coef[1]/modl$coef[2]^2)
se = sqrt(dr %*% summary(modl)$cov.un %*% dr)

# so the 95%CI are:
c(2 - 1.96*se,2 + 1.96*se)

## gosh... need to learn about Taylor Series, jacobian matrix and matrix multiplication

# use MASS package, for simplicity:
dose.p(modl, p=c(0.5,0.9))
# using the same principle as the manual approach


#### overdispersion ####

# if we do not include the correct term, variable or structure then
# binomial GLM may be overdispersed

data(troutegg)

bmod = glm(cbind(survive, total-survive) ~ location + period, family = binomial, troutegg)
summary(bmod)

# check for outliers
halfnorm(residuals(bmod))

# check for interaction
elogits = log((troutegg$survive+0.5)/(troutegg$total - troutegg$survive+0.5))
with(troutegg, interaction.plot(period, location, elogits))
# interaction plot is always hard to interpret, but from the graph
# we can conclude that there isn't any significant interaction at play

# now we have excluded the possibility of outliers & interaction
# we can say that overdispersion is in play here


# we can estimate the dispersion parameter:
sigma2 = sum(residuals(bmod, type="pearson")^2 /12)
# the square root of dispersion = RSS in gaussian model (linear model)
# the dispersion parameter calculated above is larget than it would be in the
# standard binomial GLM (will be discussed later)

drop1(bmod, scale = sigma2, test = "F") # quasi-binomial family

# no goodness of fit test is possible due to free dispersion parameter.
# we can use the dispersion parameter to scale up the estimates of the SE:
summary(bmod, dispersion = sigma2)

#

#### matched case-control studies ####

# matched case-control will eliminate the effect of the matched covariates
# which means, we cannot evaluate the effect of those matched variables

# test data: x-rays with childhood acute myeloid leukemia
data(amlxray) # matched by age only

# there are only 7 people with down syndrome
# best to exclude these cases
ii = which(amlxray$down == "yes")
ramlxray = amlxray[-c(ii, ii+1),]


require(survival)
cmod = clogit(disease ~ Sex + Mray + Fray + CnRay + strata(ID), data = ramlxray)
# the strata function indicates matched design

summary(cmod) 
# the CnRay indicate x-ray received on the child (ordinal, 4 levels)
# the summary spews linear, quadratic and cubic function (only linear is significant)


# let's see what happen if we conver the CnRay into numerical values + drop insignificant predictors:
cmodr = clogit(disease ~ Fray + unclass(CnRay) + strata(ID), ramlxray)
summary(cmodr)

#

##### CH3: Count Regression ####
# page 61



