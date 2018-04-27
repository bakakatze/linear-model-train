## Tutorial for Linear Model and its extension by Faraway (2006)

library(tidyverse)
library(faraway)
library(ggfortify)

library(MASS) # glm

library(splines)

library(survival)
library(factoextra) # to do correspondence analysis/cluster analysis

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


#### Exercise ####


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
##### CH2: Binomial Data ####

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

data(gala) 
# number of species of tortoise found on each Galapagos Island
# and the number of the endemics in that region (endemic/species)


# let's throw out the endemic variable first
gala = gala[,-2]

# plot linear model
modl = lm(Species ~ ., gala) #pick all columns as predictors
plot(predict(modl), residuals(modl), xlab = "Fitted", ylab = "Residuals")

# non-constant variance ==> play with box-cox method ==> choose square-root transformation
modt = lm(sqrt(Species) ~ ., gala)
plot(predict(modt), residuals(modt), xlab ="fitted", ylab="residuals")
summary(modt)
# the model is fairly good, but the existence of small response values (single digits)
# is worrying

#### Poisson model ####
modp = glm(Species ~ ., family = poisson, gala)
summary(modp)

# see if large deviance can be explained by an outlier
halfnorm(residuals(modp))

# plot the estimated variance against the mean
# for Poisson distribution, mean = variance
plot(log(fitted(modp)), log((gala$Species - fitted(modp))^2), 
     xlab = expression(hat(mu)), ylab = expression((y-hat(mu))^2))
abline(0,1)

# Poisson distribution has only one parameter and is not very flexible for empirical fitting purposes
# If we know specific mechanism, we can model the response as negative binomial
# If we don't know the specific mechanism, we can introduce a dispersion parameter
# and, the dispersion parameter can be estimated:

dp = sum(residuals(modp, type = "pearson")^2) / modp$df.residual

# then let's adjust the standard errors as follows:
summary(modp, dispersion = dp)


# when comparing Poisson models with overdispersion, an F-test rather than
# X2 test should be used
drop1(modp, test ="F")
# the z-statistics from the summary() are less reliable
# so we used the F-test instead

#### Rate models ####

# sometimes we can model counts data as binomial, for example:
# number of burglaries depends on the number of households, so >>
# each household will be assigned 0 for not burglared OR 1 for being burglared
# BUT, if the proportion is small, poisson approaximaton is effective
# and binomial cannot be used if one household got burglared twice

# sometimes we can model the ratio instead
# however often there are normality problem and unequal variance when taking this approach

# load data: effect of gamma radiation on the number of abnormal chromosomes (ca)
# the number of cells exposed differ everytime (cells)
data(dicentric)

round(xtabs(ca/cells ~ doseamt+doserate, dicentric), 2)

# let's plot the ratio
with(dicentric, interaction.plot(doseamt, doserate, ca/cells))


# we can try to model the rate directly
# the plot tells us that the dose is multiplicative, so we log this variable
lmod = lm(ca/cells ~ log(doserate)*factor(doseamt), dicentric)
summary(lmod)
# the model fits really well

# but, this is the residual plot
plot(residuals(lmod) ~ fitted(lmod), xlab="fitted", ylab="residuals")
abline(h=0)


## we may prefer to model the count response instead.
# we need to log the number of cells because we expect this to have a multiplicative effect
dicentric$dosef = factor(dicentric$doseamt)

pmod = glm(ca ~ log(cells) + log(doserate) * dosef, family = poisson, dicentric)
summary(pmod)

# we can relate this poisson model with a log link back to linear model for ratio response
# rate model (fix the coefficient as one by using an offset):
rmod = glm(ca ~ offset(log(cells)) + log(doserate) * dosef, family = poisson, dicentric)
summary(rmod)

# not much difference

autoplot(rmod, which = 1:6)

#
#### Negative Binomial ####

# load the data: five factors on the number of skips on a solder plate
data(solder)
modp = glm(skips ~ ., family=poisson, data=solder)
deviance(modp)
df.residual(modp) # gosh, bad...

# adding interaction
modp2 = glm(skips ~ (Opening + Solder + Mask + PadType + Panel)^2, family=poisson, solder)
deviance(modp2)
autoplot(modp2, which=1:6)

# at this point we may see that this is not going to help by 
# adding more complex transformation or interaction terms (difficult interpretation)
modn = glm(skips~., family=negative.binomial(1), solder)
modn

# we can play around with different values of the k, in the neg. binomial code
# but, it can be estimated using max likelihood
modn = glm.nb(skips~., solder)
summary(modn)

#
#### CH4: Contingency Tables ####

# create the data: it is about the existence of particle in wafers
# and its association with good or bad quality
y = c(320,14,80,36)
particle = gl(2,1,4, labels = c("no","yes"))
quality = gl(2,2, labels = c("good", "bad"))
wafer = data.frame(y, particle, quality)

# we need the data to be in that form for analysis
# but to observe it we prefer:
ov = xtabs(y ~ quality + particle)

## Poisson Model
# suppose we assume that the process is observed for some period of time
# and count the number of occurence of possible outcomes
# then, it is natural to view these outcomes ocurring at different rates and
# we could use poisson model
modl = glm(y ~ particle + quality, poisson, wafer)
summary(modl)

# to test the significance of individual predictors use the likelihood ratio test
# based on the differences in the deviance (does not matter in this dataset)
drop1(modl, test="Chi")

# the model coefficients are closely related to the marginal totals in the table
# the maximum likelihood satisfy the XTy = XTu
(t(model.matrix(modl)) %*% y)[,]

# if we add interaction term it will saturate the model
# so it will have 0 deviance and degrees of freedom
modls = glm(y ~ (particle+quality)^2, poisson, wafer)
summary(modls)

# 
#### Multinomial Model ####
pp = prop.table(xtabs(y ~ particle))

qp = prop.table(xtabs(y ~ quality))


# fitted values
fv = outer(qp, pp) * 450

# to test the fit, we compare this model against saturated model and calculate the deviance:
2 * sum(ov * log(ov/fv))
# note the deviance is the same as in Poisson Model
# so the test for independence in multinomial model = test for no interaction in Poisson model

# alternative to deviance is the Pearson X2 statistics:
sum ((ov-fv)^2/fv)

# to use Yate's continuity correction:
prop.test(ov) # deviance based test is preferable still

# Fisher's exact test:
fisher.test(ov) 
# this is good if you have small sample size and calculate exact confidence interval
# However, for larger tahles, X2 approximation tend to be very accurate anyway (deviance and X2 pearson test)

#### Large two-ways tables ####
data(haireye)

ct = xtabs(y ~ hair + eye, haireye)

# Pearson's X2 test for independence:
summary(ct)

# visualise them with dot plot and mosaic plot
dotchart(ct)
mosaicplot(ct, color = TRUE, main = NULL, las = 1)

# now try the Poisson model:
modc = glm(y ~ hair+eye, family = "poisson", haireye)
summary(modc)  
# based on the deviance of 146.44 we know that those two variables are dependent
# but if we want to know how dependent they are, we can use correspondence analysis


## correspondence analysis
z = xtabs(residuals(modc, type ="pearson") ~ hair + eye, haireye)

# calculating singular value decomposition
svdz = svd(z, 2, 2)

leftsv = svdz$u %*% diag(sqrt(svdz$d[1:2])) 
rightsv = svdz$v %*% diag(sqrt(svdz$d[1:2]))
ll = 1.1*max(abs(rightsv), abs(leftsv))

# plot them
plot(rbind(leftsv, rightsv), asp =1, xlim=c(-ll, ll), ylim=c(-ll,ll), xlab = "SV1", ylab = "SV2", type = "n")
abline(h=0, v = 0)
text(leftsv,dimnames(z)[[1]])
text(rightsv,dimnames(z)[[2]])


temp = corresp(ct, nf=3)
get_eigenvalue(temp)
fviz_ca_biplot(temp)


#
#### Matched Pairs ####

data(eyegrade)

ct = xtabs(y ~ right + left, eyegrade)

# check for independence
summary(ct)

# checking for symmetry
symfac = factor(apply(eyegrade[,2:3],1, function(x) paste (sort(x), collapse="-")))

mods = glm(y ~ symfac, eyegrade, family = poisson)
pchisq(deviance(mods), df.residual(mods), lower = F)

# check the residuals
round(xtabs(residuals(mods) ~ right + left, eyegrade), 3)

# calculate marginal counts
margin.table(ct, 1)
margin.table(ct, 2)
# we see that left eyes tend to be worse than the paired right eyes
# so perhaps marginal homogeniety does not hold here
# The assumption of symmetry implies marginal homogeniety (but the reverse is not necessarily true)

# we may observe data with different rows and columns frequencies, but still interested in symmetry:
# so we can fit a model (quasi-symmetry): log(np_ij) = log n + log a_i + log b_j + log gamma_ij

modq = glm(y ~ right + left + symfac, eyegrade, family = poisson)
pchisq(deviance(modq), df.residual(modq), lower = F) 
# we see this model fit. It can be shown that marginal homogeniety together with quasi-symmetry implies symmetry.

# we can test for marginal homogeniety by comparing the symmetry and quasi-symmetry models:
anova(mods, modq, test = "Chi") # so we find evidence of a lack of marginal homogeniety. This test is only appropriate if quasi-symmetry already holds


# when we examine the data, many people do have symmetric vision.
# we may ask whether there is independence between right and left vision among those who do not have symmetric vision.

# so, we test for quasi-independence hypothesis:
modqi = glm(y ~ right + left, eyegrade, family = poisson, subset = -c(1,6,11,16))

pchisq(deviance(modqi), df.residual(modqi), lower = F) # does not fit


#
#### Three-Way Contingency Tables ####

# Data: 20-year follow up study on the effects of smoking
# variable: alive/dead, smoker/non-smoker, age at 20-year mark
data(femsmoke)

# ignore the age-group first and see
ct = xtabs(y ~ smoker + dead, femsmoke)

prop.table(ct,1) # we see more smoker survive 20-year mark

# we can test for this significance
summary(ct)

# it is significant. But let's see if we only look at 55-64
cta = xtabs(y ~ smoker + dead, femsmoke, subset = (age == "55-64"))
prop.table(cta,1) # now we got more deaths among smoker

# this is an example of Simpson's Paradox
# let's see what happened here:
prop.table(xtabs(y ~ smoker + age, femsmoke), 2) # smokers are younger

# now let's incorporate everything
ct3 = xtabs(y ~ smoker + dead + age, femsmoke)
apply(ct3, 3, function(x) (x[1,1]*x[2,2])/(x[1,2]*x[2,1]))

# we can compute Cochran-Mantel-Haenszel test (2x2 tables across k strata)
mantelhaen.test(ct3, exact = TRUE)

# we can test for independence using linear model approach, using Pearson's X2 test:
summary(ct3)

# we can also fit the appropriate linear model:
modi = glm(y ~ smoker + dead + age, femsmoke, family = poisson)
c(deviance(modi), df.residual(modi))

# we can show the coefficient of this model correspond to the marginal proportion:
coefsmoke = exp(c(0, coef(modi)[2]))
coefsmoke/sum(coefsmoke)

# we see that these are just marginal proportions for smokers and nonsmokers in the data
prop.table(xtabs(y ~ smoker, femsmoke))
# the main effects of the model just convey information that we already know and is not the main interest of the study


## Join independence:
# let P_ij be the marginal probability that the observation falls into a (i,j, .) cells.
# now suppose that the first and second variable are dependent, but jointly independent of the third
# then: P_ijk = P_ij * P_k

modj = glm(y ~ smoker * dead + age, femsmoke, family = poisson)
c(deviance(modj), df.residual(modj)) # although there is an improvement compared with the mutual independence model, the deviance is still very high


## Conditional independence:
# the nature of the conditional independence can be determined by observing which of one of the three possible two-way interactions does not appear in the model
# the most plausible conditional independence model for our data is: [P_ijk = P_ik * P_jk | P_k]
modc = glm(y ~ smoker*age + age*dead, femsmoke, family = poisson)
c(deviance(modc), df.residual(modc)) # we see that the deviance is only slightly larger than the df, which indicate a fairly good fit

# however, we have smoe zeroes and other small numbers. So, there is some doubt of the accuracy of the X2 approximation here.
# it is better to compare model rather than assess the goodness of fit.

## Uniform association:
# we may consider a model with all three-way interactions
modu = glm(y ~ (smoker + age + dead)^2, femsmoke, family = poisson)

# now we compute the fitted values and determine the oods ratios for each age group based on the fitted values:
ctf = xtabs(fitted(modu) ~ smoker + dead + age, femsmoke)
apply(ctf, 3, function(x) (x[1,1]*x[2,2])/(x[1,2]*x[2,1]))

# this will be precisely the coefficient for the smoking and life-status term.
exp(coef(modu)['smokerno:deadno'])


## Model Selection
# log-linear models are hierarchical, so it makes sense to start with the most complex and see how far it can be reduced
modsat = glm(y ~ smoker * age * dead, femsmoke, family = poisson)
drop1(modsat, test = "Chi")
# three-way interactions suck

drop1(modu, test = "Chi") # the smoker:dead is the test for conditional independence in the Cochran-Mantel-Haenszel test


## Binomial Model
# for some three-way tables, it may be reasonable to regard one variable as the response.

# let's construct a binomial response
ybin = matrix(femsmoke$y, ncol = 2)
modbin = glm(ybin ~ smoker * age, femsmoke[1:14,], family = binomial)

drop1(modbin, test = "Chi") # can drop the interaction term

modbinr = glm(ybin ~ smoker+age, femsmoke[1:14,], family = binomial)
drop1(modbinr, test = "Chi")
# we see that both terms are significant, further simplification is impossible

# this model is equivalent to the uniform association model above:
deviance(modu)
deviance(modbinr)

# we can extract the same odds ratio:
exp(-coef(modbinr)[2]) # the sign is only for changing the reference (in Binomial model)

## NOTE:
# we prefer Binomial GLM where one factor can be identified as the response
# we prefer Poisson GLM when the relationship between variables is more symmetric

# However, the null model of Binomial GLM = two-way interaction of Poisson GLM
modbinnull = glm(ybin ~ 1, femsmoke[1:14,], family = binomial)
deviance(modbinnull)

modj = glm(y ~ smoker * age + dead, femsmoke, family = poisson)
deviance(modj)

# So, the binomial model implicitly assumes an association between smoker and age

#### Ordinal Variables ####
data(nes96)

# let's see just the party affiliation and level of education
xtabs(~ PID + educ, nes96) # both variables are ordinal here

# let's transform the data into data frame
partyed = as.data.frame(xtabs(~ PID + educ, nes96))

# fit nominal - by - nominal model (ignoring the order)
nomod = glm(Freq ~ PID + educ, partyed, family = poisson)
pchisq(deviance(nomod), df.residual((nomod), lower = F)) # no evidence against independence

# let's restructure the data so that it includes the order
partyed$oPID = unclass(partyed$PID)
partyed$oeduc = unclass(partyed$educ)

ormod = glm(Freq ~ PID + educ + I(oPID*oeduc), partyed, family = poisson)
anova(nomod, ormod, test = "Chi") # using ordinal information gives us more power to detect association
# this chisquare model base test is robust in trying to find association

# Let's examine the gamma in our latent continuous variables (i.e. I(oPID*oeduc))
summary(ormod)$coef['I(oPID * oeduc)',]
# the gamma here is positive (0.0287) which means the higher the education the greater the probability of tending to the Republican end of spectrum


## to check the robustness of the scores, we can assign different scores and test. 
# For example we assign score that can distinguish democrat and republic
# and we'd say people who complete High school or less are not different
apid = c(1,2,5,6,7,10,11)
aedu = c(1,1,1,2,2,3,3)

ormoda = glm(Freq ~ PID + educ + I(apid[oPID]*aedu[oeduc]), partyed, family = poisson)
anova(nomod, ormoda, test = "Chi") # the numerical outcome is different but the association is still significant

## for evenly spaced scores, the log-odds ratio of the 2x2 cells subset are all equal = gamma
round(xtabs(predict(ormod, type = 'response') ~ PID + educ, partyed), 2)
# now calculate log odds ratio of the lower right 2x2 cells
log(39.28*28.85/(47.49*23.19))

# it is worth examining the residuals to check if there is more struture than the model suggests
round(xtabs(residuals(ormod, type = 'response') ~ PID + educ, partyed),2)
# we see there are more weakRep with college degree than expected and with less Master degree than expected.
# this indicate the relationship is not monotone

# we can investigate this effect by considering an odrinal-by-nominal model where we treat education as nomial (i.e. column effect)
cmod = glm(Freq ~ PID + educ + educ:oPID, partyed, family = poisson)
anova(nomod, cmod, test = "Chi") # column model seems to do better

summary(cmod)
# if the relationship is monotone, than the coefficient in the interaction terms would show so
# but this is clearly not the case

# But, let's compare the linear-by-linear association:
anova(ormod, cmod, test = "Chi") # we see the simpler linear-by-libear association is preferred than the more complex model

# we see the coefficient for high school ~ Master degree is no different
# this may suggest we need a different scoring
aedu = c(1,1,2,2,2,2,2)
ormodb = glm(Freq ~ PID + educ + I(oPID*aedu[oeduc]), partyed, family = poisson)
deviance(ormodb)
deviance(ormod)
# better! this means we may be able to reduce/simplify the education categories

#### CH 5: Multinomial Data ####
# page 106
