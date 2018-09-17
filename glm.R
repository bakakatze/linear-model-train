## Tutorial for Linear Model and its extension by Faraway (2006)

library(tidyverse)
library(faraway)
library(ggfortify)

library(MASS) # glm
library(nnet) # multinomial
library(lme4) # multilevel model

library(splines)

library(survival)
library(factoextra) # to do correspondence analysis / cluster analysis


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
data(nes96)

sPID = nes96$PID

# simplify the levels into three
levels(sPID) = c("Democrat", "Democrat", "Independent", "Independent", "Independent", "Republican", "Republican")
summary(sPID)

inca = c(1.5,4,6,8,9.5,10.5,11.5,12.5,13.5,14.5,16,18.5,21,23.5,27.5,32.5,37.5,42.5,47.5,55,67.5,82.5,97.5,115)
nincome = inca[unclass(nes96$income)] # take mid point of the income
summary(nincome)

table(nes96$educ)


# plot them
# Democrat = solid, republican = dashed(green), independent = dotted(red)
matplot(prop.table(table(nes96$educ, sPID), 1), type ="l", xlab = "Education", ylab = "Proportion", lty = c(1,2,5))
cutinc = cut(nincome, 7)
il = c(8,26,42,58,74,90,107)

matplot(il,prop.table(table(cutinc,sPID),1),lty=c(1,2,5),type="l",ylab="Proportion",xlab="Income")

cutage = cut(nes96$age, 7)
al = c(24,34,44,54,65,75,85)

matplot(al,prop.table(table(cutage,sPID),1),lty=c(1,2,5),type="l",ylab="Proportion",xlab="Age")


# build the multinomial model
mmod = multinom(sPID ~ age + educ + nincome, nes96)

mmodi = step(mmod) 
# we see that removing education at the first step reduce the AIC
# age is removed on the second step

## we can also use the likelihood methods to derive a test to compare nested models.
mmode = multinom(sPID ~ age + nincome, nes96)
deviance(mmode) - deviance(mmod)

pchisq(16.206, mmod$edf - mmode$edf, lower = F)
# we see that education is not significant relative to the full model


# we can predict using midpoints income
predict(mmodi, data.frame(nincome = il), type = "probs") # the probability of being republican or independent increases with income


## summary of the relationship between the predictors and the response:
summary(mmodi)

# the intercept term model the probability of the party identification for an income of zero
# we can see:
cc = c(0, -1.17493, -0.95036)
exp(cc)/sum(exp(cc))

predict(mmodi, data.frame(nincome = 0), type = "probs")


# we can fit a multinomial logit model using a Poisson GLM.
# we can exploit this by declaring a factor that has a level for each multinomial observation (we call this response factor)
cm = diag(3)[unclass(sPID), ]

# let's see the firsr 4 individual
cm[1:4,] # republican, democrat 3x

y = as.numeric(t(cm))
resp.factor = gl(944,3)

# let's label them
cat.factor = gl(3,1,3*944, labels = c("D", "I", "R"))

# replicate the predictor
rnincome = rep(nincome, each = 3)

head(data.frame(y, resp.factor, cat.factor, rnincome))


## null model
nullmod = glm(y ~ resp.factor + cat.factor, family = poisson)

## income is modeled with an interaction with party affiliation
glmod = glm(y ~ resp.factor + cat.factor + cat.factor:rnincome, family = poisson)

deviance(glmod)
deviance(mmodi) # the deviance is the same as the multinomial model above

coef(glmod)[c(1,945:949)]
coef(mmodi)
# the point is that the multinomial logit can be viewed as a GLM type model


#### Hierarchical or Nested Response ####

data(cns)
cns

# because we got a huge number of normal birth
# it's better to model it based on 2 hierarchies
# Y = Binom(no cns malform, cns malform)
#                           cns malform = Multinom(anencephalus, spina bifida, other)


# let's start with the binom model
cns$CNS = cns$An+cns$Sp+cns$Other

plot(log(CNS/NoCNS) ~Water, cns, pch=as.character(Work))


## because Water and Area is confouded we cannot put these predictors in one model
# try both and compare

binmodw = glm(cbind(CNS,NoCNS) ~ Water + Work, cns, family=binomial)
binmoda = glm(cbind(CNS,NoCNS) ~ Area + Work, cns, family=binomial)

anova(binmodw, binmoda, test = "Chi") # no significant difference

halfnorm(residuals(binmodw))
# Newport is an outlier

summary(binmodw)
exp(coef(binmodw)[3]) # non manual worker have a 29% lower chance of CNS malform
exp(coef(binmodw)[2]*100) # 100 lower difference in water hardness havea 27% lower chance


## now let's move on to a multinomial model for the three malformation types
cmmod = multinom(cbind(An, Sp, Other) ~ Water + Work, cns)

nmod = step(cmmod)
nmod # which leaves us with a null final model

# fitted proportion
cc = c(0,0.28963, -0.98083)
names(cc) = c("An", "Sp", "Other")
exp(cc)/sum(exp(cc))
# so we see that water hardness and parents' profession are related
# to the probability of malformed birth but not the types of malformation

# let's see if we add all four categories in one model
multinom(cbind(NoCNS, An, Sp, Other) ~ Water + Work, cns)
# Water and WOrk variables are significant here, but we cannot distinguish
# the effect on the type of malformation easily in this model


#### Ordinal Multinomial Responses ####

# proportional odds logistic regression
pomod = polr(sPID ~age + educ + nincome, nes96)

c(deviance(pomod), pomod$edf)

# can be compared to the corresponding multinomial logit model
c(deviance(mmod), mmod$edf)

# the proportional odds model uses fewer parameters, but does not fit quite as well.
# we can then use an AIC-based variable selection method:

pomodi = step(pomod)
# the finished model only include income, just like the multinomial model

# we can also use a likelihood ratio test
deviance(pomodi) - deviance(pomod)
pchisq(11.151, pomod$edf - pomodi$edf, lower = F)
# no significant difference, which means simplification of the model is justifiable


## we can check the proportional odds assumption by computing the observed odds proportions
# with respect to income levels

pim = prop.table(table(nincome, sPID), 1)

logit(pim[,1]) - logit(pim[,2] + pim[,3])
# It is questionable whether these can be considered sufficiently constant,
# but at least there is no trend.

# now consider the interpretation of the fitted coefficients:
summary(pomodi)

# we can say the odds of moving from democrat to independent/republican
# or moving from democrat/independent to republican increase by a factor of
# exp(0.013120) = 1.0132 as income increases by one unit ($1000).

# notice the log-odds are similar to those obtained in the multinomial logit model.
# the intercepts correspond to the theta j.
# so for income $0, the predicted probability of being a democrat is:
ilogit(0.2091)

# while being an independent is:
ilogit(1.292) - ilogit(0.209)

# we can compute predicted values:
predict(pomodi, data.frame(nincome = il, row.names = il),
        type = "probs")

# notice how the probs of being democrat uniformly decreases whit income
# while the inverse is true of the probs of being republican
# but, the middle category, increases then decreases
# this type of behavior can be expected from the latent variable representation of the model

# we can illustrate the latent variable interpretation of proportional odds
# by computing the cutpoints for incomes of $0, $50,000 and $100,000:
x = seq(-4, 4, by = 0.05)
plot(x, dlogis(x), type = "l")
abline(v=c(0.209, 1.292))
abline(v=c(0.209,1.292) - 50*0.013120, lty=2)
abline(v=c(0.209,1.292) - 100*0.013120, lty=3) # i see....


# applying the ordered probit model to the nes96 data, we find:
opmod = polr(sPID ~ nincome, method = "probit")
summary(opmod)

# the deviance is similar to the logit model
# but the coefficient appear to be different. However, if we compute the same predictions:
dems = pnorm(0.128 - il * 0.008182)
demind = pnorm(0.798 - il*0.008182)
cbind(dems, demind-dems, 1-demind)
# the predicted values are very similar to those seen for the logit.


# to use the cloglog method, just put "cloglog" under method

#### CH6: Generalised Linear Models ####

# let's see how the fitting algorithm works (Iteratively Reweighted Least Squares [IRWLS])
data(bliss)
head(bliss)

modl = glm(cbind(dead,alive) ~ conc, family = binomial, bliss)
summary(modl)$coef

y = bliss$dead/30 # we use y as our initial guess, cause there is no 0 / 1
mu = y
eta = logit(mu)
z = eta + (y-mu)/(mu*(1-mu))
w = 30*mu*(1-mu)

lmod = lm(z ~ conc, weights = w, bliss)
coef(lmod)

# it is very close to the converged values from modl
# this is not uncommon, but to get more precise we need to do more iterations

for (i in 1:5) {
  eta = lmod$fitted.values
  mu = ilogit(eta)
  z = eta + (y-mu)/(mu*(1-mu))
  w = 30*mu*(1-mu)
  
  lmod = lm(z ~ bliss$conc, weights = w)
  cat(i, coef(lmod), "\n")
}

# this converged quite fast, if it is not fast there are usually some problems with the model
summary(lmod)

#
#### Hypothesis Test for GLM ####

summary(modl)
# we can do a goodness of fit test by examining the size of the residual deviance compared to its degrees of freedom
1- pchisq(deviance(modl), df.residual(modl))

# p value is large indicating no evidence of a lack of fit.

# we can also test for singificance of the linear concentration term by comparing the current model with the null
anova(modl, test = "Chi")

# we can also test a more complex model
modl2 = glm(cbind(dead,alive) ~ conc + I(conc^2), family = binomial, bliss)
anova(modl, modl2, test = "Chi")
anova(modl2, test = "Chi")

# so, there is no point in adding quadratic term

# differences in deviance is the preffered method. Due to problems noted by Hauck and Donner (1977)

#### GLM diagnostics ####

# Pearson residuals can be skewed for nonnormal responses
# let's see different type of residuals in GLM

residuals(modl) # deviance residuals (degault choices)

residuals(modl, "pearson") # pearson residuals

residuals(modl, "response") # response residuals
# which are the same as
bliss$dead/30 - fitted(modl)

residuals(modl, "working") # working residuals
# which is the same as
modl$residuals
# warning!! this is usually not needed for diagnostic purposes

residuals(lmod) # the working residuals are by-product of IRWLS fitting procedure


## Leverage
# for linear model:: y_hat = H*y
# H here is the hat matrix that projects the data onto the fitted values (y_hat)
# the leverages h_i are given by the diagonal of H and represent the potential of the point to influence the fit

# Leverages are somewhat different for GLMs.
# H (in GLMs) = W^(1/2) %*% X %*% inv(t(X) %*% W %*% X) %*% t(X) %*% W^(1/2)

# we extract the diagonal elements of H to get the leverages h_i. A large value of h_i
# indicates that the fit may be sensitive to the response at case _i. Large leverages typically mean that
# the predictors values are unusual in some way. Leverages in GLMs are no longer just a function of X, but
# they depend on the response through the weights W.

# calculate the leverages for GLMs
influence(modl)$hat

# jacknife approximation to scale the residuals
rstudent(modl)

# outliers may be detected by observing particularly large jacknife residuals

# Leverage only measures the potential to affect the fit, whereas measures of influence more directly assess 
# the effect of each case on the fit. We can examine the change in the fit from omitting a case 
# by looking at the changes in the coefficients:
influence(modl)$coef

# alternatively, we can examine the Cook statistics
cooks.distance(modl)
# this shows that the biggest change would occur by omitting the first observation.


#### Two types of Model Diagnostics ####

# one is to detect single cases or small groups of cases that do not fit the pattern of the rest of the data
# the other one is to check assumptions of the model (2 == structural or stochastic)

# for linear models, the plot of residuals against the fitted values is probably the single most valuable graphic

# for GLMs, we must decide on the appropriate scale for the fitted values.
# usually, it is better to plot the linear predictors rather than the predicted responses.

# let's revisit the Galapagos data
data(gala)
gala = gala[,-2]

modp = glm(Species ~ ., family = poisson, gala)
plot(residuals(modp) ~ predict(modp, type = "response"), xlab = expression(hat(mu)), ylab = "Deviance Residuals")

# there are just a few islands with a large predicted number of species while mose predicted response values are small
# this makes it difficult to see the relationship between the residuals and the fitted values because most of the points
# are compressed on the left

# now, let's try plotting the predictors instead (eta)
plot(residuals(modp) ~ predict(modp, type = "link"), xlab = expression(hat(eta)), ylab="Deviance Residuals")

# In GLMs, it is better to transform the predictors because it causes the least disruption to the GLM
# for this particular plot, there is no evidence of non-linearity

## For all GLMs but the Gaussian, we have a nonconstant variance function.
# however, by using the deviance residuals we have already scaled out the variance function
# so, provided the variance function is correct, we do excpect to see constant variance in the plot

# if we use the response residuals instead the deviance:
plot(residuals(modp, type = "response") ~ predict(modp, type = "link"), xlab=expression(hat(eta)),
     ylab = "Response Residual")
# we see a pattern of increasing variation consistent with Poisson

# plot of the residuals are not helpful for binary responses and small sample binomial

# investigating the nature of the relationship between predictors and the response is another 
# primary objective of diagnostic plots

plot(Species ~ Area, gala)
plot(Species ~ log(Area), gala)

# we see a curvilinear relationship between the predictor and the response
# However, the default GLM uses a log link which we need to take into account
# to allow for the choice of link function, we can plot the linearised response:

mu = predict(modp, type = "response")
z = predict(modp) + (gala$Species - mu)/mu

plot(z ~ log(Area), gala, ylab = "linearised response")

# we now see a linear relationship suggesting that no further transformation of area is necessary.
modpl = glm(Species ~ log(Area) + log(Elevation) + log(Nearest) + log(Scruz + 0.1) + log(Adjacent),
            family = poisson, gala)
c(deviance(modp), deviance(modpl))
# the log transformation produces smaller deviance


# the disadvantage of examining the raw relationship between response and the predictors is that it fails
# to take into account the effect of other predictors

## Partial residual plots are used for linear models to make allowance for the effect of the other predictors
# while focusing on the relationship of interest.
mu = predict(modpl, type ="response")
u = (gala$Species - mu)/mu + coef(modpl)[2]*log(gala$Area)
plot(u ~ log(Area), gala, ylab = "Partial Residual")
abline(0, coef(modpl)[2])

## CHecking unsual points in GLM
# use half-normal plot that compares the sorted absolute residuals and the quantiles of the half-normal distribution
# the residuals are not expected to be normally distributed
# so we are only looking at outliers
halfnorm(rstudent(modpl))

# no sign of outliers

# The half normal plot is also useful for positive-value diagnostics such as leverages and Cook statistics
gali = influence(modpl)
halfnorm(gali$hat)

# let's see the influence
halfnorm(cooks.distance(modpl))

# again, Santa Cruz island is influential
# we can examine the change in the fitted coefficients. For example,
# consider the change in the Scruz coefficient as shown:
plot(gali$coef[, 5], ylab = "Change in Scruz coef", xlab = "Case no.")


# let's see what happened if we exclude Santa Cruz from the model
modplr = glm(Species ~ log(Area) + log(Elevation) + log(Nearest) + log(Scruz+0.1) + log(Adjacent), gala,
             family = poisson, subset = -25)
cbind(coef(modpl), coef(modplr))

# Scruz variable changed sign.. another solution for the full model is to add a larger amount, say 0.5 instead of 0.1

# other than user-introduced anomaly, we find no difficulty.
# Using our earlier discovery of log transformation, some variable selection and allowing for remaining overdispersion
# our final model is:
modpla = glm(Species ~ log(Area) + log(Adjacent), family = poisson, gala)
dp = sum(residuals(modpla, type = "pearson")^2) / modpla$df.residual
summary(modpla, dispersion = dp)

#
#### Exercises CH8 ####

# fit the orings data with quasibinomial
data(orings)

orings$pdmg = orings$damage/6

qbinom = glm(pdmg ~ temp, family = quasibinomial, orings)
summary(qbinom)

plot(residuals(qbinom) ~ predict(qbinom, type = "link"), xlab = expression(hat(eta)), ylab = "Deviance Residuals")

halfnorm(influence(qbinom)$hat)

predict(qbinom, data.frame(temp = 31))

# fit the orings data with a binomial response and a logit link
logitmod = glm(cbind(damage, 6-damage) ~ temp, family = binomial, orings)
summary(logitmod)


#
#### CH7: Other Generalised Linear Models (GLMs) ####


# Gamma GLM
x = seq(0, 8 , by = 0.1)

plot(x, dgamma(x, 3), type = "l", ylab = "", xlab = "", ylim = c(0,1.25),
     xaxs = "i", yaxs = "i")


# revisit data from a step in the manufacturing process for semiconductors
# four factors are believed to influence the resistivity of the wafer and so
# full factorial experiment was run
# Previous experience led to the expectation that resistivity would have a 
# skewed distribution, so need to transform this
data(wafer)
summary(wafer)


# we can use log transform (Box-Cox method)
llmdl = lm(log(resist) ~.^2, wafer)
rlmdl = step(llmdl)
summary(rlmdl)

# we find a model with three two-way interactions, all with x3

# Now.. we fit corresponding gamma GLM
gmdl = glm(resist ~.^2, family = Gamma(link=log), wafer)
rgmdl = step(gmdl)
summary(rgmdl)

# the results are similar to the linear model

# the maximum likelihood estimate of phi (dispersion parameter):
require(MASS)
gamma.dispersion(rgmdl)


## more example
# Data on payment insurance claims for various areas of Sweden in 1977
# the data is subdivided by mileage driven, bonys of not having made previous claim,
# type of car

data(motorins)
motori = motorins[motorins$Zone == 1,]
gl = glm(Payment ~ offset(log(Insured)) + as.numeric(Kilometres) +
           Make + Bonus, family = Gamma(link = log), motori)
summary(gl)


# in comparison to the lognormal model:
llg = glm(log(Payment) ~ offset(log(Insured)) + as.numeric(Kilometres) +
            Make + Bonus, family = gaussian, motori)
summary(llg)


## let's compare the shapes of the distributions for the response using the dispersion
# estimates from the two models
x = seq(0,5, by = 0.05)

plot(x, dgamma(x, 1/0.55597, scale = 0.55597), type = "l", yaxs = "i", ylim = c(0,1))

plot(x, dlnorm(x, meanlog = -0.30551, sdlog = sqrt(0.55597)), type = "l",
     yaxs = "i", ylim = c(0,1))
# the lognormal model has higher kurtosis (sharper peak)

# we may also make predictions from both models
x0 = data.frame(Make = "1", Kilometres = 1, Bonus = 1, Insured = 100)

predict(gl, new = x0, se = TRUE, type = "response")
predict(llg, new = x0, se = TRUE, type = "response")
c(exp(10.998), exp(10.998)*0.16145)

#
#### Inverse Gaussian GLM ####

require(SuppDists)

x = seq(0, 8, by = 0.1)

plot(x, dinvGauss(x, 2, 10), type = "l", ylab = "", xlab = "", ylim = c(0,1.5), xaxs = "i", yaxs = "i")

# this inverse gaussian distribution is useful for modeling of lifetime distributions with non-monotone failure rates

# let's see an example of projected and actual sales data
data(cpd)

# use linear model first
lmod = lm(actual ~ projected-1, cpd)
summary(lmod)
plot(actual ~ projected, cpd)
abline(lmod)

# now consider the inverse Gaussian GLM, we must specify the identity link because we have y_i = beta * x_i
igmod = glm(actual ~ projected-1, family = inverse.gaussian(link="identity"), cpd)
summary(igmod)

plot(residuals(igmod) ~ log(fitted(igmod)), ylab ="Deviance residuals", xlab = expression(log(hat(mu))))
abline(h = 0)

#
#### Joint Modeling of the Mean and Dispersion ####

# In industrial experiments.
# sometimes we wisht to manufacture an item with a target mean or optimised response
# we set the predictors to produce items as close as possible to the mean. this requires a model for the mean
# We would also prefer that the variance of the response be small at the chosen value of the predictors for production

## example: we want to make a cake mix so that everyone can produce mediocre cake with less variance

# now look into a welding-strength experiment data
data(weldstrength)
lmod = lm(Strength ~ Drying + Material + Preheating, weldstrength)
summary(lmod)

# let's model the dispersion
# squared studentised residuals as the response in the dispersion with a gamma GLM using a log-link and weights of 1-h_i
h = influence(lmod)$hat
d = residuals(lmod)^2/(1-h)
gmod = glm(d ~ Method + Preheating, family = Gamma(link=log), weldstrength, weights = 1 - h)

# now feedbackt the estimated weights to the linear model
w = 1/fitted(gmod)
lmod = lm(Strength ~ Drying + Material + Preheating, weldstrength, weights = w)

# Iterate this until it converges
summary(gmod)
summary(lmod)

#
#### Quasi-Likelihood ####

# the benefit of quasi-binomial and quasi-poisson GLMs is that it allows for the dispersion parameter
# to be a free parameter, which is useful in modeling overdispersion
# you now can add a scale parameter (usually denoted by a sigma squared symbol) to the variance function (cause you have binomial, poisson)


# let's see data on sleep behaviour of 62 mammals
data(mammalsleep)
mammalsleep$pdr = with(mammalsleep, dream/sleep)
summary(mammalsleep$pdr)

# the dream/sleep proportion ranges from 0 to 0.46
# logit link seems sensible here
# Further, we expect the variance to be greater for moderate values of the proportion mu
# and less as mu approaches zero or one.
# This suggest a variance function of the approximate from mu(1-mu)
# this corresponds to the binomial GLM with the canonical logit link, yet the response is not binomial.
# we propose quasi-binomial:
modl = glm(pdr ~ log(body) + log(brain) + log(lifespan) + log(gestation) + predation + exposure + danger,
           family = quasibinomial, mammalsleep)
drop1(modl, test = "F") # since we use a free dispersion parameter, we must use F-tests

# backward elimination result in:
modl = glm(pdr ~ log(body) + log(lifespan) + danger, family = quasibinomial, mammalsleep)
summary(modl)

# diagnostics:
ll = row.names(na.omit(mammalsleep[,c(1,6,10,11)]))
halfnorm(cooks.distance(modl), labs = ll)
plot(predict(modl), residuals(modl, type = "pearson"), xlab = "Linear Predictor", ylab = "Pearson Residuals")

#
#### CH8: Random Effects ####

# let's start with the simplest possible random effects model
# one-way ANOVA with 'a' factor at 'a' levels
# y_ij = mu + alpha_i + E_ij
# variances: sigma^2_alpha & sigma^2_E
# Notice that this includes a correlation between observations at the same level =
# rho = sigma^2_alpha / (sigma^2_alpha + sigma^2_E)
# this is also known as Intraclass Correlation Coefficient (ICC)

data(pulp)
op = options(contrasts = c("contr.sum", "contr.poly"))
lmod = aov(bright ~ operator, pulp)
summary(lmod)
coef(lmod)
options(op)

# turning to random effects model, we can compute the variance of the operator effects:
(0.447-0.106)/5

require(lme4)
mmod = lmer(bright ~ 1+(1|operator), pulp) # this use Restricted Maximum Likelihood Estimator (REML)
summary(mmod)

# this model got a fixed effect = the intercept (denoted by the first '1')
# random effect is the (1|operator), the '1' here indicates the random effect is constant within each group
# we see here the REML method gives identical result with ANOVA. This may not be the case for unbalance design

# we can also compute the maximum likelihood estimates:
smod = lmer(bright ~ 1 + (1|operator), pulp, REML = FALSE)
summary(smod)


#
#### Inference ####

### Testing the Fixed Effect
# if you plan to use likelihood ratio test to compare two nested models that
# differ only in their fixed effects, you cannot use REML estimation.
# REML estimates the random effects by considering linear combinations of the data
# that removes the fixed effects.
# If the fixed effects are changed, the likelihoods of the two models
# will not be comparable.
# OPTIONS: 
# >> ordinary maximum likelihood
# >> parametric bootstrap
# >> F-test/t-test assuming the random part is equal to its estimated value

### Testing the Random Effect
# OPTIONS:
# chi-square approximation (conservative)
# bootsrap methods

### Expected Mean Squares
# this is a hypothesis test based on the SS found in ANOVA decompositions
# however, this requires adjustment for each model
# and the test cannot be used if the experiment is unbalanced

## ok let's see using the pulp data:
# NULL model
nullmod = lm(bright ~ 1, pulp)

# then compare with 
# smod (with random effect on the operator level) + using MLE
tempL = as.numeric(2* (logLik(smod) - logLik(nullmod)))
pchisq(tempL, 1, lower = FALSE)

# we can use parametric bootstrap to obtain a more accurate p-value
# In practice, we do not know the true value of mu and var, but
# we can use estimated values (this distinguishes parametric bootstrap from the simulation approach)

# we can simulate responses under the null:
y = simulate(nullmod)

# now taking the data we generate, we fit both the null and alt model
lrstat = numeric(1000)

for(i in 1:1000) {
  y = unlist(simulate(nullmod))
  bnull = lm(y ~ 1)
  balt = lmer(y ~ 1 + (1|operator), pulp, REML = FALSE)
  
  lrstat[i] = as.numeric(2*(logLik(balt)-logLik(bnull)))
  
}

# we may examine the distribution of the bootstrapped LRTs
# then compute the proportion that are close to zero
mean(lrstat < 0.00001)

# our estimated p-value
mean(lrstat > tempL)

# standard error of the estimate
sqrt(mean(lrstat > tempL) * 0.98/1000)

# there is evidence that the random effect is significantly explain the data
# relative to the null model

#### Predicting Random Effects ####

# in a fixed effect model, the parameters can be estimated
# BUT, in a random effects model, there ara no longer parameters
# it is just random variables under normal assumption: param ~ N(0,var)

# let's use bayesian approach
# Let f represent density then the posterio density for alpha (parameter)
# is given by:
# f(a_i|y) = f(y|a_i)f(a_i)

# we can then find posterior mean of a_i by integrating them, or..
# for general case, this work out to be:
# a_hat = D T(Z) (V)^-1 (y - XB)

# we take an empircal Bayes point of view and substitute the MLEs into
# D, V, and B to obtain the predicted random effects

ranef(mmod)$operator
cc = model.tables(lmod)
cc

# then compute the ratio to the random effect
cc[[1]]$operator / ranef(mmod)$operator
# typically the predicted random effects are smaller

# we can predict new operators by using the best linear unbiased predictors (BLUPs):
fixef(mmod) + ranef(mmod)$operator
# if we do not know the operator, we can just use the mu_hat = 60.4

# diagnostic plot for one-way random effects model:
qqnorm(resid(mmod), main = "")
plot(fitted(mmod), resid(mmod), xlab = "fitted", ylab = "residuals")
abline(0,0)

# random effects models are sensitive to outliers
# the residual fitted plot is also important cause we assumed that the
# error variance was constant

#### Blocks as Random Effects ####

# it is important to treat blocks as random effects
# we illustrate with an experiment to compare four processes: A, B, C, D
# for the production of penicillin
# the raw material: corn step liquor, is quite variable and can only be made
# in blends sufficient for four runs.
# thus a randomised complete block design is suggested by the nature of the experimental units

data(penicillin)

# we start with fixed effect analysis
op = options(contrasts=c("contr.sum", "contr.poly"))
lmod = aov(yield ~ blend + treat, penicillin)
summary(lmod)
coef(lmod)

# no significant effect between treatments but there is between blends

# now, let's see using a mixed model:
# fixed treatment effects and random blend effects
# this seems natural since the blends we use can be viewe as having been
# selected from some notional population of blends

mmod = lmer(yield ~ treat + (1|blend), penicillin)
summary(mmod)

# the residual variance is the same in both cases: 18.8
# this is true because we have a balanced design and so REML is
# equivalent to the ANOVA estimator
# The treatment effects are also the same as is the overall mean.

# The best linear unbiased predictors for the random effects are:
ranef(mmod)$blend

# we can test the significance of the fixed effects in two ways
# we can use ANOVA method, where we assume that the random effect parameters
# take their estimated values:
anova(mmod)

# the result is identical to the fixed effects analysis above.

# we can also test for a treatment effect using the maximum likelihood ratio method:
amod = lmer(yield ~ treat + (1|blend), penicillin, REML = FALSE)
nmod = lmer(yield ~ 1 + (1|blend), penicillin, REML = FALSE)

anova(amod, nmod)
# remember, we can only compare models with different fixed effects using ML method
# This is because in REML, the likelihood of linear combination not involving the
# fixed parameters is maximised.

# we can improve accuracy with the parametric bootstrap approach.
# we can generate a response from the null model and use this to compute LRT.
# We repeat this 1000 times saving LRT each time:

lrstat = numeric(1000)

for(i in 1:1000) {
  
  ryield = unlist(simulate(nmod))
  
  nmodr = lmer(ryield ~ 1 + (1|blend), penicillin, REML = FALSE)
  
  amodr = lmer(ryield ~ treat + (1|blend), penicillin, REML = FALSE)
  
  lrstat[i] = 2*(logLik(amodr) - logLik(nmodr))
  
}

# Under the standard likelihood theory, the LRT 
plot(qchisq((1:1000)/1001,3), sort(lrstat), xlab = expression(chi[3]^2),
     ylab = "Simulated LRT")
abline(0,1)

# as we can see from the plot, the approximation is not good
# we can compute aour estimated p-value as:
mean(lrstat > 4.05)

# which is much closer to the F-test result than the X_3 based approximation

## We can also test for the blends. As wth the fixed analysis, we are 
# no directly interested in the size of the blocking effects.
# However, we may wish the examine the blocking effects for information
# useful for future experiments

# we can compute the LRT:
rmod = lmer(yield ~ treat + (1|blend), penicillin)
nlmod = lm(yield ~ treat, penicillin)
2*(logLik(rmod) - logLik(nlmod, REML = TRUE))

# now we perform parametric bootstrap as before:
lrstatf = numeric(1000)

for(i in 1:1000) {
  
  ryield = unlist(simulate(nlmod))
  
  nlmodr = lm(ryield ~ treat, penicillin)
  
  rmodr = lmer(ryield ~ treat + (1|blend), penicillin)
  
  lrstatf[i] = 2*(logLik(rmodr) - logLik(nlmodr, REML = TRUE))
  
}

# again the distribution is far from X^2_1
mean(lrstatf < 0.00001)

mean(lrstatf > 2.7629)
# we found a significant blend effect, the p-value is close to 5% though
# we might wish to increase the number of bootstrap samples to increase
# our confidence

## In this example, we saw no major advantage in modeling the blocks as random effects
# so we might prefer to use the fixed effect analysis as it is simpler

#### Split Plots ####

# the difference between split plots and blocks design is that the earlier is used when there is one
# hard-to-change factor. Like, irrigation system in one plot or land.

# example: 8 fields, 2 types of crops, 4 irrigation methods.

data(irrigation)
summary(irrigation)

# Irrigation and variety are fixed effects
# Field is clearly a random effect
# interaction between field and crop is also random, cause one is random

# We consider this full model: Y_ijk = mu + I_i + V_j + (IV)_ij + F_k + (VF)_jk + e_ijk
# I = irrigation, V = variety, F = field

# We did not include (IF) term because we only have one type of irrigation used in a given field.
# It's impossible to estimate such.

lmodr = lmer(yield ~ irrigation * variety + (1|field), data = irrigation)
logLik(lmodr)

summary(lmodr)
# the largest variance component is that due to the field effect: 4.02 with residual variance 1.45

# we check the fixed effects for significance:
anova(lmodr) # no evidence for fixed effect for either irrigation or variety or their interaction

plot(fitted(lmodr), resid(lmodr), xlab = "Fitted", ylab = "Res")
qqnorm(resid(lmodr), main = "")
# no problem with non-constat variance
# so the result is valid


#### Nested Effects ####

# example: 4 samples (2 labelled G, 2 labelled H), 6 laboratories, gave 2 samples to 2 different technician,
# 1 sample divide by half > measure

data(eggs)
summary(eggs)

# we want to test consistency across labs
# so we will treat labs as random effects

cmod = lmer(Fat ~ 1 + (1|Lab) + (1|Lab:Technician) + (1|Lab:Technician:Sample), data = eggs)
summary(cmod)

# let us see if we can simplify by removing the lowest level of variance component
cmodr = lmer(Fat ~ 1 + (1|Lab) + (1|Lab:Technician), data = eggs)
summary(cmodr)

anova(cmod, cmodr) # no significant, but the p-value is conservative

VarCorr(cmodr) # the variance from Samples have been absorbed by the other components

# let's check the accuracy of these p-values
lrstat = numeric(1000)

for(i in 1:1000) {
  
  rFat = unlist(simulate(cmodr))
  
  nmod = lmer(rFat ~ 1 + (1|Lab) + (1|Lab:Technician), data = eggs)
  
  amod = lmer(rFat ~ 1 + (1|Lab) + (1|Lab:Technician) + (1|Lab:Technician:Sample), data = eggs)
  
  lrstat[i] = 2*(logLik(amod) - logLik(nmod))
  
}

mean(lrstat < 0.00001)

# we can estimate the p-value as:
2*(logLik(cmod) - logLik(cmodr))
mean(lrstat > 1.6034) # we can now say that the variation between samples can be ignore
# we can use the same method to test the variation between technicians


#### Crossed Effects ####

# full factorial design = all factors are completely crossed
# if only partially crossed = Crossed effects
# special case of Crossed effects = latin square design

# Example: 4 materials were fed into wear-testing machine.
# THe machine can take 4 samples at a time but the position is known to affect the result.

data(abrasion)
matrix(abrasion$material, 4, 4) # Latin Square!!

# fixed effect analysis
lmod = aov(wear ~ material + run + position, abrasion)
summary(lmod)

# we might regard the run and position as random effects
mmod = lmer(wear ~ material + (1|run) + (1|position), abrasion)
anova(mmod)
summary(mmod)

# lmer function is able to recognise that the run and position effects are crossed and fits the model appropriately.
# the F-test for the fixed effects is almost the same as the corresponding fixed effects analysis.
# The only differnce is that the fixed effects analysis uses a denominator degrees of freedom sixe while the random effects analysis
# is made conditional on the esitemated random effects ~ 12 degrees of freedom. The difference is not crucial here.

#### Multilevel Models ####

# Example:
# Classes within a school (up to four)
# social class of the father (1~9:: 7 = long term unemployed, 8 = not currently employed, 9 = father absent)
# raven's test score in year 1, student id number, english test score, math test score, school year (0 = year 1, 2 = year 3)

data(jsp)
# ignore the data from the first two years
jspr = jsp[jsp$year == 2,]


plot(jitter(math) ~ jitter(raven), data = jspr, xlab = "Raven score", ylab = "Math score")
boxplot(math ~ social, data = jspr, xlab = "Social class", ylab = "Math score")

# one possible approach: multiple linear regression
glin = lm(math ~ raven * gender * social, jspr)
anova(glin)

# remove gender
glin = lm(math ~ raven * social, jspr)
anova(glin)

# remove raven*social interaction even though it is significant, we have a large dataset
glin = lm(math ~ raven + social, jspr)
summary(glin)

## The above analysis assumes the students are from independent observations, but they are not!
# they come from 50 different schools. If we ignore this, we will overstate the significance.
table(jspr$school)

# let's try: random effects of school and the social class is nested within school
mmod = lmer(math ~ raven * social * gender + (1|school) + (1|school:class), jspr)
anova(mmod)

# gender is still not significant, drop it
# center the raven score, so that we can interpret social factor at the mean raven score instead of raven score of 0
jspr$craven = jspr$raven - mean(jspr$raven)
mmod = lmer(math ~ craven * social + (1|school) + (1|school:class), jspr)
anova(mmod)
summary(mmod)

# standard diagnostics:
qqnorm(resid(mmod), main = "")
plot(fitted(mmod), resid(mmod), xlab = "Fitted", ylab = "Residuals")
# hmmm... lower variance at the tail of predicted math score

# check assumption of normally distributed random effects
qqnorm(ranef(mmod)$school[[1]], main = "school effects")
qqnorm(ranef(mmod)$"school:class"[[1]], main = "class effects")
# approximately normal

## see sorted school effects
adjscores = ranef(mmod)$school[[1]]

# compare this with an adjusted ranking that simply takes the average score achieved by the school, centered by the overall average:
rawscores = coef(lm(math ~ school - 1, jspr)) #without intercept
rawscores = rawschores-mean(rawscores)

# plot this
plot(rawscores, adjscores)
sint = c(1,9,14,29)
text(rawscores[sint], adjscores[sint] + 0.2, c("1","9","15","30"))

## Compositional effects
# fixed effect predictors in this example so far have been at the lowest level, the student, but it is not improbable that factors
# at the school or class level might be important predictors of success in math test.
# We can construct such predictors from the individual-level information; such factors are called compositional effects.
# For example, the average entering score for a school might be an important predictor.

schraven = lm(raven ~ school, jspr)$fit

mmodc = lmer(math ~ craven * social + schraven * social + (1|school) + (1|school:class), jspr)
anova(mmodc)
# not significant, meh

#### Exercise ####

# page 199




