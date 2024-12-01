rm(list = ls())
library(rstanarm)

#' First illustrate IPTW with a very simplistic simulated example
n=100000
x1=rnorm(n,2,2)
x2=rnorm(n,5,2)
logitps = -6 + 0.5*x1+0.5*x2
ps = exp(logitps)/(1+exp(logitps))
a = rbinom(n,1,ps)
y=rnorm(n, x1+x2+2*a, 1) #true ATE is 2

dat_use = as.data.frame(cbind(y,a,x1,x2))

#'Estimate the propensity score
dat_use$ps = glm(a~x1+x2, family= binomial(link = "logit"), data=dat_use)$fitted.values

#' Let's look at the propensity score distribution across treatment groups.
with(dat_use, hist(ps[a==0], breaks=50, col = rgb(0,0,1, alpha=0.4),
                   main = "Histograms of PS", xlab = "ps"))
with(dat_use, hist(ps[a==1], breaks=50, col = rgb(1,0,0,alpha=0.4), add=TRUE))
legend("topright", c("Control", "Treated"), 
       fill = c(rgb(0,0,1, alpha=0.4), rgb(1,0,0, alpha=0.4)))

#' Now use the ps to define inverse probability weights
dat_use$wt[dat_use$a==1] = 1/dat_use$ps[dat_use$a==1]
dat_use$wt[dat_use$a==0] = 1/(1-dat_use$ps[dat_use$a==0])
with(dat_use, tapply(wt, a, max))

with(dat_use, hist(wt[a==1], breaks=500, col = rgb(1,0,0, alpha=0.4),xlim = c(0,20),
                   main = "Histogram of Weights", xlab = "IPW"))
with(dat_use, hist(wt[a==0], breaks=500, col = rgb(0,0,1,alpha=0.4), add=TRUE))
legend("topright", c("Control", "Treated"), 
       fill = c(rgb(0,0,1, alpha=0.4), rgb(1,0,0, alpha=0.4)))


#' Show that the propensity score distribution is equivalent across the treatment groups in the pseudopopulation
with(dat_use, tapply(ps, a, mean))
with(subset(dat_use, a==0), weighted.mean(ps, wt))
with(subset(dat_use, a==1), weighted.mean(ps, wt))

#' Show that weighted covariate means are more balanced than raw data
with(dat_use, tapply(x1, a, mean))
with(subset(dat_use, a==0), weighted.mean(x1, wt))
with(subset(dat_use, a==1), weighted.mean(x1, wt))

#'
#' Analyze the Outcome
#' 
#' 
#' Unadjusted mean comparison
with(dat_use, tapply(y, a, mean))
with(dat_use, mean(y[a==1]) - mean(y[a==0]))
summary(lm(y~a, data=dat_use))

#' Adjusted comparison
#' In this case we know the true model, so simplre regression adjustment should work fine
summary(lm(y~a+x1+x2, data=dat_use))

#' Weighted mean comparison
with(subset(dat_use, a==1), weighted.mean(y, wt)) - with(subset(dat_use, a==0), weighted.mean(y, wt))

#' Or with the survey package
library(survey)
d_wt<- svydesign(~1, weights=dat_use$wt, data=dat_use)
wtmod<- svyglm(y~a, design=d_wt)
summary(wtmod)


#'
#' Now let's see what this looks like with a real data set.  We'll use the data from
#' the Zigler et al (2018) Impact of Nonattainment Designations paper.  Note that the 
#' example below is probably not a defensible way to analyze these data, but is just 
#' being used to illustrate how IPTW might work in reality.
#'
load('pm_withps_nomed.RData') #loads in a data frame called 'dat'
head(dat)
dim(dat)

#' Subset only to the variables that will be used for the illustration
covariates<- c("pmbase2002_2004","avgdewpt", "avgtemp", "avgrelhum", "CompletelyRural", "ps", "Longitude", "Latitude",
"TotPop", "PctUrban", "PctBlack", "PctHisp", "PctHighSchool", "MedianHHInc", "PctPoor", "PctFemale", "PctOccupied",
"PctMovedIn5", "MedianHValue", "smokerate2000")
treatment<- "a"
response<- "pmfu"

dat_use = dat[!is.na(dat$pmfu), c(response, treatment,  covariates)]
dim(dat_use)

#' These data come with an estimate of the propensity score ('ps').  Let's look at its distribution
#' across treatment groups.
with(dat_use, hist(ps[a==0], breaks=50, col = rgb(0,0,1, alpha=0.4), 
                   main = "Histograms of PS", xlab = "ps"))
with(dat_use, hist(ps[a==1], breaks=50, col = rgb(1,0,0,alpha=0.4), add=TRUE))
legend("topright", c("Control", "Treated"), 
       fill = c(rgb(0,0,1, alpha=0.4), rgb(1,0,0, alpha=0.4)))

#' Now use the ps to to define inverse probability weights
dat_use$wt[dat_use$a==1] = 1/dat_use$ps[dat_use$a==1]
dat_use$wt[dat_use$a==0] = 1/(1-dat_use$ps[dat_use$a==0])
with(dat_use, tapply(wt, a, max))

with(dat_use, hist(wt[a==1], breaks=50, col = rgb(1,0,0, alpha=0.4),xlim = c(0,35),
                   main = "Histogram of Weights", xlab = "IPW"))
with(dat_use, hist(wt[a==0], breaks=50, col = rgb(0,0,1,alpha=0.4), add=TRUE))
legend("topright", c("Control", "Treated"), 
       fill = c(rgb(0,0,1, alpha=0.4), rgb(1,0,0, alpha=0.4)))

#' Show that the propensity score distribution is equivalent across the treatment groups in the pseudopopulation
with(dat_use, tapply(ps, a, mean))
with(subset(dat_use, a==0), weighted.mean(ps, wt))
with(subset(dat_use, a==1), weighted.mean(ps, wt))

#' Show that weighted covariate means are more balanced than raw data
#' (should actually check this for all covariates)
with(dat_use, tapply(pmbase2002_2004, a, mean))
with(subset(dat_use, a==0), weighted.mean(pmbase2002_2004, wt))
with(subset(dat_use, a==1), weighted.mean(pmbase2002_2004, wt))

#'
#' Analyze the Outcome
#' 
#' Unadjusted mean comparison
with(dat_use, tapply(pmfu, a, mean))
with(dat_use, mean(pmfu[a==1]) - mean(pmfu[a==0]))
summary(lm(pmfu~a, data=dat_use))

#' Adjusted comparison
#' This puts a lot of faith in our regression model....
summary(lm(pmfu~., data=dat_use[, !(names(dat_use) %in% c("wt", "ps"))]))

#' Weighted mean comparison
with(subset(dat_use, a==1), weighted.mean(pmfu, wt)) - with(subset(dat_use, a==0), weighted.mean(pmfu, wt))

#' Or with the survey package
d_wt<- svydesign(~1, weights=dat_use$wt, data=dat_use)
wtmod<- svyglm(pmfu~a, design=d_wt)
summary(wtmod)

#' We could also fit a weighted model that *also* adjusts for some (or all) covariates, which
#' might be able to use the parametric regression specification to adjust for any 
#' residual imbalances that remain after the weighting
d_wt<- svydesign(~1, weights=dat_use$wt, data=dat_use[, !(names(dat_use) %in% c("wt", "ps"))])
wtmod2<- svyglm(pmfu~a + pmbase2002_2004 + PctUrban + MedianHHInc, design=d_wt)
summary(wtmod2)

