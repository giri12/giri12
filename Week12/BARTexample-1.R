rm(list = ls())
library(BayesTree)
library(rstanarm)

#'
#' We'll use the data from the Zigler et al (2018) Impact of Nonattainment Designations 
#' paper.  Note that the example below may not be a defensible way to analyze these data, but is just 
#' being used to illustrate how BART might work in reality
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

#' Let's do some simple analyses for comparison
#' 
#' First, let's check overlap by looking at the propensity score histograms
#' across the treatment groups 
#' (note: the propensity score has already been estimated and appears in the data set)
with(dat_use, hist(ps[a==0], breaks=50, col = rgb(0,0,1, alpha=0.4), 
                   main = "Histograms of PS", xlab = "ps"))
with(dat_use, hist(ps[a==1], breaks=50, col = rgb(1,0,0,alpha=0.4), add=TRUE))
legend("topright", c("Control", "Treated"), 
       fill = c(rgb(0,0,1, alpha=0.4), rgb(1,0,0, alpha=0.4)))


#' Identify whether each unit has an estimated propensity score that 'overlaps' with
#' with the distribution in the other treatment group
minps_a0 = with(dat_use[dat_use$a==0,], min(ps))
maxps_a0 = with(dat_use[dat_use$a==0,], max(ps))
minps_a1 = with(dat_use[dat_use$a==1,], min(ps))
maxps_a1 = with(dat_use[dat_use$a==1,], max(ps))

#' Just determine overlap based on observed range of PS in the 'other' group
overlap = rep(0, dim(dat_use)[1])
overlap[dat_use$a==0 & dat_use$ps >= minps_a1 & dat_use$ps <= maxps_a1] = 1
overlap[dat_use$a==1 & dat_use$ps >= minps_a0 & dat_use$ps <= maxps_a0] = 1
table(dat_use$a,overlap)

#' Fit a linear regression model that includes all the covariates as linear terms
mod1 <- stan_glm(pmfu~., data=dat_use[, !(names(dat_use) %in% c("ps"))])
#' The coefficient corresponding to the treatment variable estimates the ATE
ATE_1 = mod1$coefficients["a"]
sdATE_1 = mod1$ses["a"]
ATE_1
sdATE_1

#' Now estimate the ATT using posterior predictions for each treated observation
#' under the case where a=0
preddat = dat_use
preddat$a=0

#' Obtain posterior predictions for every unit under the case where a=0
ypred1 = posterior_predict(mod1, newdata=preddat)
#' Calculate a matrix of estimates of the ATT by subtracting the observed outcome
#' for treated units from the posterior predictions under a=0 for treated units
ATTmat = dat_use$pmfu[dat_use$a==1] - ypred1[, dat_use$a==1] # a 4000 x 218 matrix for the 218 treated units
ATT_1 = mean(rowMeans(ATTmat))
sdATT_1 = sd(rowMeans(ATTmat))
ATT_1
sdATT_1

#' 
#' Now fit the same model, but only use the overlapping units (where overlap is determined
#' by the value of the estimated propensity score)
#' 
dat_use_overlap = subset(dat_use, overlap==1)
mod2 <- stan_glm(pmfu~., data=dat_use_overlap[, !(names(dat_use) %in% c("ps"))])
ATE_2 = mod2$coefficients["a"]
sdATE_2 = mod2$ses["a"]
ATE_2
sdATE_2

preddat = dat_use_overlap
preddat$a=0
ypred2 = posterior_predict(mod2, newdata=preddat)
ATTmat = dat_use_overlap$pmfu[dat_use_overlap$a==1] - ypred2[,dat_use_overlap$a==1] # a 4000 x 131 matrix for the 131 treated units
ATT_2 = mean(rowMeans(ATTmat))
sdATT_2 = sd(rowMeans(ATTmat))
ATT_2
sdATT_2


c(ATE_1, ATE_2)
c(ATT_1, ATT_2)


#' 
#' Now analyze the data with  BART
#' 


####GIRI
names(dat_use)
!(names(dat_use) %in% c("pmfu", "ps"))
#### GIRI


#' Create training data that uses all the covariates (not the outcome)
xt=as.matrix(dat_use[,!(names(dat_use) %in% c("pmfu", "ps"))])  

#' Create test data that includes the covariates of all the treated units
#' but sets the treatment variable = 0
#' This will be used for predicting the *other* potential outcome for the treated units only
#' i.e., for estimating the ATT
xp=as.matrix(dat_use[dat_use$a==1,!(names(dat_use) %in% c("pmfu", "ps"))]) 
xp[,1]=0

y=as.numeric(dat_use[,1])

#' Fit the BART model
bart.tot <- bart(x.train=xt,   y.train=y,  x.test=xp)

# check convergence
library(coda)
plot(as.mcmc(bart.tot$sigma))

#' Use MCMC samples to calculate individual and average treatment effects

#' First calculate MCMC simulations of individual treatment effects by subtracting
#' the observed outcome among treated units from the predicted values had they been
#' untreated (i.e., the "test predictions")
diffs=bart.tot$yhat.train[,dat_use$a==1]-bart.tot$yhat.test 
head(diffs) # a matrix with 1000 MCMC samples (rows) for each of 218 treated units (columns)
dim(diffs)

#' Row means correspond to the SATE for each MCMC iteration
mndiffs=apply(diffs,1,mean)
length(mndiffs) #A vector of 1000 simulations of the SATE
ATT_bart = mean(mndiffs) # Posterior mean SATE
ATT_bart

sdATT_bart = sd(mndiffs) # Posterior standard deviation of the SATE
sdATT_bart

c(ATT_1, ATT_2, ATT_bart)
c(sdATT_1, sdATT_2, sdATT_bart)

###############################################################################


#' Now look at the Individual Treatment Effects, (ITEs), noting that estimation may have
#' wide uncertainty with this sample size
#' 
#' The posterior ITEs are the columns of the 1000x218 matrix 'diffs'
ite_means<- apply(diffs, 2, mean)
ite_sds<- apply(diffs, 2, sd)
ite_ql = apply(diffs, 2, quantile, .025)
ite_qu = apply(diffs, 2, quantile, .975)

#' Just get a sense of the treatment effect heterogeneity by looking at a histogram
#' of the posterior mean ITEs across the sample
hist(ite_means, breaks=100)


#' Now plot the ITEs (posterior means and 95% intervals) across the values of a covariate
#' Do this for every covariate just for illustration
for (cov in covariates){
  covplot = dat_use[, cov]
  plot(covplot[dat_use$a==1], ite_means, pch=16, cex=0.75, col="red", ylim = c(-1,0.5), 
       main = paste("ITEs as a funciton of:", cov), xlab = cov, ylab = "ITE")
  arrows(covplot[dat_use$a==1], ite_ql, covplot[dat_use$a==1], ite_qu, col = rgb(0.5,0,0, alpha=0.5), angle=90, length=0.01, lwd=0.5)
}

###############################################################################

#' Now look at BART performance but only in the overlap sample
#' Create training data that uses all the covariates (not the outcome)
xt=as.matrix(dat_use_overlap[,!(names(dat_use) %in% c("pmfu", "ps"))])  
xp=as.matrix(dat_use_overlap[dat_use_overlap$a==1,!(names(dat_use) %in% c("pmfu", "ps"))]) 
xp[,1]=0

y=as.numeric(dat_use_overlap[,1])

#' Fit the BART model
bart.tot.overlap <- bart(x.train=xt,   y.train=y,  x.test=xp)

diffs=bart.tot.overlap$yhat.train[,dat_use_overlap$a==1]-bart.tot.overlap$yhat.test 
mndiffs=apply(diffs,1,mean)
ATT_bart2 = mean(mndiffs) # Posterior mean SATE
ATT_bart2

sdATT_bart2 = sd(mndiffs) # Posterior standard deviation of the SATE
sdATT_bart2

c(ATT_1, ATT_2, ATT_bart, ATT_bart2)
c(sdATT_1, sdATT_2, sdATT_bart, sdATT_bart2)




