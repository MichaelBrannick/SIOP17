######################### Call needed packages###############################
library(MASS)           # need this library for multivariate sampling
library(psych)          # need this for descriptive statistics routine
library(metafor)        # need this for meta-regression and for Morris
rm(list=ls()) #Cleans (empties) global environment
dev.off()     # cleans plots/graphics
cat("\014")   #cleans out the console
#########################################################
################## PART 1 Parameters Input for Simulation ###################
rhodistribution = "normal"
numbersamples = 30     # k, the number of studies per meta-analysis
rho = .26                # population correlation = rhoTxTy
sigmarho = .05          # population standard deviation, SDrho
numbermetas = 100      # number of metas, each based on numbersamples - replications

popCR.L <- rho-1.28*sigmarho
popCR.U <- rho+1.28*sigmarho
pop.width.CR <- popCR.U-popCR.L
pop.CR <- cbind(popCR.L, popCR.U, pop.width.CR)
########################################################
#  DESIGN ELEMENTS
########################################################
#  rho, SDrho and distribution shape are fixed factors
#  estimators to be evaluated (S&H, S&Hk, Hedges, Morris, Unit)
#  outputs for evaluation - Bias, RMSE, coverage of overall mean 
########################################################
# rho = .14, .26, .42 (based on Paterson, 2016)
# SDrho = 0, .08, .13, .20 (based on Paterson, 2016)
# number of samples (number of effect sizes, k) = 5, 10, 30, and 100
# reliability of x and y are random variables (uses Le & Schmidt 2006)
# sample size is a random variable (uses Gamma distribution, median close to Paterson)
#########################################################
# Let Tx signify true scores on x, let x represent observed scores
# Let Ty signify true scores on y, let y represent observed scores
# Let rhoTxTy signify the correlation between true scores
# Let rhoxy signify the correlation between observed scores on x and y
# Let rhoTxx signify the correlation between true and observed score on x, = sqrt(rxxp)
# Let rhoTyy signify the correlation between true and observed scores on y, = sqrt(ryyp)
# Then rhoxy = rhoTxTy*rhoTxx*rhoTyy
#
#reliability values for independent and dependent variables from Le & Schmidt 2006
rxxp <- c(.5, .5, .6, .6, .6, .6, .7, .7, .7, .7, .75, .75, .75, .75, .75,
          .75, .75, .75, .75, .75, .75, .75, .75, .75, .75, .75, .75, .75, .75, .75,
          .8, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8,
          .8, .8, .8, .8, .8, .8, .8, .85, .85, .85, .85, .85, .85, .85, .85, .85,
          .85, .85, .85, .85, .85, .85, .85, .85, .85, .85, .85, .85, .85, .85, .85,
          .85, .85, .85, .85, .85, .85, .9, .9, .9, .9, .9, .9, .9, .9, .9, .9,
          .9, .9, .9, .9, .9)
#
ryyp <- c(.3, .3, .3, .35, .35, .35, .35, .4, .4, .4, .4, .4, .4, .45, .45, .45, .45,
         .45, .45, .45, .45, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .55, .55, .55,
         .55, .55, .55, .55, .55, .55, .55, .55, .55, .6, .6, .6, .6, .6, .6, .6, 
         .6, .6, .6, .6, .6, .6, .6, .65, .65, .65, .65, .65, .65, .65, .65, .65, 
         .65, .65, .65, .7, .7, .7, .7, .7, .7, .7, .7, .7, .7, .75, .75, .75, .75,
         .75, .75, .75, .75, .8, .8, .8, .8, .8, .8, .85, .85, .85, .85, 
         .9, .9, .9)
#################################################################################
# to simulate without reliability artifacts, set the values of rxx and ryy to 1.0
################## PART 2 Generating data #######################################

###################### Revised Sample Sizes - Min(N) = 30 #####################
# large samples on average, distribution is right skewed to mimic
# observed distributions of sample sizes in meta-analyses.
###############################################################################
set.seed(1234)
bigN2 = rgamma(1000000, shape = .57, scale = 500)
bigN2 = ceiling(bigN2)

set.seed(1234)
# Keep only samples greater than 29.  
bigN = sample(bigN2[bigN2 > 29], size = 100000, replace = T)

# Use the following vectors for the actual simulation.
s.big = matrix(nrow = numbermetas, ncol = numbersamples)
for (i in 1:numbermetas){
  s.big[i, ] = sample(bigN, size = numbersamples, replace = T)
}

############################placeholders for output ##################
rhobox <- matrix(0,4,4)+diag(4) # placeholder for correlation matrix for sampling
mu <- rep(0,4)          # multivariate means of zero
names.sim <- cbind("Tx", "Ty", "X", "Y")
colnames(rhobox)<- names.sim
rownames(rhobox)<- names.sim
out1 <- 1:numbersamples # for meta-analysis sampled ri  
out2 <- out1            # sampled rxxp - population value
out3 <- out1            # sampled ryyp - population value
out4 <- out1            # sample from rxxi - sample estimate of reliability of x
out5 <- out1            # sample from ryyi - sample estimate of reliability of y
Results <- matrix(0,numbermetas,47) # for results of meta-analyses

#######################################################################
#begin outer loop - do the following for each meta-analysis until numbermetas is reached
for (j in 1:numbermetas){
  # begin inner loop - create data for one meta-analysis #
  for (i in 1:numbersamples){                       # do the folowing for number of samples per meta
    rxxi <-  sample(rxxp, size = 1, replace = T)     # sample from reliability of x
    ryyi <-  sample(ryyp, size = 1, replace = T)     # sample from reliability of y
    rho.r <- rnorm(1,rho,sigmarho)                  # sample a value of rho from the parameters specified
    if (rho.r > .99){rho.r <-rnorm(1,rho,sigmarho)} # check for admissibility & resample if needed
    if (rho.r > .99){rho.r <-.99}                   # set to boundary on second try
    if (rho.r < -.99) {rho.r <- -.99}               # set to boundary if necessary - didn't use the second sample
    att.rho.r <- rho.r*sqrt(rxxi*ryyi)              # attenuate rho for reliability
    # find the population correlation matrix from which to sample observations
    rhobox[2,1] <- rho.r                            # correlation between true scores
    rhobox[1,2] <- rho.r
    rhobox[3,1] <- sqrt(rxxi)                       # true x and observed X
    rhobox[1,3] <- sqrt(rxxi)
    rhobox[3,2] <- rho.r*sqrt(rxxi)                 # true Y and observed X               
    rhobox[2,3] <- rho.r*sqrt(rxxi)
    rhobox[4,2] <- sqrt(ryyi)                       # true Y and observed Y
    rhobox[2,4] <- sqrt(ryyi)
    rhobox[4,1] <- rho.r*sqrt(ryyi)                 # true X and observed Y
    rhobox[1,4] <- rho.r*sqrt(ryyi)
    rhobox[4,3] <- att.rho.r                        # observed X and Y at pop level
    rhobox[3,4] <- att.rho.r
    
    sample.rs <-mvrnorm(n=s.big[j,i],mu=mu,Sigma=rhobox) # sample ni data for obs r from attenuated rho
    cor1 <-cor(sample.rs)                           # compute correlation from sampled data
    out1[i] <- cor1[4,3]                            # save the sampled correlation for output
    out2[i] <- rxxi                                 # save the pop reliability of x
    out3[i] <- ryyi                                 # save the pop reliability of y
    out4[i] <- cor1[3,1]                            # save sample sqrt(rel.x)
    out5[i] <- cor1[4,2]                            # save sample sqrt(rel.y)
  } 
  # end inner loop
  
  ################## PART 3 Define values for meta-analysis, run meta-analysis #########
  ri <-out1                       # correlations for one meta-analysis
  rxx.s <- out4                     # associated estimated reliabilities of x
  ryy.s <- out5                     # associated estimated reliabilities of y
  ni <- s.big[j, ]                # sample sizes for each correlation
  k <- length(ri)                 # number of correlations
  Higgins.t <- abs(qt(.1, k-2))
  ai <- rxx.s*ryy.s                 # compound attenuation factor
  checker <- data.frame(ri,ni,ai) # create dataset for use by metafor (meta-regression)

  ########################################################################
  # Hedges 
  # Borenstein et al. (2009) method - use meta-regression
  # uses metafor (Viechtbauer, 2010) to estimate a slope for the combined
  # attenuation factor (sqrt(rxx*ryy))
  ########################################################################
  resi <- rma(ri=ri, ni=ni, measure = 'ZCOR', method='REML', data=checker,
              mods=~ai, control=list(maxiter=1000, stepadj=.5))
  preds <- predict(resi, newmods=1)     # preds contains the predicted value if no attenuation
  Hedges.M <-tanh(preds$pred)           # backtranslated value of disattenuated mean
  Hedges.CI95.U <- tanh(preds$ci.ub)    # backtranslated value of upper bound CI
  Hedges.CI95.L <- tanh(preds$ci.lb)    # backtranslated value of lower bound CI
  Hedges.CR.80.U <- tanh(preds$cr.ub)   # backtranslated value of upper bound CR
  Hedges.CR.80.L <- tanh(preds$cr.lb)   # backtranslated value of lower bound CR
  Hedges.coverage <-0
  if(Hedges.CI95.L<= rho & rho<= Hedges.CI95.U) {Hedges.coverage<-1}
  Hedges.tau.est <- sqrt(resi$tau2)
  Hedges.CR.Ldiff <- abs(Hedges.CR.80.L-popCR.L)
  Hedges.CR.Udiff <- abs(Hedges.CR.80.U-popCR.U)
  Hedges.CR.width <- Hedges.CR.80.U-Hedges.CR.80.L
  Hedges.CR.alldiff <- Hedges.CR.Ldiff+Hedges.CR.Udiff
  ########################################################################
  # Schmidt & Hunter (2015)
  # formulas for correcting for unreliability
  # attenuation factor (sqrt(rxx*ryy))
  ########################################################################
  rbar <- sum(ri*ni)/sum(ni)              # mean weighted observed r
  ##########################################
  # Correct for reliability only
  ##########################################
  A.compound <- ai                          # compound attenuation factor: sqrt(rxx*ryy))
  rC3 <- ri/A.compound                      # find corrected (disattenuated) correlations
  wi <- A.compound^2*ni                     # find the weights for the meta (S&H 2015 p.147, Eq 3.29)
  rbarC.3 <- sum(wi*rC3)/sum(wi)            # find the mean of corrected correlations
  V.rC.3 <- sum(wi*(rC3-rbarC.3)^2)/sum(wi) # find the wtd variance of the corrected correlations
  V.rC.3b <- k*V.rC.3/(k-1)                 # k correction for small number of studies
  V.eo.3 <- (1-rbar^2)^2/(ni-1)             # simple observed error variance for each study
  V.ec3 <- V.eo.3/A.compound^2              # simple error variance of corrected correlations
  V.ve3 <- sum(wi*V.ec3)/sum(wi)            # find error variance of corrected correlations
  # do not adjust sampling error for effects of range restriction
  V.rho.3 <- V.rC.3-V.ve3                   # find the random-effects variance component
  V.rho.3b <- V.rC.3b-V.ve3                 # k correction find the random-effects variance compoent 
  if (V.rho.3 < 0) {V.rho.3 <- 0}           # if REVC is less than zero, set to zero
  if (V.rho.3b < 0) {V.rho.3b <- 0}           # if REVC is less than zero, set to zero
  SD.rho.3 <- sqrt(V.rho.3)                 # find SD rho
  SD.rho.3b <- sqrt(V.rho.3b)                 # k correction find SD rho
  SEM3 <- sqrt(V.rC.3/k)                    # S&H 2015 Equation 5.2; find standard error of the mean of corrected correlations
  SEM4 <- sqrt(V.rC.3b/k)                   # k corrected standard error
  HS.CI95.U <- rbarC.3 + 1.96*SEM3           # find the confidence interval for the mean of corrected corrs
  HS.CI95.L <- rbarC.3 - 1.96*SEM3
  HS.CR.80.U <- rbarC.3 + 1.28*SD.rho.3      # credibility Interval
  HS.CR.80.L <- rbarC.3 - 1.28*SD.rho.3
  HSk.CI95.U <- rbarC.3 + 1.96*SEM4          # k corrected confidence interval
  HSk.CI95.L <- rbarC.3 - 1.96*SEM4
  HSk.CR.80.U <- rbarC.3 + 1.28*SD.rho.3b    # k corrected credibilty interval
  HSk.CR.80.L <- rbarC.3 - 1.28*SD.rho.3b
  HS.coverage <-0
  if(HS.CI95.L<= rho & rho<= HS.CI95.U) {HS.coverage<-1}
  HSk.coverage <-0
  if(HSk.CI95.L<= rho & rho<= HSk.CI95.U) {HSk.coverage<-1}
  HS.CR.Ldiff <- abs(HS.CR.80.L-popCR.L)
  HS.CR.Udiff <- abs(HS.CR.80.U-popCR.U)
  HS.CR.width <- HS.CR.80.U-HS.CR.80.L
  HS.CR.alldiff <- HS.CR.Ldiff+HS.CR.Udiff
  #
  HSk.CR.Ldiff <- abs(HSk.CR.80.L-popCR.L)
  HSk.CR.Udiff <- abs(HSk.CR.80.U-popCR.U)
  HSk.CR.width <- HSk.CR.80.U-HSk.CR.80.L
  HSk.CR.alldiff <- HSk.CR.Ldiff+HSk.CR.Udiff
  HS.tau.est <- SD.rho.3
  HSk.tau.est <- SD.rho.3b
  ########################################################################
     # Unit weights (Bonnet)
     # unweighted mean attenutation (sqrt(rxx*ryy)) computed over samples 
     # mean attenuation applied to overall mean and CI endpoints             
  ########################################################################
  Unit.ri <- ri/ai
  Unit.Mean <- sum(Unit.ri)/k
  Unit.V <- sum((Unit.ri-Unit.Mean)^2)/(k-1)
  Unit.SEM <- sqrt(Unit.V/k)
  Un.CI95.U <-Unit.Mean + 1.96*Unit.SEM               #confidence interval
  Un.CI95.L <- Unit.Mean - 1.96*Unit.SEM
  Un.CR.80.U <-Unit.Mean + 1.28*sqrt(Unit.V)          #Attempt at credibility interval to provide comporable values to other methods        
  Un.CR.80.L <-Unit.Mean - 1.28*sqrt(Unit.V)
  Unit.coverage <-0
  if(Un.CI95.L<= rho & rho<= Un.CI95.U) {Unit.coverage<-1}
  Unit.CR.Ldiff <- abs(Un.CR.80.L-popCR.L)
  Unit.CR.Udiff <- abs(Un.CR.80.U-popCR.U)
  Unit.CR.width <- Un.CR.80.U-Un.CR.80.L
  Unit.CR.alldiff <- Unit.CR.Ldiff+Unit.CR.Udiff
  Unit.tau.est <- sqrt(Unit.V)
  ########################################################################
  # Morris
  # This simulation uses metafor and the REML estimator
  # Attenuation factor applied to each study's effect size and sampling error.             
  ########################################################################
  # Computation of Morris Weights estimates, method 1.
  M1 <- sum(ni*ri)/(k*ni)                       # weighted mean, uncorrected
  var.i <-((1-M1^2)^2)/(ni-1)                     # individual study sampling variance
  var.i2 <- var.i/ai^2                            # corrected individual study sampling variance
  ri.m <- ri/ai                                   # corrected individual study ES
  morris.dat <- data.frame(cbind(ri.m,var.i2))    # collect corrected estimates
  morris1 <- rma(yi=ri.m,vi=var.i,data=morris.dat,
                 control=list(maxiter=1000, stepadj=.5))        # run the random-effects meta with REML
  Morris.M.rho <- morris1$b                      # output the mean
  Morris.CI95.L <- morris1$ci.lb                 # output the lower CI bound
  Morris.CI95.U <- morris1$ci.ub                 # output the upper CI bound
  Morris.SD.rho <- sqrt(morris1$tau2)            # output for Morris RHO sd
  Morris.coverage <-0                            # figure the covarage
  if(Morris.CI95.L<= rho & rho<= Morris.CI95.U) {Morris.coverage<-1}
  Morris.CR80.U <- Morris.M.rho + 1.28*Morris.SD.rho   # lower CR Bound                      
  Morris.CR80.L <- Morris.M.rho - 1.28*Morris.SD.rho   # Upper CR Bound
  Morris.CR.Ldiff <- abs(Morris.CR80.L-popCR.L)
  Morris.CR.Udiff <- abs(Morris.CR80.U-popCR.U)
  Morris.CR.width <- Morris.CR80.U-Morris.CR80.L
  Morris.CR.alldiff <- Morris.CR.Ldiff+Morris.CR.Udiff
  Morris.tau.est <- Morris.SD.rho
  ##########################################
  # Output results
  ##########################################
  Results[j, ] <- cbind(rho, sigmarho, 
                        Hedges.M, Hedges.CI95.L, Hedges.CI95.U, Hedges.coverage, Hedges.CR.80.L, Hedges.CR.80.U,
                        rbarC.3, HS.CI95.L, HS.CI95.U, HS.coverage, HS.CR.80.L, HS.CR.80.U,
                        rbarC.3, HSk.CI95.L, HSk.CI95.U, HSk.coverage, HSk.CR.80.L, HSk.CR.80.U,
                        Unit.Mean, Un.CI95.L, Un.CI95.U, Unit.coverage, Un.CR.80.L, Un.CR.80.U,
                        Morris.M.rho, Morris.CI95.L, Morris.CI95.U, Morris.coverage, Morris.CR80.L, Morris.CR80.U,
                        Hedges.CR.alldiff, Hedges.CR.width, HS.CR.alldiff, HS.CR.width,
                        HSk.CR.alldiff, HSk.CR.width, Morris.CR.alldiff, Morris.CR.width,
                        Hedges.tau.est, HS.tau.est, HSk.tau.est, Morris.tau.est,
                        Unit.CR.alldiff, Unit.CR.width, Unit.tau.est)
  # Results[j, ] is the results for one meta. 'j' is the subscript for numbermetas.
} # end outer loop################################################################


# 
# Rename the columns of the Results  
resultsnames = c("rho", "sigmarho", 
                 "Hedges.M", "Hedges.CI95.L", "Hedges.CI95.U", "Hedges.coverage", "Hedges.CR.80.L", "Hedges.CR.80.U",
                 "S&H.M", "S&H.CI95.L", "S&H.CI95.U", "S&H.coverage", "S&H.CR.80.L", "S&H.CR.80.U",
                 "S&Hk.M", "S&Hk.CI95.L", "S&Hk.CI95.U", "S&Hk.coverage", "S&Hk.CR80.L", "S&Hk.CR80.U",
                 "Unit.M", "Unit.CI95.L", "Unit.CI95.U", "Unit.coverage", "Unit.CR80.L", "Unit.CR80.U",
                 "Morris.M", "Morris.CI95.L", "Morris.CI95.U", "Morris.coverage", "Morris.CR80.L", "Morris.CR80.U",
                 "HedgesDiff", "Hedgeswidth", "HSDiff", "HSwidth",
                 "HSkDiff", "HSkwidth", "MorrisDiff", "Morrswidth",
                 "Hedgestau", "HStau", "HSktau", "Morristau",
                 "UnitDiff", "Unitwidth", "Unittau")
colnames(Results) = resultsnames
# resultsnames
# Results # remove the comment hashtag if you want to see each meta's results
# Calculate the biases for the 4 types of estimates. 
################## PART 4 Compare the Estimators ###################################################
# Estimator bias

HedgesBias <- sum(Results[,3]-Results[,1])/numbermetas
SnHBias <- sum(Results[,9]-Results[,1])/numbermetas
SnHkBias <- sum(Results[,15]-Results[,1])/numbermetas
UnitBias <- sum(Results[,21]-Results[,1])/numbermetas
MorrisBias <- sum(Results[,27]-Results[,1])/numbermetas
# Estimator means
HedgesMean <- mean(Results[,3])
SnHMean <- mean(Results[,9])
SnHkMean <- mean(Results[,15])
UnitMean <- mean(Results[,21])
MorrisMean <- mean(Results[,27])
# Estimator RMSE 
HedgesRMSE <- sqrt(sum((Results[,3]-Results[,1])^2)/numbermetas)
SnHRMSE <- sqrt(sum((Results[,9]-Results[,1])^2)/numbermetas)
SnHkRMSE <- sqrt(sum((Results[,15]-Results[,1])^2)/numbermetas)
UnitRMSE <- sqrt(sum((Results[,21]-Results[,1])^2)/numbermetas)
MorrisRMSE <- sqrt(sum((Results[,27]-Results[,1])^2)/numbermetas)
# Estimator Coverage
HedgesCoverage <-mean(Results[,6])
SnHCoverage <-mean(Results[,12])
SnHkCoverage <- mean(Results[,18])
UnitCoverage <- mean(Results[,24])
MorrisCoverage <- mean(Results[,30])
# Calculate the average upper bounds for the confidence intervals
# Lower CI bounds
HedgesLower <-mean(Results[,4])
SnHLower <-mean(Results[,10])
SnHkLower <- mean(Results[,16])
UnitLower <- mean(Results[,22])
MorrisLower <- mean(Results[,28])
# Upper CI bounds
HedgesUpper <-mean(Results[,5])
SnHUpper <-mean(Results[,11])
SnHkUpper <- mean(Results[,17])
UnitUpper <- mean(Results[,23])
MorrisUpper <- mean(Results[,29])
# Lower CR bounds
HedgesCRLower <-mean(Results[,7])
SnHCRLower <-mean(Results[,13])
SnHkCRLower <- mean(Results[,19])
UnitCRLower <- mean(Results[,25])
MorrisCRLower <- mean(Results[,31])
# Upper CR bounds
HedgesCRUpper <-mean(Results[,8])
SnHCRUpper <-mean(Results[,14])
SnHkCRUpper <- mean(Results[,20])
UnitCRUpper <- mean(Results[,26])
MorrisCRUpper <- mean(Results[,32])
# Credibility Interal Discrepancy
HedgesDiff <- mean(Results[,33])
HSDiff <- mean(Results[,35])
HSkDiff <- mean(Results[,37])
MorrisDiff <- mean(Results[,39])
UnitDiff <- mean(Results[45])
# Credibility Interval Width
HedgesWide <- mean(Results[,34])
HSWide <- mean(Results[,36])
HSkWide <- mean(Results[,38])
MorrisWide <- mean(Results[,40])
UnitWide <- mean(Results[,46])
#Tau Estimates
Hedgestau <- mean(Results[,41])
HStau <- mean(Results[,42])
HSktau <- mean(Results[,43])
Morristau <- mean(Results[,44])
Unittau <- mean(Results[,47])


# Output
Output <- data.frame(numbersamples, rho, sigmarho, rhodistribution,
                     HedgesBias, HedgesCoverage, HedgesRMSE, HedgesLower, HedgesUpper, HedgesCRLower, HedgesCRUpper,
                     SnHBias, SnHCoverage, SnHRMSE, SnHLower, SnHUpper, SnHCRLower, SnHCRUpper,
                     SnHkBias, SnHkCoverage, SnHkRMSE, SnHkLower, SnHkUpper, SnHkCRLower, SnHkCRUpper,
                     UnitBias, UnitCoverage, UnitRMSE, UnitLower, UnitUpper, UnitCRLower, UnitCRUpper,
                     MorrisBias, MorrisCoverage, MorrisRMSE, MorrisLower, MorrisUpper, MorrisCRLower, MorrisCRUpper, 
                     HedgesDiff, HSDiff, HSkDiff, MorrisDiff, HedgesWide, HSWide, HSkWide, MorrisWide,
                     UnitDiff, UnitWide, Hedgestau, HStau, HSktau, Morristau, Unittau)
Output


# THIS WILL EXPORT RESULTS TO YOUR DESKTOP 
# CHANGE seanpotter TO CORRECT USERNAME 
# Make sure to set col.names to FALSE after running first cell to prevent headers from being repeatedly added
# write.table(Output, file = "specific to your computer", row.names=F, append=T, col.names=F, sep=",")

