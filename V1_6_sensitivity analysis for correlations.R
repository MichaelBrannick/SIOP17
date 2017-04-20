###########################################################################################################
###########################################################################################################
### Version 1.6
### Program name: V1_6_sensitivity analysis for correlations
### Authors: Sven Kepes, Frank Bosco, Jamie Field, Mike McDaniel
### Version 1.6 (January 26, 2017)
### - Use separate syntax file with source command to run this syntax file
### - Now anticipate having 3 titles passed:
###   Output.Title1, Output.Title2, Output.Title3
### - Re-organized all analyses (moved soem analyses to a different syntax file)
###   - Naive random-effects meta-analysis (line ~101)
###   - One-sample removed analysis (line ~126)
###   - Trim and fill (line ~139)
###     - Fixed-effects trim and fill (line ~144)
###     - Random-effects trim and fill (line ~161)
###   - Contour-enhanced funnel plot with FE t&f imputation (line ~180)
###   - A priori selection models (line ~194)
###   - PET-PEESE (line ~592)
###   - Cumulative meta-analysis by precision (line ~608)
###   - Probability of excess significance (line ~643)
###   - Comprehensive results table (line ~907)
###########################################################################################################
### Version 1.5 (July 26, 2016)
### Included analyses:
### - Naive random-effects meta-analysis
### - One-sample removed analysis
### - Trim and fill
###   - Fixed-effects trim and fill
###   - Random-effects trim and fill
### - Contour-enhanced funnel plot with FE t&f imputation
### - A priori selection models from Vevea and Woods (2005)
### - PET-PEESE (line ~568)
### - Cumulative meta-analysis by sample size
### - Viechtbauer and Cheung's (2010) influence diagnostics
### - Probability of Excess Significance  (New in Version 1.4)
### - Creating a new dataframe without the identified outliers
### - Comprehensive results table
### - Creating a dataset without the identified outliers
###########################################################################################################
###########################################################################################################
### Summary of version 1.5 changes
###  1) Permit a title at the top of the output to identify the distribition analyzed.
###  2) Install and load library if needed for metafor and pwr packages
###  3) Remove outliers on the belief that such analyses are often itertaive
###     and would benefit from human attention.
###  4) Add note at end that results are in clipboard and may be pasted into spreadsheet.
###  5) Results to clipboard are only those rounded to two decimals
###  6) Used sprintf("%4.2f", varname) to get two decimal places rather than rounding
###     which truncated number to one decimal point when second digit is zero
###########################################################################################################
###########################################################################################################
### Summary of version 1.4 changes 
###  1) Reading of data file is no longer part of the code.  A datafarme named dat must exist and contain
###     the variables r and n.
###  2) P-TES (probability of excess statistics has been added)
###  3) Names for the output are made more user friendly
###  4) Not yet but soon we an read the code with a source command
###########################################################################################################
###########################################################################################################
### Hope for the future:
###  1) tune program into a series of small functions called by a main function
###  2) N with outlier(s) removed 
###  3) need significance test for Q
###########################################################################################################
###########################################################################################################


###########################################################################################################
###########################################################################################################
### START OF SYNTAX
###########################################################################################################
###########################################################################################################

#####################################################################
### make sure needed libraries are loaded; load, or install and load, needed libraries
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(metafor,pwr)

###########################################################################################################
###########################################################################################################
### Some calculations to set-up the dataset
###########################################################################################################

Output.Title1 <- ifelse(is.na(Output.Title1) == TRUE, "No Output.Title1 Provided", Output.Title1)
Output.Title2 <- ifelse(is.na(Output.Title2) == TRUE, "No Output.Title2 Provided", Output.Title2)
Output.Title3 <- ifelse(is.na(Output.Title3) == TRUE, "No Output.Title3 Provided", Output.Title3)

### Add z (yi) and the variance of z (vi) to the dataset (needed for calculations; e.g., for the selection models)
dat <- escalc(measure="ZCOR", ri=r, ni=N, data = dat)
dat$sei <- sqrt(dat$vi)

### Sorting the dataset by N (descending)
dat <- dat[order(-dat$N),] 
### Calculating cumulative N values
dat$cumsumN <- cumsum(dat$N)
dat$cumsumN <- paste0(dat$N," (Ncum = ", dat$cumsumN,")")


###########################################################################################################
###########################################################################################################
### Meta-analysis
###########################################################################################################

### Runing a naive meta on z scores. Back-transform the results to r.
res <- rma(measure="ZCOR", ri=r, ni=N, data=dat, level=95, method="DL", slab = dat$cumsumN)
res_raw_r <- rma(measure="COR", ri=r, ni=N, data=dat, level=95, method="DL")
res
OUTPUT.naive.RE.k <- as.numeric(nrow(dat))
OUTPUT.naive.RE.N <- sum(dat$N)
res.ztor.for.ci <- predict(res, level=95, transf=transf.ztor)
OUTPUT.naive.RE.mean.r <- res.ztor.for.ci$pred
OUTPUT.naive.RE.95CI.lower <- as.numeric(res.ztor.for.ci$ci.lb)
OUTPUT.naive.RE.95CI.upper <- as.numeric(res.ztor.for.ci$ci.ub)
res.ztor.for.pi <- predict(res, level=90, transf=transf.ztor)
OUTPUT.naive.RE.90PI.lower <- as.numeric(res.ztor.for.pi$cr.lb)
OUTPUT.naive.RE.90PI.upper <- as.numeric(res.ztor.for.pi$cr.ub)
OUTPUT.naive.RE.QE <- res$QE
OUTPUT.naive.RE.I2 <- res$I2
OUTPUT.naive.RE.tau <- res$tau2^.5

###End of the naive meta-analysis syntax###############################


###########################################################################################################
###########################################################################################################
### One-sample removed analysis (osr; aka leave-one-out analysis or loo)
###########################################################################################################

re.osr <- leave1out(res)
OUTPUT.osr.min <- transf.ztor(min(re.osr$estimate))
OUTPUT.osr.max <- transf.ztor(max(re.osr$estimate))
OUTPUT.osr.median <- transf.ztor(median(re.osr$estimate))

###End of one-sample removed syntax###############################


###########################################################################################################
###########################################################################################################
### Trim and fill analysis
###########################################################################################################

### Fixed-effects trim and fill
res.tf.fe <- rma(measure="ZCOR", ri=r, ni=N, data=dat, method="FE")     ### fit FE model
res.tf.fe <- trimfill(res.tf.fe, estimator="L0")                        ### trim and fill based on FE model
fill <- res.tf.fe$fill
OUTPUT.TF.FE.side <- res.tf.fe$side
OUTPUT.TF.FE.k.removed <- as.numeric(res.tf.fe$k0)
res.tf.fe <- rma(res.tf.fe$yi, res.tf.fe$vi, method="DL")             ### fit RE model based on imputed data
res.tf.fe$fill <- fill
class(res.tf.fe) <- c("rma.uni.trimfill", "rma")

res.tf.fe
res.tf.fe.ztor <- predict(res.tf.fe, level=95, transf=transf.ztor)
OUTPUT.TF.FE.mean.r <- as.numeric(res.tf.fe.ztor$pred)
OUTPUT.TF.FE.95CI.lower <- as.numeric(res.tf.fe.ztor$ci.lb)
OUTPUT.TF.FE.95CI.upper <- as.numeric(res.tf.fe.ztor$ci.ub)


### Random-effects trim and fill
res.tf.re <- rma(measure="ZCOR", ri=r, ni=N, data=dat, level=95, method = "DL")
res.tf.re <- trimfill(res.tf.re, estimator="L0")
fill <- res.tf.re$fill
OUTPUT.TF.RE.side <- res.tf.re$side
OUTPUT.TF.RE.k.removed <- as.numeric(res.tf.re$k0)
res.tf.re <- rma(res.tf.re$yi, res.tf.re$vi, method="DL")             ### fit RE model based on imputed data
res.tf.re$fill <- fill
class(res.tf.re) <- c("rma.uni.trimfill", "rma")

res.tf.re
res.tf.re.ztor <- predict(res.tf.re, level=95, transf=transf.ztor)
OUTPUT.TF.RE.mean.r <- as.numeric(res.tf.re.ztor$pred)
OUTPUT.TF.RE.95CI.lower <- as.numeric(res.tf.re.ztor$ci.lb)
OUTPUT.TF.RE.95CI.upper <- as.numeric(res.tf.re.ztor$ci.ub)

###End of the trim and fill syntax###############################


###########################################################################################################
###########################################################################################################
### Contour-enhanced funnel plot (with trim and fill imputations)
###########################################################################################################

### Contour-enhanced funnel plot with fixed-effects trim and fill imputation
funnel(res.tf.fe, yaxis="seinv", level=c(90, 95), shade=c("white","darkgray"), refline = 0)

### Contour-enhanced funnel plot with random-effects trim and fill imputation
#funnel(res.tf.fe, yaxis="seinv", level=c(90, 95), shade=c("white","darkgray"), refline = 0)

###End of the contour-enhanced funnel plot syntax###############################


#########################################################################################################
###########################################################################################################
### A priori selection models (with moderate and severe assumptions of publication bias
### Following Vevea and Woods (2005) a priory selection models
### Reference: Vevea, J. L., & Woods, C. M. (2005). Publication bias in research synthesis: Sensitivity analysis using a priori weight functions. Psychological Methods, 10, 428-443. doi: 10.1037/1082-989X.10.4.428
### Syntax is adapted from Andy Field's code, available at http://www.statisticshell.com/meta_analysis/pub_bias_r.R
###########################################################################################################

options(warn=-1)
require(foreign)
myT <- dat$yi
n <- length(myT)
v <- dat$vi
sezr<-dat$sei

npred <- 0
X <- matrix(c(rep(1,times=n)),n,npred+1,byrow=T)

if(length(v) != n) stop("Number of conditional variances must equal number of effects.")
for(i in 1:length(v)) if(v[i] <0) stop("Variances must be non-negative.")
s <- sqrt(v)

if(dim(X)[1] != n) stop("Number of rows in predictor matrix must match number of effects.")

#
# Set p = number of predictors
p <- dim(X)[2]-1

#
# In this example, cutpoints and weights are for the condition denoted 
# "severe two-tailed selection" in the paper.

# Enter cutpoints for p-value intervals
a <- c(.005,.010,.050,.100,.250,.350,.500,.650,.750,.900,.950,.990,.995,1.00)
#
# Enter fixed weight function
w1 <- matrix(
  c(1.0,.99,.95,.90,.80,.75,.65,.60,.55,.50,.50,.50,.50,.50),
  ncol=1)
w2 <- matrix(
  c(1.0,.99,.90,.75,.60,.50,.40,.35,.30,.25,.10,.10,.10,.10),
  ncol=1)
w3 <- matrix(
  c(1.0,.99,.95,.90,.80,.75,.60,.60,.75,.80,.90,.95,.99,1.0),
  ncol=1)
w4 <- matrix(
  c(1.0,.99,.90,.75,.60,.50,.25,.25,.50,.60,.75,.90,.99,1.0),
  ncol=1)

#
# Do not modify below this point
#
# Set k = number of intervals
k <- length(a)

if(length(w1) != k) stop("Number of weights must match number of intervals")
if(length(w2) != k) stop("Number of weights must match number of intervals")
if(length(w3) != k) stop("Number of weights must match number of intervals")
if(length(w4) != k) stop("Number of weights must match number of intervals")

#
# Do weighted least squares. (DonÃ???t do this at home; the algorithm is unstable
# for general regression problems, but itÃ???s easy and unlikely to cause problems
# in this context.)
varmat <- diag(v)
tmat <- matrix(myT,nrow=n)
beta <- solve( t(X)%*%solve(varmat)%*%X)%*%t(X)%*%solve(varmat)%*%tmat

#
# Set starting value for parameter vector
params <- beta

#
# Here we set up Bij and bij
#
# First, we need predicted values
xb <- X%*%beta
sig <- sqrt(v)

#
# Now bij = -si probit(aj)
b <- matrix(rep(0,n*(k-1)),nrow=n,ncol=k-1)
for(i in 1:n) 
  for(j in 1:k-1) 
    b[i,j] <- -s[i] * qnorm(a[j])

#
# And now Bij (Equation 5)
Bij <- function(params) {
  B <- matrix(rep(0,n*k),nrow=n,ncol=k)
  beta <- matrix(params[1:(p+1)],ncol=1)
  xb <- X%*%beta
  for(i in 1:n)
    for(j in 1:k) {
      if(j==1) B[i,j] <- 1.0 - pnorm( (b[i,1]-xb[i]) / sig[i])
      else if(j<k) B[i,j] <- pnorm( (b[i,j-1]-xb[i]) / sig[i]) - pnorm( (b[i,j]-xb[i]) / sig[i])
      else B[i,j] <- pnorm( (b[i,k-1]-xb[i]) / sig[i])
    }
  return(B)
}

#
# Unadjusted likelihood function
like1 <- function(params) {
  beta <- matrix(params[1:(p+1)],ncol=1)
  xb <- X%*%beta
  -2*(-1/2 * (sum( (myT-xb)^2/v ) + sum(log(v)))   )
}

#
# adjusted likelihood function Moderate One-Tailed Selection (Equation 6)
like2 <- function(params) {
  beta <- matrix(params[1:(p+1)],ncol=1)
  xb <- X%*%beta
  ll <- -1/2 * (sum( (myT-xb)^2/v ) + sum(log(v)))   
  B <- Bij(params)
  Bw <- B%*%w1
  ll <- ll - sum(log(Bw))
  return(-2*ll)
}

#
# adjusted likelihood function Severe One-Tailed Selection (Equation 6)
like3 <- function(params) {
  beta <- matrix(params[1:(p+1)],ncol=1)
  xb <- X%*%beta
  ll <- -1/2 * (sum( (myT-xb)^2/v ) + sum(log(v)))   
  B <- Bij(params)
  Bw <- B%*%w2
  ll <- ll - sum(log(Bw))
  return
  
}

#
# adjusted likelihood function Moderate Two-Tailed Selection (Equation 6)
like4 <- function(params) {
  beta <- matrix(params[1:(p+1)],ncol=1)
  xb <- X%*%beta
  ll <- -1/2 * (sum( (myT-xb)^2/v ) + sum(log(v)))   
  B <- Bij(params)
  Bw <- B%*%w3
  ll <- ll - sum(log(Bw))
  return(-2*ll)
}

#
# adjusted likelihood function Severe Two-Tailed Selection (Equation 6)
like5 <- function(params) {
  beta <- matrix(params[1:(p+1)],ncol=1)
  xb <- X%*%beta
  ll <- -1/2 * (sum( (myT-xb)^2/v ) + sum(log(v)))   
  B <- Bij(params)
  Bw <- B%*%w4
  ll <- ll - sum(log(Bw))
  return(-2*ll)
}

#
# Second derivatives, unadjusted model 
fese <- function(params) {
  beta <- matrix(params[1:(p+1)],ncol=1)
  xb <- X%*%beta
  hessian <- matrix(rep(0,(p+1)*(p+1)),p+1,p+1)
  for(i in 1:n) {
    for(j in 1:(p+1))
      for(m in 1:(p+1)) 
        hessian[j,m] <- hessian[j,m] - X[i,j]*X[i,m]/v[i]
  }
  temp <- solve(hessian)
  retval <- rep(0,(p+1))
  for(i in 1:(p+1)) retval[i] <- sqrt(-temp[i,i])
  return(retval)
}

# Fit unweighted model by maximum likelihood
fe1 <- nlminb(objective=like1, start=params)
#se1 <- fese(fe1$par)
fe1$par[2] <- (exp(2*fe1$par[1])-1)/(exp(2*fe1$par[1])+1)

#
# Fit adjusted model Moderate one-Tailed Selection
fe2 <- nlminb(objective=like2, start=fe1$par)
fe2$par[2] <- (exp(2*fe2$par[1])-1)/(exp(2*fe2$par[1])+1)

#
# Fit adjusted model Severe one-Tailed Selection
fe3 <- nlminb(objective=like3, start=fe1$par)
fe3$par[2] <- (exp(2*fe3$par[1])-1)/(exp(2*fe3$par[1])+1)

#
# Fit adjusted model Moderate Two-Tailed Selection
fe4 <- nlminb(objective=like4, start=fe1$par)
fe4$par[2] <- (exp(2*fe4$par[1])-1)/(exp(2*fe4$par[1])+1)

#
# Fit adjusted model Severe Two-Tailed Selection
fe5 <- nlminb(objective=like5, start=fe1$par)
fe5$par[2] <- (exp(2*fe5$par[1])-1)/(exp(2*fe5$par[1])+1)

#
# Also, set variance component to balanced method of moments start
rss <- sum((myT-X%*%beta)^2)
vc <- rss/(n-p-1) - mean(v)
if (vc < 0.0) vc <- 0.0

#
# Set starting value for parameter vector
params <- c(beta,vc)

#
# Here we set up Bij and bij
#
# First, we need predicted values
xb <- X%*%beta
sig <- sqrt(v + vc)

#
# Now bij = -si probit(aj)
b <- matrix(rep(0,n*(k-1)),nrow=n,ncol=k-1)
for(i in 1:n) 
  for(j in 1:k-1) 
    b[i,j] <- -s[i] * qnorm(a[j])

#
# And now Bij (Equation 5)
Bij <- function(params) {
  B <- matrix(rep(0,n*k),nrow=n,ncol=k)
  beta <- matrix(params[1:(p+1)],ncol=1)
  xb <- X%*%beta
  for(i in 1:n)
    for(j in 1:k) {
      if(j==1) B[i,j] <- 1.0 - pnorm( (b[i,1]-xb[i]) / sig[i])
      else if(j<k) B[i,j] <- pnorm( (b[i,j-1]-xb[i]) / sig[i]) - pnorm( (b[i,j]-xb[i]) / sig[i])
      else B[i,j] <- pnorm( (b[i,k-1]-xb[i]) / sig[i])
    }
  return(B)
}

#
# Unadjusted likelihood function
like1 <- function(params) {
  vc <- params[p+2]
  beta <- matrix(params[1:(p+1)],ncol=1)
  xb <- X%*%beta
  vall <- v+vc
  -2*(-1/2 * (sum( (myT-xb)^2/vall ) + sum(log(vall)))   )
}

#
# Adjusted likelihood function Moderate One-Tailed Selection (Equation 6)
like2 <- function(params) {
  vc <- params[p+2]
  beta <- matrix(params[1:(p+1)],ncol=1)
  xb <- X%*%beta
  vall <- v+vc
  ll <- -1/2 * (sum( (myT-xb)^2/vall ) + sum(log(vall)))   
  B <- Bij(params)
  Bw <- B%*%w1
  ll <- ll - sum(log(Bw))
  return(-2*ll)
}

#
# Adjusted likelihood function Severe One-Tailed Selection (Equation 6)
like3 <- function(params) {
  vc <- params[p+2]
  beta <- matrix(params[1:(p+1)],ncol=1)
  xb <- X%*%beta
  vall <- v+vc
  ll <- -1/2 * (sum( (myT-xb)^2/vall ) + sum(log(vall)))   
  B <- Bij(params)
  Bw <- B%*%w2
  ll <- ll - sum(log(Bw))
  return(-2*ll)
}

#
# Adjusted likelihood function Moderate two-Tailed Selection (Equation 6)
like4 <- function(params) {
  vc <- params[p+2]
  beta <- matrix(params[1:(p+1)],ncol=1)
  xb <- X%*%beta
  vall <- v+vc
  ll <- -1/2 * (sum( (myT-xb)^2/vall ) + sum(log(vall)))   
  B <- Bij(params)
  Bw <- B%*%w3
  ll <- ll - sum(log(Bw))
  return(-2*ll)
}

#
# Adjusted likelihood function Severe Two-Tailed Selection (Equation 6)
like5 <- function(params) {
  vc <- params[p+2]
  beta <- matrix(params[1:(p+1)],ncol=1)
  xb <- X%*%beta
  vall <- v+vc
  ll <- -1/2 * (sum( (myT-xb)^2/vall ) + sum(log(vall)))   
  B <- Bij(params)
  Bw <- B%*%w4
  ll <- ll - sum(log(Bw))
  return(-2*ll)
}

#
# Second derivatives, unadjusted model 
rese <- function(params) {
  vc <- params[p+2]
  beta <- matrix(params[1:(p+1)],ncol=1)
  xb <- X%*%beta
  vall <- v+vc
  hessian <- matrix(rep(0,(p+2)*(p+2)),p+2,p+2)
  for(i in 1:n) {
    for(j in 1:(p+1))
      for(m in 1:(p+1)) 
        hessian[j,m] <- hessian[j,m] - X[i,j]*X[i,m]/vall[i]
    for(j in 1:(p+1)) {
      hessian[j,p+2] <- hessian[j,p+2] - (myT[i]-xb[i])/vall[i]/vall[i]
      hessian[p+2,j] <- hessian[p+2,j] - (myT[i]-xb[i])/vall[i]/vall[i]
    }
    hessian[p+2,p+2] <- hessian[p+2,p+2] + 1/vall[i]/vall[i]/2 - 
      (myT[i]-xb[i])*(myT[i]-xb[i])/vall[i]/vall[i]/vall[i]
  }
  temp <- solve(hessian)
  retval <- rep(0,(p+2))
  for(i in 1:(p+2)) retval[i] <- sqrt(-temp[i,i])
  return(retval)
}

#
# Fit unadjusted model by maximum likelihood
re1 <- nlminb(objective=like1, start=params,lower=c(rep(-Inf,p+1),0.0))
#se1 <- rese(re1$par)
re1$par[3] <- (exp(2*re1$par[1])-1)/(exp(2*re1$par[1])+1)

#
# Fit adjusted model Moderate one-Tailed Selection
re2 <- nlminb(objective=like2, start=re1$par,lower=c(rep(-Inf,p+1),0.0))
re2$par[3] <- (exp(2*re2$par[1])-1)/(exp(2*re2$par[1])+1)

#
# Fit adjusted model Severe one-Tailed Selection
re3 <- nlminb(objective=like3, start=re1$par,lower=c(rep(-Inf,p+1),0.0))
re3$par[3] <- (exp(2*re3$par[1])-1)/(exp(2*re3$par[1])+1)

#
# Fit adjusted model Moderate Two-Tailed Selection
re4 <- nlminb(objective=like4, start=re1$par,lower=c(rep(-Inf,p+1),0.0))
re4$par[3] <- (exp(2*re4$par[1])-1)/(exp(2*re4$par[1])+1)

#
# Fit adjusted model Severe Two-Tailed Selection
re5 <- nlminb(objective=like5, start=re1$par,lower=c(rep(-Inf,p+1),0.0))
re5$par[3] <- (exp(2*re5$par[1])-1)/(exp(2*re5$par[1])+1)

unadjfe <- rbind(unadjfe1 <- fe1$par, unadjfe2 <- fese(fe1$par))
unadjfe[2,2]=NA
rownames(unadjfe) <- c("Parameter Estimates", "Standard errors                  ")
colnames(unadjfe) <- c("zr", "r")

unadjre <- rbind(unadjre1 <- re1$par, unadjre2 <- c(rese(re1$par),NA))
rownames(unadjre) <- c("Parameter Estimates", "Standard errors                  ")
colnames(unadjre) <- c("zr", "v", "r")

adjfe <- rbind(adj1 <- fe2$par, adj2 <- fe3$par, adj3 <- fe4$par, adj4 <- fe5$par)
rownames(adjfe) <- c("Moderate One-Tailed Selection    ", "Severe One-Tailed Selection", "Moderate Two-Tailed Selection", "Severe Two-Tailed Selection")
colnames(adjfe) <- c("zr", "r")

adjre <- rbind(adj1 <- re2$par, adj2 <- re3$par, adj3 <- re4$par, adj4 <- re5$par)
rownames(adjre) <- c("Moderate One-Tailed Selection    ", "Severe One-Tailed Selection", "Moderate Two-Tailed Selection", "Severe Two-Tailed Selection")
colnames(adjre) <- c("zr", "v", "r")

#funnel(myT, sezr, xlim=NULL, ylim=NULL, xlab="Effect Size (Zr)", ylab="Standard Error",comb.f=FALSE, axes=TRUE, pch=1, text=NULL, cex=1.5, col=1,log="", yaxis="se", sm=NULL,level=.95)

out1 <- paste("**********  SENSITIVITY OF EFFECT-SIZE ESTIMATES TO PUBLICATION BIAS  **********")
out1 <- paste(out1, "EFFECT-SIZE PARAMETER:   Correlation")
out1 <- paste(out1, "EFFECT-SIZE ESTIMATOR:   r ")     
out1 <- paste(out1, "FIXED-EFFECTS Publication Bias Model: Vevea & Woods (2005), Psychological Methods") 
out1 <- paste(out1, "Unadjusted Parameter Estimate")
out2 <- paste("Adjusted Parameter Estimate")
out3 <- paste("RANDOM-EFFECTS Publication Bias Model: Vevea & Woods (2005), Psychological Methods")
out3 <- paste(out3, "In this model v estimates population effect-size variance")
out3 <- paste(out3, "Unadjusted Parameter Estimates")
out4 <- paste("Adjusted Parameter Estimates")

cat(out1); unadjfe; cat(out2); adjfe; cat(out3); unadjre; cat(out4); adjre

OUTPUT.SelMod.1tmod.est.zr <- adjre[1, 1]
OUTPUT.SelMod.1tmod.est.var <- adjre[1, 2]
OUTPUT.SelMod.1tmod.est.r <- adjre[1, 3]
OUTPUT.SelMod.1tsev.est.zr <- adjre[2, 1]
OUTPUT.SelMod.1tsev.est.var <- adjre[2, 2]
OUTPUT.SelMod.1tsev.est.r <- adjre[2, 3]

###End of the selection model syntax###############################


###########################################################################################################
###########################################################################################################
### PET-PEESE analysis
### Based on Stanley and Doucouliagos (2014)
### Reference: Stanley, T. D., & Doucouliagos, H. (2014). Meta-regression approximations to reduce publication selection bias. Research Synthesis Methods, 5, 60-78. doi: 10.1002/jrsm.1095
###########################################################################################################

### Calculate se of correlation, variance, and precisison, and precision_sq
dat$FisherZ <-  0.5 * log((1 + dat$r) / (1 - dat$r))
dat$FisherZSE <- 1 / (sqrt(dat$N - 3))
dat$r_se <- (1 - (dat$r^ 2)) * dat$FisherZSE

dat$r_variance <- dat$r_se^2
dat$r_precision <- 1/dat$r_se
dat$r_precision_sq <- dat$r_precision^2

reg1 <- lm(r ~ r_se, weights = r_precision_sq, data=dat)
summary(reg1)
OUTPUT_PET <- reg1$coefficients[1]

### Intercept
OUTPUT_PET_pval <- (coef(summary(reg1))[1, 4]) / 2

reg2 <- lm(r ~ r_variance, weights = r_precision_sq, data=dat)

OUTPUT_PEESE <- reg2$coefficients[1]

PETPEESEFINAL <- ifelse(OUTPUT_PET_pval < .05, OUTPUT_PEESE, OUTPUT_PET)

###End of the PET-PEESE syntax###############################


###########################################################################################################
###########################################################################################################
### Cumulative meta-analysis by precision, sorting by var (~sample size)
###########################################################################################################

res.cma <- cumul(res, order=order(dat$vi))
res.cma
forestOUTPUT <- forest(res.cma, transf=transf.ztor)

### Estimate of the five most precise (i.e., largest) samples (in the r-metric)
res.cma <- cumul(res, order=order(dat$vi), transf=transf.ztor)
OUTPUT.CMA.5thEstimateTemp <- res.cma[5,0]
OUTPUT.CMA.5thEstimate <- OUTPUT.CMA.5thEstimateTemp$estimate

class(OUTPUT.CMA.5thEstimate)

###End of the cumulative meta-analysis syntax###############################


###########################################################################################################
###########################################################################################################
### Test of excess significance (for correlation effect sizes)
### Based on Ioannidis and Trikalinos (2007) and Francis (2014)
### References: Ioannidis, J. P., & Trikalinos, T. A. (2007). An exploratory test for an excess of significant findings. Clinical Trials, 4, 245-253. doi: 10.1177/1740774507079441 
###             Francis, G. (2014). The frequency of excess success for articles in Psychological Science. Psychonomic Bulletin & Review, 21, 1180-1187. doi: 10.3758/s13423-014-0601-x
###########################################################################################################

dat$n <- dat$N
rho <- OUTPUT.naive.RE.mean.r
rho
#table(dat$n)
#str(dat)

Statistical.Sig.r <- function(r = dat$r, n = dat$n)
{     
  # calculating statistical significance with both a t and a z test. reporting the z test
  # get t and and its significance value
  t <- abs(r)/ sqrt((1 - (r*r))/ (n - 2))
  pval.t <- ((1 - (pt(t, 20)))) * 2   
  # convert R to Fisher z, get the z test statistic, get the significance value  
  FisherZ = 0.5 * log((1 + r) / (1 - r))
  FisherZ
  FisherZSE = 1 / (sqrt(n - 3))
  z.test.statistic <- FisherZ/FisherZSE
  z.test.statistic
  pval.z <- 2*(1- pnorm(abs(z.test.statistic)))
  pval.z
  
  return(pval.z)
}

# call statistical significance function to add the statistical significance for each r to the dataframe
dat$stat.sig <- Statistical.Sig.r(r = dat$r, n = dat$n)
dat$stat.sig


# create a dichotomous variable in the dataframe that equals 1 if the significance value <- to .05
dat$significant <- 0  #initialize
dat$significant[dat$stat.sig <= .05] <- 1
dat$significant

# count the number of statistically significant correlations
observed.significant = sum(dat$significant)
observed.significant

# calculate the post-hoc power of each study. rho is the point estimate from the meta-analysis. p = .05
power <-pwr.r.test(n = dat$n, r = rho, sig.level = .05, power = NULL, 
                   alternative = c("two.sided", "less","greater"))
str(power)
#add the power of each study to the dataframe
dat$power <- power$power 
dat$power


# calculate the expected number of significant studies
expected.significant = sum(dat$power)
expected.significant

# define k as the number of rows in the dataframe
k = nrow(dat)

#print the k, observed.significant, and  expected.significant
#k
#observed.significant 
#expected.significant


# calculate chi-square p 246 Ioannidis & Trikalinos (2007)
chi.square <-   (((observed.significant - expected.significant)^2) / expected.significant) + 
  (((observed.significant - expected.significant)^2) / (k - expected.significant))
chi.square

# Chi-square statistical significance
p.sig.chi.square <-  1 - pchisq(chi.square,1) 
p.sig.chi.square


# binomial test for 10 or greater
mean.power = mean(dat$power)
mean.power
# probability vector contains the probabilities of obtaining a significant value for the
#      observed through the total k
prob.vector <- dbinom(x = observed.significant:k, size = k, prob = mean.power)
prob.vector

# a loop to sum the probability vector
binomial.probability <- 0  #initialize the counter
binomial.probability
# get length of vector
length.of.prob.vector <- length (prob.vector)
length.of.prob.vector

for (i in 1:length.of.prob.vector)
{
  print(i)
  print(binomial.probability)
  binomial.probability <- binomial.probability + prob.vector[i]
}
binomial.probability


#fetch Greg Francis Code for Fisher z test
#function is called TestSetProbability()
#install.packages("combinat")
#source("C:/Users/Mike McDaniel/Dropbox/Active McDaniel Research/Excess significance/TestSetProbabilty.R")
#library(combinat)
#This function runs Fisher's exact test to compute the probability of getting the observed or more rejections
# of the null hypothesis from a set of experiments with the indicated power values.

# Greg Francis (gfrancis@purdue.edu)
# 17 November 2012

TestSetProbabilty <- function(Observed, Power){
  
  numExperiments = sum(!is.na(Power)) 
  
  sumProb = 0
  
  # Use binomial distribution if power values are identical for all experiments in the set
  if(var(Power)==0){
    test <- binom.test(Observed, numExperiments, p = Power[1], alternative="greater")
    sumProb = test$p.value
  }
  else {  # have to do a Fisher's exact test or a Chi Square test for varying power values
    
    # Fisher's exact test for small enough number of experiments
    if(Observed >0 && Observed <numExperiments && numExperiments <=15)
    {  
      # for Observed and up  
      if(Observed>= (numExperiments/2)){
        BetaValues = 1-Power
        prodBeta = prod(BetaValues)
        sumProb = 0
        for(chkObsv in Observed:numExperiments){  		
          # get all possible combinations of rejections
          patterns = combn(numExperiments, chkObsv)
          numPatterns = length(patterns)/chkObsv
          
          # go through each pattern and compute probability of that pattern
          sumchkObsvProb = 0
          if(chkObsv>0){
            for(i in 1:numPatterns){
              probRejectsY = 1
              probRejectsN =1
              probNReject =1
              # probability
              for(j in 1:chkObsv) {
                # product of powers for those that reject
                probRejectsY = probRejectsY * Power[patterns[(i-1)* chkObsv +j]]
                #		cat(patterns[(i-1)* chkObsv +j], "\t", Power[patterns[(i-1)* chkObsv +j]], "\n")
                # product of Beta for those that do not reject
                probNReject = probNReject * BetaValues[patterns[(i-1)* chkObsv +j]]
              }
              #	cat("------\n")
              # probability for entire set pattern
              if(probNReject==0){probNReject=1} # in case all experiments reject
              sumchkObsvProb = sumchkObsvProb + prodBeta*probRejectsY/probNReject		
            }
          }
          sumProb = sumProb + sumchkObsvProb
        }
      }
      else {  # compute observed down and then subtract from one
        BetaValues = 1-Power
        prodBeta = prod(BetaValues)
        sumProb = 0
        for(chkObsv in 0:(Observed-1)){	
          
          if(chkObsv>0){  # 0 is a special case, handled below		
            # get all possible combinations of rejections
            patterns = combn(numExperiments, chkObsv)
            numPatterns = length(patterns)/chkObsv
            
            # go through each pattern and compute probability of that pattern
            sumchkObsvProb = 0
            
            for(i in 1:numPatterns){
              probRejectsY = 1
              probRejectsN =1
              probNReject =1
              # probability
              for(j in 1:chkObsv) {
                # product of powers for those that reject
                probRejectsY = probRejectsY * Power[patterns[(i-1)* chkObsv +j]]
                #		cat(patterns[(i-1)* chkObsv +j], "\t", Power[patterns[(i-1)* chkObsv +j]], "\n")
                # product of Beta for those that do not reject
                probNReject = probNReject * BetaValues[patterns[(i-1)* chkObsv +j]]
              }
              #	cat("------\n")
              # probability for entire set pattern
              if(probNReject==0){probNReject=1} # in case all experiments reject
              sumchkObsvProb = sumchkObsvProb + prodBeta*probRejectsY/probNReject		
            }
          }
          else {  # special case when none reject
            
            sumchkObsvProb = prod(BetaValues)
          }
          sumProb = sumProb + sumchkObsvProb
        }		
        sumProb = 1- sumProb		
      }
    }
    else if (Observed == numExperiments)
    {
      sumProb = prod(Power)
    }
    else if (Observed == 0)
    {
      sumProb = 1
    }
    else # Chi Square test for more than 15 experiments
    {
      Expected = sum(Power)
      
      A = ((Observed - Expected)^2)/(Expected) + ( (Observed-Expected)^2)/(numExperiments-Expected)
      
      
      sumProb <- 1-pchisq(A,df=1)
      
      cat("O=", Observed, " E=", Expected, " A=", A, " p=", sumProb, "\n")
      
      # Only looking for over-publication, not under-publication
      if(Observed < Expected){
        sumProb = 1
      }
    }
  }
  
  
  return(sumProb)
}
# End of function of TestSetProbabilty 

# This function runs Fisher's exact test to compute the probability of getting the observed or more rejections
# of the null hypothesis from a set of experiments with the indicated power values.
# str(power)
power2 <-dat$power
#str(power2)
fisher.exact.test.probability <- TestSetProbabilty(observed.significant, power2)
fisher.exact.test.probability

# print to console
cat("             Summary of Excess Significance Analysis", "\n")
cat("Estimate of Rho", rho, "\n")
cat("Number of Correlations: ", k, "\n")
cat("Number of statistically significant r: ", observed.significant, "\n")
cat("Expected number of statistically significant r: ", expected.significant, "\n")
cat("Chi-square test probability: ", p.sig.chi.square, "\n")
cat("Binomial test probability: ", binomial.probability, "\n")
cat("Fisher exact test probability: ", fisher.exact.test.probability, "\n")

# OUTPUT.SelMod.1tsev.est.var <- adjre[2, 2]
# OUTPUT.SelMod.1tsev.est.r <- adjre[2, 3]
OUTPUT.PTES.Number.of.Observed.Significant <- observed.significant
OUTPUT.PTES.Number.of.Expected.Significant <- expected.significant
OUTPUT.PTES.Chi.square.test.probability <- p.sig.chi.square
OUTPUT.PTES.Binomial.test.probability <- binomial.probability
OUTPUT.PTES.Fisher.exact.test.probability <- fisher.exact.test.probability

##########  End of excess significance analysis ##########################


###########################################################################################################
###########################################################################################################
### Creating the results table
###########################################################################################################


output.table.names <- c("Title1:","Title2:","Title3:", 
                        "k (number of effect sizes)", "N (Cumulative samples size)", 
                        "Random effects mean effect size estimate", 
                        "95% Confidence interval lower value", "95% Confidence interval lower value ", 
                        "90% Prediction interval lower value", "90% Prediction interval upper value",
                        "Q statistic", "I squared", "Tau",
                        "\n",
                        "**Leave One Out (LOO) results **",
                        "Minimum random effects LOO mean", "Maximum random effects LOO mean", 
                        "Median random effects LOO mean",
                        "\n",
                        "**Trim & Fill Results**", 
                        "FE Model: Side of distribution in which data were imputed", 
                        "FE Model: Number of effects imputed", 
                        "FE Model: Trim and fill adjusted mean r", 
                        "FE Model: 95% Confidence interval lower value of T&F adjusted distribution",
                        "FE Model: 95% Confidence interval upper value of T&F adjusted distribution",
                        "RE Model: Side of distribution in Which data were imputed (with FE imputation)", 
                        "RE Model: Number of effects imputed (with FE imputation)", 
                        "RE Model: Trim and fill adjusted mean r  (with FE imputation)", 
                        "RE Model: 95% Confidence interval lower value of T&F adjusted distribution",
                        "RE Model: 95% Confidence interval upper value of T&F adjusted distribution", 
                        "\n",
                        "** Cumulative Meta-Analysis by Descending N Results",
                        "RE Model estimated mean of five largest samples",
                        "\n",
                        "**Selection Model Results**",
                        "Estimated mean r with moderate bias", "Estimated mean Fisher z with moderate bias",
                        "Estimated variance of Fisher z with moderate bias)",  
                        "Estimated mean r with severe bias", "Estimated mean Fisher z with severe bias",
                        "Estimated variance of Fisher z with severe bias",  
                        "\n",
                        "**PET-PEESE Results**",
                        "PET value", "PET p value", "PEESE value", "PET-PEESE (use this value)",
                        "\n",
                        "**Probability of Excess Significance Results**",
                        "Number of Observed Statistically Significant Results",
                        "Number of Expected Statistically Significant Results",
                        "Chi-square test probability of these results",
                        "Binomial test probability of these results",
                        "Fisher exact test probability of these results")



output.table.data <- c(Output.Title1, Output.Title2,Output.Title3, 
                       OUTPUT.naive.RE.k, OUTPUT.naive.RE.N, OUTPUT.naive.RE.mean.r, OUTPUT.naive.RE.95CI.lower, OUTPUT.naive.RE.95CI.upper, 
                       OUTPUT.naive.RE.90PI.lower, OUTPUT.naive.RE.90PI.upper, OUTPUT.naive.RE.QE, OUTPUT.naive.RE.I2, OUTPUT.naive.RE.tau,
                       "\n",
                       "---------- Leave one out Results----------",
                       OUTPUT.osr.min, OUTPUT.osr.max, OUTPUT.osr.median, 
                       "\n",
                       "---------- Trim & Fill Results----------", 
                       OUTPUT.TF.FE.side, OUTPUT.TF.FE.k.removed, OUTPUT.TF.FE.mean.r, OUTPUT.TF.FE.95CI.lower, OUTPUT.TF.FE.95CI.upper, 
                       OUTPUT.TF.RE.side, OUTPUT.TF.RE.k.removed, OUTPUT.TF.RE.mean.r, OUTPUT.TF.RE.95CI.lower,OUTPUT.TF.RE.95CI.upper,
                       "\n",
                       "---------- Cumulative Meta-Analysis by Descending N Results ----------",
                       OUTPUT.CMA.5thEstimate, 
                       "\n",
                       "---------- Selection Model Results----------",
                       OUTPUT.SelMod.1tmod.est.r, OUTPUT.SelMod.1tmod.est.zr, OUTPUT.SelMod.1tmod.est.var, 
                       OUTPUT.SelMod.1tsev.est.r, OUTPUT.SelMod.1tsev.est.zr, OUTPUT.SelMod.1tsev.est.var,  
                       "\n",
                       "---------- PET-PEESE Results----------",
                       OUTPUT_PET, OUTPUT_PET_pval, OUTPUT_PEESE, PETPEESEFINAL, 
                       "\n",
                       "---------- Probability of Excess Significance Results----------",
                       OUTPUT.PTES.Number.of.Observed.Significant, 
                       OUTPUT.PTES.Number.of.Expected.Significant,
                       OUTPUT.PTES.Chi.square.test.probability,
                       OUTPUT.PTES.Binomial.test.probability,
                       OUTPUT.PTES.Fisher.exact.test.probability)

output.table.data.round2 <- c(Output.Title1, Output.Title2,Output.Title3,
                              OUTPUT.naive.RE.k, OUTPUT.naive.RE.N, 
                              sprintf("%4.2f", OUTPUT.naive.RE.mean.r), 
                              sprintf("%4.2f", OUTPUT.naive.RE.95CI.lower), 
                              sprintf("%4.2f", OUTPUT.naive.RE.95CI.upper), 
                              sprintf("%4.2f", OUTPUT.naive.RE.90PI.lower), 
                              sprintf("%4.2f", OUTPUT.naive.RE.90PI.upper), 
                              sprintf("%4.2f", OUTPUT.naive.RE.QE), 
                              sprintf("%4.2f", OUTPUT.naive.RE.I2), 
                              sprintf("%4.2f", OUTPUT.naive.RE.tau), 
                              "\n",
                              "-- Leave one out Results--",
                              sprintf("%4.2f", OUTPUT.osr.min), 
                              sprintf("%4.2f", OUTPUT.osr.max), 
                              sprintf("%4.2f", OUTPUT.osr.median),
                              "\n",
                              "-- Trim & Fill Results--", 
                              OUTPUT.TF.FE.side, OUTPUT.TF.FE.k.removed, 
                              sprintf("%4.2f", OUTPUT.TF.FE.mean.r), 
                              sprintf("%4.2f", OUTPUT.TF.FE.95CI.lower),
                              sprintf("%4.2f", OUTPUT.TF.FE.95CI.upper),
                              OUTPUT.TF.RE.side, OUTPUT.TF.RE.k.removed, 
                              sprintf("%4.2f", OUTPUT.TF.RE.mean.r), 
                              sprintf("%4.2f", OUTPUT.TF.RE.95CI.lower), 
                              sprintf("%4.2f", OUTPUT.TF.RE.95CI.upper), 
                              "\n",
                              "-- Cumulative Meta-Analysis by Descending N Results --",
                              sprintf("%4.2f", OUTPUT.CMA.5thEstimate),
                              "\n",
                              "-- Selection Model Results--",
                              sprintf("%4.2f", OUTPUT.SelMod.1tmod.est.r), 
                              sprintf("%4.2f", OUTPUT.SelMod.1tmod.est.zr),
                              sprintf("%4.2f", OUTPUT.SelMod.1tmod.est.var),  
                              sprintf("%4.2f", OUTPUT.SelMod.1tsev.est.r),
                              sprintf("%4.2f", OUTPUT.SelMod.1tsev.est.zr), 
                              sprintf("%4.2f", OUTPUT.SelMod.1tsev.est.var), 
                              "\n",
                              "-- PET-PEESE Results--",
                              sprintf("%4.2f", OUTPUT_PET), 
                              sprintf("%4.2f", OUTPUT_PET_pval), 
                              sprintf("%4.2f", OUTPUT_PEESE), 
                              sprintf("%4.2f", PETPEESEFINAL), 
                              "\n",
                              "-- Probability of Excess Significance Results--",
                              sprintf("%4.2f", OUTPUT.PTES.Number.of.Observed.Significant), 
                              sprintf("%4.2f", OUTPUT.PTES.Number.of.Expected.Significant),
                              sprintf("%4.2f", OUTPUT.PTES.Chi.square.test.probability),
                              sprintf("%4.2f", OUTPUT.PTES.Binomial.test.probability),
                              sprintf("%4.2f", OUTPUT.PTES.Fisher.exact.test.probability))

output.dataframe <- data.frame(output.table.names, output.table.data.round2)
output.dataframe

#write table to clipboard than to Excel then can move it into Word
write.table(output.dataframe,file="clipboard",sep="\t",col.names=NA)
paste("**************************************************************************")
paste("Results are copied to clipboard. Open Spreadsheet (e.g., Excel) and paste.")
paste("**************************************************************************")

###End of results table syntax###############################