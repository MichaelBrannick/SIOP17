#Sensitivity analysis of correlations.R

#Meta-Analysis of Brain Volume and Intelligence

# This program illustrates how to conduct a meta-analysis of correlation coefficients 
# in R using metafor. It does the following:
#  1. Using the pacman::p_load, loads packages, installing packages only when 
#     not already installed.
#  2. After the meta-analysis is conducted, the metafor functions for outlier detection
#     are used, outliers are deleted, and the meta-analysis is rerun.
#  3. Publication bias analyses are run on the data set containing outliers and then
#     on the data set without outliers.  These results can be copied into an Excel file.



### This next block of code is not a required step, but reflects a preference to
### clean out the global environment, plots and graphs, and the console.
### Clearing global environment, plots/graphics, and the console.
### one does not need to do this, but I like to start with an empty environment
rm(list=ls()) #Cleans (empties) global environment
dev.off()     #Cleans plots/graphics
cat("\014")   #Cleans out the console


### Use pacman to install and load packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
### the next line of code installs packages if not alreday installed and also loads the packages
### listed in the parentheses if not already loaded
pacman::p_load(readxl,metafor,meta,foreign,beepr)

### set directory so data files can be located (note the direction of the slashes)
setwd("C:/Users/skepes/Desktop/Meta-Analysis Symposium SIOP 2017")
dir()

      ####################################################################
      ####################################################################
      ####################################################################
      ### Meta-analysis of brain volume and intelligence
      ###################################################################

### Read in data from an Excel file.  Data are from: McDaniel, M.A. (2005).
### Big-brained people are smarter: A meta-analysis of the relationship between 
### in vivo brain volume and intelligence. Intelligence, 33, 337-346. 
BigBrain <- read_excel("Illustrative data.xls", sheet = "r", col_names = TRUE, col_types = NULL, na = "", skip = 0)

### Look at the data
str(BigBrain)
nrow(BigBrain)
ncol(BigBrain)

### Look at the top 6 cases
head(BigBrain)

### Create an ID number that corresponds to the data row. This ID is handy when dealing with outliers
for (i in 1:nrow(BigBrain))
  {
  BigBrain$ID[i] <- i
  }
head(BigBrain)    
    


### Metafor needs the sample size to be at least 4 every observation.  This data set
### includes an observation in whihc the sample size is 3.  Below we drop the observation.
nrow(BigBrain)  #37 observations
BigBrain <- BigBrain[BigBrain$Sample_Size > 3,]
nrow(BigBrain)  #36 observations

### Metafor analyzes correlations in a Fisher z metric
### Create effect sizes in which r is expressed as a Fisher z
### The escalc function creates an effect size yi is the Fisher z transformation of the r
### vi is created and is the sampling error variance of the Fisher z
data <- escalc(measure="ZCOR", ri=r, ni=Sample_Size, data = BigBrain)
head(data)

### Run the random effects meta-analysis without any "adjustments" (i.e., the naive meta-analysis)
### meta-analysis of the transformed correlations using a random-effects model
ma.result <- rma(yi, vi, data=data, level=95, method="DL")
ma.result
### The results shown in ma.result are in a Fisher z metric

### The predict function converts the Fisher z results into correlation coefficients
ma.result.in.r.metric <- predict(ma.result, level=95, transf=transf.ztor)
ma.result.in.r.metric

### now we look for outliers
### get outlier and influential case diagonstics
inf <- influence(ma.result)
inf   # note how observations 27 and 29 have an asterisk next to them. This indicates
      # that they are influential cases (i.e., outliers)

### Plot the diagnostics
plot(inf)

### Display the two outliers
BigBrain[c(27,29),]  # the ID numbers are 28 and 30

### Drop the two outliers and save the data file as BigBrain2
nrow(BigBrain)  #36 observations
BigBrain2 <- BigBrain[BigBrain$ID != 28 & BigBrain$ID != 30,]
nrow(BigBrain2) #34 observations

### The outlier identifier will sometimes find additional outliers if one runs
### it repeatedly.  Below we run it again on the data set BigBrain2

BigBrain2 <- escalc(measure="ZCOR", ri=r, ni=Sample_Size, data = BigBrain2)

ma.result2 <- rma(yi, vi, data=BigBrain2, level=95, method="DL")
ma.result2

inf <- influence(ma.result2)
inf   #row 1 is an outlier

plot(inf)

### Display the outlier
BigBrain2[1,]  #display the 1st row. Note that the value of ID is 1

### Drop the outlier and save the data file as BigBrain3
nrow(BigBrain2)  #34 observations
BigBrain3 <- BigBrain2[BigBrain2$ID != 1,]
nrow(BigBrain3) #33 observations

### The outlier identifier will sometimes find additional outliers if one runs
### it repeatedly.  Below we run it again on the data set BigBrain3

BigBrain3 <- escalc(measure="ZCOR", ri=r, ni=Sample_Size, data = BigBrain3)

ma.result3 <- rma(yi, vi, data=BigBrain3, level=95, method="DL")
ma.result3

inf <- influence(ma.result3)
inf 
plot(inf) #no more outliers

### The predict function converts the Fisher z results into correlation coefficients
ma.result.in.r.metric3 <- predict(ma.result3, level=95, transf=transf.ztor)
ma.result.in.r.metric3


### There are no more outliers. It is recommended that one reports the meta-analytic 
### results with and without outliers (Kepes, & McDaniel, 2015). 
### Reference: Kepes, S. & McDaniel, M.A. (2015). The validity of conscientiousness is 
### overestimated in the prediction of job performance. PLoS ONE 10(10): e0141468. doi:10.1371/journal.pone.0141468 


### The meta-analysis results ma.result is for the data set that contains outliers (BigBrain).
### ma.result3 has the results for the data set without outliers (BigBrain3)

### show ma.result for the data with outliers
ma.result
ma.result.in.r.metric

### show ma.result for the data without outliers
ma.result3
ma.result.in.r.metric3

### note that the mean changed from .3105 to .3291


######################################################################################
## The term "sensitivity analysis" can be used to refer to both outlier and 
## publication bias analyses. At this point in the analysis, we have a data 
## frame (BigBrain) that contains outliers and another data frame (BigBrain3)
## which does not include outliers.
##
## The file V1_6_sensitivity analysis for correlations.R contains a set of 
## publication bias code that can be brought into the program with a source command.
## This code is not yet a function and one needs to do a few things before 
## calling the code with a source command:
##   1) The data file needs to be called "dat"
##   2) The correlation needs to be called "r"
##   3) The sample size needs to be called "N"
##   4) Assign a character string to the variable Output.Title. That string will
##      appear on top of the output.
######################################################################################

### First, we run publication bias analyses on the data set with outliers

Output.Title1 <- "Sensitivity analysis results"
Output.Title2 <- "BigBrain data"
Output.Title3 <- "Original data; incl. outliers"
dat <- BigBrain

names(dat)
dat$N <- dat$Sample_Size

source("V1_6_sensitivity analysis for correlations.R")
beep(2)

### Source code commentary
### Note that there are two plots. The first plot is a countour-enhanced 
### funnel plot. The second plot is a cumulative meta-analysis by precision.
### One may wish to export them.
###
### Results have been written to a text file; open Excel and paste.


### Second, we run publication bias analyses on the data set with outliers removed

Output.Title1 <- "Sensitivity analysis results"
Output.Title2 <- "BigBrain data"
Output.Title3 <- "Data excl. outliers"
dat <- BigBrain3

names(dat)
dat$N <- dat$Sample_Size

source("V1_6_sensitivity analysis for correlations.R")
beep(2)

### Source code commentary
### Note that there are two plots. The first plot is a countour-enhanced 
### funnel plot. The second plot is a cumulative meta-analysis by sample size.
### One may wish to export them.
###
#### Results have been written to a text file; open Excel and paste.
