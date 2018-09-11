# GOAL
# The goal of this simulation is to illustrate how the NPS rebinning computation
# into promoters and detractors reduces the likelihood of detecting an effect when one exists.
# In this example, two datasets of N=10000 length are created on a 10-point scale: 
# (1) a dataset with a uniform distribution of 1000 cases per point and (2) a dataset where all [2:6] 
# values are rebinned to a value of 1. According to NPS computations, the second dataset is identical 
# to the first since all values ranging from 1:6 are considered 'detractors'. 

# PROCEDURE
# This simulation creates 500 samples from each dataset ('a' and 'b'). 
# Each sample consists of 10% of the values of 'a' and 'b'.
# 100% of the cases from each sample are resampled 200 times with replacement (i.e., boostrapped). 
# A statistic is applied to summarize each resample: both an NPS computation and an arthmetic mean.
# Summarized resamples are compare by calculating the percent difference. 
# The middle 95% of values from the resamples are computed to construct a confidence interval. 
# All 500 confidence intervals (per statistic) are evaluated as to whether they contain 
# the true mean percent difference between 'a' and 'b'.

# HYPOTHESIS
# If NPS performs as well as the mean, it should contain the true difference as often as the mean. 

source("helpers.R")
options(scipen=999)


#Generate two populations, a uniform population 'a' and a low-value population 'b'
a <- generateSimplePopulation(1000, 1,2,3,4,5,6,7,8,9,10)
b <- generateSimplePopulation(1000, 1,1,1,1,1,1,7,8,9,10) #more 1s, No 2-6s

mean_pctdiff <- pctDiff(mean(a),mean(b)) #-27.27% difference (5.5 --> 4.0)
nps_pdiff <- pctDiff(convertNPS10(a), convertNPS10(b)) #Zero% difference (-40 --> -40)

#Generate 500 samples using 10% of the cases from 'a' and 'b'.
#Resample each sample 500 times using 100% of the cases.
#Calculate CI of percent differences. 
#Evaluate whether the CI contains the true mean percent difference. 
simdata <- compare2SampleSets(a, b, 
                              summaryStat1 = mean,
                              summaryStat2 = convertNPS10,
                              numOfGeneratedSamples = 500, 
                              pctSizeOfGeneratedSample = .1, 
                              numOfResamples = 200, 
                              pctSizeOfResample = 1)

#ANALYSIS
#Count number of misses for NPS vs. Mean
#NPS FAILED to detect a 27% actual difference about 93% of the time. Mean failed about 5% of the time. 
t <- table(simdata$source, simdata$containsTrueDiff) 
prop.table(t, 1) * 100 

#Calculate average CI width for both
simdata$ci_range <- simdata$`97.5%`-simdata$`2.5%`
#CI range for mean = 9.64; NPS = 35.78; NPS is 4x larger
aggregate(simdata$ci_range ~ simdata$source,
          FUN = mean,
          data = simdata)


#GRAPHS
#Graph simulation CIs
plot.new()
plot(x = c(0, 500), 
     y = c(min(simdata$`2.5%`), max(simdata$`97.5%`)), 
     xlab = "population #", 
     ylab = "CI range", type = "n")
abline(h = unique(simdata$truediff), col = "black", lwd = 2)

#GRAPH 1
#Plot Mean CIs and midpoints
segments(x0 = seq(0, 500), 
         y0 = simdata$`2.5%`[simdata$source == "meanCI"], 
         x1 = seq(0, 500), 
         y1 = simdata$`97.5%`[simdata$source == "meanCI"], 
         col = "#ff5252", 
         lwd = 2)
points(x = seq(0, nrow(simdata)/2-1),
       y = (((simdata$`97.5%`+simdata$`2.5%`)/2)[simdata$source == "meanCI"]),
       pch = 19,
       col = "#c50e29",
       bg = "#c50e29")
#Plot NPS CIs and midpoints
segments(x0 = seq(0, 500), 
         y0 = simdata$`2.5%`[simdata$source == "npsCI"], 
         x1 = seq(0, 500), 
         y1 = simdata$`97.5%`[simdata$source == "npsCI"], 
         col = "#536dfe",
         lwd = 2)
points(x = seq(0, nrow(simdata)/2-1),
       y = (((simdata$`97.5%`+simdata$`2.5%`)/2)[simdata$source == "npsCI"]),
       pch = 19,
       col = "#0043ca",
       bg = "#0043ca")

#GRAPH 2
#Graph histograms of CI midpoints
hist(((simdata$`97.5%`+simdata$`2.5%`)/2)[simdata$source == "meanCI"], 
     breaks = 30, 
     main = "",
     xlab = "Midpoint of each samples mean-percent-difference",
     col = "#ff5252",
     lty = "blank",
     xlim = c(-40, 40))
hist(((simdata$`97.5%`+simdata$`2.5%`)/2)[simdata$source == "npsCI"], 
     breaks = 30, 
     col = "#536dfe", 
     lty = "blank",
     add = TRUE)
abline(v = unique(simdata$truediff), col = "black", lwd = 4)







