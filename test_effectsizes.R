#GOAL:
#For a given set of random population differences, what % of the time will the mean capture the true diff vs. NPS?

#PROCEDURE:
#step 1. Generate 2 populations (populatin pairs)
#step 2. Calc their diff, save it
#step 3. resample X times. calc stat (mean or NPS). Create 95% CI around this stat. save this CI range.
#step 4. go to step 1 Y times
#step 5. count total misses per population-pair diff

#Setup
source("helpers.R")
options(scipen=999)

#Libraries
library(ggplot2)

#Define simulation function
diffTest <- function(totalPopulations, numOfSamples) {
  popDiff <- c()
  a_mean <- c()
  b_mean <- c()
  simulate_mean_num_missed <- c()
  simulate_nps_num_missed <- c()
  simulate_mean_ci_range <- c()
  simulate_nps_ci_range <- c()
  mean_misses <- c()
  nps_misses <- c()
  for (eachpop in 1:totalPopulations) {
    #Generate populations
    a <- generateNormalPopulation(sizeofPopulation = 10000, scaleMin = 1, scaleMax = 10)
    b <- generateNormalPopulation(sizeofPopulation = 10000, scaleMin = 1, scaleMax = 10)
    
    #Record population parameters
    popDiff[eachpop] <- pctDiff(mean(a), mean(b))
    a_mean[eachpop] <- mean(a)
    b_mean[eachpop] <- mean(b)
    
    #Simulate bootstrapping and count number of times CI does not contain true diff (the CI "missed")
    simdata <- compare2SampleSets(a, b, numOfGeneratedSamples = numOfSamples, numOfResamples = 200)
    simulate_mean_num_missed[eachpop] <- nrow(simdata[simdata$source == "mean_CI" & simdata$containsTrueDiff == FALSE,]) #/ nrow(simulate[simulate$source == "meanCI",])
    simulate_nps_num_missed[eachpop] <- nrow(simdata[simdata$source == "convertNPS10_CI" & simdata$containsTrueDiff == FALSE,]) #/ nrow(simulate[simulate$source == "npsCI",])
    
    #Calculate CI range and compute average range pre source(NPS vs. Mean)
    #simdata$ci_range <- simdata$`97.5%` - simdata$`2.5%`
    #simulate_mean_ci_range[eachpop] <- mean(simdata$ci_range[simdata$source == "meanCI"])
    #simulate_nps_ci_range[eachpop] <- mean(simdata$ci_range[simdata$source == "npsCI"])
  }
  
  results <- data.frame(samplenum = 1:totalPopulations,
                        aMean = a_mean,
                        bMean = b_mean,
                        pdiffs = popDiff,
                        #mean_ci_range = simulate_mean_ci_range,
                        #nps_ci_range = simulate_nps_ci_range,
                        mean_misses = simulate_mean_num_missed, 
                        nps_misses = simulate_nps_num_missed)
  return(results)
}

#Simulate data and sort by population differences
output <- diffTest(1000, 200)
sout <- output[order(output$pdiffs),]
#Calculate absolute value of population differences and sort them
sout$abs_pdiff <- abs(sout$pdiffs)
ab <- sout[order(sout$abs_pdiff),]

#ANALYSIS
#NPS misses 9.84 times more often compared to mean (Mean = 11k, NPS = 111k)
time_nps_misses <- sum(sout$nps_misses) / sum(sout$mean_misses)

#Divide pdiffs into three equal parts. Count number of misses for NPS vs. Mean for each section
#TODO

#Measure the width of the CI for both
#TODO

#Calculcate the correlation of the original dataset with the average correlation of each sample of Mean vs. NPS
#TODO

#What is the reliability of each sample? 
#TODO

#As sample size increases, how large to CIs become?
#TODO

#As SD decreases, how often does each miss?

#GRAPHS
#Plot the size of the population differences by the number of misses for Mean and NPS
ggplot(data = sout, aes(x = pdiffs)) + 
  geom_point(aes(y = mean_misses, color = "mean", alpha = 1/10)) +
  geom_point(aes(y = nps_misses, color = "nps", alpha = 1/10))

#Same as above but pop differences are absolute values
ggplot() +
  geom_point(data = ab, aes(x = abs_pdiff, y = nps_misses, color = "nps", alpha = 1/10)) +
  geom_point(data = ab, aes(x = abs_pdiff, y = mean_misses, color = "mean", alpha = 1/10))

#Scatterplot of population means, color coded by number of NPS misses
ggplot(sout, aes(x = aMean, y = bMean, color = nps_misses, size = pdiffs, alpha = 1/100)) + 
  geom_point() + 
  scale_color_gradientn(colors = rainbow(10))

#Plot cumulative distribution of misses by size of population difference
ggplot(data = ab, aes(x = ab$abs_pdiff)) + 
  geom_point(aes(y = cumsum(nps_misses), color = "nps")) +
  geom_point(aes(y = cumsum(mean_misses), color = "mean")) + 
  ylab("Cumulative # of Misses") + xlab("Abs Population Difference. Mean(a) - Mean(b)")





