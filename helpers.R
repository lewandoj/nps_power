setwd("C:/Users/Palmer/Google Drive/R projects/nps_power")
options(scipen=999)

#Libraries
library(truncnorm)


#bootstrapping functions
#Takes a vector and resample it with replacement X times
#Args:
# data = vector of numerics 
# numOfSamples = integer indicating # of resamples from the data, with replacement
# pctOfCasesFromPopulation = numeric from 0 - 1
# statistic = an aggregation function that will summarize the data into a single value (such as mean, median)
# alpha = the percentile of cases that will be included within the sample (default = .05)
# Returns the 2.5 and 97.5 percentiles of values that represent the lower and upper CI values  
#Example: bootstrap(x, 100, .5, mean) eats a vector of numbers, resamples 100 times, returns top/bottom 2.5% values
bootstrap <- function(data, numOfSamples, pctOfCasesFromPopulation, statistic, alpha = .05) {
  
  savedStat <- c()
  #Pick statistical aggregation function like mean, median, or NPS. Should only have 1 arg (x, ...)
  applyStatistic <- match.fun(statistic)
  
  for (i in 1:numOfSamples) {
    #Generate a sample from data
    bs_sample <- sample(data, (length(data) * pctOfCasesFromPopulation), replace = TRUE)
    #Apply that statistic to the sample
    calculate_statistic <- applyStatistic(bs_sample)
    #Create object savedStat to hold result from each resample
    savedStat[i] <- calculate_statistic
  }
  #Order all results
  ordered <- sort(savedStat)
  #Caculuate upper and lower alpha values (divided by 2 bc it's a 2-tail test), save results as a list
  CI <- quantile(ordered, c(1-(1-(alpha/2)), (1-(alpha/2))), type = 4)
  #lowerCI <- quantile(ordered, 1-(1-(alpha/2)))
  result <- list(savedStat, CI)
  names(result) <- c("savedStat", "CIs")
  
  return(result)
}


bootstrapDiffs <- function(sample1, sample2, maxTimesResampled, pctOfCasesFromSamples, statistic, alpha = .05) {
  
  savedStat <- c()
  #Pick statistical aggregation function like mean, median, or NPS. Should only have 1 arg (x, ...)
  applyStatistic <- match.fun(statistic)
  
  for (i in 1:maxTimesResampled) {
    #Generate a sample from data
    sample1_resample <- sample(sample1, (length(sample1) * pctOfCasesFromSamples), replace = TRUE)
    sample2_resample <- sample(sample2, (length(sample2) * pctOfCasesFromSamples), replace = TRUE)
    #Apply that statistic to the sample
    calculate_statistic <- pctDiff(applyStatistic(sample1_resample), applyStatistic(sample2_resample))
    #calculate_statistic <- applyStatistic(sample1_resample) - applyStatistic(sample2_resample) #OLD
    
    #Create object savedStat to hold result from each resample
    savedStat[i] <- calculate_statistic
  }
  #Order all results
  sortedStat <- sort(savedStat)
  #Caculuate upper and lower alpha values (divided by 2 bc it's a 2-tail test), save results as a list
  CI <- quantile(sortedStat, c((1-(1-(alpha/2))), (1-(alpha/2))), type = 5, na.rm = TRUE)
  #lowerCI <- quantile(ordered, 1-(1-(alpha/2)))
  result <- list(savedStat, CI)
  names(result) <- c("savedStat", "CIs")
  
  return(result)
}




#NPS conversion functions
#Args: data = a vector of integers with a max of 10, 7, or 5
#Returns: NPS score of data ranging from -1 to 1
convertNPS10 <- function(data) {
  #Recode: {9,10} are promoters; {7,8} are neutrals; {1,2,3,4,5,6} are detractors    
  total <- length(data)
  pctPromoters <- length(data[data > 8])/total
  pctNeutrals <- length(data[data == 7 | data == 8])/total
  pctDetractors <- length(data[data < 7])/total
    
  nps <- ((1 * pctPromoters) + (0 * pctNeutrals) + (-1 * pctDetractors)) * 100
  return(nps)
} 

convertNPS7 <- function(data){
  #Recode: {6,7} are promoters; {5,6} are neutrals; {1,2,3,4} are detractors
  total <- length(data)
  pctPromoters <- length(data[data > 5])/total
  pctNeutrals <- length(data[data == 5 | data == 6])/total
  pctDetractors <- length(data[data < 5])/total
  
  nps <- ((1 * pctPromoters) + (0 * pctNeutrals) + (-1 * pctDetractors)) * 100
  return (nps)
}

convertNPS5 <- function(data){
    #Recode: {4,5} are promoters; {3} is neutral; {1,2} are detractors
    total <- length(data)
    pctPromoters <- length(data[data > 3])/total
    pctNeutrals <- length(data[data == 3])/total
    pctDetractors <- length(data[data < 3])/total
    
    nps <- ((1 * pctPromoters) + (0 * pctNeutrals) + (-1 * pctDetractors))
    return(nps)
    
  }

# Compare two sets of samples to population. 
# Each row is the difference between two bootstrapped samples and whether that boostrap CI contains the pop diff
#Args: 
# p1/p2 = vectors of population raw values to be compared
# numOfGeneratedSamples = integer indicating number of samples to generate from each pop
# numOfResamples = integer indicating number of times to boostrap each sample to generate the CI
#Returns: data frame with 5 rows: 
# Lower and Upper CIs; 
# statistic that was calculated; 
# true population diff value
# bool of whether the CI contains the pop diff
compare2SampleSets <- function(p1, p2, 
                               numOfGeneratedSamples = 50, pctSizeOfGeneratedSample = .01, 
                               numOfResamples = 100, pctSizeOfResample = 1.0) {

  #Sample each population 'numOfGeneratedSamples' times using .1% of population
  first <- resampler(p1, numOfGeneratedSamples, pctSizeOfGeneratedSample)
  second <- resampler(p2, numOfGeneratedSamples, pctSizeOfGeneratedSample)
  
  #Calculate difference in mean/NPS for each sample and boostrap those differences
  bsMeanDiffs <- c()
  bsnpsDiffs <- c()
  for (i in 1:numOfGeneratedSamples) {
    bsMeanDiffs[[i]] <- bootstrapDiffs(first[[i]], second[[i]], numOfResamples, pctSizeOfResample, mean)$CIs 
    bsnpsDiffs[[i]] <- bootstrapDiffs(first[[i]], second[[i]], numOfResamples, pctSizeOfResample, convertNPS10)$CIs
  }
  result <- as.data.frame(do.call(rbind, c(bsMeanDiffs, bsnpsDiffs)))
  
  #rename source column
  result$source[1:length(bsMeanDiffs)] <- "meanCI"
  result$source[(length(bsnpsDiffs)+1):nrow(result)] <- "npsCI"
  
  result$truediff <- pctDiff(mean(p1), mean(p2))
    #result$truediff[result$source == "meanCI"] <- (mean(p1) - mean(p2))
    #result$truediff[result$source == "npsCI"] <- (convertNPS10(p1) - convertNPS10(p2))
  #test if CIs contain truediff
  result$containsTrueDiff <- (result$truediff > result$`2.5%`) & (result$truediff < result$`97.5%`)
  return(result)
}

#Generate a population of user responses on NPS scale
generateNormalPopulation <- function(sizeofPopulation, scaleMin, scaleMax, scaleMean = NULL, scaleSD = NULL) {
  round(rtruncnorm(n = sizeofPopulation,
                   a = scaleMin, 
                   b = scaleMax,
                   mean = 
                     if(is.null(scaleMean)) {
                       runif(1, scaleMin, scaleMax) #pick random number between scaleMin and scaleMax
                       } else {
                         scaleMean
                         }, 
                   sd = 
                     if(is.null(scaleSD)) {
                       runif(1, 1.0, 10.0) #pick a random number between 1 and 10
                     } else {
                       scaleSD
                     }), 0)
}

#Generate x numbers of each [a, j]. x times 10 numbers total
generateSimplePopulation <- function(x,a,b,c,d,e,f,g,h,i,j) {
  c(rep.int(a, x),
    rep.int(b, x),
    rep.int(c, x),
    rep.int(d, x),
    rep.int(e, x),
    rep.int(f, x),
    rep.int(g, x),
    rep.int(h, x),
    rep.int(i, x),
    rep.int(j, x))
}

#Calculate the percent difference between two numeric values
#Return the percent difference
pctDiff <- function(first, second) {
  if(first == 0) {
    return(NA)
    } else { 
      return(((second - first) / first) * 100)
    }
  }


#Generate a set of random samples from a given population 
#Args:
# population = vector from which you'll resample 
# numOfResamples = numeric from 1 - infinity, number of times to resample the populatin
# pctSampledFromPop = numeric from 0 - 1. What % of the population do you want to sample?
# Returns: List of vectors from each sample
resampler <- function(population, numOfResamples, pctSampledFromPop) {
  data <- c()
  for (i in 1:numOfResamples) {
    data[[i]] <- sample(population, (length(population) * pctSampledFromPop), replace = TRUE)
  }
  return(data)
}


  
  

