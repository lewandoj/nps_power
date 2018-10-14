source("helpers.R")
options(scipen=999)

#Libraries
library(ggplot2)

#Define simulation function
sampleSizeTest <- function(maxPopSize, numOfSamples) {
  popDiff <- c()
  a_mean <- c()
  b_mean <- c()
  misses <- c()
  pos <- c()
  popLength <- c()
  missesList <- list(c(), c(), c(), c())
  iter <- seq(from = .001, to = maxPopSize, by = .001)
  
  #Generate populations
  a <- generateSimplePopulation(1000, 1,2,3,4,5,6,7,8,9,10)
  b <- generateSimplePopulation(1000, 1,1,1,1,1,1,7,8,9,10) #more 1s, No 2-6s
  
  for (eachpop in iter) {
    pos <- eachpop*1000
    #Record population parameters
    popDiff[pos] <- pctDiff(mean(a), mean(b))
    a_mean[pos] <- mean(a)
    b_mean[pos] <- mean(b)
    popLength[pos] <- length(a)*eachpop
    
    
    #Simulate bootstrapping and count number of times CI does not contain true diff (the CI "missed")
    simdata <- compare2SampleSetskWay(a, b, 
                                      summaryStats = list(mean = mean, median = median, NPS10 = convertNPS10, top2 = convertTopTwo10),
                                      numOfGeneratedSamples = numOfSamples, 
                                      numOfResamples = 100, 
                                      pctSizeOfGeneratedSample = eachpop)
    #From each item in simdata, select col 5, then count num of FALSE values. FALSE means number of times CI didn't capture true diff
    misses <- lapply(lapply(simdata, function(x) x[,5]), function(x) length(x[x == FALSE]))
    missesList <- mapply(append, missesList, misses, SIMPLIFY = F)
  }
  
  results <- data.frame(sampleSize = popLength,
                        aMean = a_mean,
                        bMean = b_mean,
                        pdiffs = popDiff,
                        meanMisses = missesList[[1]],
                        medianMisses = missesList[[2]],
                        npsMisses = missesList[[3]],
                        top2Misses = missesList[[4]])
  return(results)
}

output <- sampleSizeTest(.10, 100)


ggplot(data = output, aes(x = sampleSize)) + 
  geom_point(aes(y = meanMisses, color = "mean", alpha = 1/100)) +
  geom_point(aes(y = npsMisses, color = "nps", alpha = 1/100)) +
  geom_point(aes(y = medianMisses, color = "median", alpha = 1/100)) +
  geom_point(aes(y = top2Misses, color = "top2", alpha = 1/100))


#Why are the # misses increasing with N?
#Each calculation is converging on its own measure of central tendency, away from the mean. 
#Note that NPS and Top2 are converging less slowly (their misses reach max at larger Ns)




