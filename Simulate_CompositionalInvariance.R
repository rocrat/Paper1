##Simulate changes in compostion induced by changes in the total.

###############################################################################
###############################################################################
## Load libraries and set up workspace
###############################################################################
###############################################################################

library(HTGPackage)
library(tidyverse)

MFtrans <- function(x){
  delta <- 0.55/sum(x)
  tdelta <- sum(x == 0) * delta
  cx <-  x/sum(x)
  cxt <- ifelse(cx == 0, delta, cx * (1-tdelta))
  return(cxt)
}

###############################################################################
###############################################################################
## Read in real data to guide simulations
###############################################################################
###############################################################################

dfraw <- readHTG("C:/Projects/CompositionalInvariance/RawData/Pfizer Sample Input OBP parsed 18-Nov-2016.xls",
                 include.probe = FALSE)
totals <- colSums(dfraw)
totals <- totals[order(totals)]

#close data to get cell probabilities using MF transformation
dfclo <- apply(dfraw, 2, function(x) x/sum(x))
#Find component probabilities by taking the average over all samples. Since 
#samples are not techinical reps we separate by type
tnbcMeans <- rowMeans(dfclo[, which(grepl("Trip", colnames(dfclo)))])
lumAMeans <- rowMeans(dfclo[, which(grepl("Breast", colnames(dfclo)))])

# tnbcMeans <- tnbcMeans[order(tnbcMeans)]
# lumAMeans <- lumAMeans[order(lumAMeans)]

plot(log(tnbcMeans), log(lumAMeans))
abline(a=0, b = 1)

###############################################################################
###############################################################################
## Simulate Data with compositional invariance
###############################################################################
###############################################################################

# The total number of reads assigned to a sample is a random variable.  Under
# normal sequencing we assume these follow a multinomial distribution with
# each sample having approximately the same probability of receiving reads.
# In practice, the sample allocation is much less even than a multinomial with 
# equal probability, however.  Therefore we generate sample cell probabilities
# using the normal distribution with mean and standard deviation estimated 
# from real data (the Pfizer sample input study data)

simtotals <- rmultinom(1, size = 500000000, prob = rep(1/69, 69))#not enough variation

set.seed(1)
simtotals <- rnorm(n = 69, mean = 7222000, sd = 1000000)
simtotals <- simtotals[order(simtotals)]

simdf <- data.frame(tots = c(totals, simtotals), type = rep(c("real", "simulated"), each = 69), rank = c(1:69, 1:69))
ggplot(simdf, aes(x = rank, y = tots)) + 
  geom_point(aes(color = type)) + 
  ylab("Totals") +
  xlab("Rank") +
  plotTheme()


#Use tnbc means as probabilities for data generation. This method under-
#represents the realized number of 0's. Therefore, I set the limit of
#detection at 15 so that all values <= 15 -> 0
simdat <- matrix( nrow = length(tnbcMeans), ncol = length(simtotals))

for (i in 1:length(simtotals)) {
  simdat[, i] <- rmultinom(1, size = simtotals[i], prob = tnbcMeans)
  simdat[, i] <- ifelse(simdat[, i] <= 15, 0, simdat[, i])
}
#this is a data set which should exhibit compositional invariance
simdat <- as.data.frame(simdat) 

###############################################################################
###############################################################################
## Simulate data WITHOUT compositional invariance
###############################################################################
###############################################################################


#Create data set without compositional invariance by multiplying the probability 
#vector of each sample by a vector of scaling factors determined by 
