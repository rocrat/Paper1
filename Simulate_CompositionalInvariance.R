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

MFtrans.clr <- function(x){#apply to a vector
  cxt <- MFtrans(x)
  ccxt <- log(cxt) - sum(log(cxt))/length(cxt)
  return(ccxt)
}

MFtrans.alr <- function(x, ivar){#apply to a vector
  cxt <- MFtrans(x)
  lcxt <- log(cxt)
  ccxt <- lcxt[-ivar] - lcxt[ivar]
  return(ccxt)
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
simdat <- cbind(probe = paste0("p", 1:nrow(simdat)), simdat)

beta1.alr <- compInvTest(data = simdat, trans = "alr")
beta1.clr <- compInvTest(data = simdat, trans = "clr")
beta1.clo <- compInvTest(data = simdat, trans = "clo")

sum(beta1.alr$betas$adjP <= 0.05)
sum(beta1.clr$betas$adjP <= 0.05)
sum(beta1.clo$betas$adjP <= 0.05, na.rm = TRUE)
sum(beta1.clr$betas$est > 0 & beta1.clr$betas$adjP <= 0.05, na.rm = TRUE)
sum(beta1.clr$betas$est < 0 & beta1.clr$betas$adjP <= 0.05, na.rm = TRUE)

clrdf <- beta1.clr$data
clrbeta <- beta1.clr$betas

probProbes <- paste0("p", beta1.clr$betas[which(beta1.clr$betas$adjP <= 0.05), ]$probe)
ranProbes <- sample(probProbes, 100)
ggplot(clrdf[which(clrdf$probe %in% ranProbes), ], aes(x = total, y = trans.count)) +
  geom_point(aes(color = probe)) + 
  geom_smooth(aes(color = probe), method = "lm") +
  plotTheme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            legend.position = "none")

###############################################################################
###############################################################################
## Simulate data WITHOUT compositional invariance
###############################################################################
###############################################################################


#Create data set without compositional invariance by multiplying the probability 
#vector of each sample by a vector of scaling factors determined by the total
# reads for a sample.  The number of multipiers (betas) != 1 should be the same 
# within a data set (same betas) but this number should vary between data sets.
#This model assumes that some probes are sensitive to the total number of reads
# allocated to a sample. I.e. some probes are more likely to be counted if there
# are a large number of total reads. (not sure what the chemical mechanism would
# be) (empirically these effects should be on the order of 1.5x10^-7)

percBeta <- c(.01, .05, .15, .25)


simdatIV <- list()
for (i in 1:length(percBeta)){
  #Get the number of non-zero multipliers
  nB <- ceiling(percBeta[i], length(tnbcMeans))
  for (j in length(simtotals)){
    newProbs <- 
  }
}


