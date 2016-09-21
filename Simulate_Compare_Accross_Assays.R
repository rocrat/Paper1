#simulate effect of comparing accross different assays
library(compositions)
library(ggplot2)
MFtrans <- function(x){
  delta <- 0.55/sum(x)# .55 times the smallest detectable value (1 read) after closure
  tdelta <- sum(x == 0) * delta
  cx <- clo(x)
  cxt <- ifelse(cx == 0, delta, cx * (1-tdelta))
  ccxt <- clr(cxt)
  return(ccxt)
}

cpm.stand <- function(x){
  # browser()
  if(any(apply(x, 2, function(y) any(is.na(y))))) stop("No NA's are allowed in the data matrix")
  totals <- apply(x, 2, sum)+1
  x2 <- x+.5
  CPM <- log2(t(t(x2)/totals) * 10^6)
  return(CPM)
}

#Start with 100 samples of the same genes with the exact same 
#level of expression
#create 100 samples of the "same" genes
set.seed(1014)
samps <- list()
for(i in 1:100){
  samps[[i]] <- rbeta(400, shape1 = 1.5, shape2 = 1.5)
}

simDat <- function(minprop, maxprop){
  #create two datasets from two different length assays
  a1 <- list()
  #simulate random allocation of reads to each sample (different number of total reads for each sample)
  a1.totals <- rmultinom(1, size = 25e6, prob = clo(runif(100, min = minprop, max = maxprop)))
  for(i in 1:length(samps)){
    #add additional values for genes which do not overlap between assays and close composition
    s <- clo(c(samps[[i]], runif(220)))
    #generate reads from cell probabilities using the multinomial
    a1[[i]] <- rmultinom(n = 1, size = a1.totals[i], prob = s)
  }
  
  #repeat for the second assay which has a ton more probes (similar to IO vs OBP)
  a2 <- list()
  a2.totals <- rmultinom(1, size = 25e6, prob = clo(runif(100, min = minprop, max = maxprop)))
  for(i in 1:length(samps)){
    #add additional values and close composition
    s <- clo(c(samps[[i]], runif(2000)))
    #generate reads from cell probabilities
    a2[[i]] <- rmultinom(n = 1, size = a2.totals[i], prob = s)
  }
  
  
  a1 <- as.data.frame(a1)
  a2 <- as.data.frame(a2)
  names(a1) <- names(a2) <- paste0("sample", 1:ncol(a1))
  rownames(a1) <- paste0("probe", 1:nrow(a1))
  rownames(a2) <- paste0("probe", 1:nrow(a2))
  
  #select only the "overlapping" probes from each assay
  a1s <- a1[1:400, ]
  a2s <- a2[1:400, ]
  
  #Calculate correlation between paired samples in the two assays
  sample.cors <- vector()
  for(i in names(a1s)){
    sample.cors <- c(sample.cors, stats::cor(a1s[, i], a2s[, i]))
  }
  
  #Calculate correlation between genes in the two assays
  probe.cors <- vector()
  for(i in rownames(a1s)){
    probe.cors <- c(probe.cors, stats::cor(as.numeric(a1s[i, ]), as.numeric(a2s[i, ])))
  }
  
  cordf <- data.frame(rho = c(sample.cors, probe.cors), cor.type = c(rep("sample", length(sample.cors)), rep("probe", length(probe.cors))))
  return(cordf)
}

cordf1 <- simDat(minprop = .01, maxprop = .01)
cordf2 <- simDat(minprop = .009, maxprop = .011) 
cordf3 <- simDat(minprop = .005, maxprop = .015) 
cordf4 <- simDat(minprop = .003, maxprop = .017) 

cordf1$variation <- "None"
cordf2$variation <- "Small"
cordf3$variation <- "Medium"
cordf4$variation <- "Large"

cordf <- rbind(cordf1, cordf2, cordf3, cordf4)
cordf$variation <- factor(cordf$variation, levels = c("None", "Small", "Medium", "Large"))

ggplot(cordf, aes(x = rho)) + geom_density(aes(fill = cor.type), alpha = .5) + facet_wrap(~variation, ncol = 2, scales = "free_y") + ggtitle("Correlation between probes and samples\nwith differing inter-sample variation in total reads") + theme_bw()

simDatCLR <- function(minprop, maxprop){
  #create two datasets from two different length assays
  a1 <- list()
  #simulate random allocation of reads to each sample
  a1.totals <- rmultinom(1, size = 25e6, prob = clo(runif(100, min = minprop, max = maxprop)))
  for(i in 1:length(samps)){
    #add additional values for genes which do not overlap between assays and close composition
    s <- clo(c(samps[[i]], runif(220)))
    #generate reads from cell probabilities using the multinomial
    a1[[i]] <- rmultinom(n = 1, size = a1.totals[i], prob = s)
  }
  
  #repeat for the second assay which has a ton more probes (similar to IO vs OBP)
  a2 <- list()
  a2.totals <- rmultinom(1, size = 25e6, prob = clo(runif(100, min = minprop, max = maxprop)))
  for(i in 1:length(samps)){
    #add additional values and close composition
    s <- clo(c(samps[[i]], runif(2000)))
    #generate reads from cell probabilities
    a2[[i]] <- rmultinom(n = 1, size = a2.totals[i], prob = s)
  }
  
  
  a1 <- as.data.frame(a1)
  a2 <- as.data.frame(a2)
  names(a1) <- names(a2) <- paste0("sample", 1:ncol(a1))
  rownames(a1) <- paste0("probe", 1:nrow(a1))
  rownames(a2) <- paste0("probe", 1:nrow(a2))
  
  #select only the "overlapping" probes from each assay
  a1s <- as.data.frame(apply(a1[1:400, ], 2, MFtrans))
  a2s <- as.data.frame(apply(a2[1:400, ], 2, MFtrans))
  
  #Calculate correlation between paired samples in the two assays
  sample.cors <- vector()
  for(i in names(a1s)){
    sample.cors <- c(sample.cors, stats::cor(a1s[, i], a2s[, i]))
  }
  
  #Calculate correlation between genes in the two assays
  probe.cors <- vector()
  for(i in rownames(a1s)){
    probe.cors <- c(probe.cors, stats::cor(as.numeric(a1s[i, ]), as.numeric(a2s[i, ])))
  }
  
  cordf <- data.frame(rho = c(sample.cors, probe.cors), cor.type = c(rep("sample", length(sample.cors)), rep("probe", length(probe.cors))))
  return(cordf)
}

cordf1.clr <- simDatCLR(minprop = .01, maxprop = .01)
cordf2.clr <- simDatCLR(minprop = .009, maxprop = .011) 
cordf3.clr <- simDatCLR(minprop = .005, maxprop = .015) 
cordf4.clr <- simDatCLR(minprop = .003, maxprop = .017) 

cordf1.clr$variation <- "None"
cordf2.clr$variation <- "Small"
cordf3.clr$variation <- "Medium"
cordf4.clr$variation <- "Large"

cordf.clr <- rbind(cordf1.clr, cordf2.clr, cordf3.clr, cordf4.clr)
cordf.clr$variation <- factor(cordf.clr$variation, levels = c("None", "Small", "Medium", "Large"))

ggplot(cordf.clr, aes(x = rho)) + geom_density(aes(fill = cor.type), alpha = .5) + facet_wrap(~variation, ncol = 2, scales = "free_y") + ggtitle("Correlation between probes and samples\nwith differing inter-sample variation in total reads") + theme_bw() + xlim(c(0.4,1))

#impact of cpm on whole sample
simDatCPM <- function(minprop, maxprop, otherGeneMax = 1){
  #create two datasets from two different length assays
  a1 <- list()
  #simulate random allocation of reads to each sample
  a1.totals <- rmultinom(1, size = 25e6, prob = clo(runif(100, min = minprop, max = maxprop)))
  for(i in 1:length(samps)){
    #add additional values for genes which do not overlap between assays and close composition
    s <- clo(c(samps[[i]], runif(220, min = 0, max = otherGeneMax)))
    #generate reads from cell probabilities using the multinomial
    a1[[i]] <- rmultinom(n = 1, size = a1.totals[i], prob = s)
  }
  
  #repeat for the second assay which has a ton more probes (similar to IO vs OBP)
  a2 <- list()
  a2.totals <- rmultinom(1, size = 25e6, prob = clo(runif(100, min = minprop, max = maxprop)))
  for(i in 1:length(samps)){
    #add additional values and close composition
    s <- clo(c(samps[[i]], runif(2000, min = 0, max = otherGeneMax)))
    #generate reads from cell probabilities
    a2[[i]] <- rmultinom(n = 1, size = a2.totals[i], prob = s)
  }
  
  
  a1 <- as.data.frame(cpm.stand(as.data.frame(a1)))
  a2 <- as.data.frame(cpm.stand(as.data.frame(a2)))
  names(a1) <- names(a2) <- paste0("sample", 1:ncol(a1))
  rownames(a1) <- paste0("probe", 1:nrow(a1))
  rownames(a2) <- paste0("probe", 1:nrow(a2))
  
  #select only the "overlapping" probes from each assay
  a1s <- as.data.frame(a1[1:400, ])
  a2s <- as.data.frame(a2[1:400, ])
  
  #Calculate correlation between paired samples in the two assays
  sample.cors <- vector()
  for(i in names(a1s)){
    sample.cors <- c(sample.cors, stats::cor(a1s[, i], a2s[, i]))
  }
  
  #Calculate correlation between genes in the two assays
  probe.cors <- vector()
  for(i in rownames(a1s)){
    probe.cors <- c(probe.cors, stats::cor(as.numeric(a1s[i, ]), as.numeric(a2s[i, ])))
  }
  
  cordf <- data.frame(rho = c(sample.cors, probe.cors), cor.type = c(rep("sample", length(sample.cors)), rep("probe", length(probe.cors))))
  return(cordf)
}

cordf1.cpm <- simDatCPM(minprop = .01, maxprop = .01)
cordf2.cpm <- simDatCPM(minprop = .009, maxprop = .011) 
cordf3.cpm <- simDatCPM(minprop = .005, maxprop = .015) 
cordf4.cpm <- simDatCPM(minprop = .003, maxprop = .017) 

cordf1.cpm$variation <- "None"
cordf2.cpm$variation <- "Small"
cordf3.cpm$variation <- "Medium"
cordf4.cpm$variation <- "Large"

cordf.cpm <- rbind(cordf1.cpm, cordf2.cpm, cordf3.cpm, cordf4.cpm)
cordf.cpm$variation <- factor(cordf.cpm$variation, levels = c("None", "Small", "Medium", "Large"))

ggplot(cordf.cpm, aes(x = rho)) + geom_density(aes(fill = cor.type), alpha = .5) + facet_wrap(~variation, ncol = 2, scales = "free_y") + ggtitle("Correlation between probes and samples\nwith differing inter-sample variation in total reads") + theme_bw() + xlim(c(0.4,1))

ggplot(cordf4.cpm, aes(x = rho)) + geom_density(aes(fill = cor.type), alpha = .5)

cordf4.cpm <- simDatCPM(minprop = .003, maxprop = .017, otherGeneMax = 5) 
