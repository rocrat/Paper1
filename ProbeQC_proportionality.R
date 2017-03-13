## Proportionality as Probe QC metric


###############################################################################
###############################################################################
## This is a trial analysis to see if we can use proportionality of POS 
## controls for evaluating no sample controls
###############################################################################
###############################################################################



###############################################################################
###############################################################################
## Read in data and load packages
###############################################################################
###############################################################################

library(tidyverse)
library(HTGPackage)

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


df <- readHTG("./Data/ProbeQC/ParserResults_vF_MercuryProbeFeasibility_06Feb2017-rev.xls",
              include.probe = FALSE)

###############################################################################
###############################################################################
## Try using POS as denominator of ALR transformation first
###############################################################################
###############################################################################


df.alr <- apply(df, 2, MFtrans.alr, ivar = which(rownames(df) == "POS1"))

df.ns <- df.alr %>%
  as.data.frame() %>%
  mutate(probe = rownames(df.alr)) %>%
  gather(sample, alr.count, 1:ncol(df)) 

df.ns$type <- gsub("166\\-(\\w+)\\s.*|\\d+.*", "\\1", df.ns$sample)
df.ns$type <- ifelse(grepl("HTG", df.ns$type), "HTG", df.ns$type)

ggplot(df.ns, aes(x = alr.count)) + 
  geom_density( aes(color = sample)) + 
  facet_wrap(~type, ncol = 2) + 
  plotTheme()
