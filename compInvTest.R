#Test for compositional Invariance using either ALR or CLR transofrmation
#Function takes a HTG data as produced by read HTG with include.probe = TRUE
#

compInvTest <- function(data, #data in wide form with first variable as probe
                        trans, #one of "alr" or "clr"
                        ivar = NULL, #The denominator if using alr trans
                        ctrlProbes = NULL, #character vector of control probes
                        alpha = 0.05) { #significance level for test (after fdr adj)
  
  #Change to long form and add identifying information
  
  if(!is.character(data[, 1]) & !is.factor(data[, 1])) stop(
    "The first column must contain probe names"
  )
  df.l <- data %>%
    gather(sample, count, 2:ncol(data)) 
  
  
  #Gather the sample totals
  totals <- data.frame(sample = names(data)[-1],
                       total = log(colSums(data[, -1])),
                       stringsAsFactors = FALSE)
  
  #Add the totals as a variable to the long data
  df.l <- dplyr::full_join(df.l, totals, by = "sample")
  
  
  if(trans == "alr"){
    
    if(is.null(ivar)){
      #Find highest expressing probe for use as ALR denominator
      avg <- df.l %>%
        group_by(probe) %>%
        summarize(avg = mean(count)) 
      den <- as.character(avg$probe[which.max(avg$avg)])
      
    }else{
      den <- ivar
    }
    
    
    
    ##Creat alr transformed count values
    df.l.s <- split(df.l, f = df.l$sample)
    alr.s <- list()
    for (i in 1:length(df.l.s)) {
      df <- df.l.s[[i]]
      #remove pos and neg
      if(!is.null(ctrlProbes)){
        df <- df[-which( df$probe %in% ctrlProbes), ]
      }
      
      alrvec <- MFtrans.alr(df$count, ivar = which(df$probe == den))
      df.alr <- data.frame(probe = df$probe[-which(df$probe == den)],
                           trans.count = alrvec,
                           sample = unique(df$sample),
                           total = unique(df$total))
      alr.s[[i]] <- df.alr
    }
    #combine back into log data
    df.l <- do.call(rbind, alr.s)
  }
  
  if(trans == "clr"){
    df.l <- df.l %>%
      group_by(sample) %>%
      mutate(trans.count = MFtrans.clr(count))
  }
  
  if(trans == "clo"){
    df.l <- df.l %>%
      group_by(sample) %>%
      mutate(trans.count = MFtrans(count))
  }
  
  betas <- testInvariance(df.l,
                          countVar = "trans.count",
                          splitVar = "probe",
                          predVar = "total")
  numSig <- sum(betas$adjP < alpha)
  
  return(list(betas = betas, data = df.l, numSig = numSig))
}



testInvariance <- function(x, countVar, splitVar, predVar) {
  #x = data frame 
  #countVar = variable with transformed read counts
  #splitVar = generally the column of probe names
  #predVar = generally the total number of reads
  
  #first split data into a list of data frames for each probe
  xs <- split(x, f = factor(x[[splitVar]])) 
  betas <- matrix(nrow = length(xs), ncol = 4) %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  names(betas) <- c(splitVar, "est", "err", "pval")
  for ( i in 1:length(xs)){
    #For each probe determine whether the sample total affects the count
    mod <- lm(as.formula(paste0(countVar, " ~ ", predVar)), data = xs[[i]])
    est <- summary(mod)$coef[2, 1]
    err <- summary(mod)$coef[2, 2]
    pval <- summary(mod)$coef[2, 4]
    betas[i, 1] <- unique(xs[[i]]$probe)
    betas[i, 2:4] <- c(est, 
                       err,
                       pval)
  }
  #Adjust for multiple testing
  betas$adjP <- p.adjust(betas$pval, method = "BH")
  return(betas)
}
