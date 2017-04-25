###############################################################################
###############################################################################
## Compositional invariance using single model
###############################################################################
###############################################################################


compInvOne <- function(x, #data in wide form with first variable as probe
                       ivar = NULL, #The denominator if using alr trans
                       ctrlProbes = NULL, #character vector of control probes
                       alpha = 0.05){
  df.l <- x %>%
    gather(sample, count, 2:ncol(x)) 
  
  
  #Gather the sample totals
  totals <- data.frame(sample = names(x)[-1],
                       total = log(colSums(x[, -1])),
                       stringsAsFactors = FALSE)
  
  #Add the totals as a variable to the long data
  df.l <- dplyr::full_join(df.l, totals, by = "sample")
  
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
  #combine back into long data
  df.l <- do.call(rbind, alr.s)
  
  mod <- lm(trans.count ~ probe + total, data = df.l)
  return(mod)
}
