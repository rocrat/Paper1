##########################################################################
##########################################################################
### Script to apply Compositional Invariance measures to various data sets
##########################################################################
##########################################################################


##These functions produce two types of compositional invariance plots
heatPlot <- function(df){
  #Create annotation for  samples
  anno <- data.frame(sample = names(df),
                     total = colSums(df),
                     stringsAsFactors = FALSE)
  #Set the order of the samples based on total read count
  anno <- anno[order(anno$total), ]
  
  clrdf <- apply(df, 2, MFtrans.clr) %>%
    as.data.frame()
  
  #find the distance between all pairwise samples using euclidean distance
  d <- stats::dist(t(clrdf), method = "euclidean") %>%
    as.matrix()
  
  #set the order of samples
  dforder <- anno$sample
  
  dl <- d %>%
    as.data.frame() %>%
    mutate(row.sample = rownames(d)) %>%
    gather(col.sample, distance, 1:ncol(d))
  
  dl$row.sample <- factor(dl$row.sample, levels = dforder)
  dl$col.sample <- factor(dl$col.sample, levels = dforder)
  
  
  #Add the total reads along the diagonal
  dltots <- dl[which(dl$distance == 0) ,] 
  dltots$total <- anno[dltots$row.sample, ]$total
  dltots$row.sample <- factor(dltots$row.sample, levels = dforder)
  dltots$col.sample <- factor(dltots$col.sample, levels = dforder)
  
  #Plot pairwise distance between all samples
  heat <- ggplot(dl, aes(x = col.sample, y = row.sample)) + 
    geom_tile(aes(fill = distance)) + 
    geom_tile(data = dltots, 
              aes(color = log(total)), 
              fill = "lightgrey", 
              lwd = 1) +
    scale_fill_gradient("Distance", na.value = "lightgrey", limits = c(0,100)) +
    scale_color_gradientn("log Total Reads", 
                          colors = c("red", "orange", "white")) +
    xlab("Sample Ordered by Total Reads") + 
    ylab("Sample Ordered by Total Reads") +
    plotTheme(axis.text.x = element_blank(),
              axis.text.y = element_blank())
  
  return(list(p = heat, df = dl))
}

distPlot <- function(df, quant){
  ##This function takes in a data frame with samples as columns
  ## and generates a plot showing the distance between each sample
  ## and the mean of the top quant samples
  totals <- apply(df, 2, sum)
  qt <- quantile(totals, prob = quant)
  
  dfclr <- apply(df, 2, MFtrans.clr) %>%
    as.data.frame()
  
  cent <- apply(dfclr[, which(totals > quant)], 1, mean)
  distances <- apply(t(dfclr), 1, function(x) dist(rbind(x, cent)))
  
  distdf <- data.frame(SampleID = names(distances),
                       distance = distances,
                       total = as.numeric(totals),
                       inMean = ifelse(totals > qt, TRUE, FALSE))
  
  
  distdf$SampleID <- factor(distdf$SampleID, 
                            levels = names(totals)[order(totals)])
  
  p <- ggplot(distdf, aes(x = SampleID, y = distance)) +
    geom_point(aes(size = log(total), color = log(total))) +
    plotTheme(axis.text.x = element_blank()) +
    ylab(paste0("Distance between sample and center\n", 
                100*quant, "% of high reads samples" )) +
    xlab("Samples ordered by total aligned reads") +
    scale_color_continuous(guide = FALSE)
  
  return(list(p = p, df = distdf))
}

