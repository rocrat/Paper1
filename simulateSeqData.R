#simulate data with batch effects
# source("http://bioconductor.org/biocLite.R")
# biocLite("polyester")
# biocLite("ballgown")
library(polyester)
library(Biostrings)
library(ballgown)

plotPCA2 <- function(data, annodata, projname, Max=80, pc1 = "PC1", pc2 = "PC2", scale = FALSE){
  x <- t(data)
  pcares<-prcomp(x, scale = scale)
  print(summary(pcares))
  pvar <- summary(pcares)$importance[3, 2]
  plot(pcares, type="l")
  pcavect<-pcares$x[,1:5]
  indanno<-intersect(row.names(pcavect), row.names(annodata))
  pcavect<-pcavect[indanno,]
  anno.pca<-annodata[indanno,]
  pcadata <- cbind(anno.pca, pcavect)
  return(list(pcadata, pvar))
}

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

create_my_reads <- function(mu, fit, p0, m = NULL, n = NULL, mod = NULL, beta = NULL, seed = NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(mod) | is.null(beta)) {
    cat("Generating data from baseline model.\n")
    if (is.null(m) | is.null(n)) {
      stop(.makepretty("create_read_numbers error: if you don't specify\n            mod and beta, you must specify m and n.\n"))
    }
    index = 1:length(mu)
    mus = mu[index]
    p0s = p0[index]
    mumat = log(mus + 0.001) %*% t(rep(1, n))
  }
  else {
    m = dim(beta)[1]
    n = dim(mod)[1]
    index = 1:length(mu)
    mus = mu[index]
    p0s = p0[index]
    ind = !apply(mod, 2, function(x) {
      all(x == 1)
    })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log(mus + 0.001) + beta %*% t(mod)
  }
  muvec = as.vector(mumat)
  sizevec = predict(fit, muvec)$y
  sizemat = matrix(sizevec, nrow = m)
  counts = sizemat * NA
  for (i in 1:m) {
    counts[i, ] = rbinom(n, prob = (1 - p0s[i]), size = 1) * 
      rnbinom(n, mu = exp(mumat[i, ]), size = exp(sizemat[i, ]))
  }
  return(counts)
}

fasta_file <- "./Paper1/mart_export_io.txt"
fasta <- readDNAStringSet(fasta_file)

small_fasta = fasta[1:20]
writeXStringSet(small_fasta, 'io_small.txt')

smallfc <- matrix(1, 20, 2)
# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx = round(20 * width(small_fasta) / 100)

# simulation call:
simulate_experiment('io_small.txt', reads_per_transcript = readspertx, 
                    num_reps=c(10,10), fold_changes = smallfc, outdir = './Paper1/simulated_reads') 

load("./Paper1/Data/IO_Repro.Rdata")

raw.w <- raw.w[-c(1:4), -1]

#remove questionable samples
suspects <- c("uRNA_25ng_P88_R2", "uRNA_25ng_P88_R3", "uRNA_D2_P02_R15", "uRNA_D2_P02_R16", "uRNA_D2_P02_R17")
raw.clean <- data.matrix(raw.w[, -which(names(raw.w) %in% suspects)])

#get parameters for simulation
params <- get_params(raw.clean)

N <- 30

simdat1 <- create_my_reads(params$mu, params$fit, params$p0, m = nrow(raw.clean), n = N, seed = 1014) 
#shift parameters for second data set
newmu <- round(params$mu - (.2 * params$mu))
newmu <- ifelse(newmu < 0, 0, newmu)
simdat2 <- create_my_reads(newmu, params$fit, params$p0, m = nrow(raw.clean), n = N, seed = 1410) 
colnames(simdat1) <- paste0("samp1_", 1:N)
colnames(simdat2) <- paste0("samp2_", 1:N)

#combine into 1 data set
simdat <- cbind(simdat1, simdat2)
rownames(simdat) <- rownames(raw.clean)

annodf <- data.frame(SampleID = c(colnames(simdat1), colnames(simdat2)), Plate = rep(c(1, 2), each = N))
rownames(annodf) <- annodf$SampleID

#See if there is a random batch effect

pcaDatl2 <- apply(simdat, c(1,2), function(x) log2(x + 1))
pcplotdat <- plotPCA(pcaDatl2, annodf, scale = FALSE)
ggpca <- ggplot(pcplotdat[[1]], aes(x= PC1, y= PC2, color=factor(Plate), label = SampleID)) + 
  geom_point( size = 2)  + 
  # geom_text(size = 1.7) + 
  ggtitle(paste0("Proportion of variance =", round(pcplotdat[[2]], 2)))
# xlim(-50, 50) + ylim(-50, 50) 
ggsave("./Paper1/Figures/polySim_PCA_Raw(log2).png", ggpca, width = 9, height = 8 )

pcaDatn <- quantile_normalisation(simdat)
pcplotdatn <- plotPCA(pcaDatn, annodf, scale = FALSE)
ggpcan <- ggplot(pcplotdatn[[1]], aes(x= PC1, y= PC2, color=factor(Plate), label = SampleID)) + 
  geom_point( size = 2) + 
  ggtitle(paste0("Proportion of variance =", round(pcplotdatn[[2]], 2)))
# geom_text(size = 1.7) 
ggsave("./Paper1/Figures/polySim_PCA_qnorm.png", ggpcan, width = 9, height = 8 )


pcaDatc <- as.data.frame(apply(simdat, 2, FUN = MFtrans ))
pcplotdatc <- plotPCA(pcaDatc, annodf, scale = FALSE)
ggpcac <- ggplot(pcplotdatc[[1]], aes(x= PC1, y= PC2, color=factor(Plate), label = SampleID)) + 
  geom_point( size = 2)   + 
  # geom_text( size = 1.7) +
  ggtitle(paste0("Proportion of variance =", round(pcplotdatc[[2]], 2)))
# xlim(-10, 50) + ylim(-10, 10)
ggsave("./Paper1/Figures/polySim_PCA_CLR.png", ggpcac, width = 9, height = 8 )
