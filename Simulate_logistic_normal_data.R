###############################################################################
###############################################################################
## Simulate logistic normal data.  As per Dean's recommendation I start with 
## simulating multivariate log-normal data as a basis and then imposing the 
## sum constraint through the closure operation
###############################################################################
###############################################################################

###############################################################################
# The multivariate log-normal  is related to the multivariate normal 
# distribution through the following identities:
# X is multi-variate normal N_D(mu, Sigma), then Y = exp(X) is multivariate 
# log-normal with:
# E[Y]_i = exp(mu_i + Sigma_ii/2) and,
# Var[Y]_ij = exp(mu_i + mu_j + 1/2(Sigma_ii + Sigma_jj)) * (exp(Sigma_ij) - 1)
# 
# these properties can be used to generate multivariate log-normal samples 
# which can then be used to generate logistic normal data through closure
###############################################################################

biocLite("limma")

rlogistic <- function(N, D, mu, sig){
  #########################################################
  ### This function generates multivariate logistic normal
  ### data given:
  ### N = the number of samples
  ### D = the number of components
  ### mu = vector of component means in units of mn norm
  ### sig = the covariance matrix in units of the multi-
  ###     variate normal (see above for conversion)
  ### Returns a matrix with samples as columns and 
  ### components as rows
  #########################################################
  
  #Create matrix of N(0,1) random variables
  Z <- matrix(rnorm(N*D), nrow = N, ncol = D)
  
  #Factor sigma into Q^TQ
  ev <- eigen(sig, symmetric = TRUE)
  
  #Get Q from eigen values and eigen vectors
  Q <- ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors)
  
  #Create matrix of component means
  mu <- matrix(mu, N, D, byrow = TRUE)
  
  X <- Z %*% Q + mu 
  
  #Convert to log-normal
  Y <- exp(X)
  
  #Convert to logistic normal through closure
  U <- apply(Y, 1, function(x) x/sum(x))
  
  return(U)
}



#Start with low dimension composition D=3


#Create covariance matrix sigma
sig <- matrix(1, D, D)
sig[2, 1] <- sig[1, 2] <- .5
sig[3, 1] <- sig[1, 3] <- .2
sig[2, 3] <- sig[3, 2] <- 0

X <- rlogistic(N = 100, D = 3, mu = c(1, 3, 6), sig)

cor(t(X))
cov(t(X))

svd(X)

#Estimate mu and sigma from real data
dfraw <- readHTG("Data/Pfizer Sample Input OBP parsed 18-Nov-2016.xls",
                 include.probe = FALSE)

dfclo <- apply(dfraw, 2, function(x) x/sum(x))

ecov <- cov(t(dfclo))
