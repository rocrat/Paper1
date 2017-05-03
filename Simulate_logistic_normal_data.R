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



#Start with low dimension composition D=3

#generate d N(0,1) vectors 
nlist <- list()

#number of components
D <- 3
#number of samples
N <- 1000


#Create matrix of N(0,1) random variables
Z <- matrix(rnorm(N*D), nrow = N, ncol = D)


#Create covariance matrix sigma
sig <- matrix(1, D, D)
sig[2, 1] <- sig[1, 2] <- .5
sig[3, 1] <- sig[1, 3] <- .2
sig[2, 3] <- sig[3, 2] <- 0

#Factor sigma into Q^TQ
ev <- eigen(sig, symmetric = TRUE)

#Get Q from eigen values and eigen vectors
Q <- ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors)

#Create matrix of component means
mu <- matrix(c(1, 4, 6), N, D, byrow = TRUE)

X <- Z %*% Q + mu 

cor(X)


sigma <- matrix(0, nrow = 3, ncol = 3)
diag(sigma) <- 1
tz <- mvrnorm(10, rep(0, 3), sigma)


