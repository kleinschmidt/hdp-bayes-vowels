### finiteNormMixtureGibbs.r
###
### Implementation of a Gibbs sampler for a finite normal mixture model
### of (hopefully) arbitrary dimensions and numbers of components.  Also
### probably grossly inefficient, but it's a start.
###
### (This implementation uses the 'bayesm' and 'mnormt' libraries for
### the Wishart and multivariate normal distributions, repsectively.
### They're available from CRAN)
###
### The variable names follow the presentation in McLachlan & Peel (200)
### Chapter 4, which if you're at UMD can be found at
###   http://files.ling.umd.edu/locker/Statistics/Books/
###
### Dave Kleinschmidt
###   Oct 20 2009


library(bayesm)
library(mnormt)


## helper functions for visualization:
# plot side-by-side histograms (good for 1-d case)
compHist <- function(pts, comps) {
  h <- NULL
  h0 <- hist(pts, plot=FALSE)
  for (ii in 1:nrow(comps)) {
    hn <- hist(pts[comps[ii,]==1], breaks=h0$breaks, plot=FALSE)
    h <- rbind(h, hn$counts)
  }
  colnames(h) <- h0$mids
  barplot(h, beside=TRUE)
}

# 2-d scatterplot, labled by component assignment
comp2dplot <- function(pts, comps) {
  plot(pts, type="n")
  for (ii in 1:nrow(comps)) {
    ppts <- pts[comps[ii,]==1, ]
    points(ppts[,1], ppts[,2], pch=as.character(ii))
  }
}

# wrapper for bayesm rwishart function (which returns a list...)
rrwishart <- function(nu, V) {
  samp <- rwishart(nu,V)
  return(samp$W)
}


## Generate data points:
d <- 2      # number of dimensions
n <- 200    # number of data points

g <- 2      # number of components

# underlying component precisions
C <- c(2, 0, 0, 2, 2, 0, 0, 2)
dim(C) <- c(d, d, g)

# underlying component means
MU <- c(-1, 1, -1, 1)
dim(MU) <- c(d,g)

# generate observations from the mixture model
X <- rbind(rmnorm(n/2, mean=MU[1,], varcov = solve(C[,,1])),
           rmnorm(n/2, mean=MU[2,], varcov = solve(C[,,2])) )



# Hyperparameters
alpha <- 1
W <- MU
K <- rep(1, g)
R <- rep(5, g)

# initial component assignments
Z <- rmultinom(n, size=1, prob=rdirichlet(rep(alpha, g)))

# Number of full sweeps through the sampler
nIter <- 100


# open a plot window
x11()

## GIBBS SAMPLER ###########################################
for (iter in 1:nIter) {
  N <- rowSums(Z)

  ## update component parameters
  for (ii in 1:g) {
    x <- X[Z[ii,] == 1, ]
    xbar <- colMeans(x)

    # sample from the conditional posterior for the precision
    v <- 0
    for (jj in 1:N[ii]) {
      v <- v + (x[jj, ]-xbar)%*%t(x[jj, ]-xbar)
    }
    Cstar <- solve(solve(C[,,ii]) +
                   v +
                   (N[ii]*R[ii]/(N[ii]+R[ii]) *
                    (xbar-W[ii,]) %*% t(xbar-W[ii,])))

    S[,,ii] <- rrwishart(N[ii]+R[ii], Cstar)
    
    
    # sample from the conditional posterior for the mean
    w <- (N[ii]*xbar + K[ii]*W[ii,]) / (N[ii]+K[ii])
    MU[ii,] <- rmnorm(1, mean=w, varcov=solve( (N[ii]+K[ii])*S[,,ii] ))
  }
  
  ## re-assign observations to categories 
  tau <- rep(0,g)
  for (jj in 1:n) {
    for (ii in 1:g) {
      tau[ii] <- (N[ii] - Z[ii,jj] + alpha/g) / (n - 1 + alpha) *
        dmnorm(X[jj,], mean=MU[ii,], varcov=solve(S[,,ii]))
    }
    #tau <- tau/sum(tau)    ###### rmultinom autmoatically normalizes probs 
    Z[,jj] <- rmultinom(n=1, size=1, prob=tau)
  }

  # visualization
  comp2dplot(X,Z)
}


    
