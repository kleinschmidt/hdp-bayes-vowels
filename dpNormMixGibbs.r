### dpNormMixGibbs.r
###
### Implementation of a Dirichlet Proccess Normal Mixture Model Gibbs sampler
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
###   Oct 20 2009 -- cloned from finiteNormMixtureGibbs.r
###   Oct 22 2009 -- change component/label representation (Z=[1, 2, 1, 1, ...])
###   Oct 23 2009 -- added ellipse/contour plotting to 2D visualization


library(bayesm)
library(mnormt)
library(ellipse)


## helper functions for visualization:
# plot side-by-side histograms (good for 1-d case)
compHist <- function(pts, comps) {
  h <- NULL
  h0 <- hist(pts, plot=FALSE)
  for (c in unique(comps)) {
    hn <- hist(pts[comps==c], breaks=h0$breaks, plot=FALSE)
    h <- rbind(h, hn$counts)
  }
  colnames(h) <- h0$mids
  barplot(h, beside=TRUE)
}

# 2-d scatterplot, labled by component assignment
comp2dplot <- function(pts, labels, comps, means, precisions) {
  plot(pts, type="n")
  for (co in comps) {
    ppts <- pts[labels==co, ]
    points(ppts[,1], ppts[,2], pch=as.character(co))
    lines(ellipse(x=solve(precisions[1:2,1:2,which(comps==co)]),
                  centre=means[which(comps==co),1:2]))
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
precision <- 1
C <- array(rep(diag(precision, d), g), dim=c(d,d,g))

# underlying component means
W <- array(data=rep(c(1,-1), d), dim=c(g,d))

# generate observations from the mixture model
X <- rbind(rmnorm(n/2, mean=W[1,], varcov = solve(C[,,1])),
           rmnorm(n/2, mean=W[2,], varcov = solve(C[,,2])) )
dim(X) <- c(n, d)


# Hyperparameters
alpha <- 1
K <- rep(1, g)
R <- rep(5, g)

# initial component assignments
Z <- apply(rmultinom(n, size=1, prob=rdirichlet(rep(alpha, g))),
           2, function(x) {which(x==1)})
comps <- 1:g

## initialize component mean and precision variables
S <- array(data=rep(0, d*d*g), dim=c(d,d,g))
MU <- array(data=rep(0, d*g), dim=c(g,d))


# Number of full sweeps through the sampler
nIter <- 100


# open a plot window
#x11()

## GIBBS SAMPLER ###########################################
for (iter in 1:nIter) {
    ## get component counts
    N <- sapply(comps, function(x) {sum(Z==x)})

    ## update component parameters
    for (ii in comps) {
        x <- array(X[Z==ii, ], c(N[ii],d))
        xbar <- colMeans(x)

        ## sample from the conditional posterior for the precision
        v <- 0
        for (jj in 1:N[ii]) {
            v <- v + (x[jj, ]-xbar)%*%t(x[jj, ]-xbar)
        }
        Cstar <- solve(solve(C[,,ii]) +
                       v +
                       (N[ii]*R[ii]/(N[ii]+R[ii]) *
                        (xbar-W[ii,]) %*% t(xbar-W[ii,])))

        S[,,ii] <- rrwishart(N[ii]+R[ii], Cstar)
        
        
        ## sample from the conditional posterior for the mean
        w <- (N[ii]*xbar + K[ii]*W[ii,]) / (N[ii]+K[ii])
        MU[ii,] <- rmnorm(1, mean=w, varcov=solve( (N[ii]+K[ii])*S[,,ii] ))
    }
    
    ## re-assign observations to categories 
    tau <- rep(0,length(comps))
    for (jj in 1:n) {
        for (ii in comps) {
            tau[ii] <- (N[ii] - (Z[jj]==ii) + alpha/g) / (n - 1 + alpha) *
                dmnorm(X[jj,], mean=MU[ii,], varcov=solve(S[,,ii]))
        }
        ##tau <- tau/sum(tau)    ###### rmultinom autmoatically normalizes probs 
        Z[jj] <- which(rmultinom(n=1, size=1, prob=tau)==1)
    }

    ## visualization
    if (d==1) {
        compHist(X,Z)
    } else {
        comp2dplot(X[,1:2], Z, comps, MU, S)
    }
}



