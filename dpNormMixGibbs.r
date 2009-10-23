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


## helper functions ############################################################
## plot side-by-side histograms (good for 1-d case)
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

## 2-d scatterplot, labled by component assignment
comp2dplot <- function(pts, labels, comps, means, precisions) {
  plot(pts, type="n")
  for (co in comps) {
    ppts <- pts[labels==co, ]
    if (!is.matrix(ppts)) {ppts <- matrix(ppts, length(ppts)/2, 2)}
    points(ppts[,1], ppts[,2], pch=as.character(co))
    lines(ellipse(x=solve(precisions[1:2,1:2,which(comps==co)]),
                  centre=means[which(comps==co),1:2]))
  }
}


## wrapper for bayesm rwishart function (which returns a list...)
rrwishart <- function(nu, V) {
  samp <- rwishart(nu,V)
  return(samp$W)
}

#### INITIALIZATION ###############################################################

## Generate data points:
d <- 2      # number of dimensions
n <- 200    # number of data points
g <- 2      # number of components (for GENERATION)

## underlying component precisions
precision <- 4
C <- array(rep(diag(precision, d), g), dim=c(d,d,g))
## underlying component means
W <- array(data=rep(c(1,-1), d), dim=c(g,d))
## generate observations from the mixture model
X <- rbind(rmnorm(n/2, mean=W[1,], varcov = solve(C[,,1])),
           rmnorm(n/2, mean=W[2,], varcov = solve(C[,,2])) )
dim(X) <- c(n, d)

## initial component assignments
Z <- rep(1,n)
comps <- 1
N <- n

## initialize component mean and precision variables (to 0)
S <- array(data=rep(0, d*d), dim=c(d,d,1))
MU <- array(data=rep(0, d), dim=c(1,d))

## open a plot window
x11()

#### CONSTANTS ETC ##################################################################
## Number of full sweeps through the sampler
nIter <- 100
## Number of samples used to estimate p(x) for new components
nPriorSamples <- 5

## Hyperparameters
alpha <- 1
K <- 1
R <- 5
W <- colMeans(X)
C <- solve(cov(X))

######################################################################################
## GIBBS SAMPLER #####################################################################
######################################################################################
for (iter in 1:nIter) {
  ## get component counts (NOW UPDATED BELOW)
  #N <- sapply(comps, function(x) {sum(Z==x)})

  ## update component parameters #################################################
  for (co in comps) {
    nco <- which(comps==co)
    x <- array(X[Z==co, ], c(N[nco],d))
    xbar <- colMeans(x)

    ## sample from the conditional posterior for the precision
    v <- 0
    for (obs in 1:N[nco]) {
      v <- v + (x[obs, ]-xbar)%*%t(x[obs, ]-xbar)
    }
    Cstar <- solve(solve(C) +
                   v +
                   (N[nco]*R/(N[nco]+R) *
                    (xbar-W) %*% t(xbar-W)))
    
    S[,,nco] <- rrwishart(N[nco]+R, Cstar)
    
    ## sample from the conditional posterior for the mean
    w <- (N[nco]*xbar + K*W) / (N[nco]+K)
    MU[nco,] <- rmnorm(1, mean=w, varcov=solve( (N[nco]+K)*S[,,nco] ))
  }
  
  ## re-assign observations to categories #########################################
  tau <- rep(0,length(comps))
  for (obs in 1:n) {
    ## get multinomial probabilities
    ## for expressed components:
    for (k in 1:length(comps)) {
      co <- comps[k]
      tau[k] <- (N[comps==co]-(Z[obs]==co)) / (n-1 +alpha) *
        dmnorm(X[obs,], mean=MU[comps==co,], varcov=solve(S[,,comps==co]))
    }
    ## for a new component:
    samps <- rep(0,nPriorSamples)
    for (j in 1:nPriorSamples) {
      sigma <- solve(rrwishart(R, C))
      mu <- as.vector(rmnorm(mean=W, varcov=sigma/K))
      samps[j] <- dmnorm(X[obs,], mean=mu, varcov=sigma)
    }
    tau[length(comps)+1] <- alpha / (n-1+alpha) * mean(samps)

    ## draw the new component from a multinomial distribution
    newComp <- which(rmultinom(n=1, size=1, prob=tau)==1)

    if (newComp > length(comps)) {
      ## new component
      ## find a new label
      comps <- append(comps, min(setdiff(1:(max(comps)+1), comps)))
      Znew <- tail(comps,1)
      ## store initial component precision and mean
      S <- array(c(S, solve(sigma)), c(d,d,length(comps)))
      MU <- rbind(MU, mu)
      ## update counts
      N[comps==Znew] <- 1
      N[comps==Z[obs]] <- N[comps==Z[obs]] - 1
    } else {
      ## old component
      ## just update counts:
      Znew <- which(comps==newComp)
      N[comps==Znew] <- N[comps==Znew] + 1
      N[comps==Z[obs]] <- N[comps==Z[obs]] - 1
    }

    ## Delete a component if its last observation is removed
    if (N[comps==Z[obs]] == 0) {
      whichComp <- which(comps==Z[obs])
      comps <- comps[-whichComp]
      N <- N[-whichComp]
      MU <- MU[-whichComp,]
      S <- S[,,-whichComp]
    }
  
    ## record new label
    Z[obs] <- Znew
  }

  
  ## visualization
  if (d==1) {
    compHist(X,Z)
  } else {
    comp2dplot(X[,1:2], Z, comps, MU, S)
  }
}



