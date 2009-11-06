
library(ars)
## import some display functions...
source("displayFcns.r")

## density function for the conditional posterior on the component
## precision shape parameter
dBetaPosterior <- function(Beta, S, w) {
    if (length(Beta) < 1) {
        return(NULL)
    } else {
        k <- dim(S)[3]
        beta <- Beta[1]
        p <- gamma(beta/2)^(-k) *
            exp(-1/(2*beta)) *
                (beta/2)^((2*beta-3)/2) *
                    prod(S*w)^(beta/2) *
                        exp(-beta*w*sum(S)/2)
        return(c(p, dBetaPosterior(Beta[-1], S, w)))
    }
}

## Sample from the conditional posterior on beta (component precision prior shape) using
## ARS
rBetaPosterior <- function(num=1, S, w) {
    k <- dim(S)[3]
    h <- function(y) {
        1.5*log(2) - 0.5*y - k*log(gamma(exp(y)/2)) - 0.5/exp(y) +
          0.5*k*exp(y) * (y-log(2)+sum(log(S*w)-S*w)/k)
    }
    hprim <- function(y) {
        -0.5 + 0.5/exp(y) +
            0.5*k*exp(y) * (y - digamma(exp(y)/2) + 1-log(2) + sum(log(S*w)-S*w)/k)
    }
    return(exp(ars(n=num, f=h, fprima=hprim)))
}

## Sample from the conditional posterior on the category spread parameter alpha
## using ARS
rAlphaPosterior <- function(num=1, N) {
    n <- sum(N)
    k <- length(N)
    h <- function(y)
        {y*(k-0.5) - 0.5*exp(-y) -
             sapply(y, function(x) {sum(log(exp(x)+seq(0,n-1)))} ) }
    hprim <- function(y)
        {k-0.5 + 0.5*exp(-y) -
             exp(y)*sapply(y, function(x) {sum(1/(exp(x)+seq(0,n-1)))} ) }
    return(exp(ars(n=num, f=h, fprima=hprim, x=c(-2,0,2))))
}



## Gibbs sampler for parameter-free DP Gaussian mixture model

## data is in X, dim(X) = c(nObs,nDim)
meanX <- mean(X)
varX <- var(X)
precisionX <- solve(varX)

lambda <- rmnorm(n=1, mean=colMeans(X), sd=sd(X))
r <- rgamma(n=1, shape=1, scale=solve(varX))
MU <- rnorm(n=1, mean=lambda, sd=solve(sqrt(r)))

beta <- solve(rgamma(n=1, shape=1, scale=1))
w <- rgamma(n=1, shape=1, scale=varX)
S <- rgamma(n=1, shape=beta, rate=w)

Z <- rep(1,nObs)
comps <- 1
N <- nObs

## initialize hyperparameters
for (iter in 1:nIters) {
    ## update hyperparameters
    lambda <- rnorm(n=1,
                    mean=(meanX*precisX + r*sum(MU))/(precisX+k*r),
                    sd=solve(sqrt(precisX + k*r)))
    r <- rgamma(n=1,
                shape=(k+1),
                rate=(varX + sum((MU-lambda)^2))/(k+1))

    w <- rgamma(n=1,
                shape=(k*beta+1),
                rate=(precisX + beta*sum(S))/(k*beta+1))
    beta <- rBetaPosterior(num=1, S, w)

    alpha <- rAlphaPosterior(num=1, N)

    ## update component parameters
    for (co in 1:k) {
        Xco <- X[Z==comps[co]]
        MU[co] <- rnorm(n=1,
                        mean=(mean(Xco)*N[co]*S[co] + lambda*r) / (N[co]*S[co] + r),
                        sd=solve(sqrt(N[co]*S[co] + r)) )
        S[co] <- rgamma(n=1,
                        shape=(beta+N[co]),
                        rate=(w*beta + sum((Xco-MU[co])^2))/(beta+N[co]) )
    }

    ## update latent variables (category/component assignments)

    for (j in 1:nObs) {
        probs <- rep(0,k+1)
        for (co in 1:k) {
            probs[co] <- dnorm(X[j], mean=MU[co], sd=solve(sqrt(S[co]))) *
                (N[co]-(Z[j]==comps[co]))
        }
        ## draw component parameter values from their priors for unrepresented comps
        mu <- rnorm(n=1,mean=lambda,sd=solve(sqrt(r)))
        s <- rgamma(n=1,shape=beta, rate=w)
        probs[k+1] <- dnorm(X[j],
                            mean=mu,
                            sd=solve(sqrt(s))) * alpha
        newComp <- which(rmultinom(n=1, size=1, prob=probs)==1)
        if (newComp > nComps) {
            ## new component
            nComps <- nComps + 1
            ## find a new label
            comps <- append(comps, min(setdiff(1:(max(comps)+1), comps)))
            Znew <- tail(comps,1)
            ## store initial component precision and mean
            S <- c(S,s) #array(c(S, solve(sigma)), c(d,d,nComps))
            MU <- c(MU,mu) #rbind(MU, mu)
            ## update counts
            N <- append(N, 1)
            N[comps==Z[obs]] <- N[comps==Z[obs]] - 1
        } else {
            ## old component
            ## just update counts:
            Znew <- comps[newComp]
            N[newComp] <- N[newComp] + 1
            N[comps==Z[obs]] <- N[comps==Z[obs]] - 1
        }

        ## Delete a component if its last observation is removed
        if (N[comps==Z[obs]] == 0) {
            whichComp <- which(comps==Z[obs])
            comps <- comps[-whichComp]
            nComps <- nComps - 1
            N <- N[-whichComp]
            MU <- array(data=MU[-whichComp,], dim=c(nComps,d))
            S <- array(data=S[,,-whichComp], dim=c(d,d,nComps))
        }

        ## Record the new component assignment for this observation
        Z[obs] <- Znew
    }

    ## visualization
    compHist(X, Z, comps, N, MU, S)
}
