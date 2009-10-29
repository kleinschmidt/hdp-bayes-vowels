
library(ars)

## density function for the conditional posterior on the component
## precision shape parameter
dBetaPosterior <- function(Beta, S, w) {
    if (length(Beta) < 1) {
        return(NULL)
    } else {
        k <- dim(S)[3]
        beta <- Beta[1]
        p <- gamma(beta/2)^(-k) *
            exp(1/(2*beta)) *
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
        -k*log(gamma(exp(y)/2)) - .5/exp(y) - .5*y - 1.5*log(2) +
            .5*exp(y) * (k*(y+log(2)) + sum(w*S + log(w*S)))
    }
    hprim <- function(y) {
        -k*exp(y)*0.5*digamma(0.5*exp(y)) + 0.5/exp(y) - 0.5 +
            0.5*exp(y)(k*(y+log(2)+1) + sum(w*S + log(w*S)))
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
