
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
