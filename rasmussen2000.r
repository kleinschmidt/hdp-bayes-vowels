
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

## density function for the conditional posterior on the concentration parameter
dAlphaPosterior <- function(alpha, N) {
    p <- alpha^(k-3/2) * exp(-1/2/alpha) * gamma(alpha) / gamma(n+alpha)
    return(p)
}
