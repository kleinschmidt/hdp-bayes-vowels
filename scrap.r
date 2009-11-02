f <- function(b, S, w) {
    if (length(b)==0) {
        return(NULL)
    } else {
        return (c(F(b[1], S, w), f(b[-1], S, w)))
    }
}
F <- function(b, S, w) {
    k <- dim(S)[3]
    return(gamma(b/2)^(-k) * exp(-1/2/b) * (b/2)^((k*b-3)/2) *
           prod( (S*w)^(b/2) * exp(-b*w/2*S) ))
}

g <- function(x) {log(x)}
gprim <- function(x) {1/x}
ginv <- function(y) {exp(y)}

fy <- function(y, fx, ginv, gprim, ...) {
    if (length(y)==0) {
        return(NULL)
    } else {
        return(c(Fy(y[1], fx, ginv, gprim, ...), fy(y[-1], fx, ginv, gprim, ...)))
    }
}

Fy <- function(y, fx, ginv, gprim, ...) {
    abs(1/gprim(ginv(y))) * fx(ginv(y), ...)
}

h <- function(y, S, w) {
    if (length(y)==0) {
        return(NULL)
    } else {
        return( c(H(y[1], S, w), h(y[-1], S, w)))
    }
}
H <- function(y, S, w) {
    k <- dim(S)[3]
    return(y+log((gamma(exp(y)/2))^(-k)) - 1/(2*exp(y)) +
           (k*exp(y)-3)/2 * (y-log(2)) +
           exp(y)/2 * sum(log(S*w) - S*w))
}


hprim <- function(y, S, w) {
    if (length(y)==0) {
        return(NULL)
    } else {
        return(c(Hprim(y[1], S, w), hprim(y[-1], S, w)))
    }
}
Hprim <- function(y, S, w) {
    k <- dim(S)[3]
    return(y - k*exp(y)/2*digamma(exp(y)/2) + 1/2/exp(y) +
           (k*exp(y)/2)*(y-log(2)) + (k*exp(y)-3)/2 +
           exp(y)/2 * sum(log(S*w) - S*w))
}
