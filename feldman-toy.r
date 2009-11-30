mu0 <- 0
s0 <- 1
nu0 <- 1.001

phonlab <- c(1,2,3,4)
numphon <- 4

phonmeans <- c(-5, -1, 1, 5)
phonvar <- c(1,1,1,1)

lexlab <- c(1,2,3,4)
numlex <- 4
lexlen <- c(2,2,3,1)
l <- list(c(1,2), c(4,3), c(1,4,1), 4)

lh <- list(matrix(c(1,1, 3,1, 3,3), byrow=TRUE, ncol=2),
           matrix(c(1,2), byrow=TRUE, ncol=2),
           matrix(c(2,2), byrow=TRUE, ncol=2),
           matrix(c(2,1, 3,2, 4,1), byrow=TRUE, ncol=2))

Nl <- c(200, 200, 100, 100)
Np <- c(400, 200, 200, 400)


## generate observations from a lexical item
lexgen <- function(n, lexeme, phonlab, mu, s, width=0) {
    nc <- max(width, length(lexeme))
    output <- matrix(rep(NA, n*nc), ncol=nc)
    for (p in 1:length(lexeme)) {
        phon <- which(phonlab==lexeme[p])
        output[,p] <- rnorm(n, mean=mu[phon], sd=sqrt(s[phon]))
    }
    return(output)
}

#X <- matrix(rep(NA, sum(Nl)*max(lexlen)), ncol=max(lexlen))
w <- NULL
z <- NULL
for(n in 1:numlex) {
    w <- rbind(X, lexgen(Nl[n], l[[n]], phonlab, phonmeans, phonvar, max(lexlen)))
    z <- c(z, rep(n, Nl[n]))
}

