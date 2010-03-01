### feldman2009.r --- implementations of the various models from Feldman,
### Griffiths, and Morgan (2009) paper on HDP learning of lexical and
### phonetic categories


################################################################################
### Toy 1d model ###############################################################
################################################################################

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

makeW <- function(l, Nl) {
    phonlab <- c(1,2,3,4)
    numphon <- 4

    phonmeans <- c(-5, -1, 1, 5)
    phonvar <- c(1,1,1,1)

    lexlen <- sapply(l, length)

    w <- NULL
    for(n in 1:length(l)) {
        w <- rbind(w, lexgen(Nl[n], l[[n]], phonlab, phonmeans, phonvar, max(lexlen)))
    }
    return(w)
}

lnoMinPairs <- list(c(1,2), c(4,3), c(1,4,1), 4)
NlnoMinPairs <- c(200, 200, 100, 100)

w_noMinPairs <- makeW(lnoMinPairs, NlnoMinPairs)


lminPairs <- list(c(1,2), c(1,3), c(4,2), c(4,3), c(1,4,1), c(4))
NlminPairs <- rep(100, 6)

w_minPairs <- makeW(lminPairs, NlminPairs)

source("lex_hdp.r")



#cat("----------------------- 1-D model with no minimum pairs -----------------------------")
#ms_noMinPairs <- lexhdp(w_noMinPairs, nIter=30)
#cat("------------------------ 1-D model with minimum pairs -------------------------------")
#ms_minPairs <- lexhdp(w_minPairs, nIter=30)




################################################################################
## test different mixtures of words...##########################################
################################################################################


for (p in seq(.05,.45,.10)) {
    cat("----------------------- 1-D model with probabilistic sparseness of", p, " percent---------------------\n")
    Nl = c(200*(1-p), 200*p, 200*p, 200*(1-p), 100, 100)
    w = makeW(lminPairs, Nl)
    ms_probSparse[(p+.05)/.1] <- list(w=w, Nl=Nl, p=p, output=lexhdp(w, nIter=30))
    
    save(list=ls(all=TRUE), file='1dsparseness.RData')
}

