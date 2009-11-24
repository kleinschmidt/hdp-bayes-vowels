## LEX_HDP.R -- functions to simulate from a lexical-distributional learner as described
## by Feldman et al. 2009


lexhdp <- function(w, alpha, beta, mu0, sigma0, nu0) {
    ## initialize labels
    z <- sapply(w, function(x) {length(x)})
    
    ## "data structures"
    ## observed words:
    ## w[[j]] [i]   value of position i of word j
    ## z[j]         lexical label of word j
    ##
    ## lexicon:
    ## lexlab[k]    ID of kth lexeme (used for mapping between z and l)
    ## lexlen
    ## numlex       total number of lexical items
    ## l[[k]] [i]   phoneme label of segment i of lexeme k (list of vectors of "h"s)
    ## Nl[k]        number of words labeled as kth lexeme
    ## wk[k] [n]    observations associated with kth lexeme (vector of "j"s)
    ##
    ## phone inventory:
    ## phonlab[h]   ID of hth phoneme (used for mapping between l and p)
    ## numphon      total number of phonemes
    ## Np[h]        number of segments associated with hth phoneme
    ## lh[[h]]      segments associated with hth phoneme: list of sum(Nl) by 2 matrices

    for (iter in 1:nIter) {
        ## first sweep: update lexicon
        for (k in 1:numlex) {
            for (i in 1:length(l[[k]])) {
                obs <- w[[
