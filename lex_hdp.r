## LEX_HDP.R -- functions to simulate from a lexical-distributional learner as described
## by Feldman et al. 2009


lexhdp <- function(w, alpha, beta, mu0, sigma0, nu0) {
    ## initialize labels
    z <- sapply(w, function(x) {length(x)})
    
    
