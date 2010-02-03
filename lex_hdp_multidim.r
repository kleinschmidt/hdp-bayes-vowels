## LEX_HDP.R -- functions to simulate from a lexical-distributional learner as described
## by Feldman et al. 2009


lexhdp <- function(w, nIter=10, r=5, Alpha=1, Beta=1, mu0=0, s0=1, nu0=1.001, OUTPUT="FULL") {
    ## "data structures"
    ## observed words:
    ## w[d,i,j]     value of dimension d at position i of word j (NULL entires for short words)
    ## numwords     length of the corpus
    ## z[j]         lexical label of word j
    ##
    ## lexicon:
    ## lexlab[k]    ID of kth lexeme (used for mapping between z and l)
    ## lexlen       length (in segments) of each lexeme
    ## numlex       total number of lexical items
    ## l[[k]] [i]   phoneme label of segment i of lexeme k (list of vectors of phone IDs)
    ## Nl[k]        number of words labeled as kth lexeme
    ## wk[[k]] [n]  observations associated with kth lexeme (vector of "j"s)
    ##
    ## phone inventory:
    ## phonlab[h]   ID of hth phoneme (used for mapping between l and p)
    ## numphon      total number of phonemes
    ## Np[h]        number of segments associated with hth phoneme
    ## lh[[h]]      segments associated with hth phoneme: list of matrices of [lexID,pos]s
    ## 
    ##
    ## other parameters:
    ## Alpha        concentration parameter for phone inventory
    ## Beta         concnetration parameter for lexicon
    ## r            number of new samples to draw to estimate prior propability

    ## initialize phone/lexeme inventories
    z <- apply(w[1,,], 2, function(x) {sum(!is.na(x))})
    numwords <- dim(w)[3]

    lexlab <- unique(z)
    lexlen <- lexlab          # okay since z is chosen to be length
    numlex <- length(lexlab)

    l <- lapply(lexlen, function(x) {rep(1, x)})
    Nl <- sapply(lexlab, function(x) {sum(z==x)})
    wk <- lapply(lexlab, function(x) {which(z==x)})

    phonlab <- 1
    numphon <- 1
    dims <- nrow(w)
    Np <- sum(lexlen)
    lh <- list(cbind(unlist(lapply(1:numlex, function(x) {rep(lexlab[x],lexlen[x])})),
                     unlist(lapply(lexlen, function(x) {1:x}))))


    modelState <- vector("list", nIter)
    

    for (iter in 1:nIter) {
        cat("ITERATION NUMBER", iter, "\n")
        ## first sweep: update word labels
        ## for each word in w:
        for (jj in 1:numwords) {
            ## hold-out w:
            ## get value for jjth word
            kkz <- decodeZ(z[jj], lexlab)
            obs <- w[,1:lexlen[kkz],jj]
            ## remove jjth word from wk, Nl
            wk[[kkz]] <- wk[[kkz]] [wk[[kkz]]!=jj]
            Nl[kkz] <- Nl[kkz] - 1

            ## if Nl is 0, set it to Beta/(r+1)
            if (Nl[kkz]==0) {
                heldOutLastObs <- TRUE
                Nl[kkz] <- Beta/(r+1)
            } else {
                heldOutLastObs <- FALSE
            }
                

            ## calculate likelihoods:
            probs <- vector("numeric", numlex+r)
            ## for each lexeme in l:
            for (kk in 1:numlex) {
                if(lexlen[kk]==lexlen[kkz]) {
                    ## get the observed words associated with k
                    ##lexObs <- w[wk[[kk]],1:lexlen[kk]]
                    ## calculate/store the likelihood(obs given lexObs)
                    probs[kk] <- wordlhood(obs, l[[kk]], phonlab, lh, w, wk, lexlab, mu0, nu0, s0)
                }
            }
            ## calculate/store the prior prob of obs
            newlexemes <- rlex(r, len=lexlen[kkz], phonlab, Np)
            #print(newlexemes)
            probs[(1:r)+numlex] <- apply(newlexemes, 1,
                                         function(l) {wordlhood(obs, l, phonlab, lh, w, wk, lexlab, mu0, nu0, s0)})
            
            ## multiply by Ns/beta
            ## sample new p from multinomial distribution
            newkk <- which(rmultinom(1, 1, probs*c(Nl,rep(Beta/(r+heldOutLastObs),r)))==1)

            ## correct for treating a totally held out component as unexpressed
            if (heldOutLastObs) {Nl[kkz] <- 0}
            
            ## update counts/labels
            if (newkk > numlex) {
                ## new component
                ## generate a new comp ID
                lexlab <- updateIDs(lexlab)
                newID <- tail(lexlab,1)
                ## store the new sampled lexeme
                newlex <- newlexemes[newkk-numlex, ]
                l[[numlex+1]] <- newlex
                cat("  Adding new lexeme with ID ", newID,
                    " (word number ", jj, ")\n", sep="")
                cat("    l[", numlex, "] =", newlex, "\n", sep=" ")
                ## update phoneme->segment mappings and counts
                
                lh <- updateLh(lh, newlex, newID, numphon, phonlab)
                Np <- updateNp(Np, newlex, phonlab)
                ## and add count/length entries
                Nl <- append(Nl, 1)
                lexlen <- append(lexlen, length(newlex))
                numlex <- numlex + 1
                ## ...as well as a lexeme->observation map entry
                wk[[numlex]] <- jj
            } else {
                ## old component
                newID <- lexlab[newkk]
                Nl[newkk] <- Nl[newkk] + 1
                wk[[newkk]] <- append(wk[[newkk]], jj)
            }

            ## record the new label
            z[jj] <- newID

            ## check to see if the old lexeme is empty, and delete it if so
            if (Nl[kkz] == 0) {
                cat("  Removing lexeme ID number ", lexlab[kkz],
                    " (Nl = ", Nl[kkz],
                    ", word number ", jj, ")\n", sep="")
                for (hh in decodeZ(unique(l[[kkz]]), phonlab)) {
                    lh[[hh]] <- matrix(lh[[hh]] [lh[[hh]][,1] != lexlab[kkz], ], ncol=2)
                }
                Np <- Np - apply(sapply(l[[kkz]], function(p) {phonlab==p}), 1, sum)
                lexlab <- lexlab[-kkz]
                lexlen <- lexlen[-kkz]
                numlex <- numlex - 1
                l <- l[-kkz]
                Nl <- Nl[-kkz]
                wk <- wk[-kkz]
            }
            
            
        } ## end: first sweep


        ## second sweep: update lexical segment labels
        ## for each segment in each lexical item:
        for (kk in 1:numlex) {
            lexobs <- getLexObs(kk, lexlen, w, wk)
            for (ii in 1:length(l[[kk]])) {
                ## get the observations associated with this segment
                obs <- lexobs[,ii]
                
                ## hold out these observations
                hhl <- decodeZ(l[[kk]][ii], phonlab)
                lh[[hhl]] <- holdoutLh(lh[[hhl]], kk, ii, lexlab)
                Np[hhl] <- Np[hhl] - 1

                ## get the probability for each phone category and the prior
                ## use log lhood because of the higher sample sizes... 
                probs <- sapply(c(1:numphon, 0),
                                function(h) phonlhood(obs, h, log=TRUE, Np, lh, w, wk, lexlab, mu0, nu0, s0))
                ## multiply by Ns/alpha
                probs <- probs - max(probs) + log(c(Np, Alpha))
                ## sample new p from multinomial distro
                newhh <- which(rmultinom(n=1, size=1, prob=exp(probs))==1)

                ## update labels/counts
                if (newhh > numphon) {
                    ## new component
                    ## generate new ID
                    phonlab <- updateIDs(phonlab)
                    newID <- tail(phonlab,1)
                    numphon <- numphon + 1
                    ## record the new ID
                    Np[newhh] <- 1
                    lh[[newhh]] <- matrix(c(lexlab[kk], ii), ncol=2)
                    cat("  Adding new phoneme with ID ", newID,
                        " (lexID = ", lexlab[kk], ", segment = ", ii, ")\n", sep="")
                } else {
                    ## old component
                    newID <- phonlab[newhh]
                    Np[newhh] <- Np[newhh] + 1
                    lh[[newhh]] <- rbind(lh[[newhh]], c(lexlab[kk], ii))
                }

                ## record new label
                l[[kk]] [ii] <- newID

                ## remove old comp if it's empty
                if (Np[hhl] == 0) {
                    ## a helpful message!
                    cat("  Removing phoneme ID number ", phonlab[hhl],
                        " (Np = ", Np[hhl],
                        ", lexID = ", lexlab[kk], ", segment = ", ii, ")\n", sep="")
                    lh <- lh[-hhl]
                    Np <- Np[-hhl]
                    numphon <- numphon - 1
                    phonlab <- phonlab[-hhl]
                }
            }
        } ## end: second sweep

        if (OUTPUT=="FULL") {
            modelState[[iter]] <-
              list(z=z, lexlab=lexlab, lexlen=lexlen, numlex=numlex, l=l, Nl=Nl, wk=wk,
                   phonlab=phonlab, numphon=numphon, Np=Np, lh=lh, Alpha=Alpha, Beta=Beta, r=r)
        } else {
            modelState[[iter]] <-
              list(z=z, l=l)
        }
    } ## end: one iteration

    
    return(modelState)
    
    
    

}

#####################################################################
###                 HELPER FUNCTIONS                              ###
#####################################################################

## helper function to retrieve observations assigned to lexeme kk
getLexObs <- function(kk, lexlen, w, wk) {
    obs <- matrix(w[wk[[kk]],1:lexlen[kk]], ncol=lexlen[kk])
    return(obs)
}

## helper function to get observations assigned to phoneme hh
getPhonObs <- function(hh, lh, w, wk, lexlab) {
    obs <- unlist(apply(lh[[hh]],   # matrix of [lexID,i] rows
                        1,          # "margainalize" over rows
                        function(x) {w[wk[[decodeZ(x[1],lexlab)]], x[2]]}))
                                        # index into w matrix, getting proper wk index
                                        # using lexlab
    
    return(obs)
}

## convert stable component labels (ID number) into vector index
decodeZ <- function(z, lexlab) {
    sapply(z, function(l) {which(lexlab==l)})
}


## calculate the likelihood of a single word observation given a class
wordlhood <- function(obs, lexeme, phonlab, lh, w, wk, lexlab, mu0, nu0, s0) {
    lhood <- 1
    for (ii in 1:length(lexeme)) {
        hh <- decodeZ(lexeme[ii], phonlab)
        phobs <- getPhonObs(hh, lh, w, wk, lexlab)
        ## function from Feldman 2009:
        #Ncm <- sum(sapply(lh[[hh]][,1], function(x) {length(wk[[decodeZ(x,lexlab)]])}))
        Ncm <- length(phobs)
        if (Ncm < 1) {meanphobs <- 0}
        else {meanphobs <- mean(phobs)}
        mun <- ((nu0 * mu0) + (Ncm * meanphobs)) / (nu0 + Ncm)
        nun <- nu0 + Ncm
        nusn <- nu0*s0 + sum((phobs-meanphobs)^2) + nu0*Ncm/(nu0+Ncm)*(meanphobs-mu0)^2
        lhoodii <- 1/(beta(nun/2,0.5)*sqrt(nusn*(1+nun)/nun)) *
          (1 + (obs[ii]-mun)^2 / (nusn*(1+nun)/nun)) ^ (-(nun+1)/2)
        lhood <- lhood*lhoodii
    }
    return(lhood)
} 

## convert summary statistics for a phonetic category for use in calculating likelihood
phonstats <- function(hh, lh, w, wk, lexlab, mu0, nu0, s0) {
  # for each expression of this phoneme

  # theoretically, this will produce a dims-by-Ncm matrix of observations for this phone
  phonObs = matrix(unlist(apply(lh[[hh]], 1, function(x) {w[ , x[2], wk[[decodeZ(x[1],lexlab)]] ]} )), nrow=dims)
  # need to compute mean, number, and variance
  pnum = ncol(phonObs)
  pmean = rowMeans(phonObs)
  pcov = matrix(rowSums(apply(phonObs, 2, function(ob) {(ob-pmean) %*% t(ob-pmean)})), dims)

  pcov = s0 + pcov + nu0*pnum / (nu0+pnum) * (pmean-mu0) %*% t(pmean-mu0)
  pmean = (nu0*mu0 + pnum*pmean) / (nu0 + pnum)
  pnum = nu0 + pnum

  return(list(nuN = pnum, muN = pmean, sigmaN = pcov))
}

## word likelihood for multidimensional data (each observation has >1 dimensions)
mwordlhood <- function(obs, lexeme, phonlab, lh, w, wk, lexlab, mu0, nu0, s0) {
  llhood = 0
  obs = as.matrix(obs)
  for (ii in 1:length(lexeme)) {
    hh = decodeZ(lexeme[ii], phonlab)
    catstats = phonstats(hh, lh, w, wk, lexlab, mu0, nu0, s0)
    llhood = llhood + with(catstats, mydmvt(obs[,ii], mu=muN, sigma=sigmaN, nu=nuN, d=dims, log=TRUE))
  }
  return(llhood)
}

mydmvt <- function(x, mu, sigma, nu, d, log=FALSE) {
  x = as.matrix(x)
  logdf = function(X) {lgamma((nu+1)/d) - lgamma((nu+1-d)/2) - log(det(pi*sigma*(nu+1)/nu))/2 -
                         (nu+1)/2 * log(1 + t(X-mu) %*% solve(sigma*(nu+1)/nu) %*% (X-mu))}
  dens = apply(x, 2, logdf)
  if (log) return(dens)
  else return(exp(dens))
}

## multivariate gamma function
mgamma <- function(x, d, log=FALSE) {
  out = 0.25*d*(d-1) * log(pi) + sapply(x, function(x) sum(lgamma(x - seq(0,d-1)/2)))
  if (log) return(out)
  else return(exp(out))
}

mlgamma <- function(x, d) {mgamma(x,d,log=TRUE)}
  
  
  


## compute the likelihood of assigning given observations to phoneme hh.  returns
## the prior probability of the observations if hh=0.  returns log-likelihood if log=TRUE
phonlhood <- function(obs, hh=0, log=FALSE, Np, lh, w, wk, lexlab, mu0, nu0, s0) {
    if (hh==0 || Np[hh]==0) {phobs <- 0; Ncm <- 0}  # if prior specified or all obs held-out
    else {phobs <- getPhonObs(hh, lh, w, wk, lexlab); Ncm <- length(phobs)}
    ## function from Feldman 2009
    mun <- ((nu0 * mu0) + (Ncm * mean(phobs))) / (nu0 + Ncm)
    nun <- nu0 + Ncm
    nusn <- nu0*s0 + sum((phobs-mean(phobs))^2) + nu0*Ncm/(nu0+Ncm)*(mean(phobs)-mu0)^2
    obsmean <- mean(obs)
    n <- length(obs)
    if(log==FALSE) {
        lhood <- sqrt(nun/(n+nun)) / beta(nun/2, n/2) / nusn^(n/2) *
          (1 + (sum((obs-obsmean)^2) + (obsmean-mun)^2*(n*nun)/(n+nun))/nusn) ^ (-(nun+n)/2)
    } else {
        lhood <- 0.5*(log(nun) - log(n+nun)) - lbeta(nun/2, n/2) - (n/2)*log(nusn) -
          (nun+n)/2 * log(1 + (sum((obs-obsmean)^2) + (obsmean-mun)^2*(n*nun)/(n+nun))/nusn)
    }
    return(lhood)
}




## draw n new lexical items randomly
rlex <- function(num, len=1, phonlab, Np) {
                                        #print(phonlab)
    matrix(phonlab[apply(rmultinom(num*len, 1, Np), 2, function(x){which(x==1)})],
           ncol=len)
}

## update the phoneme-->lexeme segment mappings, and return the list of new mappings
updateLh <- function(Lh, newlexeme, lexID, numphon, phonlab) {
    for (hh in 1:numphon) {
        for (n in 1:length(newlexeme)) {
            if (newlexeme[n] == phonlab[hh]) {Lh[[hh]] <- rbind(Lh[[hh]], c(lexID,n))}
        }
    }
    return(Lh)
}

## NOT SURE THAT THIS NEEDS TO BE DONE HERE....
updateNp <- function(Np, newlexeme, phonlab) {
    Np + sapply(phonlab, function(p) {sum(newlexeme==p)})
}

## add a new, unique ID to a list of component labels, and return the list
updateIDs <- function(comps) {
    return(append(comps, min(setdiff(1:(max(comps)+1), comps))))
}


## takes a single entry from lh (aka lh[h]) and removes the entry corresponding
## to position ii of lexical item kk
holdoutLh <- function(lhh, kk, ii, lexlab) {
    matrix(lhh[(lhh[,1]!=lexlab[kk] | lhh[,2]!=ii), ], ncol=2)
}

ltolh <- function(l, lh, phonlab, lexlab, lexlen) {
    LH <- vector("list", numphon)
    for (hh in 1:numphon) {
        for (kk in 1:numlex) {
            matches <- which(l[[kk]] == phonlab[hh])
            LH[[hh]] <- rbind(LH[[hh]], cbind(rep(lexlab[kk],length(matches)), matches))
        }
    }
}


