## LEX_HDP.R -- functions to simulate from a lexical-distributional learner as described
## by Feldman et al. 2009


lexhdp <- function(w, alpha, beta, mu0, sigma0, nu0) {
    ## initialize labels
    z <- sapply(w, function(x) {length(x)})
    
    ## "data structures"
    ## observed words:
    ## w[j,i]       value of position i of word j (NULL entires for short words)
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

    for (iter in 1:nIter) {
        ## first sweep: update word labels
        ## for each word in w:
        for (jj in 1:numwords) {
            ## hold-out w:
            ## get value for jjth word
            kkz <- decodeZ(z[jj], lexlab)
            obs <- w[jj,1:lexlen[kkz]]
            ## remove jjth word from wk, Nl
            wk[[kkz]] <- wk[[kkz]][-which(wk[[kkz]]==jj)]
            Nl[kkz] <- Nl[kkz] - 1
            ## NEED TO CHECK IF Nl[kkz] == 0 NOW...

            ## calculate likelihoods:
            probs <- vector("numeric", numlex+r)
            ## for each lexeme in l:
            for (kk in 1:numlex) {
                if(lexlen[kk]==lexlen[kkz]) {
                    ## get the observed words associated with k
                    ##lexObs <- w[wk[[kk]],1:lexlen[kk]]
                    ## calculate/store the likelihood(obs given lexObs)
                    probs[kk] <- wordlhood(obs, l[[kk]])
                }
            }
            ## calculate/store the prior prob of obs
            newlexemes <- rlex(r, len=lexlen[kkz])
            probs[(1:r)+numlex] <- apply(newlexemes, 1, function(l) {wordlhood(obs, l)})
            
            ## multiply by Ns/beta
            ## sample new p from multinomial distribution
            newkk <- which(rmultinom(1, 1, probs*c(Nl,rep(Beta/r,r)))==1)
            
            ## update counts/labels
            if (newkk > numlex) {
                ## new component
                ## generate a new comp ID
                lexlab <- updateIDs(lexlab)
                newID <- tail(lexlab,1)
                ## store the new sampled lexeme
                newlex <- newlexemes[newkk-numlex, ]
                l[[numlex+1]] <- newlex
                ## update phoneme->segment mappings and counts
                lh <- updateLh(lh, newlex, newID)
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
                for (hh in decodeZ(unique(l[[kkz]]), phonlab)) {
                    lh[[hh]] <- lh[[hh]] [lh[[hh]][,1] != lexlab[kkz], ]
                }
                lexlab <- lexlab[-kkz]
                lexlen <- lexlen[-kkz]
                numlex <- numlex - 1
                l <- l[-kkz]
                Nl <- Nl[-kkz]
                wk <- wk[-kkz]
            }
        }

        ## second sweep: update lexical segment labels
        ## for each segment in each lexical item:
        for (kk in 1:numlex) {
            lexobs <- getLexObs(kk)
            for (ii in 1:length(l[[k]])) {
                ## get the observations associated with this segment
                obs <- lexobs[,ii]
                
                ## hold out these observations
                hhl <- decodeZ(l[[kk]][ii], phonlab)
                lh[[hhl]] <- holdoutLh(lh[[hhl]], kk, ii)
                Np[hhl] <- Np[hhl] - 1

                probs <- vector("numeric", numphon+1)
                ## for each phoneme c:
                for (hh in 1:numphon) {
                    ## get the [k,i] tuples associated with this phone
                    ## get the observations associated with each [k,i] tuple
                    ## calculate/store likelihood(obs given phonObs)
                    probs[hh] <- phonlhood(obs, hh)
                }
                ## calculate/store prior probability of obs
                ## multiply by Ns/alpha
                ## sample new p from multinomial distro
                ## update labels/counts
            }
        }

    }

    ## helper function to retrieve observations assigned to lexeme kk
    getLexObs <- function(kk) {
        obs <- w[wk[[kk]],1:lexlen[kk]]
        return(obs)
    }

    ## helper function to get observations assigned to phoneme hh
    getPhonObs <- function(hh) {
        obs <- unlist(apply(lh[[hh]],   # matrix of [nk,i] rows
                            1,          # "margainalize" over rows
                            function(x) {w[wk[[which(lexlab==x[1])]], x[2]]}))
                                        # index into w matrix, getting proper wk index
                                        # using lexlab
        return(obs)
    }

    ## convert stable component labels (ID number) into vector index
    decodeZ <- function(z, lexlab) {
        sapply(z, function(l) {which(lexlab==l)})
    }

    ## calculate the likelihood of a single word observation given a class
    wordlhood <- function(obs, lexeme) {
        lhood <- 1
        for (ii in 1:length(lexeme)) {
            hh <- decodeZ(lexeme[ii], phonlab)
            phobs <- getPhonObs(hh)
            ## function from Feldman 2009:
            #Ncm <- sum(sapply(lh[[hh]][,1], function(x) {length(wk[[decodeZ(x,lexlab)]])}))
            Ncm <- length(phobs)
            mun <- ((nu0 * mu0) + (Ncm * mean(phobs))) / (nu0 + Ncm)
            nun <- nu0 + Ncm
            nusn <- nu0*s0 + sum((phobs-mean(phobs))^2) + nu0*Ncm/(nu0+Ncm)*(mean(phobs)-mu0)^2
            lhoodii <- 1/(beta(nun/2,0.5)*sqrt(nusn*(1+nun)/nun)) *
                (1 + (obs[ii]-mun)^2 / (nusn*(1+nun)/nun)) ^ (-(nun+1)/2)
            lhood <- lhood*lhoodii
        }
        return(lhood)
    }

    ## compute the likelihood of assigning given observations to phoneme hh.  returns
    ## the prior probability of the observations if hh=0.  returns log-likelihood if log=TRUE
    phonlhood <- function(obs, hh=0, log=FALSE) {
        if (hh==0) {phobs <- 0; Ncm <- 0}
        else {phobs <- getPhonObs(hh); Ncm <- length(phobs)}
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
    rlex <- function(num, len=1) {
    matrix(phonlab[apply(rmultinom(num*len, 1, Np), 2, function(x){which(x==1)})],
               ncol=len)
    }

    ## update the phoneme-->lexeme segment mappings, and return the list of new mappings
    updateLh <- function(Lh, newlexeme, lexID) {
        for (hh in 1:numphon) {
            for (n in 1:length(newlexeme)) {
                if (newlexeme[n] == phonlab[hh]) {Lh[[hh]] <- rbind(Lh[[hh]], c(lexID,n))}
            }
        }
        return(Lh)
    }

    ## NOT SURE THAT THIS NEEDS TO BE DONE HERE....
#    updateNp <- function(Np, newlexeme, oldlexeme) {
#        Np + sapply(phonlab, function(p) {sum(newlexeme==p)}) - 
#             sapply(phonlab, function(p) {sum(oldlexeme==p)})
#    }

    ## add a new, unique ID to a list of component labels, and return the list
    updateIDs <- function(comps) {
        return(append(comps, min(setdiff(1:(max(comps)+1), comps))))
    }


    ## takes a single entry from lh (aka lh[h]) and removes the entry corresponding
    ## to position ii of lexical item kk
    holdoutLh <- function(lhh, kk, ii) {
        lhh[(lhh[,1]!=kk | lhh[,2]!=ii), ]
    }
