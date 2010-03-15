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

test_sparseness <- function() {
    minp = 0.05
    startp = 0.05
    maxp = 0.45
    pstep = 0.10

    ms_probSparse = vector("list", (maxp-minp)/pstep + 1)

    for (p in seq(startp,maxp,pstep)) {
        seed = .Random.seed
        cat("----------------------- 1-D model with probabilistic sparseness of", p, " percent---------------------\n")
        Nl = round(c(200*(1-p), 200*p, 200*p, 200*(1-p), 100, 100))
        w = makeW(lminPairs, Nl)
        output = lexhdp(w, nIter=30)
        ms_probSparse[[round((p-minp)/pstep + 1)]] <- list(w=w, Nl=Nl, p=p, seed=seed, output=output)
        
        save(list="ms_probSparse", file='1dsparseness.RData')
    }
}

################################################################################
## ANALYSIS   ##################################################################
################################################################################


# extract a list of segment labels based on z and l for precision/recall analysis
getposlabs <- function(z, l, lexlab) {
    unlist(lapply(z, function(zz) l[[which(lexlab==zz)]]))
}



this = ms_probSparse[[3]]

getHitsFAs <- function(df, ns=length(df$output), masks=NA) {
    cat('analyzing p = ', df$p, ', n = ', sep='')

    # get actual segment labels
    baseNl = df$Nl
    realposlabs = getposlabs( unlist(lapply(1:length(baseNl), function(N) rep(N, baseNl[N]))),
        lminPairs,
        1:length(baseNl) )
    len = length(realposlabs)

    # compute pairwise same-different
    same = unlist(lapply(1:(len-1), function(N) {realposlabs[N]==realposlabs[(N+1):len]}))
    if (!any(is.na(masks))) {
        mask = unlist(lapply(1:(len-1),
                      function(N) {any(realposlabs[N]==masks) |
                                   sapply(realposlabs[(N+1):len], function(l) {any(l==masks)})}))
    }

    if (any(ns<0)) {ns = 1:length(df$output)}
    out = NULL
    out.masked = NULL

    for (n in ns) {
        cat(n, ', ')
        obsposlabs = with(df$output[[n]], getposlabs(z, l, lexlab))
        resp = unlist(lapply(1:(len-1), function(N) {obsposlabs[N]==obsposlabs[(N+1):len]}))
        
        tab = xtabs(~same+resp, data.frame(same=same, resp=resp))
        #print(tab)
        out = rbind(out, data.frame(p=df$p, iter=n,
                                    hits=tab['TRUE', 'TRUE'],
                                    miss=tab['TRUE', 'FALSE'],
                                    FAs =tab['FALSE', 'TRUE'],
                                    rej =tab['FALSE', 'FALSE']))
        if (!any(is.na(masks))) {
            tab.masked = xtabs(~same+resp, data.frame(same=same[mask], resp=resp[mask]))
            out.masked = rbind(out.masked, data.frame(p=df$p, iter=n,
                hits=tab.masked['TRUE', 'TRUE'],
                miss=tab.masked['TRUE', 'FALSE'],
                FAs =tab.masked['FALSE', 'TRUE'],
                rej =tab.masked['FALSE', 'FALSE']))
        }
    }
    cat('\n')

    confusions = xtabs(~real+obs, data.frame(real=realposlabs, obs=obsposlabs))
    
    return(list(full=out, masked=out.masked, confusions=confusions))
}

analyze_probSparse <- function(dfs) {
    out = list(full=NULL, masked=NULL)
    for (df in ms_probSparse) {
        res = getHitsFAs(df, 1:length(df$output), c(2,3))
        out$full = rbind(out$full, res$full)
        out$masked = rbind(out$masked, res$masked)
    }
    return(out)
}

crunchHitsFAs <- function(res) {
    res$acc = res$hits / (res$hits + res$FAs)
    res$compl = res$hits / (res$hits + res$miss)
    res$p = factor(res$p)

    res.m = melt(res[,c('p', 'iter', 'acc', 'compl')], id=c('p', 'iter'))
    p = ggplot(data=res.m, aes(x=iter, y=value))
    p + geom_line() + facet_grid(variable ~ p)
    return(res.m)
}

confusion_tables <- function(model) {
    baseNl = model$Nl
    realposlabs = getposlabs( unlist(lapply(1:length(baseNl), function(N) rep(N, baseNl[N]))),
        lminPairs,
        1:length(Nl) )
    obsposlabs = with(model$output[[length(model$output)]], getposlabs(z, l, lexlab))
    return(xtabs(~real+obs, data.frame(real=realposlabs, obs=obsposlabs)))
}

#full_analysis <- function() {
#    analyze_probSparse(ms_probSparse)
#    
