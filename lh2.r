running.var <- function(X) {
    n = 0
    mold = 0
    sold = 0

    for (x in X) {
        n = n+1
        mnew = mold + (x-mold)/n
        snew = ((n-1)*sold + (x-mnew)*(x-mold))/n
        mold = mnew
        sold = snew
    }
    return (list(n=n, m=mnew, s=snew))
}

var.increment <- function(xnew, v) {
    mold = v$m
    sold = v$s
    n = v$n+1

    if (sold < 0) {stop(paste('bad sample variance:', sold))}
    else if (n < 1) {stop(paste('bad sample count:', n))}

    mnew = mold + (xnew-mold)/n
    snew = ((n-1)*sold + (xnew-mnew)*(xnew-mold))/n
    v$n = n
    v$m = mnew
    v$s = snew
    return(v)
}

var.decrement <- function(xold, v) {
    mnew = v$m
    snew = v$s
    n = v$n - 1

    if (n<1) {
        warning('Removed last item from running stats')
        n = 0
        mold = 0
        sold = 0
    } else {
        mold = (mnew*(n+1) - xold) / n
        sold = (snew*(n+1) - (xold-mnew)*(xold-mold)) / n
    }
    v$n = n
    v$s = sold
    v$m = mold
    return(v)
}

# comput summary stats from pooling observations corresponding to x and y.  other fields
# of x are preserved in output z (y's are ignored).
var.increment.pooled <- function(x, y) {
    mx = x$m
    nx = x$n
    sx = x$s

    my = y$m
    ny = y$n
    sy = y$s

    z = x
    z$n = nx + ny
    z$m = (nx*mx + ny*my)/nz
    z$s = (nx*sx + ny*sy)/nz + nx*ny/(nx+ny)^2 * (mx-my)^2
    return(z)
}

# compute summary stats (n, mean, var) for group z with y removed.  other fields of z are
# preserved.
var.decrement.pooled <- function(z, y) {
    nz = z$n
    mz = z$m
    sz = z$s

    ny = y$n
    my = y$m
    sy = y$s

    x = z
    x$n = nz - ny
    x$m = (nz*mz - ny*my)/nx
    x$s = sz + ny/nx*(sz-sy) - ny/nz*(mx-my)^2
    return(x)
}

# compute the likelihood of the group of observations corresponding to obs belonging
# to the component described by phon.  both obs and phon are vectors of summary statistics
# e.g. names(.) = 'n', 'm', 's' for count, mean, variance (sample variance -- divided by n)
lhood <- function(obs, phon, params, log=TRUE) {
    n = obs$n
    nuN = with(params, nu + phon$n)
    muN = with(params, mu * nu/nuN + phon$m * phon$n/nuN)
    sN = with(params, s + phon$n*phon$s + (phon$m - mu) * nu*phon$n/(nu+phon$n))

    # normalization constants
    L = lgamma((nuN + n)/2) + 0.5*log(nuN/(n+nuN)) - lgamma(nuN/2) - n/2*log(pi*sN)
    # observations difference from phon
    L = L - (nuN + n)/2 * log(1 + (n * obs$s + (obs$m-muN)^2 * n*nuN/(n+nuN)) / sN)

    if (!log) L = exp(L)
    return(L)
}

lhood1 <- function(x, phon, params, log=TRUE) {lhood(list(n=1,m=x,s=0), phon, params, log)}


word.lhood <- function(word, lexeme, phons, params, log=TRUE) {
    L = 0
    for (i in 1:length(word)) {
        phon = phon.get(phons, lexeme[[i]]$p)
        L = L + lhood1(word[[i]], phon, params)
    }
}

# add an observation or group of observations to a phonetic category
phon.push <- function(phon, obs) {
    for (o in obs) {
        phon = var.increment(o, phon)
    }
    return(phon)
}

# remove an observation of group of observations from a phonetic cateogry
# (and update summary statistics)
phon.pull <- function(phon, obs) {
    for (o in obs) {
        phon = var.decrement(o, phon)
    }
    return(phon)
}

# add an observation or group of observations to a lexical category
lex.push <- function(lex, obs) {}
lex.pull <- function(lex, obs) {}

phon.get <- function(ps, name) {
    n = which(ps$names == name)
    if (length(n) < 1) {stop(paste('No phon with name', name))}
    return(ps$vals[[n]])
}

phons.make <- function(obs, labs) {
    names = unique(labs)
    vals = lapply(names, function(l) {running_var(obs[labs==l])})
    return(list(names=names, vals=vals))
}




lex.hdp <- function(words, params=list(nu=1.001, mu=0, s=1, alpha=1, r=5), n.iter=5) {
    for (i in 1:n.iter) {
        # first sweep: update word labels
        for (word in words) {
            lxcn = lex.holdout(lxcn, word)
            lhoods = sapply(
                c(lxcn, lex.prior.samp(word, phons, params)),
                function(l) {word.lhood(word, l, phons, params)})
            

        # end one interation
    }
}



######################################################################################################
### Some code for testing...
######################################################################################################

makeSomeWords <- function() {
    l = list(c(1,2), c(2,1), 1)
    means = c(-2,2)
    s = 1
    Nl = c(100,100,100)
    words = NULL
    for (i in 1:length(l)) {
        thisWord = NULL
        for (j in 1:Nl[i]) {
            w = NULL
            for (pos in l[[i]]) {
                p = l[[i]][pos]
                o = rnorm(n=1, mean=means[p])
                w = c(w, o)
            }
            thisWord = c(thisWord, list(w))
        }
        words = c(words, thisWord)
    }

    
    
    return(words)
}
               
