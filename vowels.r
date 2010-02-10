library(mvtnorm)

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

## generate observations from a lexical item
## returns a sub-array of dim()=c(d,i,n) (one sheet per repetition) with columns
## padded with NAs for short words
##
## mu should be a be a d-by-nphon matrix, each column is a category mean
## s is a d-by-d-by-nphon 3-d array with sheets as covariance matrices
lexgen <- function(n, lexeme, phonlab, mu, s, width=0) {
  d=nrow(mu)
  # make a (padded) dim()=c(n,d,i) array first,
  num = max(width,length(lexeme))
  output = array(NA, c(n,d,num))
  for (pos in 1:length(lexeme)) {
    phon = which(phonlab==lexeme[pos])
    output[,,pos] = rmvnorm(n=n, mean=mu[,phon], sigma=s[,,phon])
  }
  # and permute
  return(aperm(output, c(2,3,1)))
}

makeEasyW <- function(l = list(c(1), c(2), c(1,2), c(2,1)) , Nl = c(100, 100, 100, 100)) {
  phonlab = c(1,2)
  numphon = length(phonlab)

  lexlen = sapply(l, length)
  
  covs = array(rep(diag(2), 2), c(2,2,2))
  means = cbind(c(2,2), c(-2,-2))

  w = NULL
  for (n in 1:length(l)) {
    w = c(w, lexgen(Nl[n], l[[n]], phonlab, means, covs, max(lexlen)))
  }

  return(array(w, c(nrow(means), max(lexlen), sum(Nl))))
}

  

