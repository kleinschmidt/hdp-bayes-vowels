### feldman2009.r --- implementations of the various models from Feldman,
### Griffiths, and Morgan (2009) paper on HDP learning of lexical and
### phonetic categories


################################################################################
### Toy 1d model ###############################################################
################################################################################

mu0 <- array(data=c(-5, -1, 1, 5), dim=c(4,1))
n0 <- c(400, 200, 200, 400)
s0 <- array(data=c(1,1,1,1), dim=c(1,1,4))

