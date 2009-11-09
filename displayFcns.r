### displayFcns.r --- functions to visualize the state of a gaussian mixture model
### with arbitrary components.  Works for 1d and 2d (or higher, but only displays
### the first two dimensions)

## plot side-by-side histograms (good for 1-d case)
compHist <- function(pts, labels, comps, counts, means, precisions) {
  h <- NULL
  h0 <- hist(pts, plot=FALSE)
  for (co in comps) {
    hn <- hist(pts[labels==co], breaks=h0$breaks, plot=FALSE)
    h <- rbind(h, hn$counts)
  }
  colnames(h) <- h0$mids
  par(mfcol=c(2,1))
  colors <- (comps-1)%%length(palette()) + 1
  barplot(h, beside=TRUE, col=colors)
  ## plot the component pdfs, too
  counts <- counts/sum(counts)
  x <- seq(min(h0$breaks), max(h0$breaks), .1)
  plot(c(min(x), max(x)), c(0,max(counts*dnorm(means, mean=means, sd=sqrt(1/precisions)))),
       type="n", ylab="", xlab="") 
  for (compn in 1:length(comps)) {
    lines(x,
          counts[compn]*dnorm(x, mean=means[compn], sd=sqrt(1/precisions[compn])),
          col=colors[compn])
  }
}

## 2-d scatterplot, labled by component assignment
palette(c("red", "green3", "blue", "cyan", "magenta", "yellow"))
plotChars <- c(1, 2, 5, 6)
comp2dplot <- function(pts, labels, comps, means, precisions) {
  plot(pts, type="n")
  for (co in comps) {
    ppts <- pts[labels==co, ]
    if (!is.matrix(ppts)) {ppts <- matrix(ppts, length(ppts)/2, 2)}
    color <- (co-1)%%length(palette()) + 1
    points(ppts[,1], ppts[,2],
           pch=plotChars[ceiling(co/length(palette()))],
           col=color)
    lines(ellipse(x=solve(precisions[1:2,1:2,which(comps==co)]),
                  centre=means[which(comps==co),1:2]),
          col=color)
  }
}
