mode<-function(x) {
  x=x[!is.na(x)]
  
  den=density(x)
  den$x[which.max(den$y)]}
