refptsEB<-function(x){
  fbar(x)=FLQuant(refpts(x)[,"harvest",drop=TRUE])
  fbar(x)=qmax(fbar(x),1e-12)
  x=brp(x)
  
  dmns=dimnames(refpts(x))
  dmns$quant=c(dimnames(refpts(x))$quant,"eb")
  
  rfpts=FLPar(NA,dimnames=dmns)
  rfpts[,-9]=computeRefpts(x)
  rfpts[, 9]=ebiomass(x)
  
  rfpts}




