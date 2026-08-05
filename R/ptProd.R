# Pella–Tomlinson surplus production at biomass B
ptFn<-function(B, param) 
  
  (param["r"] %/% param["p"]) %*% B * (1 - exp(log(B %/% param["k"])%*%param["p"]))

# Yield at biomass B under equilibrium (Y(B) = surplus production)
ptY<-function(x) 
  ptFn(stock(x), params(x))

ptSp<-function(x){
  res=model.frame(FLQuants(x,yield=FLCore::catch,biomass=FLCore::stock,pf=ptY),drop=TRUE)
  res=data.frame(biomass=head(res$biomass,-1),
                 sp     =head(res$biomass,-1)-tail(res$biomass,-1)+head(res$yield,-1))
  res}