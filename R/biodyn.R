require(plyr)
require(mpb)

gtBd<-function(x) 1/params(x)["r"]

rebuildBiodyn<-function(r,p,k=1e3,b0=1,nyrs=50,niters=101){
  
  object=biodyn(params=FLPar(r=r,p=p,k=k,B0=b0))
  
  shape=c(refpts(object)["bmsy"]%/%params(object)["k"])
  
  object=window(object,end=nyrs)
  object@stock[]=refpts(object)["bmsy"]
  object@catch[]=0

  target=c(refpts(object)["bmsy",1])

  object=propagate(object,niters)   

  object@stock=object@stock%*%FLQuant(rep(seq(0,1,length.out=niters),each=nyrs),dimnames=dimnames(stock(object)))
  object=fwd(object,harvest=stock(object)[,-1]%=%0)

  dat=cbind(target=target,as.data.frame(stock(object),drop=TRUE))
  dat=ddply(dat,.(iter), with, {
  
    rtn=try(data.frame(year=year[(data-target)^2==min((data-target)^2)][1]))
    
    if ("try-error"%in%is(rtn)) return(NULL)
    
    return(rtn)
    })
  
  dat=cbind(dat,initial=c(stock(object)[,1,,,,dat$iter])/target)
  
  dat=dat[order(dat$initial),c("year","initial")][-1,]
  cbind(shape=shape,dat[dat$year<nyrs,])}

# dat=mdply(expand.grid(r=seq(0.1,   1,length.out=10),
#                       p=seq(2.0001,5,length.out=10)),rebuild)
# 
# ggplot(dat)+
#   geom_line(aes(initial,year,col=ac(r),group=paste(p,r)))+
#   #geom_hline(aes(yintercept=mean(gt(object))),col="red")+
#   theme_minimal()+
#   theme(legend.position="bottom")+
#   labs(y="Year", x=expression(Biomass/B[MSY]))



