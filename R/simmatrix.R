
############ funcitons
mafsim<-function(p, type='uniform',rate=1, cutoff=1){
  stopifnot(type %in% c('uniform','exponential'))

  if(type=='uniform'){
    runif(p,0,0.5)
  }else if(type=='exponential'){
    es<-rexp(p,rate = 1)
    # es/max(es)
    0.5* es/ (cutoff*max(es))
  }
}

Xsim<-function(n=100,p=1000, maf){
  # X<-cbind(sapply(1:p,function(i) rbinom(n = n,size=1,prob = maf[i])))
#  X<-cbind(sapply(1:p,function(i) sample(c(-1,+1),size=n,replace=T,prob = c(maf[i],1-maf[i] ) )))
  X<-cbind(sapply(1:p,function(i) sample(c(-1,+1),size=n,replace=T,prob = c(1-maf[i],maf[i] ) )))
return(X)
}