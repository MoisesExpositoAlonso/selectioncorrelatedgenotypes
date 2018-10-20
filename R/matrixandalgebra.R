

Xc<-function( X){
  Xnew<-apply(X,2, function(i){ (i-mean(i))/sd(i) } )
  return(Xnew);
}

Kcalc<-function(X){
  K=(1/nrow(X))*(X_ %*% t(X_))
  return(K);
}

Vcalc<-function(X){
  V=t(X) %*% X
  return(V)
}

gwamarg<-function(y,X){
  b<-sapply(1:ncol(X), function(i){
    r<-lm(y~X[,i])
    return(coefficients(r)[2])
  })
}
  

gwapca<-function(y,X,PCs){
  bpca<-sapply(1:ncol(X), function(i){
    r<-lm(y~X[,i]+PCs)
    return(coefficients(r)[2])
  })
}
  