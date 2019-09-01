set.seed(181225)
#set.seed(800607)

#library('ggplot2')
#library('RColorBrewer')
#library('cowplot')
#library('png')
#library('plot3D')
library('doParallel')
library('foreach')
library('LaplacesDemon')
library('mvtnorm')
#library('MASS')
#library('dplyr')


## update hyperparameter for given data

updatepar <- function(given,parms){
    ngiven <- ncol(given)
    if(ngiven>0){
        mugiven <- rowMeans(given)
        ncovgiven <- crossprod(t(given-mugiven))
        k <- parms$k + ngiven
        return(list(k=k,
        mu=(parms$k*parms$mu + ngiven*mugiven)/k,
        nu=parms$nu + ngiven,
        s=parms$s + ncovgiven + parms$k * ngiven * crossprod(t(mugiven-parms$mu))/k
        ))
    }else{return(parms)}
}

logpredpost <- function(want,given,parms){
    d <- nrow(want)
    nwant <- ncol(want)

    parms <- updatepar(given,parms)
    scol <- diag(rep(1,nwant)) + 1/parms$k

    dmvt(want,delta=)
    
}
