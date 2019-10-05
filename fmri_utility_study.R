## Author: PGL  Porta Mana
## Created: 2017-11-29T12:03:29+0100
## Last-Updated: 2019-09-05T15:39:54+0200
##
library('fBasics')
library('ggplot2')
library('RColorBrewer')
library('cowplot')
library('doParallel')
library('foreach')
library('LaplacesDemon')
library('mvtnorm')

## Colour-blind friendly palette, from http://www.sron.nl/~pault/
mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
mycolours <- c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black')
palette(mycolours)
dev.off()

dd <- 40 ## number of graph quantities used in the paper

hdata <- t(unname(read.csv('weights_healthy_40cons.csv',header=FALSE,sep=',')))
sdata <- t(unname(read.csv('weights_schizo_40cons.csv',header=FALSE,sep=',')))

## data transformations and their initial parameters
## modified logit + Jacobian
logit2 <- function(x){log(1+x)-log(1-x)}
ilogit2 <- function(x){ex <- exp(x)
    (ex-1)/(ex+1)}
jlogit <- function(x){2/((1+x)*(1-x))}
parlogit2 <- list(k=1, mu=rep(0,dd), nu=dd+1, cov=diag(rep(dd+2,dd)))


## modified tangent + Jacobian
tg2 <- function(x){tan(x*pi/2)}
itg2 <- function(x){2*atan(x)/pi}
jtg <- function(x){pi/(1+cos(pi*x))}
partg2 <- list(k=1, mu=rep(0,dd), nu=dd+1, cov=diag(rep((dd+2)/4,dd)))

## normal + Jacobian
id2 <- function(x){x}
iid2 <- function(x){x}
jid <- function(x){1}
parid2 <- list(k=1, mu=rep(0,dd), nu=dd+1, cov=10*diag(rep(1,dd)))

## Jacobian of ln
jlog <- function(x){1/x}

## update parameters given new data
## parameters: k, mu, nu, sigma
updateparameters <- function(data,parameters){
    n <- nrow(data)
    k0 <- parameters$k
    mu0 <- parameters$mu
    
    k <- k0 + n
    mud <- colMeans(data)
    covd <- cov(data)*(n-1)

    list(k=k,
         mu=(k0*mu0+n*mud)/k,
         nu=parameters$nu + n,
         cov=parameters$cov + covd + k0*n*crossprod(t(mud-mu0))/k)
}


loglikelihood <- function(x,parameters){
    k <- parameters$k
    dof <- parameters$nu - dd + 1
    mvtnorm::dmvt(x, delta=parameters$mu, sigma=parameters$cov*(k+1)/k/dof,
         df=dof, log=T, type='shifted')
}

## (u(h|h) - u(s|h), u(s|s) - u(h|s)
logutility <- log(c(1,1))


## train modified logit
hdataparlogit2 <- updateparameters(logit2(hdata),parlogit2)
sdataparlogit2 <- updateparameters(logit2(sdata),parlogit2)


## train modified tangent
hdatapartg2 <- updateparameters(tg2(hdata),partg2)
sdatapartg2 <- updateparameters(tg2(sdata),partg2)

## train normal
hdataparid2 <- updateparameters(hdata,parid2)
sdataparid2 <- updateparameters(sdata,parid2)


generatef <- function(parameters){
    k <- parameters$k
    dof <- parameters$nu - dd + 1
    mvtnorm::rmvt(1, delta=parameters$mu,sigma=parameters$cov*(k+1)/k/dof,
         df=dof,type='shifted')
}

generatesamples <- function(seqname,
                            transf,gtransf,
                            hgenparams,sgenparams,
                            hparams,sparams,
                            logutility=log(c(1,1)),
                            nsamples,nsubsamples,
                            cores=1,
                            rseed=999
                            ){
    set.seed(rseed)
    if(cores>1){
        cl <- makeCluster(cores, outfile="")
        registerDoParallel(cl)}
    
    hsamples <- foreach(i=1:nsamples, .combine=c, .export=ls(globalenv()), .packages=)%dopar%{
        sample <- transf(gtransf(generatef(hgenparams)))
        logutility[1]+loglikelihood(sample,hparams) > logutility[2]+loglikelihood(sample,sparams)
    }
    ssamples <- foreach(i=1:nsamples, .combine=c, .export=ls(globalenv()))%dopar%{
        sample <- transf(gtransf(generatef(sgenparams)))
        logutility[1]+loglikelihood(sample,hparams) < logutility[2]+loglikelihood(sample,sparams)
    }

    sequtil <- foreach(i=round(seq(1,nsamples,length.out=nsubsamples)), .combine=c, .export=ls(globalenv()))%dopar%{mean(hsamples[1:i]/2+ssamples[1:i]/2)}
    write.table(sequtil,file=paste0(seqname,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')

    if(cores>1){stopCluster(cl)}
}

nsamples <- 1e5
nsubsamples <- 100
cores <- 2
## normal
generatesamples(seqname='sequtil_id',
                transf=id2,gtransf=iid2,
                hgenparams=hdataparid2,sgenparams=sdataparid2,
                hparams=hdataparid2,sparams=sdataparid2,
                logutility=log(c(1,1)),
                nsamples=nsamples,nsubsamples=nsubsamples,
                cores=cores)

## tangent
generatesamples(seqname='sequtil_tg',
                transf=tg2,gtransf=itg2,
                hgenparams=hdatapartg2,sgenparams=sdatapartg2,
                hparams=hdatapartg2,sparams=sdatapartg2,
                logutility=log(c(1,1)),
                nsamples=nsamples,nsubsamples=nsubsamples,
                cores=cores)

## logit
generatesamples(seqname='sequtil_logit',
                transf=logit2,gtransf=ilogit2,
                hgenparams=hdataparlogit2,sgenparams=sdataparlogit2,
                hparams=hdataparlogit2,sparams=sdataparlogit2,
                logutility=log(c(1,1)),
                nsamples=nsamples,nsubsamples=nsubsamples,
                cores=cores)
