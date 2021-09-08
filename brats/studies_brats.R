## Author: Battistin, Gonzalo Cogno, Porta Mana
## Last-Updated: 2021-08-21T20:35:56+0200
################
## Script for:
## - outputting samples of prior & posterior distributions
## - calculating posteriors
## Uses Dirichlet prior
################
## library('parallel')
## mycluster <- makeCluster(3)
#### Custom setup ####
## Colour-blind friendly palettes, from https://personal.sron.nl/~pault/
## (consider using khroma package instead)
## library('RColorBrewer')
## mypurpleblue <- '#4477AA'
## myblue <- '#66CCEE'
## mygreen <- '#228833'
## myyellow <- '#CCBB44'
## myred <- '#EE6677'
## myredpurple <- '#AA3377'
## mygrey <- '#BBBBBB'
## mypalette <- c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black')
## palette(mypalette)
## barpalette <- colorRampPalette(c(mypurpleblue,'white',myredpurple),space='Lab')
## barpalettepos <- colorRampPalette(c('white','black'),space='Lab')
library('data.table')
library('khroma')
palette(colour('bright')())
## palette(colour('muted')())
library('ggplot2')
library('ggthemes')
theme_set(theme_bw(base_size=18))
scale_colour_discrete <- scale_colour_bright
#library('cowplot')
library('png')
library('foreach')
##library('ash')
#library('extraDistr')
library('PReMiuM')
library('mvtnorm')
options(bitmapType='cairo')
pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5)} # to output in png format
##
library('foreach')
library('doFuture')
registerDoFuture()
library('doRNG')

##scorefiles <- c('dynUnetValScoresV0.csv', 'multi_encodnnUValScoresLossV0.csv')
scorefiles <- c('Preds_DynUMultiEncodAtten_Brats21DiceCELogRanger21_FullSet_2ModelsWithNLoss_Ep20_20_Logit_Ensemble_modifiedNoTCProc.csv', 'seg_val_m1_m2_withoutHistogramStandardized.csv')
##
models <- scorefiles
modelsdata <- data.table()
for(i in 1:length(scorefiles)){
                 modelsdata <- rbind(modelsdata,
                                    cbind(
                                        as.data.table(read.csv(scorefiles[i],header=T,sep=',')), data.table(model=models[i])))
                                    }
covariates <- setdiff(colnames(modelsdata), c('scan_id', 'model'))

##
pdff('histograms_models')
for(cov in covariates){
    print(
        ggplot(modelsdata) +
        geom_histogram(aes(x=get(cov), fill=model), color='black', bins=10, alpha=0.5, position='identity') + scale_fill_bright() + theme(legend.position='top') +
        xlab(cov) +
        geom_vline(data=modelsdata[, .(mean=mean(get(cov))), by=model], aes(xintercept=mean, color=model, linetype=model), alpha=0.5, size=2)  + theme(legend.position='top')
        )
}
dev.off()
##
covariates <- setdiff(colnames(modelsdata), c('scan_id', 'model'))
##
pdff('pairplots_models')
for(icov1 in 1:(length(covariates)-1)){
    for(icov2 in (icov1+1):length(covariates)){
        cov1 <- covariates[icov1]
        cov2 <- covariates[icov2]
    print(
        ggplot(modelsdata) +
        geom_point(aes(x=get(cov1), y=get(cov2), color=model, shape=model), size=2, alpha=0.5) +
        xlab(cov1) + ylab(cov2)  + theme(legend.position='top')
        )
}}
dev.off()







covariates <- setdiff(colnames(modelsdata), c('scan_id', 'model'))
datas <- foreach(modname=unique(modelsdata$model))%do%{
    as.matrix(modelsdata[model==modname, -c('scan_id','model')])
}
names(datas) <- unique(modelsdata$model)

covars <- foreach(i=1:length(datas))%do%{
    tcrossprod(t(datas[[i]])-rowMeans(datas[[i]]))/nrow(datas[[i]])
}
names(covars) <- names(datas)

corrs <- lapply(covars,cov2cor)

#### scale bounded covariates
modelsdatascaled <- modelsdata
scalelogit <- function(x){ qlogis(x)/2 }#-log(1/x-1)/2 } ## faster than qlogis(x, scale=1/2)
##
scalelog <- function(x){ log(x) }
##
## transform 0-Inf features to log-scale
indx <- foreach(nn=c('Haus'), .combine='|')%do%{grepl(nn, colnames(modelsdata))}
for(elem in colnames(modelsdata)[indx]){
    datum <- scalelog(modelsdata[[elem]])
    eps <- max(diff(sort(unique(datum[abs(datum)!=Inf]))))
    datum[datum==-Inf] <- min(datum[abs(datum)!=Inf])-eps
    modelsdatascaled[, elem] <- datum
}
names(modelsdatascaled)[indx] <- paste0('log-',names(modelsdata)[indx])
## transform 0-1 features to logit scale
indx <- foreach(nn=c('Dice','Sen','Spe'), .combine='|')%do%{grepl(nn, colnames(modelsdata))}
for(elem in colnames(modelsdata)[indx]){
    datum <- scalelogit(modelsdata[[elem]])
    eps <- max(diff(sort(unique(datum[abs(datum)!=Inf]))))
    datum[datum==Inf] <- max(datum[abs(datum)!=Inf])+eps
    datum[datum==-Inf] <- min(datum[abs(datum)!=Inf])-eps
    modelsdatascaled[, elem] <- datum
}
names(modelsdatascaled)[indx] <- paste0('logit-',names(modelsdata)[indx])
##
covariatesscaled <- setdiff(colnames(modelsdatascaled), c('scan_id', 'model'))

##
pdff('pairplots_models_scaled')
for(icov1 in 1:(length(covariatesscaled)-1)){
    for(icov2 in (icov1+1):length(covariatesscaled)){
        cov1 <- covariatesscaled[icov1]
        cov2 <- covariatesscaled[icov2]
    print(
        ggplot(modelsdatascaled) +
        geom_point(aes(x=get(cov1), y=get(cov2), color=model, shape=model), size=2, alpha=0.5) +
        xlab(cov1) + ylab(cov2) + theme(legend.position='top')
        )
}}
dev.off()
##
pdff('histograms_models_scaled')
for(cov in covariatesscaled){
    print(
        ggplot(modelsdatascaled) +
        geom_histogram(aes(x=get(cov), fill=model), color='black', bins=10, alpha=0.5, position='identity') + scale_fill_bright() + theme(legend.position='top') +
        xlab(cov) + theme(legend.position='top') +
        geom_vline(data=modelsdatascaled[, .(mean=mean(get(cov))), by=model], aes(xintercept=mean, color=model, linetype=model), alpha=0.5, size=2) 
        )
}
dev.off()

######################################################################
#### prediction via MCMC
mdata <- modelsdatascaled
covNames <- setdiff(colnames(mdata), c('scan_id', 'model'))
discreteCovs <- covNames[sapply(covNames, function(x){is.integer(mdata[[x]])})]
continuousCovs <- covNames[sapply(covNames, function(x){is.double(mdata[[x]])})]
##
dimsC <- length(continuousCovs)
ilogit <- grepl('logit-', continuousCovs)
ilog <- grepl('log-', continuousCovs)
## Hyperparameters
mulog <- 1
mulogit <- 0
mu0 <- mulog * ilog + mulogit * ilogit
names(mu0) <- continuousCovs
varmulog <- 2^2
varmulogit <- 2^2
sigmu0 <- varmulog * ilog + varmulogit * ilogit
names(sigmu0) <- continuousCovs
expvarlog <- (1/2)^2
expvarlogit <- (1/2)^2
expvar0 <- expvarlog * ilog + expvarlogit * ilogit
names(expvar0) <- continuousCovs
## diagvar0 <- (2)^2
## df0 <- (2*expvarlog^2)/diagvar0 + 4
df0 <- 30
hmu0 <- mu0 # sapply(continuousCovs, function(x){mu0[sapply(names(mu0),function(y){grepl(y,x)})]})
hnu0 <- df0 + dimsC - 1
hkappa0 <- expvarlog/varmulog
## hkappa0 <- 0.01
hDelta0 <- solve(diag(expvar0))/(df0-2) # solve(diag(sapply(continuousCovs, function(x){expvar0[sapply(names(expvar0),function(y){grepl(y,x)})]})))/(df0-2)
colnames(hDelta0) <- rownames(hDelta0) <- names(hmu0)
##
testhp <- setHyperparams(mu0=mu0, kappa0=hnu0, R0=hDelta0, nu0=hkappa0)
##
models <- unique(mdata$model)
##
rm(mcmcrun)
gc()
starttime <- Sys.time()
plan(sequential)
plan(multisession, workers = 3L)
mcmcrun <- foreach(case=models, .inorder=FALSE, .packages=c('data.table','PReMiuM','mvtnorm'))%dopar%{
    outfile <- paste0('_mcoutput_',case)
    datamcr <- mdata[model==case, covNames, with=F]
    ##
    ## str(datamcr)
        c(model=case, profRegr(excludeY=TRUE, xModel='Normal', nSweeps=2000e2, nBurn=3000e2, nFilter=2e2, data=as.data.frame(datamcr), nClusInit=80, covNames=c(discreteCovs,continuousCovs), continuousCovs=continuousCovs, nProgress=100e2, seed=147, output=outfile, useHyperpriorR1=FALSE, useNormInvWishPrior=TRUE, hyper=testhp, alpha=4))
}
plan(sequential)
names(mcmcrun) <- models
elapsedtime <- Sys.time() - starttime
elapsedtime
## 500: 1.22 min
## 5000: 9.14 min
## 50000: 53.22 min
## 500000: 355.8 min
## 5000000: 2.116277 days
## Save MCMC samples
MCMCdata <- as.list(rep(NA,length(mcmcrun)))
names(MCMCdata) <- names(mcmcrun)
##
for(case in models){
    outfile <- paste0('_mcoutput_',case)
    testmc <- mcmcrun[[case]]
    ## log-likelihood and log-posteriors
    fc <- file(paste0(outfile,"_logPost.txt"))
    logPost <- sapply(strsplit(readLines(fc), " +"), as.numeric)
    rownames(logPost) <- c('log-post','log-likelihood','log-prior')
    close(fc)
    ## Samples of numbers of clusters
    fc <- file(paste0(outfile,'_nClusters.txt'))
    nList <- sapply(strsplit(readLines(fc), " +"), as.integer)
    close(fc)
    ## Samples of Dirichlet-process alpha
    fc <- file(paste0(outfile,'_alpha.txt'))
    alphaList <- sapply(strsplit(readLines(fc), " +"), as.numeric)
    close(fc)
    ## Samples of cluster weights
    fc <- file(paste0(outfile,'_psi.txt'))
    psiList <- lapply(strsplit(readLines(fc), " +"), function(x){x <- as.numeric(x); x[x>=0]})
    close(fc)
    ## ##  Samples of Dirichlet-distribution phis (discrete covariates)
    ## fc <- file(paste0(outfile,'_phi.txt'))
    ## phiList <- lapply(strsplit(readLines(fc), " +"), function(x){x <- as.numeric(x); x[x>=0]})
    ## close(fc)
    ##  Samples of normal-distribution means (continuous covariates)
    fc <- file(paste0(outfile,'_mu.txt'))
    muList <- lapply(strsplit(readLines(fc), " +"), as.numeric)
    close(fc)
    ##  Samples of normal-distribution covariances (continuous covariates)
    fc <- file(paste0(outfile,'_Sigma.txt'))
    sigmaList <- lapply(strsplit(readLines(fc), " +"), as.numeric)
    close(fc)
    ##
    nContCov <- length(testmc$covNames)
##    nDiscrCov <- testmc$nDiscreteCovs
    ## nCat <- testmc$nCategories
    ## cumcats <- c(0,cumsum(nCat))
    for(i in 1:length(nList)){
        ## catgroups <- cumcats*nList[i]
        ## datum <- phiList[[i]]
        ## phiList[[i]] <- lapply(1:length(testmc$nCategories), function(j){
        ##     y <- datum[(catgroups[j]+1):catgroups[j+1]]
        ##     dim(y) <- c(nList[i], testmc$nCategories[j])
        ##     aperm(y)
        ## })
        ## names(phiList[[i]]) <- testmc$discreteCovs
        dim(sigmaList[[i]]) <- c(nList[i], nContCov, nContCov)
        sigmaList[[i]] <- aperm(sigmaList[[i]])
        rownames(sigmaList[[i]]) <- testmc$covNames
        colnames(sigmaList[[i]]) <- testmc$covNames
        dim(muList[[i]]) <- c(nList[i], nContCov)
        muList[[i]] <- aperm(muList[[i]])
        rownames(muList[[i]]) <- testmc$covNames
    }
    ## Save the samples above
    MCMCdata[[case]] <- list(model=case, nList=nList, alphaList=alphaList, psiList=psiList, muList=muList, sigmaList=sigmaList, logPost=logPost)
}
##
## save.image(file=paste0('_directmodel_N',ndata,'_',length(covNums),'covs.RData'))
##
## Diagnostic plots
pdff('mcsummary')
for(j in 1:length(MCMCdata)){
    sd <- MCMCdata[[j]]
    matplot(sd$logPost[2,], type='l',ylim=range(sd$logPost[2,],na.rm=T,finite=T),ylab='log-likelihood',col=palette()[j], main=paste0('freqs: ',models[j]))
    }
for(j in 1:length(MCMCdata)){
    sd <- MCMCdata[[j]]
matplot(sd$nList,type='l',ylim=range(sd$nList,na.rm=T,finite=T),ylab='no. clusters',col=palette()[j], main=paste0('freqs: ',models[j]))
}
for(j in 1:length(MCMCdata)){
    sd <- MCMCdata[[j]]
    matplot(sd$alphaList,type='l',ylim=range(sd$alphaList,na.rm=T,finite=T),ylab='alpha',col=palette()[j], main=paste0('freqs: ',models[j]))
    }
for(j in 1:length(MCMCdata)){
    sd <- MCMCdata[[j]]
    for(i in c(1,3)){matplot(sd$logPost[i,],type='l',ylim=range(sd$logPost[i,],na.rm=T,finite=T),col=palette()[j], main=paste0('freqs: ',models[j]))}
    }
dev.off()
##
save.image(file=paste0('_MCMCresults.RData'))


## function to normalize absolute frequencies
normalize <- function(freqs){freqs/sum(freqs)}
## for rows of frequency distributions
normalizem <- function(freqs){freqs/rowSums(freqs)}
##
meanssamples <- list()
for(case in models){
    print(case)
    plan(sequential)
    plan(multisession, workers = 6L)
    meanssamples[[case]] <- foreach(cov=covNames, .combine=cbind)%:%foreach(asample=1:length(MCMCdata[[case]]$nList), .combine=rbind, .inorder=FALSE)%dopar%{
        mulist <- MCMCdata[[case]]$muList[[asample]][cov,]
        sdlist <- sqrt(MCMCdata[[case]]$sigmaList[[asample]][cov,cov,])
        if(grepl('log-', cov)){
            newmulist <- exp(mulist + sdlist^2/2)
        } else if(grepl('logit-', cov)){
            newmulist <- sapply(1:length(mulist), function(alist){
                intf <- function(t){dnorm(-log(1/t-1)/2,mulist[alist],sdlist[alist])/(1-t)}
                intres <- integrate(intf, 0, 1, subdivisions=10000L, stop.on.error=F)
                if(intres$message=='OK'){return(intres$value/2)} else {return(NA)}
            })
        }
        W <- MCMCdata[[case]]$psiList[[asample]]
        sum(newmulist * W, na.rm=T)/sum(W[!is.na(newmulist)])
    }
    plan(sequential)
    colnames(meanssamples[[case]]) <- covariates
}

datameans <- data.table()
for(i in 1:length(meanssamples)){
    datameans <- rbind(datameans,
                       cbind(as.data.table(meanssamples[[i]]), data.table(model=names(meanssamples)[i]))
)
}

pdff('uncertainties_meanscores')
for(cov in covariates){
    print(
        ggplot(datameans) +
        geom_histogram(aes(x=get(cov), y=..density.., fill=model), color='black', bins=20, alpha=0.5, position='identity') + scale_fill_bright() + theme(legend.position='top') +
        xlab(paste0('predicted long-run mean score for ',cov)) + ylab('probability density') +
        geom_vline(data=datameans[, .(mean=mean(get(cov))), by=model], aes(xintercept=mean, color=model, linetype=model), alpha=0.75, size=2) 
        )
}
dev.off()


diffmeans <- meanssamples[[2]]-meanssamples[[1]]
probpositive <- apply(diffmeans, 2, function(x){sum(x>0)/length(x)})


save.image(file=paste0('_MCM_means_Cresults.RData'))



diffmeans <- meanssamples[[2]][sample(nrow(meanssamples[[1]])),]-meanssamples[[1]][sample(nrow(meanssamples[[1]])),]

probpositive2 <- apply(diffmeans, 2, function(x){sum(x>0)/length(x)})



library('foreach')
library('doFuture')
registerDoFuture()
library('rgl')
library('plot3D')


score <- function(p,v,t,q=1){
    ins <- seq_len(v+1) - 1
    outs <- seq_len(t-v+1) - 1
    pins <- dbinom(x=ins,size=v,prob=p)
    pouts <- dpois(x=outs,lambda=q)
    vals <- 2*ins/(v+outer(ins,outs,'+'))
    vals[is.na(vals)] <- 1
    sum(t(vals * pins) * pouts)
}

gc()
t <- 128^3
xg <- seq(0.7,1,length.out=10)
yg <- round(seq(0,3^3,length.out=10))
me <- mesh(xg,yg)

gc()
plan(sequential)
plan(multisession, workers = 6L)
zz <- foreach(i=seq_along(xg), .combine=rbind)%:%foreach(j=seq_along(yg), .combine=cbind)%do%{
    score(me$x[i,j],me$y[i,j],t) }
plan(sequential)
##

plot3d(me$x,me$y,zz, xlab='p', ylab='tumour volume/total', zlab='Dice score', type='s',size=0.5, zlim=c(0,1))

surface3d(me$x,me$y,zz, col=palette()[1])#, xlab='p', ylab='tumour volume', zlab='Dice score')


datafiles <- c('Brats21PlusBrainROI_4CVFrom1Split_MetaDF.csv',
               'ValOfTrainFiles_DynUnet_Brats21DiceCELogRanger21_CV4of1Split_fold3_0.9041_Fold3_epoch124.csv',
               'ValOfTrainFiles_DynUnet_Brats21DiceCELogRanger21_CV4of1Split_fold2_0.8814_Fold2_epoch60.csv',
               'ValOfTrainFiles_DynUnet_Brats21DiceCELogRanger21_CV4of1Split_0.8799_Fold0_epoch82.csv',
               'ValOfTrainFiles_DynUnet_Brats21DiceCELogRanger21_CV4of1Split_fold1_0.9195_Fold1_epoch66.csv')
##
voldata <- foreach(file=datafiles)%do%{as.data.table(read.csv(file,header=T,sep=','))}
names(voldata) <- c('vols',paste0('f',1:4))

vdata <- voldata$vols
volcols <- colnames(vdata)[sapply(colnames(vdata), function(x){grepl('Vol',x)})]
vdata <- vdata[ ,c('BraTS21ID',volcols,'cfold'), with=F]
##
tregions <- unique(tdata$Tumor.regions)

pdff('vol_vs_dice')
for(regi in tregions){
    for(vset in 1:4){
        sdata <- voldata[[vset+1]]
        posids <- sapply(sdata$BraTS21ID,function(x){which(vdata$BraTS21ID == x)})
        tdata <- cbind(vdata[posids], sdata[,-1])
        ##
        vcol <- paste0(regi,'Vol')
        print(
            ggplot(tdata[Tumor.regions==regi,c('Dice.score','Model',vcol),with=F]) +
            geom_point(aes(x=(Dice.score), y=(get(vcol)), color=Model, shape=Model), size=3, alpha=0.9) +
            xlab(paste0('dice, fold ',vset-1)) + ylab(vcol) + theme(legend.position='top')
        )
    }
}
dev.off()

pdff('log-vol_vs_dice')
for(regi in tregions){
    for(vset in 1:4){
        sdata <- voldata[[vset+1]]
        posids <- sapply(sdata$BraTS21ID,function(x){which(vdata$BraTS21ID == x)})
        tdata <- cbind(vdata[posids], sdata[,-1])
        ##
        vcol <- paste0(regi,'Vol')
        print(
            ggplot(tdata[Tumor.regions==regi,c('Dice.score','Model',vcol),with=F]) +
            geom_point(aes(x=qlogis(Dice.score), y=log(get(vcol)), color=Model, shape=Model), size=3, alpha=0.9) +
            xlab(paste0('logit-dice, fold ',vset-1)) + ylab(paste0('log-',vcol)) + theme(legend.position='top')
        )
    }
}
dev.off()


commonids <- vdata$BraTS21ID %in% sdata$BraTS21ID

rvdata <- vdata[BraTS21ID %in% sdata$BraTS21ID]

commonids <- intersect(sdata$BraTS21ID, vdata$BraTS21ID)

tdiff <- setdiff(sdata$BraTS21ID, commonids)

ids <- foreach(x=voldata)%do%{x$BraTS21ID}






datafile <- c('MergedScoreWithVolumesLessthan60Percent4ValOfTrain_DynUMEncodAtten_Brats21DiceCELogRng21_1Split_NLossV0_0.9031_Fold0_epoch80.csv')
##
datat <- as.data.table(read.csv(datafile,header=T,sep=','))

tregions <- unique(datat$Tumor.regions)
##
pdff('vol_vs_dice2')
for(regi in tregions){
        vcol <- paste0(regi,'Vol')
        print(
            ggplot(datat[Tumor.regions==regi,c('Dice.score','Model',vcol),with=F]) +
            geom_point(aes(x=qlogis(Dice.score), y=log(get(vcol)), color=Model, shape=Model), size=3, alpha=0.9) +
            xlab(paste0('dice score')) + ylab(vcol) + theme(legend.position='top')
        )
    }
dev.off()




#######################################################################




xg <- seq(0,1,0.01); yg <- 0.6*xg+0.2*(1-xg); yg2 <- 0.7*xg+0.3*(1-xg)
dt <- data.table(x=xg,y=c(yg,yg2),treat=rep(c('drug','placebo'),each=length(xg)))
png('drup_plac.png',height=11.7,width=16.5,units='cm',res=600)
print(
qplot(x=x,y=y,data=dt,color=treat,geom='line',lty=treat, lwd=I(1.2), xlab=TeX('$P(N_{fem}\\; |\\; A)$'), ylab=TeX('$P(N_{heal}\\;|\\; N_{treat}\\, &\\, D\\, & \\, A)$'),ylim=c(0,1)) + scale_color_manual(values=palette()[c(4,3)]) + theme(legend.pos='top')
)
dev.off()

matplot(xgrid <- seq(0,1,0.1), ygrid <- matrix(c(0.2*xgrid + 0.6*(1-xgrid), 0.3*xgrid + 0.7*(1-xgrid)),ncol=2), type='l', ylim=c(0,max(ygrid)), ylab=('p. healthy')); grid()


covariates <- setdiff(colnames(modelsdatascaled), c('scan_id', 'model'))
datas <- foreach(modname=unique(modelsdatascaled$model))%do%{
    as.matrix(modelsdatascaled[model==modname, -c('scan_id','model')])
}
names(datas) <- unique(modelsdatascaled$model)
##
covars <- foreach(i=1:length(datas))%do%{
    tcrossprod(t(datas[[i]])-rowMeans(datas[[i]]))/nrow(datas[[i]])
}
names(covars) <- names(datas)
##
corrs <- lapply(covars,cov2cor)



qt <- seq(0.1,99.9,by=0.1)/100
matplot(x=matrix(pcauchy(qnorm(qt,lower.tail=T),lower.tail=T),nrow=1), y=0, type='p', pch=16, cex=0.5, col='black')

hist(plogis(qnorm(qt,lower.tail=T),lower.tail=T),breaks=seq(-0.001,1.001,length.out=50))


matplot(xgrid <- seq(-1,1,length.out=1000)*5, (cbind(dnorm(xgrid),dcauchy(xgrid),dlogis(xgrid))), type='l')


####################################
#### Old checks
####################################


filenames <- c(
    'ValOfTrainFiles_DynUnet_Brats21DiceCELogRanger21_CV4of1Split_0.8799_Fold0_epoch82.csv',
    'ValOfTrainFiles_DynUnet_Brats21DiceCELogRanger21_CV4of1Split_fold1_0.9195_Fold1_epoch66.csv',
    'ValOfTrainFiles_DynUnet_Brats21DiceCELogRanger21_CV4of1Split_fold2_0.8814_Fold2_epoch60.csv',
    'ValOfTrainFiles_DynUnet_Brats21DiceCELogRanger21_CV4of1Split_fold3_0.9041_Fold3_epoch124.csv')

dataall <- foreach(i=1:4)%do%{as.data.table(read.csv(filenames[i],header=T,sep=','))}

ttypes <- c('TC','ET','WT')
    pdff(paste0('hist_fold'))
    for(ttype in ttypes){
for(fold in 1:4){
        hist(dataall[[fold]][Tumor.regions==ttype,Dice.score],
             breaks=seq(0,1,length.out=50),
             main=paste0('fold ',fold,'  ',ttype),xlab='dice score')
    }
}
    dev.off()

    pdff(paste0('points_fold'))
    for(ttype in ttypes){
for(fold in 1:4){
        matplot(x=matrix(0:1,nrow=2),y=matrix(rep(dataall[[fold]][Tumor.regions==ttype,Dice.score],2),nrow=2,byrow=T),
             type='l',lty=1,lwd=0.1,col='black',
             main=paste0('fold ',fold,'  ',ttype),ylab='dice score',xlab='')
    }
}
    dev.off()




















source('funcmodel3randomized.R')
##
plan(sequential)
plan(multisession, workers = 3L)
resu <- foreach(chunk=0:2, .inorder=F, .packages=c('data.table','foreach'))%dorng%{ postsamples(chunk)}







tsamples <- mcsamples
dim(tsamples) <- c(nrow(mcsamples),2,(maxS1+1))
matplot(x=0:maxS,y=t(tsamples[1:10,1,1:maxS1]),type='l',lty=1,col=paste0(mygrey,'88'))

tmeansamples <- c(apply(tsamples[,1,1:maxS],1,function(x){sum(x * (0:maxS))}))
hist(tmeansamples)

hist(tsamples[,1,-(1:(maxS1+1))])

tmeansamples2 <- mcsamples[,]
dim(tsamples) <- c(nrow(mcsamples),2,(maxS1+2))
matplot(x=0:maxS,y=t(tsamples[1:10,1,1:maxS1]),type='l',lty=1,col=paste0(mygrey,'88'))



gfactor <- function(mean,pstrength,sfreq){
    tgeo <- pstrength*geomd(x=mean,y=0:maxS)
    exp(sum(lgamma(tgeo+sfreq)-lgamma(tgeo))+lgamma(sum(tgeo))-lgamma(sum(tgeo+sfreq)))
}

postggamma <- function(mean,pstrength,sfreq,pmean,psd){
    sgamma <- (pmean/psd)^2
    rgamma <- sqrt(sgamma)/psd
    gfactor(mean,pstrength,sfreq)*dgamma(x=mean,shape=sgamma,rate=rgamma)
}

normalizerows(sampleFreqs) %*% (0:maxS)

matplot(xgrid <- seq(0,3,length.out=100),sapply(xgrid,function(x){gfactor(x,100,sampleFreqs[2,])}),type='l')

matplot(xgrid <- seq(0,3,length.out=100),sapply(xgrid,function(x){postggamma(x,100,sampleFreqs[1,]*0,0.6,0.5)}),type='l',add=F)
matplot(xgrid <- seq(0,3,length.out=100),sapply(xgrid,function(x){postggamma(x,100,sampleFreqs[1,]/80,0.6,0.5)}),type='l',add=F)
matplot(xgrid <- seq(0,3,length.out=100),sapply(xgrid,function(x){postggamma(x,100,sampleFreqs[1,],0.6,0.5)}),type='l',add=F)



## test with poisson instead of geometric
gfactor <- function(mean,pstrength,sfreq){
    rt <- 2
    tgeo <- pstrength*normalize(dgamma(x=0:maxS,shape=rt,rate=rt/mean))
    exp(sum(lgamma(tgeo+sfreq)-lgamma(tgeo))+lgamma(sum(tgeo))-lgamma(sum(tgeo+sfreq)))
}
postggamma <- function(mean,pstrength,sfreq,pmean,psd){
    sgamma <- (pmean/psd)^2
    rgamma <- sqrt(sgamma)/psd
    gfactor(mean,pstrength,sfreq)*dgamma(x=mean,shape=sgamma,rate=rgamma)
}
matplot(xgrid <- seq(0.001,3,length.out=100),sapply(xgrid,function(x){gfactor(x,100,sampleFreqs[2,])}),type='l')

normalizerows(sampleFreqs) %*% (0:maxS)


matplot(xgrid <- seq(0,3,length.out=100),sapply(xgrid,function(x){postggamma(x,100,sampleFreqs[1,]*0,0.6,0.5)}),type='l',add=F)
matplot(xgrid <- seq(0,3,length.out=100),sapply(xgrid,function(x){postggamma(x,100,sampleFreqs[1,]/80,0.6,0.5)}),type='l',add=F)
matplot(xgrid <- seq(0,3,length.out=100),sapply(xgrid,function(x){postggamma(x,100,sampleFreqs[1,],0.6,0.5)}),type='l',add=F)





plan(sequential)
plan(multisession, workers = 3L)
tme <- 0.6
tsd <- 0.5
tshape <- (tme/tsd)^2
trate <- sqrt(tshape)/tsd
tmeans <- rgamma(n=1000, shape=tshape, rate=trate)
hist(tmeans)
tsamples <- foreach(i=1:1000, .inorder=F, .combine=rbind)%dorng%{
    rdirch(n=1,alpha=100*geomd(x=tmeans[i], y=0:maxS))
}
matplot(x=0:maxS,y=t(mcsamples[1:10,]),type='l',lty=1,col=paste0(mygrey,'88'))

tmeansamples <- c(apply(tsamples,1,function(x){sum(x * (0:maxS))}))
hist(tmeansamples)

hist(tmeans)

matplot(x=0:maxS,y=t(tsamples[1:5,]),type='l',lty=1)

## library('parallel')
## mycluster <- makeCluster(3)
library('nimble')
    dlogsmoothmean <- nimbleFunction(
        run = function(x=double(1), alphas=double(1), powerexp=double(0), shapegamma=double(0), rategamma=double(0), smatrix=double(2), normstrength=double(0, default=1000), log=integer(0, default=0)){
            returnType(double(0))
            tx <- sum(x)
            f <- exp(x)/sum(exp(x))
            dmean <- inprod(f,0:(length(f)-1))
            logp <- sum((alphas+1) * log(f)) + (shapegamma-1)*log(dmean) - (rategamma*dmean)^powerexp - sum((log(f) %*% smatrix)^2) - normstrength  * tx^2 
            if(log) return(logp)
            else return(exp(logp))
        })
    assign('dlogsmoothmean', dlogsmoothmean, envir = .GlobalEnv)
##
runcode <- function(chunk,seed=147){
#### Custom setup ####
## Colour-blind friendly palettes, from https://personal.sron.nl/~pault/
## (consider using khroma package instead)
library('RColorBrewer')
mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
mypalette <- c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black')
palette(mypalette)
barpalette <- colorRampPalette(c(mypurpleblue,'white',myredpurple),space='Lab')
barpalettepos <- colorRampPalette(c('white','black'),space='Lab')
#dev.off()
####
library('data.table')
library('khroma')
library('ggplot2')
library('ggthemes')
theme_set(theme_bw(base_size=18))
scale_colour_discrete <- scale_colour_bright
#library('cowplot')
library('png')
library('foreach')
## library('doFuture')
## registerDoFuture()
## library('doRNG')
library('ash')
#library('nimble')
#library('extraDistr')
options(bitmapType='cairo')
pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in png format
#### End custom setup ####
##
#######################################
#### FUNCTION TO CALCULATE MUTUAL INFO FROM JOINT DISTRIBUTION
## freqs[S,B] = freq spike count B and stimulus S (one ROW per stimulus)
## The function calculates the conditional frequencies of B|S
## and constructs a new joint distribution with equal marginals for S
## Note: don't need to normalize input to mutualinfo
mutualinfo <- function(jointFreqs,base=2){##in bits by default
    stimulusFreqs <- 1/nrow(jointFreqs)
    ## (conditional freqs B|S) * (new freq S)
    jointFreqs <- (jointFreqs/rowSums(jointFreqs)) * stimulusFreqs
    sum(jointFreqs *
        log2(jointFreqs/outer(rowSums(jointFreqs), colSums(jointFreqs))),
        na.rm=TRUE)/log2(base)
}
## function to normalize absolute frequencies
normalize <- function(freqs){freqs/sum(freqs)}
normalizerows <- function(freqs){freqs/rowSums(freqs)}
normalizecols <- function(freqs){t(t(freqs)/colSums(freqs))}
##
t2f <- function(t){exp(t)/sum(exp(t))}
    library('nimble')
longrunDataFile  <- 'SpikeCounts_and_direction.csv'
sampleIndexFile  <- 'index_mat_160.csv'
#plan(sequential)
maxS <- 12
maxS1 <- maxS + 1
maxS2 <- 2 * maxS1
##
## load full recording
longrunData  <- as.data.table(t(read.csv(longrunDataFile,header=FALSE,sep=',')))
colnames(longrunData) <- c('nspikes', 'stimulus')
stimulusVals <- unique(longrunData[,stimulus])
nStimuli <- length(stimulusVals)
nspikesVals <- 0:maxS
## longrunData <- longrunData[nspikes<=maxS] # for debug, REMEMBER TO REMOVE
## frequencies of full recording
longrunFreqs <- foreach(stim=stimulusVals, .combine=rbind)%do%{
    tabulate(longrunData[stimulus==stim,nspikes]+1, nbins=maxS1)
}
rownames(longrunFreqs) <- nameStimulus <- paste0('stimulus',stimulusVals)
colnames(longrunFreqs) <- nameNspikes <- paste0('nspikes',nspikesVals)
longrunMI <- c(bit=mutualinfo(longrunFreqs))
##
#chunk <- 1
    ## Functions definitions
    ##
    ## Normalize rows
    Nnormrows <- nimbleFunction(
        run = function(x=double(2)){
            newx <- matrix(value=0,init=FALSE,nrow=dim(x)[1],ncol=dim(x)[2])
            for(i in 1:(dim(x)[1])){ newx[i,] <- x[i,]/sum(x[i,]) }
            return(newx)
            returnType(double(2))
        })
    assign('Nnormrows', Nnormrows, envir = .GlobalEnv)
    ## Cross-entropy
    Ncrossentropy <- nimbleFunction(
        run = function(x=double(1), y=double(1, default=x), base=double(0, default=2)){
            nzero <- which(x>0)
            return(sum(x[nzero] * log(y[nzero])/log(base)))
            returnType(double(0))
        })
    assign('Ncrossentropy', Ncrossentropy, envir = .GlobalEnv)
    ##Ccentropy <- compileNimble(Ncentropy)    
    ##
    ## Mutual info
    Nmutualinfo <- nimbleFunction(
        run = function(x=double(2), base=double(0, default=2)){
            newx <- Nnormrows(x)/(dim(x)[1])
            marg <- numeric(value=0, length=dim(x)[2])
            for(i in 1:(dim(x)[1])){marg <- marg + newx[i,]}
            return(Ncrossentropy(x=c(newx), y=c(newx), base=base) - Ncrossentropy(x=marg, y=marg, base=base) + log(dim(x)[1])/log(base))
            returnType(double(0))
        })
    assign('Nmutualinfo', Nmutualinfo, envir = .GlobalEnv)
    ## Cmutualinfo <- compileNimble(Nmutualinfo)
    ##
    ## Transform log-frequencies to frequencies
    Nx2f <- nimbleFunction(
        run = function(x=double(2)){
            return(Nnormrows(exp(x)))
            returnType(double(2))
        })
    assign('Nx2f', Nx2f, envir = .GlobalEnv)
    ##
    ##
    ## stimulus0 0.08165385
    ## stimulus1 1.23955008
    ## > summary(lmf)
    ##        V1         
    ##  Min.   :0.08165  
    ##  1st Qu.:0.37113  
    ##  Median :0.66060  
    ##  Mean   :0.66060  
    ##  3rd Qu.:0.95008  
    ##  Max.   :1.23955
    print('HEY')
    priorMeanSpikes <- 1.6 # 0.2 = 5Hz * (40Hz/1000s)
    priorSdSpikes <- 0.15 # 0.2 = 5Hz * (40Hz/1000s)
    shapegamma <- (priorMeanSpikes/priorSdSpikes)^2
    rategamma <- sqrt(shapegamma)/priorSdSpikes
    ## shapegamma <- 1
    ## rategamma <- 1
    powergamma <- 1
    print('data')
    print(c(chunk,shapegamma, rategamma, powergamma))
    prioralphas <- rep(0,maxS1)
    smoothness <- 2
    smoothm <- t(diff(diag(maxS1),differences=2))
    ##
    chunkIndices <- as.matrix(read.csv(sampleIndexFile,header=FALSE,sep=','))[chunk+(chunk==0),]
    sampleData <- longrunData[chunkIndices,]
    ##print(str(sampleData))
    sampleFreqs <- foreach(stim=stimulusVals, .combine=rbind)%do%{
        tabulate(sampleData[stimulus==stim,nspikes]+1, nbins=maxS1)
    } 
    dimnames(sampleFreqs) <- dimnames(longrunFreqs)
    nSamples <- sum(sampleFreqs)
    sampleMI <- c(bit=(chunk>0)*mutualinfo(sampleFreqs) - 2*(chunk==0))
    sampleFreqs <- sampleFreqs * (chunk>0)
    ##
    ##
    ## MONTE CARLO sampling for prior and posterior
    ##
    ##
    ## Probability density
    dlogsmoothmean <- nimbleFunction(
        run = function(x=double(1), alphas=double(1), powerexp=double(0), shapegamma=double(0), rategamma=double(0), smatrix=double(2), normstrength=double(0, default=1000), log=integer(0, default=0)){
            returnType(double(0))
            tx <- sum(x)
            f <- exp(x)/sum(exp(x))
            dmean <- inprod(f,0:(length(f)-1))
            logp <- sum((alphas+1) * log(f)) + (shapegamma-1)*log(dmean) - (rategamma*dmean)^powerexp - sum((log(f) %*% smatrix)^2) - normstrength  * tx^2 
            if(log) return(logp)
            else return(exp(logp))
        })
    assign('dlogsmoothmean', dlogsmoothmean, envir = .GlobalEnv)
    #Cdlogsmoothmean <- compileNimble(dlogsmoothmean)
    lnprob <- nimbleCode({
        for(i in 1:nStimuli){
            X[i,1:maxS1] ~ dlogsmoothmean(alphas=postalphas[i,1:maxS1], powerexp=powergammac, shapegamma=shapegammac, rategamma=rategammac, smatrix=smatrixc[1:maxS1,1:smoothdim], normstrength=1000)
        }
    })
    ##
    constants <- list(postalphas=sampleFreqs+prioralphas, powergammac=powergamma, shapegammac=shapegamma, rategammac=rategamma, smoothdim=ncol(smoothm), smatrixc=smoothness*smoothm, nStimuli=nStimuli, maxS1=maxS1, maxS=maxS)
    ##
    initX <- normalizerows(sampleFreqs+prioralphas+1)
    initX <- log(initX) - rowSums(log(initX))/maxS1
    inits <- list(X=initX+rnorm(length(initX)))
    ##
    ##
    model2 <- nimbleModel(code=lnprob, name='model2', constants=constants, inits=inits, data=list())
    Cmodel2 <- compileNimble(model2, showCompilerOutput = TRUE, resetFunctions = TRUE)
    confmodel2 <- configureMCMC(Cmodel2, nodes=NULL)
    ## confmodel2$addSampler(target='X', type='AF_slice', control=list(sliceAdaptFactorMaxIter=20000, sliceAdaptFactorInterval=1000, sliceAdaptWidthMaxIter=1000, sliceMaxSteps=100, maxContractions=100))
    for(i in 1:nStimuli){
        confmodel2$addSampler(target=paste0('X[',i,',]'), type='AF_slice', control=list(sliceAdaptFactorMaxIter=10000, sliceAdaptFactorInterval=500, sliceAdaptWidthMaxIter=500, sliceMaxSteps=100, maxContractions=100))
    }
    confmodel2$addMonitors('logProb_X')
    confmodel2
    mcmcsampler <- buildMCMC(confmodel2)
    Cmcmcsampler <- compileNimble(mcmcsampler, resetFunctions = TRUE)
    mcsamples <- runMCMC(Cmcmcsampler, nburnin=10000, niter=20000, thin=10, setSeed=seed)
    nDraws <- nrow(mcsamples)
    llsamples <- mcsamples[,-(1:(maxS1*2)),drop=F]
    ##
    condfreqSamples <- t(apply(mcsamples[,1:(maxS1*2)],1,function(x){
        dim(x) <- c(2,maxS1)
        Nx2f(x)}))
    dim(condfreqSamples) <- c(nrow(mcsamples),nStimuli,maxS1)
    dimnames(condfreqSamples) <- list(NULL, nameStimulus, nameNspikes)
    ##
    MIsamples <- apply(condfreqSamples,1,Nmutualinfo)
    MIDistr <- hist(MIsamples, breaks=seq(0,1,by=0.02), plot=F)
    MIQuantiles <- quantile(x=MIsamples, probs=c(0.025,0.5,0.975))
    ##
    meanSsamples <- t(apply(condfreqSamples,1,function(x){x %*% (0:maxS)}))
    ##
    ## PLOTS
    nPlotSamples <- 100
    maxX <- 8
    if(chunk==0){psign <- 1}else{psign <- -1}
    pdff(paste0('testsummaryN',chunk))
    matplot(x=nspikesVals, y=t(condfreqSamples[round(seq(1,nDraws,length.out=nPlotSamples)),1,]),
            type='l', lty=1, lwd=2, col=paste0(mygrey,'66'), ylim=c(min(0,psign),1),  xlim=c(0,maxX), xlab='spikes/bin', ylab='freq', cex.lab=2, cex.axis=2)
        matplot(x=nspikesVals, y=psign*t(condfreqSamples[round(seq(1,nDraws,length.out=nPlotSamples)),2,]),
                type='l', lty=1, lwd=1, col=paste0(mygrey,'66'), add=TRUE)
    ##
    if(pflag==0){matplot(x=nspikesVals, y=t(normalizerows(sampleFreqs)*c(1,psign)),
                         type='l', lty=2, lwd=5, col=myyellow, add=TRUE)}
    ##
    matplot(x=nspikesVals, y=t(normalizerows(longrunFreqs)*c(1,psign)),
            type='l', lty=4, lwd=4, col='black', add=TRUE)
    ##
    title(paste0(nSamples,' data samples,',
                 ' chunk =', chunk,
                 ', prior weight = ', sum(prioralphas),
                 '\n superdistr ',chunk), cex.main=2)
    legend('topright',c('long-run freqs','sample freqs'),lty=c(1,2),lwd=c(2,5),col=c('black',myyellow),cex=1.5)
    ##
    ##
    matplot(x=MIDistr$mids, y=MIDistr$density,
            type='h', lty=1, lwd=15, col=paste0(mypurpleblue,'88'), xlim=c(0,1),
            xlab='MI/bit', ylab='prob dens', cex.lab=2, cex.axis=2)
    for(q in MIQuantiles){
        matlines(x=rep(q,2),y=c(-1,1/2)*max(MIDistr$density), lty=2, lwd=6, col=mygreen)
    }
    matlines(x=rep(sampleMI,2),y=c(-1,2/3)*max(MIDistr$density), lty=4, lwd=6, col=myyellow)
    matlines(x=rep(longrunMI,2),y=c(-1,2/3)*max(MIDistr$density), lty=1, lwd=6, col=myredpurple)
    title('predicted MI distr', cex.main=2)
    legend('topright',c('sample MI', 'long-run MI'),lty=1,col=c(myyellow,myredpurple),lwd=4,cex=1.5)
    ##
    ## Diagnostics
    hist(meanSsamples[,1], xlim=c(0,3),ylab='mean spikecounts')
    hist(meanSsamples[,2], xlim=c(0,3),ylab='mean spikecounts')
    matplot((MIsamples),type='l', lty=1,ylab='MI samples')
    matplot((llsamples[,]),type='l', lty=1,ylab='log-posterior')
    matplot((mcsamples[,1]),type='l', lty=1,ylab='samples of first freq')
    dev.off()
    NULL
}
##
##
## plan(sequential)
## plan(multisession, workers = 3L)
## clusterExport(cl=mycluster, c('runcode'))
## alloutput <- parLapply(cl = mycluster, X = 0:2, fun = function(chunk){runcode(chunk)})
## stopCluster(mycluster)

corrp <- matrix(0,11,4)
dimnames(corrp) <- list(paste0('count',0:10), paste0('stimulus_',0:1,'_then_',rep(0:1,each=2)))
batches <- matrix(foreach(i1=0:1,.combine=c)%:%foreach(i2=0:1,.combine=c)%do%{c(i1,i2)},2)
for(i in 1:ncol(batches)){colnames(corr3p)[1+sum(c(1,2)*batches[,i])] <- paste0('seq_',paste0(batches[,i],collapse="_"))}
for(i in 2:nrow(longrunData)){
    batch <- as.matrix(longrunData[(i-1):i])
    group <- 1+sum((1:2)*batch[,2])
    corrp[batch[2,1]+1,group] <- corrp[batch[2,1]+1,group] + 1
}
corrf <- t(t(corrp)/colSums(corrp))
##
matplot(0:10,corrf,type='l',lty=c(1,2), lwd=3,col=c(myred,myredpurple,myblue,mypurpleblue),xlab='spike count',ylab='long-run relative frequency')
legend('topright',colnames(corrf),lty=c(1,2),lwd=3,col=c(myred,myredpurple,myblue,mypurpleblue))
##
##
rcorrp <- matrix(0,11,4)
dimnames(rcorrp) <- list(paste0('count',0:10), paste0('stimulus_',0:1,'_then_',rep(0:1,each=2)))
batches <- matrix(foreach(i1=0:1,.combine=c)%:%foreach(i2=0:1,.combine=c)%do%{c(i1,i2)},2)
for(i in 1:ncol(batches)){colnames(corr3p)[1+sum(c(1,2)*batches[,i])] <- paste0('seq_',paste0(batches[,i],collapse="_"))}
##
randomlrd <- longrunData
for(i in 0:1){
whichstim <- which(longrunData$stimulus==i)
randomlrd[whichstim] <- longrunData[sample(whichstim)]
}
for(i in 2:nrow(longrunData)){
    batch <- as.matrix(randomlrd[(i-1):i])
    group <- 1+sum((1:2)*batch[,2])
    rcorrp[batch[2,1]+1,group] <- rcorrp[batch[2,1]+1,group] + 1
}
rcorrp <- rcorrp
rcorrf <- t(t(rcorrp)/colSums(rcorrp))
##
pdff('2step_condfreqs')
matplot(0:10,corrf,type='l',lty=c(1,2), lwd=3,col=c(myred,myredpurple,myblue,mypurpleblue),xlab='spike count',ylab='long-run relative frequency',main='data')
legend('topright',paste0(colnames(corrf),': ',colSums(corrp)),lty=c(1,2),lwd=3,col=c(myred,myredpurple,myblue,mypurpleblue))
matplot(0:10,rcorrf,type='l',lty=c(1,2), lwd=1,col=c(myred,myredpurple,myblue,mypurpleblue),xlab='spike count',ylab='long-run relative frequency',main='randomized within each stimulus')
for(sample in 1:10){
    rcorrp <- matrix(0,11,4)
    dimnames(rcorrp) <- list(paste0('count',0:10), paste0('stimulus_',0:1,'_then_',rep(0:1,each=2)))
    batches <- matrix(foreach(i1=0:1,.combine=c)%:%foreach(i2=0:1,.combine=c)%do%{c(i1,i2)},2)
    for(i in 1:ncol(batches)){colnames(corr3p)[1+sum(c(1,2)*batches[,i])] <- paste0('seq_',paste0(batches[,i],collapse="_"))}
    ##
    randomlrd <- longrunData
    for(i in 0:1){
        whichstim <- which(longrunData$stimulus==i)
        randomlrd[whichstim] <- longrunData[sample(whichstim)]
    }
    for(i in 2:nrow(longrunData)){
        batch <- as.matrix(randomlrd[(i-1):i])
        group <- 1+sum((1:2)*batch[,2])
        rcorrp[batch[2,1]+1,group] <- rcorrp[batch[2,1]+1,group] + 1
    }
    rcorrp <- rcorrp
    rcorrf <- t(t(rcorrp)/colSums(rcorrp))
    matplot(0:10,rcorrf,type='l',lty=c(1,2), lwd=1,col=c(myred,myredpurple,myblue,mypurpleblue),xlab='spike count',ylab='long-run relative frequency',main='randomized within each stimulus',add=T)
}
legend('topright',colnames(rcorrf),lty=c(1,2),lwd=3,col=c(myred,myredpurple,myblue,mypurpleblue))
dev.off()



#### corr now & future
fcorrp <- matrix(0,11,4)
dimnames(fcorrp) <- list(paste0('count',0:10), paste0('stimulus_now',0:1,'_then_',rep(0:1,each=2)))
batches <- matrix(foreach(i1=0:1,.combine=c)%:%foreach(i2=0:1,.combine=c)%do%{c(i1,i2)},2)
for(i in 1:ncol(batches)){colnames(fcorrp)[1+sum(c(1,2)*batches[,i])] <- paste0('seq_',paste0(batches[,i],collapse="_"))}
for(i in 1:(nrow(longrunData)-1)){
    batch <- as.matrix(longrunData[i:(i+1)])
    group <- 1+sum((1:2)*batch[,2])
    fcorrp[batch[1,1]+1,group] <- fcorrp[batch[1,1]+1,group] + 1
}
fcorrf <- t(t(fcorrp)/colSums(fcorrp))
##
matplot(0:10,fcorrf,type='l',lty=c(1,2), lwd=3,col=c(myred,myredpurple,myblue,mypurpleblue),xlab='spike count',ylab='long-run relative frequency')
legend('topright',colnames(fcorrf),lty=c(1,2),lwd=3,col=c(myred,myredpurple,myblue,mypurpleblue))
##
##
rfcorrp <- matrix(0,11,4)
dimnames(rfcorrp) <- list(paste0('count',0:10), paste0('stimulus_',0:1,'_then_',rep(0:1,each=2)))
batches <- matrix(foreach(i1=0:1,.combine=c)%:%foreach(i2=0:1,.combine=c)%do%{c(i1,i2)},2)
for(i in 1:ncol(batches)){colnames(fcorrp)[1+sum(c(1,2)*batches[,i])] <- paste0('seq_',paste0(batches[,i],collapse="_"))}
##
randomlrd <- longrunData
for(i in 0:1){
whichstim <- which(longrunData$stimulus==i)
randomlrd[whichstim] <- longrunData[sample(whichstim)]
}
for(i in 1:(nrow(longrunData)-1)){
    batch <- as.matrix(randomlrd[i:(i+1)])
    group <- 1+sum((1:2)*batch[,2])
    rfcorrp[batch[1,1]+1,group] <- rfcorrp[batch[1,1]+1,group] + 1
}
rfcorrp <- rfcorrp
rfcorrf <- t(t(rfcorrp)/colSums(rfcorrp))
##
pdff('2step_condfreqs_future')
matplot(0:10,fcorrf,type='l',lty=c(1,2), lwd=3,col=c(myred,myredpurple,myblue,mypurpleblue),xlab='spike count',ylab='long-run relative frequency',main='data')
legend('topright',paste0(colnames(fcorrf),': ',colSums(fcorrp)),lty=c(1,2),lwd=3,col=c(myred,myredpurple,myblue,mypurpleblue))
matplot(0:10,rfcorrf,type='l',lty=c(1,2), lwd=1,col=c(myred,myredpurple,myblue,mypurpleblue),xlab='spike count',ylab='long-run relative frequency',main='randomized within each stimulus')
for(sample in 1:10){
    rfcorrp <- matrix(0,11,4)
    dimnames(rfcorrp) <- list(paste0('count',0:10), paste0('stimulus_',0:1,'_then_',rep(0:1,each=2)))
    batches <- matrix(foreach(i1=0:1,.combine=c)%:%foreach(i2=0:1,.combine=c)%do%{c(i1,i2)},2)
    for(i in 1:ncol(batches)){colnames(fcorrp)[1+sum(c(1,2)*batches[,i])] <- paste0('seq_',paste0(batches[,i],collapse="_"))}
    ##
    randomlrd <- longrunData
    for(i in 0:1){
        whichstim <- which(longrunData$stimulus==i)
        randomlrd[whichstim] <- longrunData[sample(whichstim)]
    }
    for(i in 1:(nrow(longrunData)-1)){
        batch <- as.matrix(randomlrd[i:(i+1)])
        group <- 1+sum((1:2)*batch[,2])
        rfcorrp[batch[1,1]+1,group] <- rfcorrp[batch[1,1]+1,group] + 1
    }
    rfcorrp <- rfcorrp
    rfcorrf <- t(t(rfcorrp)/colSums(rfcorrp))
    matplot(0:10,rfcorrf,type='l',lty=c(1,2), lwd=1,col=c(myred,myredpurple,myblue,mypurpleblue),xlab='spike count',ylab='long-run relative frequency',main='randomized within each stimulus',add=T)
}
legend('topright',colnames(rfcorrf),lty=c(1,2),lwd=3,col=c(myred,myredpurple,myblue,mypurpleblue))
dev.off()







corr3p <- matrix(0,11,8)
dimnames(corr3p) <- list(paste0('count',0:10), paste0('stimulus_',1:8))
batches <- matrix(foreach(i1=0:1,.combine=c)%:%foreach(i2=0:1,.combine=c)%:%foreach(i3=0:1,.combine=c)%do%{c(i1,i2,i3)},3)
for(i in 1:ncol(batches)){colnames(corr3p)[1+sum(c(1,2,4)*batches[,i])] <- paste0('seq_',paste0(batches[,i],collapse="_"))}
for(i in 3:nrow(longrunData)){
    batch <- as.matrix(longrunData[(i-2):i])
    group <- 1+sum(c(1,2,4)*batch[,2])
    corr3p[batch[3,1]+1,group] <- corr3p[batch[3,1]+1,group] + 1
}
corr3f <- t(t(corr3p)/colSums(corr3p))
##
matplot(0:10,corr3f,type='l',lty=c(1,2,3,4), lwd=3,col=c(myredpurple,mypurpleblue),xlab='spike count',ylab='long-run relative frequency')
legend('topright',colnames(corr3f),lty=c(1,2,3,4),lwd=3,col=c(myredpurple,mypurpleblue))
##
##
rcorr3p <- matrix(0,11,8)
dimnames(rcorr3p) <- list(paste0('count',0:10), paste0('stimulus_',1:8))
batches <- matrix(foreach(i1=0:1,.combine=c)%:%foreach(i2=0:1,.combine=c)%:%foreach(i3=0:1,.combine=c)%do%{c(i1,i2,i3)},3)
for(i in 1:ncol(batches)){colnames(rcorr3p)[1+sum(c(1,2,4)*batches[,i])] <- paste0('seq_',paste0(batches[,i],collapse="_"))}
##
randomlrd <- longrunData
for(i in 0:1){
whichstim <- which(longrunData$stimulus==i)
randomlrd[whichstim] <- longrunData[sample(whichstim)]
}
for(i in 3:nrow(longrunData)){
    batch <- as.matrix(randomlrd[(i-2):i])
    group <- 1+sum(c(1,2,4)*batch[,2])
    rcorr3p[batch[3,1]+1,group] <- rcorr3p[batch[3,1]+1,group] + 1
}
rcorr3p <- rcorr3p
rcorr3f <- t(t(rcorr3p)/colSums(rcorr3p))
##
pdff('3step_condfreqs')
matplot(0:10,corr3f,type='l',lty=c(1,2,3,4), lwd=3,col=rep(c(myredpurple,mypurpleblue),each=4),xlab='spike count',ylab='long-run relative frequency',main='data')
legend('topright',paste0(colnames(corr3f),': ',colSums(corr3p)),lty=c(1,2,3,4),lwd=3,col=rep(c(myredpurple,mypurpleblue),each=4))
matplot(0:10,rcorr3f,type='l',lty=c(1,2,3,4), lwd=1,col=rep(c(myredpurple,mypurpleblue),each=4),xlab='spike count',ylab='long-run relative frequency',main='randomized within each stimulus')
for(sample in 1:10){
    rcorr3p <- matrix(0,11,8)
    ## dimnames(rcorr3p) <- list(paste0('count',0:10), paste0('stimulus_',1:8))
    ## batches <- matrix(foreach(i1=0:1,.combine=c)%:%foreach(i2=0:1,.combine=c)%:%foreach(i3=0:1,.combine=c)%do%{c(i1,i2,i3)},3)
    ## for(i in 1:ncol(batches)){colnames(corr3p)[1+sum(c(1,2,4)*batches[,i])] <- paste0('seq_',paste0(batches[,i],collapse="_"))}
    ##
    randomlrd <- longrunData
    for(i in 0:1){
        whichstim <- which(longrunData$stimulus==i)
        randomlrd[whichstim] <- longrunData[sample(whichstim)]
    }
    for(i in 3:nrow(longrunData)){
        batch <- as.matrix(randomlrd[(i-2):i])
        group <- 1+sum(c(1,2,4)*batch[,2])
        rcorr3p[batch[3,1]+1,group] <- rcorr3p[batch[3,1]+1,group] + 1
    }
    rcorr3p <- rcorr3p
    rcorr3f <- t(t(rcorr3p)/colSums(rcorr3p))
    matplot(0:10,rcorr3f,type='l',lty=c(1,2,3,4), lwd=1,col=rep(c(myredpurple,mypurpleblue),each=4),xlab='spike count',ylab='long-run relative frequency',main='randomized within each stimulus',add=T)
}
legend('topright',colnames(corr3f),lty=c(1,2,3,4),lwd=3,col=rep(c(myredpurple,mypurpleblue),each=4))
dev.off()






















paste0('stimulusseq_',0:1,rep(0:1,each=2),rep(0:1,each=4))

corr3p <- matrix(0,11,8)
dimnames(corr3p) <- list(paste0('count',0:10), paste0('stimulusseq_',0:1,rep(0:1,each=2),rep(0:1,each=4)))
batches <- matrix(foreach(i1=0:1,.combine=c)%:%foreach(i2=0:1,.combine=c)%:%foreach(i3=0:1,.combine=c)%do%{c(i1,i2,i3)},3)
for(i in 1:ncol(batches)){colnames(corr3p)[1+sum(c(1,2,4)*batches[,i])] <- paste0('seq_',paste0(batches[,i],collapse="_"))}
for(i in 3:nrow(longrunData)){
    batch <- as.matrix(longrunData[(i-2):i])
    group <- 1+sum(c(1,2,4)*batch[,2])
    corr3p[batch[3,1]+1,group] <- corr3p[batch[3,1]+1,group] + 1
}

corr3f <- t(t(corr3p)/colSums(corr3p))

matplot(0:10,corr3f,type='l',lty=c(1,2,3,4), lwd=2,col=coluse <- rep(c(myred,myblue),each=4),xlab='spike count',ylab='long-run relative frequency')
legend('topright',colnames(corr3f),lty=c(1,2,3,4),lwd=3,col=coluse)














set.seed(333)
nn <- 20000
points <- matrix(runif(2*nn),ncol=2)
points <- points[order(points[,1]),]
f <- 0.88 #(nn)/nn
alpp <- '88'
alpm <- '88'
cexp <- 0.5
cexm <- 0.5
pointsp <- if(f>0){points[1:round(nn*f),,drop=F]}else{NULL}
pointsm <-  if(f<1){points[(round(nn*f)+1):nn,,drop=F]}else{NULL}
png(filename=paste0('f','pred2','.png'),width=10,height=10,units='cm',res=400)
par(pty='s',mar=rep(0,4))
if(f>=0.5){
matplot(pointsp[,1],pointsp[,2],type='p',pch=16,cex=cexp,xaxt='n',yaxt='n',xlab=NA,ylab=NA,col=paste0(palette()[1],alpp),xlim=c(0,1),ylim=c(0,1),pty='s')
matplot(pointsm[,1],pointsm[,2],type='p',pch=16,cex=cexm,xaxt='n',yaxt='n',ylab=NA,col=paste0(palette()[6],alpm),add=T)
} else {
matplot(pointsm[,1],pointsm[,2],type='p',pch=16,cex=cexm,xaxt='n',yaxt='n',xlab=NA,ylab=NA,col=paste0(palette()[2],alpm),xlim=c(0,1),ylim=c(0,1),pty='s')
matplot(pointsp[,1],pointsp[,2],type='p',pch=16,cex=cexp,xaxt='n',yaxt='n',xlab=NA,ylab=NA,col=paste0(palette()[1],alpp),xlim=c(0,1),ylim=c(0,1),add=T)    
}
dev.off()

library('ggplot2')
library('ggthemes')
theme_set(theme_bw(base_size=18))
scale_colour_discrete <- scale_colour_bright
library('data.table')

positive <- 13
negative <- 4

mean <- 0.5
sd <- 0.2
betaShape1 <- ((1 - mean) * mean/sd^2 - 1) * mean
betaShape2 <- betaShape1*(1-mean)/mean
predata <- function(f){dbeta(f, shape1=betaShape1, shape2=betaShape2)}
numerator <- function(f){
    choose(positive+negative, negative) * f^positive * (1-f)^negative *
        predata(f)
}
denominator <- integrate(numerator, lower=0, upper=1)$value
xgrid <- seq(0, 1, length.out=100)
toPlot <- rbind(data.frame(f=xgrid, probability=predata(xgrid), given='initial assumption'),
                data.frame(f=xgrid, probability=numerator(xgrid)/denominator, given='data')
                )
## matplot(x=xgrid, y=cbind( predata(xgrid), numerator(xgrid)/denominator),
##       type='l', xlab='f', ylab='probability')
## grid()
## legend(x='topleft',legend=c('pre-data', 'post-data'),col=)
qplot(f, probability, data=toPlot, color=given, geom='line', lty=given, lwd=I(2)) + theme(legend.pos='top')



      y=cbind(
                     predata(xgrid)/normalizationPredata,
                     numerator(xgrid)/denominator),
        type='l', xlab='f', ylab='probability')



qplot(f, p, data=toPlot, color=which y=cbind(
                     predata(xgrid)/normalizationPredata,
                     numerator(xgrid)/denominator),
        type='l', xlab='f', ylab='probability')




qplot(x=xgrid, y=cbind(
                     predata(xgrid)/normalizationPredata,
                     numerator(xgrid)/denominator),
        type='l', xlab='f', ylab='probability')
grid()
legend(x='topleft',legend=c('pre-data', 'post-data'))

print('probability for f > 0.5')
print(integrate(numerator, lower=0.5, upper=1)$value/denominator)



library('ggplot2')

## Data:
positive <- 13
negative <- 4

## Parameters for pre-data distribution (mean and standard deviation):
mean <- 0.5
sd <- 0.2

betaShape1 <- ((1 - mean) * mean/sd^2 - 1) * mean # shape-parameters of beta distribution
betaShape2 <- betaShape1 * (1 - mean)/mean

## Pre-data distribution (represented by a beta distribution, https://mathworld.wolfram.com/BetaDistribution.html)
predata <- function(f){dbeta(f, shape1=betaShape1, shape2=betaShape2)}

## Numerator and denominator of pre-data distribution:
numerator <- function(f){
    choose(positive+negative, negative) * f^positive * (1-f)^negative * predata(f)
}
denominator <- integrate(numerator, lower=0, upper=1)$value # integral approximates sum

## Plot the two distributions:
fgrid <- seq(0, 1, length.out=1000) # create a grid of f-coordinates
toPlot <- rbind(data.table(f=fgrid, probability=predata(fgrid),
                           given='initial assumption'),
                data.table(f=fgrid, probability=numerator(fgrid)/denominator, given='data'))

qplot(f, probability, data=toPlot, geom='line',
      color=given, lty=given, lwd=I(1.5)) + theme(legend.pos='top')

## Print probability for f > 0.5, given the data
print('probability for f > 0.5:')
print(integrate(numerator, lower=0.5, upper=1)$value/denominator)
