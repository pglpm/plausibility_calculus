### definitions.R 
## Author: Porta Mana, Bachmann
## Created: 2017-11-29T11:48:44+0100
## Last-Updated: 2019-08-30T13:59:59+0200
##
## Definitions of functions for files 'sample1.R' and 'plotresults.R'

library('fBasics')
library('LaplacesDemon') ## provides logit and other useful functions
library('ggplot2')
library('RColorBrewer')
library('mvtnorm')
library('cowplot')

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


## functions used by the models defined in sect 2.4.4
## modified logit + Jacobian
logit2 <- function(x){log(1+x)-log(1-x)}
jlogit <- function(x){2/((1+x)*(1-x))}

## modified tangent + Jacobian
tg2 <- function(x){tan(x*pi/2)}
jtg <- function(x){pi/(1+cos(pi*x))}

## normal + Jacobian
id2 <- function(x){x}
jid <- function(x){1}

## Jacobian of ln
jlog <- function(x){1/x}

## multivariate normal (not used in scripts)
## mvn <- function(x,mu,sigma){
##    exp(-c((x - mu) %*% solve(sigma) %*% (x - mu))/2)/sqrt(det(2*pi*sigma))
## }

# Used in function 'logevidence': multivariate t distribution as in eqn (36)
likelihood <- function(x,parameters){
    l <- length(x)
    kappa <- parameters[[2]]
    nu <- parameters[[1]] - l + 1
    dmvt(x, delta=parameters[[3]],
         sigma=parameters[[4]]*(kappa+1)/kappa/nu, df=nu, log=F)
}

## Used in function 'logevidence': log of multivariate t above
loglikelihood <- function(x,parameters){
    l <- length(x)
    kappa <- parameters[[2]]
    nu <- parameters[[1]] - l + 1
    dmvt(x, delta=parameters[[3]],
         sigma=parameters[[4]]*(kappa+1)/kappa/nu, df=nu, log=T)
}

## Used in function 'logevidence': update prior coefficients with new
## datum, as in eqn (34) with n=1
updateparameters1 <- function(datum,parameters){
    nu <- parameters[[1]]
    kappa <- parameters[[2]]
    newkappa <- kappa + 1
    mu <- parameters[[3]]
    sigma <- parameters[[4]]
    list(nu + 1, newkappa,
    (kappa*mu + datum)/newkappa,
    sigma + kappa * ((datum-mu) %*% t(datum-mu))/newkappa )
} # checked against same alg in mathematica

## (Not used) update prior coefficients with several data
## updateparameters <- function(data,parameters){
##     N <- dim(data)[2]
##     means <- colMeans(t(data))
##     scatm <- cov(t(data))*(N-1)
##     nu <- parameters[[1]]
##     kappa <- parameters[[2]]
##     newkappa <- kappa + N
##     mu <- parameters[[3]]
##     sigma <- parameters[[4]]
##     list(nu + N, newkappa,
##     (kappa*mu + N * means)/newkappa,
##     sigma + scatm + kappa * N * ((means-mu) %*% t(means-mu))/newkappa )
## }

## Prior coefficients used by the three models of sect 2.4.4
dd <- 40 ## number of graph quantities used in the paper
priorparameters0d <- function(d){list(d+1, 1, rep(0,d), (2+d)*diag(rep(1,d)))}
priorparameters0 <- priorparameters0d(dd)
priorparametersf <- function(d,kappa,mu){list(d+1, kappa, rep(0,d), mu*(2+d)*diag(rep(1,d)))}
priorparametersk <- function(kappa,mu,sigma){list(dd+1, kappa, rep(mu,dd), diag(rep(sigma,dd)))}
## cupprior <- list(dd + 1, 0.1, rep(0,dd), (2*dd+3)*diag(rep(1,dd))) ## not used

## Multivariate t distribution (not used)
## mvt <- function(x,nu,mu,sigma){
##     l <- length(mu)
##    Gamma[(nu + l)/2]/Gamma[nu/2]/
##      Sqrt@Det[
##        nu*Pi*sigma]/(1 + (x - mu).Inverse[nu*sigma].(x - mu))^((nu + 
##          l)/2)];
## }

## Names of connectivities (not used)
## graphprop <- c("f[73, 90]", "f[9, 78]", "f[54, 91]", "f[22, 50]", "f[30, 31]", "f[75, 78]", "f[4, 31]", "f[9, 51]", "f[31, 32]", "f[24, 31]", "f[23, 36]", "f[72, 78]", "f[36, 92]", "f[32, 58]", "f[15, 35]", "f[29, 84]", "f[15, 85]", "f[25, 70]", "f[54, 85]", "f[9, 50]", "f[49, 68]", "f[19, 90]", "f[51, 91]", "f[63, 91]", "f[35, 60]", "f[9, 49]", "f[27, 62]", "f[36, 85]", "f[2, 14]", "f[31, 51]", "f[51, 68]", "f[50, 68]", "f[36, 91]", "f[52, 68]", "f[32, 74]", "f[2, 15]", "f[37, 58]", "f[9, 38]", "f[51, 75]", "f[25, 65]")

## Other function, not used
## transfd <- function(data){mapply(function(f,x) f(x),transf,data)}


##
## This is the main function for the calculation of utility and surprise sequences
##
logevidence <- function(nsamples, ## number of samples
                        categories, ## list of health conditions (strings)
                        quantities, ## names of graph quantities (not used)
                        datafiles, ## list of files with the data
                        properties, ## vector of properties to choose from data
                        transformations, ## list of functions "l" for generalized normal
                        scales, ## graph-quantity rescalings
                        jacobians, ## Jacobians of transformations
                        prior, ## coefficients of prior
                        pretestprob, ## pre-test probabilities
                        utilitymatrix, ## utility matrix
                        rawmoments=F, ## whether moments of quantities are raw (not used)
                        seed=666 ## seed for reproducible calculations
                        ){
    n <- length(categories) # num. health categories
    ##d <- dim(read.matrix(datafiles[1]))[1] # num. graph quantities
    d <- length(properties) # num. graph quantities
    print(paste0(d,' graph quantities'))
    if(length(quantities)==1){quantities <- sprintf("q[%d]",1:d)}
    ## If "prior" has 4 elements, it means same prior for all health categories
    ## otherwise we have a list of priors, one for each health category
    ## (this was written to deal with 2 or 3 health categories)
    if(length(prior)==4){oldprior <- prior
        prior <- list()
        for(i in 1:n){prior[[i]] <- oldprior}
        }
    datao <- list() # original data
    datas <- list() # scaled data
    datat <- list() # logit/log-transformed scaled data
    N <- rep(NA,n) # num. patients per health category
    for(i in 1:n){
        categ <- categories[i]
        ## read data
        datao[[categ]] <- if(length(properties)>1){read.matrix(datafiles[i])[properties,]}else{t(as.matrix(read.matrix(datafiles[i])[properties,]))}
        N[i] <- dim(datao[[categ]])[2]
        print(paste0(N[i], " ", categ, ' patients'))
        if(rawmoments == T){print('Calculating raw moments...')
        ## replace central moments with raw moments
        tempdata <- datao[[categ]]
        datanew <- datao[[categ]]
        ## transform 2nd moments to raw ones
        for(j in seq(1,21,4)){
            datanew[j + 1,] <- sqrt(tempdata[j,]^2 + tempdata[j + 1,])
            ## transform 3nd moments to raw ones
            datanew[j + 2,] <- (tempdata[j,]^3 + 3*tempdata[j,]*tempdata[j + 1,] + tempdata[j + 2,])^(1/3)
            ## transform 4nd moments to raw ones
            tempdata[j + 3,] <- (tempdata[j,]^4 + 6*tempdata[j,]^2*tempdata[j + 1,] + 4*tempdata[j,]*tempdata[j + 2,] + tempdata[j + 3,])^(1/4)
        }
        datao[[categ]] <- datanew
        }
        rownames(datao[[categ]]) <- quantities
        colnames(datao[[categ]]) <- sprintf("id[%d]",1:(N[i]))
        print('Rescaling the data...')
        datas[[categ]] <- datao[[categ]]/scales # scaled data
        ## create logit/log-transformed scaled data
        datat[[categ]] <- datas[[categ]]
        for(j in 1:(N[i])){
           datat[[categ]][,j] <- mapply(function(f,x) f(x),transformations,datas[[categ]][,j])
        }
    }
    NN <- sum(N) # total number of patients
    print(paste0(NN, ' patients in total.'))
    probsamples1 <- matrix(NA,nsamples,NN) # probabilities for graph quantities
    probsamples2 <- probsamples1 # probabilities for health conditions
    ## hits <- probsamples1 # utilities
    conditionsequences <- probsamples1 # sequences of sampled health conditions
    probdistsamples2 <- array(NA,dim=c(nsamples,NN,n)) # prob. distr. for health conditions
    probdist1 <- rep(NA,n) # empty container for prob. distribution
    ##
    print('Starting sampling...')
    ## progress bar
    pb <- txtProgressBar(min = 1, max = nsamples, style = 3)
    for(k in 1:nsamples){
        set.seed(seed+k) # if we want the same patient sequence
        train <- sample(NN) # randomize order of patients
        ## set parameters to initial prior for each health category
        newprior <- prior
        
        for(i in 1:NN){ # take in randomized training data in succession
            patient <- train[i] # current patient
            ## loop to find out health category of current patient
            categ <- 1
            while(patient > N[categ]){
                patient <- patient - N[categ]
                categ <- categ + 1}
            datum <- datat[[categ]][,patient] # graph data for current patient
            ## probability of graph data given each health condition
            for(j in 1:n){ probdist1[j] <- likelihood(datum,newprior[[j]]) }
            ## probability distribution of health condition
            probdist2 <- (probdist1 * pretestprob)/c(probdist1 %*% pretestprob)
            jacobian <- prod(mapply(function(f,x) f(x),jacobians,datas[[categ]][,patient]))
            ## Now we save all results for this patient:
            ##
            ## probability of current datum
            probsamples1[k,i] <- probdist1[categ]*jacobian
            ## probability of true health condition
            probsamples2[k,i] <- probdist2[categ]
            ## full probability distribution for health conditions
            probdistsamples2[k,i,] <- probdist2
            ## true health condition
            conditionsequences[k,i] <- categ
            ## utility gained
            #hits[k,i] <- utilitymatrix[which.max(c(utilitymatrix %*% probdist2)),categ]
            ##
            ## Update normal-inv.Wishart coefficients with current patient
            ## (this step is unnecessary for the last patient; script could
            ## be made faster by taking last patient out of for-loop and
            ## eliminating this calculation)
            newprior[[categ]] <- updateparameters1(datum,newprior[[categ]])
            setTxtProgressBar(pb, k)
        }
    }
    close(pb)
    return(list(probsamples1=probsamples1,
                probsamples2=probsamples2,
                probdistsamples2=probdistsamples2,
                conditionsequences=conditionsequences#, hits=hits
                ))
}


##
## Functions to calculate averages from sampled data
##
expect <- function(mean,vari){
    mean + vari/2+
        pnorm((mean+vari)/sqrt(vari),lower=F,log.p=T) -
        pnorm(mean/sqrt(vari),lower=F,log.p=T)
}

expectlist <- function(list){
    mean <- mean(list)
    vari <- var(list)
    mean + vari/2+
        pnorm((mean+vari)/sqrt(vari),lower=F,log.p=T) -
        pnorm(mean/sqrt(vari),lower=F,log.p=T)
}


hitsplot <- function(samples,pdfname,yrange=NULL){
    x <- 1:(dim(samples[[1]][[2]])[2])
    g <- ggplot()
    for(i in 1:length(samples)){
    df <- data.frame(x, y=colMeans(samples[[i]][[5]]))
    g <- g + geom_point(data=df, aes(x,y), colour=mycolours[i]) +
        geom_line(data=df, aes(x,y), colour=mycolours[i])
    }
    if(is.vector(yrange)){ g <- g + ylim(yrange[1],yrange[2])}
    pdf(pdfname)
    print(g)
    dev.off()
}

cumhitsplot <- function(samples,pdfname,yrange=NULL){
    x <- 1:(dim(samples[[1]][[2]])[2])
    g <- ggplot()
    for(i in 1:length(samples)){
    df <- data.frame(x, y=cumsum(colMeans(samples[[i]][[5]])))
    g <- g + geom_point(data=df, aes(x,y), colour=mycolours[i]) +
        geom_line(data=df, aes(x,y), colour=mycolours[i])
    }
    if(is.vector(yrange)){ g <- g + ylim(yrange[1],yrange[2])}
    pdf(pdfname)
    print(g)
    dev.off()
}



prob2logev <- function(probsamples){t(apply(log(probsamples),1,cumsum))}
prob2meanlogev <- function(probsamples){colMeans(t(apply(log(probsamples[[2]]),1,cumsum)))}

prob2expsurprise <- function(probsamples,whichsample=2){apply(log(probsamples[[whichsample]]),2,expectlist)}

prob2expcumsurprise <- function(probsamples,whichsample=2){apply(t(apply(log(probsamples[[whichsample]]),1,cumsum)),2,expectlist)}

prob2expsurpriseemp <- function(probsamples,whichsample=2){log(colMeans(probsamples[[whichsample]]))}

prob2expcumsurpriseemp <- function(probsamples,whichsample=2){log(colMeans(exp(t(apply(log(probsamples[[whichsample]]),1,cumsum)))))}


##
## Main plot function
##
genplot <- function(samples,pdfname,pnames,yrange=NULL,xrange=NULL,xytext=NULL,legenda=NULL){
    numgraphs <- dim(samples)[1]-1
    x <- 1:length(samples[1,])
    alldata <- data.frame(x, y=samples[1,])
    for(i in 2:(numgraphs+1)){
        alldata <- rbind(alldata,data.frame(x, y=samples[i,]))}
    alldata$grp <- rep(factor(1:(numgraphs+1), labels=pnames),times=rep(length(x),numgraphs+1))
    cols <- c(mycolours[1:numgraphs],'black')
    g <- ggplot() + scale_color_manual(values=cols) + theme_classic() +
        scale_shape_manual(values=c(16,17,15,-1)) +
        scale_linetype_manual(values=c('solid','solid','solid','longdash'))
    g <- g + geom_point(data=alldata, aes(x,y, colour=grp, shape=grp)) +
        geom_line(data=alldata, aes(x,y, colour=grp, linetype=grp))
##        geom_line(data=df, aes(x,y, colour=name), linetype=2)
   if(is.vector(legenda)){ g <- g + theme(plot.title=element_text(size=10, margin = margin(b = -10),hjust=0.1),
                   legend.title=element_blank(),
                   legend.background=element_blank(),
                   legend.justification=legenda,
                   legend.position=legenda,
                   aspect.ratio=0.5
                   #panel.margin.y = unit(-200, "lines")
                  #,plot.margins=unit(c(0,0,0,0),"mm")
                   )} else {g <- g + theme(plot.title=element_text(size=10, margin = margin(b = -10),hjust=0.1),
                   ##legend.title=element_blank(),
                   ##legend.background=element_blank(),
                   ##legend.justification=legenda,
                   legend.position='none',
                   aspect.ratio=0.5)}
        #coord_cartesian() +
        #labs(x=NULL,y=NULL) 
        #scale_x_continuous(expand=c(0,0)) +
        #scale_y_continuous(expand=c(0,0))
    if(is.vector(yrange)){ g <- g + ylim(yrange[1],yrange[2])}
    if(is.vector(xrange)){ g <- g + xlim(xrange[1],xrange[2])}
    if(is.vector(xytext)){ g <- g + labs(x=xytext[1], y=xytext[2], title=xytext[3])}
    ggsave(pdfname, width = 148, height = 148*0.5, units='mm', dpi = 300)
 #   pdf(pdfname)
 #   print(g)
 #   dev.off()
}




## Other plot functions (not used)
logevplot <- function(k=2,samples,pdfname,yrange=NULL){
    x <- 1:(dim(samples[[1]][[k]])[2])
    g <- ggplot()
    for(i in 1:length(samples)){
    df <- data.frame(x, y=colMeans(prob2logev(samples[[i]][[k]])))
    g <- g + geom_line(data=df, aes(x,y), colour=mycolours[i])
    }
    if(is.vector(yrange)){ g <- g + ylim(yrange[1],yrange[2])}
    pdf(pdfname)
    print(g)
    dev.off()
}

meanlogevplot <- function(k=2,ncond,samples,pdfname,yrange=NULL){
    x <- 1:length(samples[1,])
    g <- ggplot()
    for(i in 1:(dim(samples)[1])){
    df <- data.frame(x, y=samples[i,])
    g <- g + geom_point(data=df, aes(x,y), colour=mycolours[i]) +
        geom_line(data=df, aes(x,y), colour=mycolours[i])
    }
    if(k==2){dfchance <- data.frame(x, y=-x*log(ncond))
    g <- g + geom_line(data=dfchance, aes(x,y), linetype=2)}
    if(is.vector(yrange)){ g <- g + ylim(yrange[1],yrange[2])}
    pdf(pdfname)
    print(g)
    dev.off()
}
probmeanplot <- function(k=2,samples,pdfname){
    x <- 1:(dim(samples[[1]][[k]])[2])
    g <- ggplot()
    for(i in 1:length(samples)){
    df <- data.frame(x, y=colMeans(samples[[i]][[k]]))
    g <- g + geom_line(data=df, aes(x,y), colour=mycolours[i])
    }
    pdf(pdfname)
    print(g)
    dev.off()
}
probgmeanplot <- function(k=2,samples,pdfname){
    x <- 1:(dim(samples[[1]][[k]])[2])
    g <- ggplot()
    for(i in 1:length(samples)){
    df <- data.frame(x, y=exp(colMeans(log(samples[[i]][[k]]))))
    g <- g + geom_line(data=df, aes(x,y), colour=mycolours[i])
    }
    pdf(pdfname)
    print(g)
    dev.off()
}
surpriseplot <- function(k=2,samples,pdfname){
    x <- 1:(dim(samples[[1]][[k]])[2])
    g <- ggplot()
    for(i in 1:length(samples)){
    df <- data.frame(x, y=colMeans(-log(samples[[i]][[k]])))
    g <- g + geom_line(data=df, aes(x,y), colour=mycolours[i])
    }
    pdf(pdfname)
    print(g)
    dev.off()
}
hitsdiffplot <- function(samples,pdfname){
    x <- 1:(dim(samples[[1]][[2]])[2])
    g <- ggplot()
    for(i in 2:length(samples)){
        df <- data.frame(x,
                         y=colMeans(prob2hits(samples[[i]])) - colMeans(prob2hits(samples[[1]]))
                         )
    g <- g + geom_line(data=df, aes(x,y), colour=mycolours[i])
    }
    pdf(pdfname)
    print(g)
    dev.off()
}


