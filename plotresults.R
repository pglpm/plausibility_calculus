### plotresults.R 
## Author: Porta Mana, Bachmann
## Created: 2017-11-29T12:12:43+0100
## Last-Updated: 2018-03-21T18:39:55+0100
##
## Plot graphs with results from sampling, as figure 3 in paper:
## https://dx.doi.org/10.17605/OSF.IO/R2HUZ
##
## Uses output from code similar to the one in 'sample1.R'

## File 'definitions.R' needs to be in same directory
source('definitions.R')

## We need:
## 1. Symbols of health conditions considered in the study:
healthconditions <- c('H', 'S')

## 2. Locations and names of files containing the samples to compare
savefilenames <- c('samples_logit.rds','samples_tg.rds','samples_normal.rds')

## 3. Names of different models that will be compared
modelnames <- c('logit-normal', 'tangent-normal','normal','chance')


## Creating empty lists to hold the results
graphlist <- list()
results <- list()
expsurpriselist <- list()
expcumsurpriselist <- list()
expsurpriseemplist <- list()
expcumsurpriseemplist <- list()
utilitieslist <- list()
cumutilitieslist <- list()
numberofanalyses <- length(savefilenames)

## Read sample data into list "results"
for(i in 1:numberofanalyses){
    print(paste0('Loading graph setup ',i))
    results[[i]] <- readRDS(savefilenames[i])
}

numpatients <- dim(results[[1]][[2]])[2]

## Calculate log-probability ('surprise'), cumulative log-probability,
## utility, cumulative utility, as explained in sect.2.5.1. In the paper
## only the log-probability and utility are compared
##
## The functions 'prob2expsurprise' etc. are defined in 'definitions.R'.
## They calculate the averages.
##
for(i in 1:numberofanalyses){
    expsurpriselist[[i]] <- prob2expsurprise(results[[i]])
    expcumsurpriselist[[i]] <- prob2expcumsurprise(results[[i]])
    expsurpriseemplist[[i]] <- prob2expsurpriseemp(results[[i]])
    expcumsurpriseemplist[[i]] <- prob2expcumsurpriseemp(results[[i]])
    utilitieslist[[i]] <- colMeans(results[[i]][[5]])
    cumutilitieslist[[i]] <- cumsum(colMeans(results[[i]][[5]]))
}

## Add an element to each list of results, for chance results
expsurpriselist[[numberofanalyses+1]] <- rep(-log(2),numpatients)
expcumsurpriselist[[numberofanalyses+1]] <- cumsum(rep(-log(2),numpatients))
expsurpriseemplist[[numberofanalyses+1]] <- rep(-log(2),numpatients)
expcumsurpriseemplist[[numberofanalyses+1]] <- cumsum(rep(-log(2),numpatients))
utilitieslist[[numberofanalyses+1]] <- rep(2.41-2,numpatients)
cumutilitieslist[[numberofanalyses+1]] <- cumsum(rep(2.41-2,numpatients))


## transforms the resulting lists to matrices
expsurprise <- matrix(unlist(expsurpriselist),
                              nrow=numberofanalyses+1,byrow=T)
expcumsurprise <- matrix(unlist(expcumsurpriselist),
                              nrow=numberofanalyses+1,byrow=T)
expsurpriseemp <- matrix(unlist(expsurpriseemplist),
                              nrow=numberofanalyses+1,byrow=T)
expcumsurpriseemp <- matrix(unlist(expcumsurpriseemplist),
                              nrow=numberofanalyses+1,byrow=T)
utilities <- matrix(unlist(utilitieslist),
                              nrow=numberofanalyses+1,byrow=T)
cumutilities <- matrix(unlist(cumutilitieslist),
                              nrow=numberofanalyses+1,byrow=T)

## Plot the evolution of log-probability and utility for all models
##
## Function 'genplot' defined in 'definitions.R'. It needs:
##
## 1. Data to plot, in matrix form. One row for each model.
## 2. Name of file where to save plot
## 3. Names of the models compared
## 4. y-axis range (may be NULL for automatic)
## 5. x-axis range (may be NULL for automatic)
## 6. text for x axis, y axis, and title (may be NULL)
## 7. vector of xy position of legend, ranging from 0 to 1
print('Plotting results...')

genplot(utilities,'utilities.pdf',modelnames,yrange=c(-0.1,1),xrange=c(1,104),xytext=c('patients','utility, asymmetric',''),legenda=NULL)

genplot(cumutilities,'cumutilities.pdf',modelnames,yrange=c(0,81),xrange=c(1,104),xytext=c('patients','cumulative utility',''),legenda=NULL)

genplot(expsurpriseemp,'expectedsurpriseemp.pdf',modelnames,yrange=c(-0.7,-0.25),xrange=c(1,104),xytext=c('patients','log-probability, post-test',''),legenda=c(0.3,0.1))

genplot(expcumsurprise,'expectedcumsurprise.pdf',modelnames,yrange=c(-72,0),xrange=c(1,104),xytext=c('patients','ln[P(health condition)]',''),legenda=NULL)

stop('End of script')
