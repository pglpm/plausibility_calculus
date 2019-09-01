### sample1.R 
## Author: Porta Mana, Bachmann
## Created: 2017-11-29T12:03:29+0100
## Last-Updated: 2019-08-30T14:00:35+0200
##
## Sampling as described in sect 2.5.1 of paper:
## https://dx.doi.org/10.17605/OSF.IO/R2HUZ

## File 'definitions.R' needs to be in same directory
source('definitions.R')

## For each sample we need:

## 1. Symbols of health conditions considered in the study:
healthconditions <- c('H', 'S')

## 2. List of files with graph values for each health condition. Each of
## these files must have one column for each graph property, one row for
## each patient
datafilenames <- c('weights_con_40cons', 'weights_Schizo_40cons')

## 2b. vector of graph properties to use
properties <- c(1,2)

## 3. Vector of the form c(function1, function2, function3,...), with one
## entry for each graph property: it's the function "l" for the generalized
## normal model, see sect. 2.4.3. The order is the same as the rows of the
## data files for the patients.
##
## logit-normal: logit2
## tangent-normal: tg2
## normal: id2
##
## Below: example with the logit transformation, sect 2.4.4
transformations <- c(tg2) ## defined in definitions.R
    
## 4. Vector of the form c(jacobian1, jacobian2, jacobian3,...): Jacobian
## determinants of the functions above, "l'" in the paper.
## logit-normal: jlogit
## tangent-normal: jtg
## normal: jid
jacobians <- c(jtg)

## 5. Vector of the form c(scale1, scale2, scale3...), with one entry for
## each graph property. These are rescalings of the graph properties, to
## make computation less prone to under- or over-flow. The rescalings are
## to be determined from a different set of patients!
scales <- rep(1,1) ## no scaling of graph quantities

## 6. list of prior coefficients for the normal-inv.Wishart distribution,
## as in sect. 2.4.4: (nu0, kappa0, delta0, Delta0). Alternatively, list of
## lists of prior coefficients, one list for each health condition (so
## different health conditions have different priors)
##
prior <- priorparametersk(1, 0, 42/4) ## for tangent-normal
## prior <- priorparametersk(1, 0, 20) ## for normal
## prior <- priorparameters0 ## for logit-normal

## 7. utility matrix, same form as eqn (37)
utilitymatrix <- matrix(c(3,1,0,4),nrow=2)-2

## 8. pre-test probabilities
pretestprob <- c(1,1)/2

## 9. number of samples
numberofsamples <- 104 #* 5e3

## 10. file where to save results of sampling
savefilename <- 'seqsamples_tg.rds'


## Now we calculate the log-evidence and utility, as in sect. 2.5.1. The function 'logevidence' is defined in definitions.R
print('Calculating graph setup ')

results <- logevidence(numberofsamples,healthconditions,c(1),
                       datafilenames,
                       properties,
                       transformations,
                       scales,
                       jacobians,
                       prior,
                       pretestprob,utilitymatrix,rawmoments=F)

print('saving data...')
saveRDS(results,savefilename)

stop('End of script')
