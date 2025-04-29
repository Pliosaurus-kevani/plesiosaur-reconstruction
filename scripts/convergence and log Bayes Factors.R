#This script was written by Ruizhe Jackevan Zhao 
#for the research article
#"Body Reconstruction and Size Estimation of Plesiosaurs：
#enlightenment on the ribcage restoration of extinct amniotes in
#2D environments".

#This script aims to check the convergence of parameters
#estimated using BayesTraits
#and compute the log Bayes Factors from marginal likelihoods

#note that the data files used here are uploaded to my Github repository as a .zip package
#you will need to uncompress the file "BayesTraits_strict.zip" 
#before executing this script


#Part of the following code is modified from 
#https://research-information.bris.ac.uk/en/studentTheses/macroevolution-of-early-tetrapods 
#and
#https://doi.org/10.1038/s42003-022-03322-y

#The code was run in R version 4.3.3.

#load packages-----------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(coda) #version 0.19-41; loaded to check convergence of parameters

#check convergence-------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#change the working pathway

setwd()

#create lists to contain results

skip_hom_size_convergence <- list()
hom_size_convergence <- list()
skip_het_size_convergence <- list()
het_size_convergence <- list()

#read data

for(i in 1:100)
{
  #Variable Rates
  skip_het_size_convergence[[i]] <- grep("Tree No", scan(file = paste("tree_size_het_",i,".log.txt.Log.txt", sep=""), what="c", quiet=T, sep="\n", blank.lines.skip=FALSE)) - 1
  het_size_convergence[[i]] = read.table(paste("tree_size_het_",i,".log.txt.Log.txt", sep=""), skip = skip_het_size_convergence[[i]], sep = "\t",  quote="\"", header = TRUE)
  het_size_convergence[[i]] = het_size_convergence[[i]][,-ncol(het_size_convergence[[i]])]
  
  #Homogeneous rates
  skip_hom_size_convergence[[i]] <- grep("Tree No", scan(file = paste("tree_size_hom_",i,".log.txt.Log.txt", sep=""), what="c", quiet=T, sep="\n", blank.lines.skip=FALSE)) - 1
  hom_size_convergence[[i]] = read.table(paste("tree_size_hom_",i,".log.txt.Log.txt", sep=""), skip = skip_hom_size_convergence[[i]], sep = "\t",  quote="\"", header = TRUE)
  hom_size_convergence[[i]] = hom_size_convergence[[i]][,-ncol(hom_size_convergence[[i]])]
  
}

#convert the results to coda format

res_hom_size <- list()
res_het_size <- list()

for(i in 1:100)
{
  res_het_size[[i]] <- mcmc(subset(het_size_convergence[[i]], select=-c(Iteration, Tree.No)),
                            start=min(het_size_convergence[[i]]$Iteration),
                            end=max(het_size_convergence[[i]]$Iteration),thin=20000)
  
  res_hom_size[[i]] <- mcmc(subset(hom_size_convergence[[i]], select=-c(Iteration, Tree.No)),
                            start=min(hom_size_convergence[[i]]$Iteration),
                            end=max(hom_size_convergence[[i]]$Iteration),thin=20000)
}

#get effective size(should be >200)

ess_hom_size_list <- list()
ess_het_size_list <- list()

for(i in 1:100)
{
  ess_hom_size_list[[i]] <- effectiveSize(res_hom_size[[i]])
  ess_het_size_list[[i]] <- effectiveSize(res_het_size[[i]])
}

#find minimum effective size

#Homogeneous rates
tmp_ess_hom_size_min <- do.call(rbind, ess_hom_size_list)
ess_hom_size_min <- min(tmp_ess_hom_size_min)
ess_hom_size_min

#Variable rates
tmp_ess_het_size_min <- do.call(rbind, ess_het_size_list)
ess_het_size_min <- min(tmp_ess_het_size_min)
ess_het_size_min

#calculate log Bayes Factor----------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#Import marginal likelihood from BayesTraits runs

varRates.size.marg.logLik <- list()
homRates.size.marg.logLik <- list()

for(i in 1:100)
{
  varRates.size.marg.logLik[[i]] <- readLines(paste("tree_size_het_",i,".log.txt.Stones.txt", sep=""))
  varRates.size.marg.logLik[[i]] <- tail(varRates.size.marg.logLik[[i]], 1)
  varRates.size.marg.logLik[[i]] <- as.numeric(gsub("[^[:digit:].-]", "", varRates.size.marg.logLik[[i]]))
  varRates.size.marg.logLik[[i]]
}

for(i in 1:100)
{
  homRates.size.marg.logLik[[i]] <- readLines(paste("tree_size_hom_",i,".log.txt.Stones.txt", sep=""))
  homRates.size.marg.logLik[[i]] <- tail(homRates.size.marg.logLik[[i]], 1)
  homRates.size.marg.logLik[[i]] <- as.numeric(gsub("[^[:digit:].-]", "", homRates.size.marg.logLik[[i]]))
  homRates.size.marg.logLik[[i]]
}

#construct a data frame containing log marginal likelihood

marg.logLik.size <- data.frame(TreeNo = seq(1:100), hom = unlist(homRates.size.marg.logLik), het = unlist(varRates.size.marg.logLik))

#calculate log Bayes Factor and import the results into a data frame

logBF.df.size <- data.frame(TreeNo = seq(1:100), hom = 2*(marg.logLik.size$hom-marg.logLik.size$hom), het = 2*(marg.logLik.size$het-marg.logLik.size$hom))

#summary bayes factors
positive.size <- NROW(logBF.df.size$het[logBF.df.size$het >= 2])
strong.size <- NROW(logBF.df.size$het[logBF.df.size$het >= 5 ])
very.strong.size <- NROW(logBF.df.size$het[logBF.df.size$het >= 10 ])
