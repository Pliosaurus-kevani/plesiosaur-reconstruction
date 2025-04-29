#This script was written by Ruizhe Jackevan Zhao 
#for the research article
#"Body Reconstruction and Size Estimation of Plesiosaurs：
#enlightenment on the ribcage restoration of extinct amniotes in
#2D environments".

#This file reads the post burn-in tree samples from the Bayesian inference
#and randomly select 100 trees
#then save them to a .RData file

#The code was run in R version 4.3.3.

#load packages-----------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(RevGadgets) #version 1.2.1; loaded to read and process tree data

#load and process the tree files-----------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#set the working pathway

setwd()

#read trees (this can take several minutes)

trees <- readTrees("results500000.trees", tree_name = "Tree", burnin = 0.25)

#read the maximum clade credibility tree (MCCT) found by RevBayes

MCCT <- readTrees("focal0.25.mcc.tre")

MCCT.tree <- MCCT[[1]][[1]]@phylo

#select 100 trees randomly

set.seed(1)

random <- sample(1:length(trees[[1]]),100)

tree <- trees[[1]][random]

TREE <- list()

for(i in 1:100)
{
  TREE[[i]] <- tree[[i]]@phylo
}

#remove the trees file to save memory

rm(trees)

#change the working pathway before saving

setwd()

save.image("plesiosaur trees.RData")
