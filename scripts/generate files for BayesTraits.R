#This script was written by Ruizhe Jackevan Zhao 
#for the research article
#"Body Reconstruction and Size Estimation of Plesiosaurs：
#enlightenment on the ribcage restoration of extinct amniotes in
#2D environments".

#This file aims to produce the files required by BayesTraits
#to investigate the tempo of plesiosaur size evolution

#The code was run in R version 4.3.3.

#Part of the following code is modified from 
#https://research-information.bris.ac.uk/en/studentTheses/macroevolution-of-early-tetrapods 
#and
#https://doi.org/10.1038/s42003-022-03322-y

#load packages-----------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(RevGadgets)  #version 1.2.1; loaded to read and process the trees
library(ape) #version 5.8; loaded to process the trees

#load the trees----------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#change the working pathway

setwd()

#load the tree files

load("plesiosaur trees.RData")

#this .RData file was produced by the R script "load and save trees.R"
#if you don't have the .RData, execute that R script first

#read and process data---------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#change the working pathway

setwd()

#read the body size data

size.data <- read.csv("size evolution.csv", header = T, row.names = 1)

#extract the volumes of osteologically mature specimens
#estimated using model, trunk, or DDV
#"strict" := s

size.s <- size.data[size.data$ontogenetic.status == "mature" & (size.data$source == "model" | size.data$source == "trunk" | size.data$source == "DDV"),]

#manually remove 

#extract the estimated volumes (V)

V <- setNames(size.data[,3], row.names(size.data))
V.s <- setNames(size.s[,3], row.names(size.s))


#log10-transform the data

V <- log10(V)
V.s <- log10(V.s)

#prune the trees since some OTUs in the tree are not present in the data set

TREE.V <- keep.tip.multiPhylo(TREE, row.names(size.data))
TREE.V.s <- keep.tip.multiPhylo(TREE, row.names(size.s))


#prune the maximum clade credibility tree (MCCT)

MCCT.V <- keep.tip(MCCT.tree, row.names(size.data))
MCCT.V.s <- keep.tip(MCCT.tree, row.names(size.s))

#plot to check

plot(ladderize(MCCT.V), cex = 0.5, no.margin = T)
plot(ladderize(MCCT.V.s), cex = 0.5, no.margin = T)
plot(ladderize(TREE.V[[100]]), cex = 0.5, no.margin = T)
plot(ladderize(TREE.V.s[[100]]), cex = 0.5, no.margin = T)

#create BayesTraits command files: 100 trees-----------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#Change the way numbers are stored

scipen_default <- getOption("scipen") #get default setting (should be 0)
options(scipen = 999) #Disable scientific notation in R

#create cmd files
#heterogeneous(het); homogeneous(hom)

Script_size_het.cmd <- list()
Script_size_hom.cmd <- list()

#Use global transformation (Lambda) here
#Use the threaded version (using Cores command to set the number of cores)

#heterogeneous

for(i in 1:100)
{
  Script_size_het.cmd[[i]] <- paste("7\n",
                                    "2\n",
                                    "VarRates\n",
                                    "iterations 22000000\n",
                                    "sample 20000\n",
                                    "burnin 2000000\n",
                                    "Lambda\n",
                                    "stones 100 1000\n",
                                    "Cores 8\n",
                                    "Logfile tree_size_het_",i,".log.txt\n", sep="",
                                    "Info\n",
                                    "Run")
}

#homogeneous

for(i in 1:100)
{
  Script_size_hom.cmd[[i]] <- paste("7\n",
                                    "2\n",
                                    "iterations 22000000\n",
                                    "sample 20000\n",
                                    "burnin 2000000\n",
                                    "Lambda\n",
                                    "stones 100 1000\n",
                                    "Cores 8\n",
                                    "Logfile tree_size_hom_",i,".log.txt\n", sep="",
                                    "Info\n",
                                    "Run")
}

#create BayesTraits command files: MCCT----------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

Script_size_MCCT_het.cmd <- paste("7\n",
                                  "2\n",
                                  "VarRates\n",
                                  "iterations 22000000\n",
                                  "sample 20000\n",
                                  "burnin 2000000\n",
                                  "Lambda\n",
                                  "stones 100 1000\n",
                                  "Cores 8\n",
                                  "Logfile MCCT_size_het.log.txt\n", sep="",
                                  "Info\n",
                                  "Run")

Script_size_MCCT_hom.cmd <- paste("7\n",
                                  "2\n",
                                  "iterations 22000000\n",
                                  "sample 20000\n",
                                  "burnin 2000000\n",
                                  "Lambda\n",
                                  "stones 100 1000\n",
                                  "Cores 8\n",
                                  "Logfile MCCT_size_hom.log.txt\n", sep="",
                                  "Info\n",
                                  "Run")

#export the files: 100 random trees--------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#change the working pathway

setwd()

#save the data

write.table(V.s, "size.txt", sep = "\t", quote = F)

#Note: 
#you will need to manually delete the first row of this .txt file
#which is an "x"
#otherwise BayesTraits.exe won't run.

#save the trees

for(i in 1:100)
{
  write.nexus(TREE.V.s[[i]], file = paste("tree_",i,".nex", sep = ""), translate = T)
}

#save the BayesTraits command files

for(i in 1:100)
{
  write(Script_size_het.cmd[[i]], file=paste("hetRates_size_Script_",i,".cmd", sep=""))
  write(Script_size_hom.cmd[[i]], file=paste("homRates_size_Script_",i,".cmd", sep=""))
}

#export the files: MCCT_all----------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#change the working pathway

setwd()

#save the data

write.table(V, "size.txt", sep = "\t", quote = F)

#Note: 
#you will need to manually delete the first row of this .txt file
#which is an "x"
#otherwise BayesTraits.exe won't run.

#save the pruned MCCT

write.nexus(MCCT.V, file = "MCCT_pruned.nex", translate = T)

#save the BayesTraits command files

write(Script_size_MCCT_het.cmd, file = "Script_size_MCCT_het.cmd")
write(Script_size_MCCT_hom.cmd, file = "Script_size_MCCT_hom.cmd")

#export the files: MCCT_strict-------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#change the working pathway

setwd()

#save the data

write.table(V.s, "size.txt", sep = "\t", quote = F)

#Note: 
#you will need to manually delete the first row of this .txt file
#which is an "x"
#otherwise BayesTraits.exe won't run.

#save the pruned MCCT

write.nexus(MCCT.V.s, file = "MCCT_pruned.nex", translate = T)

#save the BayesTraits command files

write(Script_size_MCCT_het.cmd, file = "Script_size_MCCT_het.cmd")
write(Script_size_MCCT_hom.cmd, file = "Script_size_MCCT_hom.cmd")


#run BayesTraits---------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#100 random trees

#run the following for-loop in cmd

#for /l %i in (1,1,100) do (BayesTraitsV4.exe tree_%i.nex size.txt < hetRates_size_Script_%i.cmd)

#for /l %i in (1,1,100) do (BayesTraitsV4.exe tree_%i.nex size.txt < homRates_size_Script_%i.cmd)

#pruned MCCT

#run the following commands in cmd

#BayesTraitsV4.exe MCCT_pruned.nex size.txt < Script_size_MCCT_het.cmd

#run the Posterior Processor---------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#100 random trees

#run the following code in cmd:

#for /l %i in (1,1,100) do (PPPostProcess.exe tree_size_het_%i.log.txt.VarRates.txt > PP_neck_tree_%i.txt)

#pruned MCCT

#run the following code in cmd:

#PPPostProcess.exe MCCT_size_het.log.txt.VarRates.txt > PP_size_MCCT.txt


