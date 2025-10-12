#This script was written by Ruizhe Jackevan Zhao 
#for the research article
#"Body Reconstruction and Size Estimation of Plesiosaurs"

#This script aims to perform time-calibration to infer the branch length
#of the strict consensus tree.

#The code was run in R version 4.3.3.


#load packages required for data analyses--------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(ape) #version 5.8; loaded for tree manipulation
library(paleotree) #version 3.4.7; loaded for tree time-calibration
library(strap) #version 1.6-1; loaded for visualization


#read data and perform time-calibration----------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#change the working pathway

#setwd()

#read the strict consensus tree

tree <- read.nexus("consensus.nex")

#read the time data

time <- read.table("time.txt", header = T, row.names = 1)

#set seed
#to ensure that the random resolution of polytomies is replicable

set.seed(1)

#perform time-calibration

tree.time <- timePaleoPhy(tree, timeData = time[tree$tip.label,], type= "mbl", vartime = 1, randres = T, plot = F)

#a warning will be returned:
#Do not interpret a single randomly-resolved tree
#this is because the strict consensus tree contains several polytomies,
#which are randomly resolved during time-calibration

#plot

geoscalePhylo(ladderize(tree.time), ages = time, boxes = "Age", units=c("Period", "Epoch", "Age"), width = 1, cex.tip = 0.4, cex.ts = 0.6, tick.scale = "no", label.offset = 0, show.tip.label=T)

#save the tree

write.nexus(tree.time, file = "time_tree.nex", translate = T)

