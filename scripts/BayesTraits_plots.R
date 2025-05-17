#This script was written by Ruizhe Jackevan Zhao 
#for the research article
#"Body Reconstruction and Size Estimation of Plesiosaurs：
#enlightenment on the ribcage restoration of extinct amniotes in
#2D environments".

#This script aims to summarize the results relevant to plesiosaur size evolution
#and creates a plot that is presented in the article

#note that the data files used here are uploaded to my Github repository as two separate .zip packages
#you will need to uncompress the files
#"BayesTraits_MCCT.zip" and "BayesTraits_MCCT_strict.zip" 
#before executing this script

#The code was run in R version 4.3.3.

#load packages-----------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(ape) #version 5.8; loaded to process the trees
library(ggplot2) #version 3.5.1; loaded for visualization
library(ggforce) #version 0.4.2; loaded for visualization
library(ggtree) #version 3.13.1; loaded for visualization
library(viridis) #version 0.6.5; loaded for generating color-blind friendly color bar
library(cowplot) #version 1.1.3; loaded to combine the plots into one figure

#read and process data---------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#change the working pathway

setwd()

#read the maximum clade credibility tree (MCCT)

MCCT.tree <- read.nexus("focal0.25.mcc.tre")

#change the working pathway

setwd()

#read the body size data

size.data <- read.csv("size evolution.csv", header = T, row.names = 1)

#extract the volumes of osteologically mature specimens
#estimated using model, trunk, or DDV
#"strict" := s

size.s <- size.data[size.data$ontogenetic.status == "mature" & (size.data$source == "model" | size.data$source == "trunk" | size.data$source == "DDV"),]

#prune the MCCT since some OTUs in the tree are not present in the data sets

MCCT.V <- keep.tip(MCCT.tree, row.names(size.data))
MCCT.s <- keep.tip(MCCT.tree, row.names(size.s))

#plot to check

plot(ladderize(MCCT.V), cex = 0.5, no.margin = T)
plot(ladderize(MCCT.s), cex = 0.5, no.margin = T)

#summarize the branch-specific evolutionary rate of body size: all-------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#note: you will need to uncompress the file "BayesTraits_MCCT.zip"
#uploaded to my Github repository
#before executing the following lines

#change the working pathway

setwd()

#read the results calculated based on the MCCT
#post process := "postproc"

bayes.postproc <- read.delim('PP_size_MCCT.txt', header = TRUE, check.names = FALSE, sep = "\t", colClasses = c(rep(NA, 22), rep("NULL", 1)))

#convert the spaces in the column names to underscores "_"

colnames(bayes.postproc) <- gsub(" ", "_", colnames(bayes.postproc))

#Split the reported strings/names, so that we can use them with getMRCA

bayes.postproc$Taxa_List <- sapply(strsplit(as.character(bayes.postproc$Taxa_List), ","), "[") 

#create new ID columns

bayes.postproc["NEW_Edge_ID"] <- NA
bayes.postproc["NEW_Node_ID"] <- NA 

for(j in 1:length(bayes.postproc$Taxa_List))
{
  if(length(bayes.postproc$Taxa_List[[j]]) == 1){
    #Single tips: find tip nodes and edges
    bayes.postproc$NEW_Node_ID[[j]]  <- which(MCCT.V$tip.label == bayes.postproc$Taxa_List[[j]]) #find tip nodes
    bayes.postproc$NEW_Edge_ID[[j]] <- which.edge(MCCT.V, bayes.postproc$Taxa_List[[j]]) #find tip edges
  }  else{
    #Multiple tips: find internal nodes and edges of MRCA
    #Notice that the root edge is not assigned an edge ID, i.e. it stays NA
    bayes.postproc$NEW_Node_ID[[j]]  <- getMRCA(MCCT.V, bayes.postproc$Taxa_List[[j]]) #find internal node
    bayes.postproc$NEW_Edge_ID[[j]] <- which.edge(MCCT.V, getMRCA(MCCT.V, bayes.postproc$Taxa_List[[j]])) #find internal edges
  }
}

#re-order the rates according to new node IDs

bayes.postproc.ordered <- bayes.postproc[order(bayes.postproc$NEW_Node_ID),]

#extract size evolutionary rates

size.rates <- log10(bayes.postproc.ordered$Mean_Scalar)

#assign NA to the root

size.rates[[113]] <- NA

#plot: log-transformed evolutionary rates of body size_all---------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#remove the down space symbol from the tip labels

MCCT.V$tip.label <- gsub("_", " ", MCCT.V$tip.label)

#change some typos of the tip labels

MCCT.V$tip.label[[3]] <- '"Aristonectes" quiriquinensis'
MCCT.V$tip.label[[65]] <- "Seeleyosaurus guilelmiimperatori"
MCCT.V$tip.label[[75]] <- '"Monquirasaurus" boyacensis'

rate.plot.all <- 
  ggtree(MCCT.V, aes(color = size.rates), size = 1.2, ladderize = F)+
  viridis::scale_color_viridis(name = expression(paste(log[10],"(rates)")), option = "D", begin = 0.2) +
  geom_tiplab(size = 2.5, color = "black", fontface="italic") +
  #add this line if you need to manually check whether the rates are projected to the correct branches
  #geom_text2(aes(label=(10^(evol_rates))),hjust=-.3,color="red", size = 1.5) +
  theme(legend.position = "none", plot.margin = margin(10, 110, 10, 10),
        legend.title = element_text(size = 15)) +
  coord_cartesian(clip = 'off')

rate.plot.all

#summarize the branch-specific evolutionary rate of body size: strict----------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#note: you will need to uncompress the file "BayesTraits_MCCT_strict.zip"
#uploaded to my Github repository
#before executing the following lines

#change the working pathway

setwd()

#read the results calculated based on the MCCT
#post process := "postproc"

bayes.postproc.s <- read.delim('PP_size_MCCT.txt', header = TRUE, check.names = FALSE, sep = "\t", colClasses = c(rep(NA, 22), rep("NULL", 1)))

#convert the spaces in the column names to underscores "_"

colnames(bayes.postproc.s) <- gsub(" ", "_", colnames(bayes.postproc.s))

#Split the reported strings/names, so that we can use them with getMRCA

bayes.postproc.s$Taxa_List <- sapply(strsplit(as.character(bayes.postproc.s$Taxa_List), ","), "[") 

#create new ID columns

bayes.postproc.s["NEW_Edge_ID"] <- NA
bayes.postproc.s["NEW_Node_ID"] <- NA 

for(j in 1:length(bayes.postproc.s$Taxa_List))
{
  if(length(bayes.postproc.s$Taxa_List[[j]]) == 1){
    #Single tips: find tip nodes and edges
    bayes.postproc.s$NEW_Node_ID[[j]]  <- which(MCCT.s$tip.label == bayes.postproc.s$Taxa_List[[j]]) #find tip nodes
    bayes.postproc.s$NEW_Edge_ID[[j]] <- which.edge(MCCT.s, bayes.postproc.s$Taxa_List[[j]]) #find tip edges
  }  else{
    #Multiple tips: find internal nodes and edges of MRCA
    #Notice that the root edge is not assigned an edge ID, i.e. it stays NA
    bayes.postproc.s$NEW_Node_ID[[j]]  <- getMRCA(MCCT.s, bayes.postproc.s$Taxa_List[[j]]) #find internal node
    bayes.postproc.s$NEW_Edge_ID[[j]] <- which.edge(MCCT.s, getMRCA(MCCT.s, bayes.postproc.s$Taxa_List[[j]])) #find internal edges
  }
}

#re-order the rates according to new node IDs

bayes.postproc.ordered.s <- bayes.postproc.s[order(bayes.postproc.s$NEW_Node_ID),]

#extract size evolutionary rates

size.rates.s <- log10(bayes.postproc.ordered.s$Mean_Scalar)

#assign NA to the root

size.rates.s[[77]] <- NA

#plot: log-transformed evolutionary rates of body size_all---------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#remove the down space symbol from the tip labels

MCCT.s$tip.label <- gsub("_", " ", MCCT.s$tip.label)

#change some typos of the tip labels

MCCT.s$tip.label[[1]] <- '"Aristonectes" quiriquinensis'
MCCT.s$tip.label[[46]] <- "Seeleyosaurus guilelmiimperatori"
MCCT.s$tip.label[[55]] <- '"Monquirasaurus" boyacensis'

rate.plot.strict <- 
  ggtree(MCCT.s, aes(color = size.rates.s), size = 1.2, ladderize = F)+
  viridis::scale_color_viridis(name = expression(paste(log[10],"(rates)")), option = "D", begin = 0.2) +
  geom_tiplab(size = 2.5, color = "black", fontface="italic") +
  #add this line if you need to manually check whether the rates are projected to the correct branches
  #geom_text2(aes(label=(10^(evol_rates))),hjust=-.3,color="red", size = 1.5) +
  theme(legend.position=c(.1, .2), plot.margin = margin(10, 90, 10, 10),
        legend.title = element_text(size = 10)) +
  coord_cartesian(clip = 'off')

rate.plot.strict

#pie chart showing the proportion of size proxies------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#create a data frame containing all types of volumetric sources
#and the corresponding number of taxa

sources <- c("DDV", "femur chord", "ischium width", "model", "skull-neck", "trunk")
N <- as.vector(table(size.data$source))

pie.data <- data.frame(source = sources, number = N)

#adjust the order of sources

pie.data$source <- factor(pie.data$source, levels = c("model", "trunk", "DDV", "ischium width", "femur chord", "skull-neck"))

#plot

pie.chart <- 
  ggplot(data = pie.data, aes(x = 2, y = number, fill = source))+
  geom_bar(stat = "identity", position = "stack", color = "white")+
  coord_polar("y")+
  geom_text(aes(label = source), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  theme_void()+
  theme(legend.position = "none", plot.margin = margin(10,10,10,10))+
  xlim(1,2.5)

pie.chart  

#combine the plots in a single figure------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

dodo <- plot_grid(pie.chart, rate.plot.strict, nrow = 2, labels = c("B", "C"), align = "hv", rel_heights = c(0.5,1))

figure <- plot_grid(rate.plot.all, dodo, nrow = 1, labels = c("A", ""), align = "h")

figure

#change the working pathway before saving

setwd()

ggsave(filename = "Figure5.pdf", figure, width = 13, height = 13)
