#This script was written by Ruizhe Jackevan Zhao 
#for the research article
#"Body Reconstruction and Size Estimation of Plesiosaurs"

#This file aims to test the performance of different skeletal measurements
#as the proxy of estimated body volume
#using ordinary least squares (OLS).

#phylogenetic generalized least squares (PGLS) is also employed
#to avoid potential type I errors.

#The code was run in R version 4.3.3.

#load packages-----------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(ape) #version 5.8; loaded to read and process tree
library(nlme) #version 3.1; loaded to perform phylogenetic generalized least squars

#read and partition the dataset------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#set the working pathway

#setwd()

#read the dataset

proxy.data <- read.csv("proxy.csv", header = T, row.names = 1, na.strings = "NA")

#define a list containing the elements

proxy.names <- c("skull_neck", "trunk", "vertebral volume", "vertebral area","humerus length", "humerus width",
           "femur length", "femur width", "coracoid length", "coracoid width",
           "pubis length", "pubis width", "ischium length", "ischium width")

#define a list for the partition performed below
#"combination between proxy and volume" := proxy.comb 

proxy.comb <- list()

#use a for-loop to partition the dataset
#a temporary variable called "dodo" is hired here

for(i in 1:(length(proxy.data)-1))
{
  #extract one proxy and the estimated masses
  
  dodo <- cbind(proxy.data[,i], proxy.data$mass)
  
  #set the row names
  
  row.names(dodo) <- row.names(proxy.data)
  
  #eliminate the rows with missing values
  
  dodo <- na.omit(dodo)
  
  #transport the data to proxy.comb
  
  proxy.comb[[i]] <- dodo
}

#set the names of combinations

names(proxy.comb) <- proxy.names

#perfromed linear regressions--------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#create a list for the regression analyses (ordinary least squares, OLS) below

OLS <- list()


#use a for-loop to perform the OLS

for(i in 1:length(proxy.comb))
{
  OLS[[i]] <- lm(log10(proxy.comb[[i]][,2]) ~ log10(proxy.comb[[i]][,1]))
}

#set the names of combinations

names(OLS) <- proxy.names

#calculate average absolute per cent prediction error (PE) of OLS models-------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#for the definition of PE, see the article and citations therein

#define a list to contain the average absolute per cent prediction error values

PE <- list()

#use a double for-loop to perform the calculation
#another temporary variables called "moa" is hired.

moa <- list()

for(i in 1:length(proxy.comb))
{
  for(j in 1:nrow(proxy.comb[[i]]))
  {
    #eliminate the selected row from the training set
    #then perform the OLS and save the result in dodo
    
    dodo <- lm(log10(V2) ~ log10(V1), data = as.data.frame(proxy.comb[[i]][-j,]))
    
    #perform prediction and save its absolute difference 
    #with body volume obtained from modelling
    #and save the result in moa
    
    moa[[j]] <- 100*abs(proxy.comb[[i]][j,2] - 10^(predict(dodo, data.frame(V1 = proxy.comb[[i]][j,1]))))/(10^(predict(dodo, data.frame(V1 = proxy.comb[[i]][j,1]))))
    
  }
  
  #unlist moa
  
  moa <- unlist(moa)
  
  #so moa is the array containing absolute per cent prediction error
  
  #compute the mean prediction error and save the result
  
  PE[[i]] <- mean(moa)
  
  #clear the memory of moa
  
  moa <- list()
}

#set the names of the PE array

names(PE) <- proxy.names

PE <- unlist(PE)

#summarize the OLS reults into one data frame----------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#define a data frame to contain the results

OLS.export <- data.frame(matrix(ncol = 6, nrow = length(proxy.comb)))

#change the row names and column names of the data frame

row.names(OLS.export) <- proxy.names
colnames(OLS.export) <- c("Slope", "Intercept", "N", "P-value", "R square", "PE")

#use a for-loop to fill OLS.export

for(i in 1:length(OLS))
{
  #slope
  OLS.export[i,1] <- OLS[[i]]$coefficients[2]
  
  #intercept
  OLS.export[i,2] <- OLS[[i]]$coefficients[1]
  
  #sample size (N)
  OLS.export[i,3] <- nrow(proxy.comb[[i]])
  
  #P-value
  if(summary(OLS[[i]])$coefficients[2,4] < 0.001)
  {
    OLS.export[i,4] <- c("< 0.001")
  }
  
  #R square
  OLS.export[i,5] <- summary(OLS[[i]])$r.squared
}

#PE

OLS.export[,6] <- PE

#change the working pathway before saving the results

#setwd()

write.csv(OLS.export, file = "proxy_OLS.csv")

#prune the dataset-------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#We now need to compare the performance of the OLS and PGLS models
#but Styxosaurus SDSM 451 is not present in the tree
#Unfortunately, we have to prune this sample from the dataset

#create a list to contain the data

proxy.pruned <- list()

for(i in 1:(length(proxy.data)-1))
{
  #extract one proxy and the estimated masses
  
  dodo <- cbind(proxy.data[,i], proxy.data$mass)
  
  #set the row names
  
  row.names(dodo) <- row.names(proxy.data)
  
  #eliminate the rows with missing values
  
  dodo <- na.omit(dodo)
  
  #prune Styxosaurus_SDSM_451
  
  dodo <- dodo[row.names(dodo) != "Styxosaurus_SDSM_451",]
  
  #transport the data to proxy.comb
  
  proxy.pruned[[i]] <- dodo
}

#set the names of combinations

names(proxy.pruned) <- proxy.names

#perform phylogenetic generalized least squares (PGLS)-------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#change the working pathway

#setwd()

#read the maximum clade credibility tree (MCCT)

tree <- read.nexus("time_tree.nex")

#plot to check

plot(ladderize(tree), no.margin = T, cex = 0.5)

#define lists to contain the regression models

PGLS <- list()
OLS.p <- list()

#use a for-loop to perform the analyses
#a new temporary variable called "auk" is employed

for(i in 1:length(proxy.pruned))
{
  #prune the tree according to the samples of each combination
  #temporary tree := "tree.t"
  
  tree.t <- keep.tip(tree, row.names(proxy.pruned[[i]]))
  
  #save the species name list in a variable called "spp"
  
  spp <- tree.t$tip.label

  #construct a matrix
  
  dodo <- vcv(tree.t, corr = T)
  
  #use the weights argument since the tree is non-ultrametric
  
  moa <- diag(dodo)
  
  #save the selected combination in auk 
  #and reorder the data according to the pruned tree
  
  auk <- as.data.frame(proxy.pruned[[i]])
  auk <- auk[tree.t$tip.label,]
  
  #fit the models
  
  PGLS[[i]] <- gls(log10(V2) ~ log10(V1), data = auk, correlation = corBrownian(phy = tree.t, form = ~spp), weights = varFixed(~moa))
  OLS.p[[i]] <- lm(log10(V2) ~ log10(V1), data = auk)
}

#define a data frame to contain the PGLS results

PGLS.export <- data.frame(matrix(ncol = 7, nrow = length(proxy.pruned)))

#change the row names and colume names

row.names(PGLS.export) <- proxy.names
colnames(PGLS.export) <- c("Slope", "Intercept", "N", "P-value", "AICc_PGLS", "AICc_OLS", "Delta_AICc")

#use a for-loop to fill PGLS.export

for(i in 1:nrow(PGLS.export))
{
  #slope
  
  PGLS.export[i,1] <- PGLS[[i]]$coefficients[2]
  
  #intercept
  
  PGLS.export[i,2] <- PGLS[[i]]$coefficients[1]
  
  #sample size (N)
  
  PGLS.export[i,3] <- nrow(proxy.pruned[[i]])
  
  #P-value
  
  if(summary(PGLS[[i]])$tTable[2,4] < 0.001)
  {
    PGLS.export[i,4] <- c("< 0.001")
  }
  
  #AICc_PGLS
  
  PGLS.export[i,5] <- AIC(PGLS[[i]]) + 2*2*3/(nrow(proxy.pruned[[i]])-2-1)
  
  #AICc_OLS
  
  PGLS.export[i,6] <- AIC(OLS.p[[i]]) + 2*2*3/(nrow(proxy.pruned[[i]])-2-1)
  
  #Difference in AICc values
  
  PGLS.export[i,7] <- PGLS.export[i,5] - PGLS.export[i,6]
}

#change the working pathway before saving the results

#setwd()

write.csv(PGLS.export, file = "proxy_PGLS.csv")
