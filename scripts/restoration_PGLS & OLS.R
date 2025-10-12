#This script was written by Ruizhe Jackevan Zhao 
#for the research article
#Body Reconstruction and Size Estimation of Plesiosaurs

#This script aims to fit regression models
#to restore the missing puzzles in plesiosaur reconstruction

#The models include:
#ordinary least squares (OLS)
#phylogenetic generalized least squares (PGLS
#non-linear regression

#The code was run in R version 4.3.3.

#load packages required for data analyses--------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(ape) #version 5.8; loaded to read and manipulate tree
library(geiger) #version 2.0.11; loaded to perform name check
library(nlme) #version 3.1-164; loaded to perform PGLS
library(drc) #version 3.0-1; loaded to perform non-linear regression
library(ggplot2) #version 3.5.1; loaded for visualization

#read the time-calibrated strict consensus tree--------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#change the working pathway

setwd("C:/Users/28299/Desktop/my papers/Body Reconstruction and Size Estimation of Plesiosaurs/Phylogeny/most parsimony_v3")

tree <- read.nexus("time_tree.nex")

#Skull-Neck-Cervical number regression model-----------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#change the working pathway

setwd("C:/Users/28299/Desktop/my papers/Body Reconstruction and Size Estimation of Plesiosaurs/Regression")

#read the data

SKL_neck <- read.csv("neck_skull.csv", row.names = 1, header = T)

#I have performed the calculation in Excel.
#The "log.x" in the dataset means log10-transformed cervical number, here regarded as a pseudo-continuous variable
#The "log.y" in the dataset means log10-transformed ratio of (skull/skull+neck)

#3 taxa in the skull-neck dataset are absent from the current phylogeny

#find them

p.names <- name.check(tree, SKL_neck)

#prune the tree

t.SKL <- drop.tip(tree, c(p.names$data_not_tree, p.names$tree_not_data))

#remove the 3 taxa from the dataset

SKL_neck.p <- SKL_neck[t.SKL$tip.label,]

#model fitting: pruned---------------------------------------------------------#

#re-order the dataset to match the tips of the pruned tree

SKL_neck.p <- SKL_neck.p[match(t.SKL$tip.label, row.names(SKL_neck.p)),]

#extract the taxon names

spp.SKL <- row.names(SKL_neck.p)

#define the Brownian motion correlation structure

cor.SKL <- corBrownian(phy = t.SKL, form = ~spp.SKL)

#create the correlation matrix 
#this step is performed since the tree isn't ultrametric

mat.SKL <- vcv(t.SKL, corr = T)

tip.heights.SKL <- diag(mat.SKL)

#fit the PGLS model

SKL.pgls.p <- gls(log.y ~ log.x, data = SKL_neck.p, correlation = cor.SKL, weights = varFixed(~tip.heights.SKL))

#fit the OLS model

SKL.ols.p <- lm(log.y ~ log.x, data = SKL_neck.p)

#fit the log-logistic model

SKL.log.p <- drm(log.y ~ log.x, fct = LL.4(), data = SKL_neck.p)

#compare the AIC values of the two models

SKL.aic.p <- setNames(c(AIC(SKL.pgls.p), AIC(SKL.ols.p), AIC(SKL.log.p)), c("PGLS_pruned", "OLS_pruned", "log-logistic_pruned"))

#compute the AICc values 
#(note that the three models have different numbers of parameters)

SKL.aicc.p <- SKL.aic.p

SKL.aicc.p[1:2] <- SKL.aicc.p[1:2] + (2*2*3/34)
SKL.aicc.p[3] <- SKL.aicc.p[3] + (2*4*5/32)

SKL.aicc.p

#view the results

summary(SKL.pgls.p)
summary(SKL.ols.p)
summary(SKL.log.p)

#calculate the mean per cent prediction error of each model: pruned dataset----#

#for the definition of PE, see the article and citations therein

#define an matrix to contain the average absolute per cent prediction error values

PE.SKL.p <- matrix(ncol = 3, nrow = nrow(SKL_neck.p))
row.names(PE.SKL.p) <- row.names(SKL_neck.p)
colnames(PE.SKL.p) <- c("PGLS", "OLS", "nonlinear")

#use a for-loop to calculate the mean PE values
#two temporary variables called "dodo" and "moa" are employed here

for(i in 1:nrow(SKL_neck.p))
{
  #remove one sample and save the rest samples in dodo
  
  dodo <- SKL_neck.p[-i,]
  
  #prune the tree
  
  t.SKL_tempo <- keep.tip(t.SKL, row.names(dodo))
  
  #there is no need to order the taxa names
  #since it has been done above
  
  #extract the taxon names
  
  spp.SKL <- row.names(dodo)
  
  #define the Brownian motion correlation structure
  
  cor.SKL <- corBrownian(phy = t.SKL_tempo, form = ~spp.SKL)
  
  #create the correlation matrix 
  #this step is performed since the tree isn't ultrametric
  
  mat.SKL <- vcv(t.SKL_tempo, corr = T)
  
  tip.heights.SKL <- diag(mat.SKL)
  
  #fit the PGLS model and save the result in moa
  
  moa <- gls(log.y ~ log.x, data = dodo, correlation = cor.SKL, weights = varFixed(~tip.heights.SKL))
  
  PE.SKL.p[i,1] <- 100*abs(10^(SKL_neck.p[i,2]) - 10^(predict(moa, data.frame(log.x = SKL_neck.p[i,1]))))/(10^(predict(moa, data.frame(log.x = SKL_neck.p[i,1]))))
  
  #fit the OLS model and save the result in moa
  
  moa <- lm(log.y ~ log.x, data = dodo)
  
  PE.SKL.p[i,2] <- 100*abs(10^(SKL_neck.p[i,2]) - 10^(predict(moa, data.frame(log.x = SKL_neck.p[i,1]))))/(10^(predict(moa, data.frame(log.x = SKL_neck.p[i,1]))))
  
  #fit the nonlinear model and save the result in moa
  
  moa <- drm(log.y ~ log.x, fct = LL.4(), data = dodo)
  
  PE.SKL.p[i,3] <- 100*abs(10^(SKL_neck.p[i,2]) - 10^(predict(moa, data.frame(log.x = SKL_neck.p[i,1]))))/(10^(predict(moa, data.frame(log.x = SKL_neck.p[i,1]))))
}

#check the mean values

setNames(c(mean(PE.SKL.p[,1]), mean(PE.SKL.p[,2]), mean(PE.SKL.p[,3])), c("PGLS_pruned", "OLS_pruned", "log-logistic_pruned"))

#model fitting: whole----------------------------------------------------------#

#fit the OLS model

SKL.ols <- lm(log.y ~ log.x, data = SKL_neck)

#fit the log-logistic model

SKL.log <- drm(log.y ~ log.x, fct = LL.4(), data = SKL_neck)

#compare the AIC values of the two models

SKL.aic <- setNames(c(AIC(SKL.ols), AIC(SKL.log)), c("OLS", "log-logistic"))

#compute the AICc values 
#(note that the two models have different numbers of parameters)

SKL.aicc <- SKL.aic

SKL.aicc[1] <- SKL.aicc[1] + (2*2*3/37)
SKL.aicc[2] <- SKL.aicc[2] + (2*4*5/35)

SKL.aicc

#view the results

summary(SKL.ols)
summary(SKL.log)

#calculate the mean per cent prediction error of each model: whole dataset-----#

#for the definition of PE, see the article and citations therein

#define an matrix to contain the average absolute per cent prediction error values

PE.SKL <- matrix(ncol = 2, nrow = nrow(SKL_neck))
row.names(PE.SKL) <- row.names(SKL_neck)
colnames(PE.SKL) <- c("OLS", "nonlinear")

#use a for-loop to calculate the mean PE values
#two temporary variables called "dodo" and "moa" are employed here

for(i in 1:nrow(SKL_neck))
{
  #remove one sample and save the rest samples in dodo
  
  dodo <- SKL_neck[-i,]
  
  #fit the OLS model and save the result in moa
  
  moa <- lm(log.y ~ log.x, data = dodo)
  
  PE.SKL[i,1] <- 100*abs(10^(SKL_neck[i,2]) - 10^(predict(moa, data.frame(log.x = SKL_neck[i,1]))))/(10^(predict(moa, data.frame(log.x = SKL_neck[i,1]))))
  
  #fit the nonlinear model and save the result in moa
  
  moa <- drm(log.y ~ log.x, fct = LL.4(), data = dodo)
  
  PE.SKL[i,2] <- 100*abs(10^(SKL_neck[i,2]) - 10^(predict(moa, data.frame(log.x = SKL_neck[i,1]))))/(10^(predict(moa, data.frame(log.x = SKL_neck[i,1]))))
}

#check the mean values

setNames(c(mean(PE.SKL[,1]), mean(PE.SKL[,2])), c("OLS", "log-logistic"))

#trunk-rib regression----------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#read the data

rib <- read.csv("trunk_rib.csv", row.names = 1, header = T)

#log10-transform the data

rib[,1] <- log10(rib[,1])
rib[,2] <- log10(rib[,2])

#model fitting-----------------------------------------------------------------#

#prune the tree first, since it contains much more taxa than the dataset

t.rib <- keep.tip(tree, row.names(rib))

#plot to check

plot(ladderize(t.rib), no.margin = T)

#re-order the dataset to match the tips of the pruned tree

rib <- rib[match(t.rib$tip.label, row.names(rib)),]

#extract the taxon names

spp.rib <- row.names(rib)

#define the Brownian autocorrelation structure

cor.rib <- corBrownian(phy = t.rib, form = ~spp.rib)

#create the correlation matrix 
#this step is performed since the tree isn't ultrametric

mat.rib <- vcv(t.rib, corr = T)

tip.heights.rib <- diag(mat.rib)

#fit the PGLS model

rib.pgls <- gls(Rib ~ Trunk, data = rib, correlation = cor.rib, weights = varFixed(~tip.heights.rib))

#fit the OLS model

rib.ols <- lm(Rib ~ Trunk, data = rib)

#compare the AIC values of the two models

rib.aic <- setNames(c(AIC(rib.pgls), AIC(rib.ols)), c("PGLS", "OLS"))

rib.aic

#calculate the AICc values 
#(note that the two models have the same number of parameters[=2] and samples[=24])

rib.aicc <- rib.aic+(2*2*3/21)

#view the results

summary(rib.pgls)
summary(rib.ols)

#calculate the mean per cent prediction error of each model--------------------#

#define an matrix to contain the average absolute per cent prediction error values

PE.rib <- matrix(ncol = 2, nrow = nrow(rib))
row.names(PE.rib) <- row.names(rib)
colnames(PE.rib) <- c("PGLS", "OLS")

for(i in 1:nrow(rib))
{
  #remove one sample and save the rest samples in dodo
  
  dodo <- rib[-i,]
  
  #prune the tree
  
  t.rib <- keep.tip(tree, row.names(dodo))
  
  #extract the taxon names
  
  spp.rib <- row.names(dodo)
  
  #define the Brownian autocorrelation structure
  
  cor.rib <- corBrownian(phy = t.rib, form = ~spp.rib)
  
  #create the correlation matrix 
  #this step is performed since the tree isn't ultrametric
  
  mat.rib <- vcv(t.rib, corr = T)
  
  tip.heights.rib <- diag(mat.rib)
  
  #fit the PGLS model and save the result in moa
  
  moa <- gls(Rib ~ Trunk, data = dodo, correlation = cor.rib, weights = varFixed(~tip.heights.rib))
  
  PE.rib[i,1] <- 100*abs(10^(rib[i,2]) - 10^(predict(moa, data.frame(Trunk = rib[i,1]))))/(10^(predict(moa, data.frame(Trunk = rib[i,1]))))
  
  #fit the OLS model and save the result in moa
  
  moa <- lm(Rib ~ Trunk, data = dodo)
  
  PE.rib[i,2] <- 100*abs(10^(rib[i,2]) - 10^(predict(moa, data.frame(Trunk = rib[i,1]))))/(10^(predict(moa, data.frame(Trunk = rib[i,1]))))
}

#check the mean values

setNames(c(mean(PE.rib[,1]), mean(PE.rib[,2])), c("PGLS", "OLS"))

#trunk-tail regression---------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#read the dataset

tail <- read.csv("tail.csv", row.names = 1, header = T)

#log10-transform the data

for(i in 1:4)
{
  tail[,i] <- log10(tail[,i])
}


#model fitting-----------------------------------------------------------------#

#prune the tree first, since it contains much more taxa than the dataset

t.tail <- keep.tip(tree, row.names(tail))

#plot to check

plot(ladderize(t.tail), no.margin = T)

#re-order the dataset to match the tips of the pruned tree

tail <- tail[match(t.tail$tip.label, row.names(tail)),]

#extract the taxon names

spp.tail <- row.names(tail)

#define the Brownian autocorrelation structure

cor.tail <- corBrownian(phy = t.tail, form = ~spp.tail)

#create the correlation matrix 
#this step is performed since the tree isn't ultrametric

mat.tail <- vcv(t.tail, corr = T)

tip.heights.tail <- diag(mat.tail)

#fit the PGLS models

tail_FL.pgls <- gls(Tail ~ FL, data = tail, correlation = cor.tail, weights = varFixed(~tip.heights.tail))
tail_FW.pgls <- gls(Tail ~ FW, data = tail, correlation = cor.tail, weights = varFixed(~tip.heights.tail))
tail_Trunk.pgls <- gls(Tail ~ Trunk, data = tail, correlation = cor.tail, weights = varFixed(~tip.heights.tail))

#fit the OLS model

tail_FL.ols <- lm(Tail ~ FL, data = tail)
tail_FW.ols <- lm(Tail ~ FW, data = tail)
tail_Trunk.ols <- lm(Tail ~ Trunk, data = tail)


#compute the AICc values of the models

tail.aic <- setNames(c(AIC(tail_FL.pgls), AIC(tail_FL.ols), AIC(tail_FW.pgls), AIC(tail_FW.ols), AIC(tail_Trunk.pgls), AIC(tail_Trunk.ols)), c("FL_PGLS", "FL_OLS", "FW_PGLS", "FW_OLS", "Trunk_PGLS", "Trunk_OLS"))

tail.aic

#compute the AICc values
#(note that the models have the same number of parameters[=2] and samples[=19])

tail.aicc <- tail.aic + (2*2*3/19)

#view the results

summary(tail_FL.pgls)
summary(tail_FL.ols)
summary(tail_FW.pgls)
summary(tail_FW.ols)
summary(tail_Trunk.pgls)
summary(tail_Trunk.ols)

#calculate the mean per cent prediction error of each model--------------------#

PE.tail <- matrix(ncol = 6, nrow = nrow(tail))
row.names(PE.tail) <- row.names(tail)
colnames(PE.tail) <- c("FL_PGLS", "FL_OLS", "FW_PGLS", "FW_OLS", "Trunk_PGLS", "Trunk_OLS")

#use a for-loop to perform the calculation 

for(i in 1:nrow(tail))
{
  #remove one sample from the data set and save the rest in dodo
  
  dodo <- tail[-i,]
  
  #prune the MCCT
  
  t.tail <- keep.tip(tree, row.names(dodo))
  
  #extract the taxon names
  
  spp.tail <- row.names(dodo)
  
  #define the Brownian autocorrelation structure
  
  cor.tail <- corBrownian(phy = t.tail, form = ~spp.tail)
  
  #create the correlation matrix 
  #this step is performed since the tree isn't ultrametric
  
  mat.tail <- vcv(t.tail, corr = T)
  
  tip.heights.tail <- diag(mat.tail)
  
  #fit the PGLS model and save the results
  
  moa <- gls(Tail ~ FL, data = dodo, correlation = cor.tail, weights = varFixed(~tip.heights.tail))
  
  PE.tail[i,1] <- 100*abs(10^(tail[i,4]) - 10^(predict(moa, data.frame(FL = tail[i,1]))))/(10^(predict(moa, data.frame(FL = tail[i,1]))))
  
  moa <- gls(Tail ~ FW, data = dodo, correlation = cor.tail, weights = varFixed(~tip.heights.tail))
  
  PE.tail[i,3] <- 100*abs(10^(tail[i,4]) - 10^(predict(moa, data.frame(FW = tail[i,2]))))/(10^(predict(moa, data.frame(FW = tail[i,2]))))
  
  moa <- gls(Tail ~ Trunk, data = dodo, correlation = cor.tail, weights = varFixed(~tip.heights.tail))
  
  PE.tail[i,5] <- 100*abs(10^(tail[i,4]) - 10^(predict(moa, data.frame(Trunk = tail[i,3]))))/(10^(predict(moa, data.frame(Trunk = tail[i,3]))))
  
  #fit the OLS models and save the results
  
  moa <- lm(Tail ~ FL, data = dodo)
  
  PE.tail[i,2] <- 100*abs(10^(tail[i,4]) - 10^(predict(moa, data.frame(FL = tail[i,1]))))/(10^(predict(moa, data.frame(FL = tail[i,1]))))
  
  moa <- lm(Tail ~ FW, data = dodo)
  
  PE.tail[i,4] <- 100*abs(10^(tail[i,4]) - 10^(predict(moa, data.frame(FW = tail[i,2]))))/(10^(predict(moa, data.frame(FW = tail[i,2]))))
  
  moa <- lm(Tail ~ Trunk, data = dodo)
  
  PE.tail[i,6] <- 100*abs(10^(tail[i,4]) - 10^(predict(moa, data.frame(Trunk = tail[i,3]))))/(10^(predict(moa, data.frame(Trunk = tail[i,3]))))
}

setNames(c(mean(PE.tail[,1]), mean(PE.tail[,2]), mean(PE.tail[,3]), mean(PE.tail[,4]), mean(PE.tail[,5]), mean(PE.tail[,6])), c("FL_PGLS", "FL_OLS", "FW_PGLS", "FW_OLS", "Trunk_PGLS", "Trunk_OLS"))

#plot the skull-neck regression model: pruned dataset--------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#manually define the log-logstic function and PGLS function

eq_LL = function(x){(1.2725)/(1+((x/1.5819)^10.3058))-1.4344}
eq_PGLS = function(x){-1.4661*x+1.5406}

SKL_plot <- 
  ggplot(data = SKL_neck.p, aes(x = log.x, y = log.y))+
  geom_point(color = "black", size = 2.5)+
  labs(x = "lg(cervical number)", y = "lg(skull/(skull+neck))")+
  stat_function(fun= eq_LL, aes(color = "log-logistic"), lwd = 1)+
  stat_function(fun= eq_PGLS, aes(color = "PGLS"), lwd = 1, linetype = "dotted")+
  geom_smooth(method = lm, se = F, aes(color = "OLS"), linetype = "dashed")+
  scale_color_manual("models", values = c("#F8766D", "#82BD4D", "#619CFF"))+
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
             axis.title = element_text(size = 15),
             axis.text = element_blank(), 
             axis.ticks = element_blank(),
             legend.text = element_text(size = 10),
             legend.title = element_text(size = 10, face = "bold"))+
  ggtitle("skull-neck regression")
  
print(SKL_plot)

#change the working pathway before saving 

setwd("C:/Users/28299/Desktop/my papers/Body Reconstruction and Size Estimation of Plesiosaurs/Figures")

ggsave("skull-neck.pdf", SKL_plot, width = 8, height = 6)
