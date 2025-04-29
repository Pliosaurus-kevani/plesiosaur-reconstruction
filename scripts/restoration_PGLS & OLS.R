#This script was written by Ruizhe Jackevan Zhao 
#for the research article
#"Body Reconstruction and Size Estimation of Plesiosaurs：
#enlightenment on the ribcage restoration of extinct amniotes in
#2D environments".

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
library(nlme) #version 3.1-164; loaded to perform PGLS
library(drc) #version 3.0-1; loaded to perform non-linear regression
library(ggplot2) #version 3.5.1; loaded for visualization

#read the Maximum Clade Credibility Tree (MCCT)--------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#change the working pathway

setwd()

tree <- read.nexus("focal0.25.mcc.tre")

#Skull-Neck-Cervical number regression model-----------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#change the working pathway

setwd()

#read the data

SKL_neck <- read.csv("neck_skull.csv", row.names = 1, header = T)

#I have performed the calculation in Excel.
#The "log.x" in the dataset means log10-transformed cervical number, here regarded as a pseudo-continuous variable
#The "log.y" in the dataset means log10-transformed ratio of (skull/skull+neck)

#remove the sample Kronosaurus_queenslandicus
#since it isn't in the tree

SKL_neck <- SKL_neck[-23,]

#model fitting-----------------------------------------------------------------#

#prune the tree first, since it contains much more taxa than the dataset

t.SKL <- keep.tip(tree, row.names(SKL_neck))

#plot to check

plot(t.SKL, no.margin = T)

#re-order the dataset to match the tips of the pruned tree

SKL_neck <- SKL_neck[match(t.SKL$tip.label, row.names(SKL_neck)),]

#extract the taxon names

spp.SKL <- row.names(SKL_neck)

#define the Brownian motion correlation structure

cor.SKL <- corBrownian(phy = t.SKL, form = ~spp.SKL)

#create the correlation matrix 
#this step is performed since the tree isn't ultrametric

mat.SKL <- vcv(t.SKL, corr = T)

tip.heights.SKL <- diag(mat.SKL)

#fit the PGLS model

SKL.pgls <- gls(log.y ~ log.x, data = SKL_neck, correlation = cor.SKL, weights = varFixed(~tip.heights.SKL))

#fit the OLS model

SKL.ols <- lm(log.y ~ log.x, data = SKL_neck)

#fit the log-logistic model

SKL.log <- drm(log.y ~ log.x, fct = LL.4(), data = SKL_neck)

#compare the AIC values of the two models

SKL.aic <- setNames(c(AIC(SKL.pgls), AIC(SKL.ols), AIC(SKL.log)), c("PGLS", "OLS", "log-logistic"))

#compute the AICc values 
#(note that the three models have different numbers of parameters)

SKL.aicc <- SKL.aic

SKL.aicc[1:2] <- SKL.aic[1:2] + (2*2*3/37)
SKL.aicc[3] <- SKL.aic[3] + (2*4*5/35)

SKL.aicc

#view the results

summary(SKL.pgls)
summary(SKL.ols)
summary(SKL.log)

#calculate the mean per cent prediction error of each model--------------------#

#for the definition of PE, see the article and citations therein

#define an matrix to contain the average absolute per cent prediction error values

PE.SKL <- matrix(ncol = 3, nrow = nrow(SKL_neck))
row.names(PE.SKL) <- row.names(SKL_neck)
colnames(PE.SKL) <- c("PGLS", "OLS", "nonlinear")

#use a for-loop to calculate the mean PE values
#two temporary variables called "dodo" and "moa" are employed here

for(i in 1:nrow(SKL_neck))
{
  #remove one sample and save the rest samples in dodo
  
  dodo <- SKL_neck[-i,]
  
  #prune the MCCT
  
  t.SKL <- keep.tip(tree, row.names(dodo))
  
  #there is no need to order the taxa names
  #since it has been done above
  
  #extract the taxon names
  
  spp.SKL <- row.names(dodo)
  
  #define the Brownian motion correlation structure
  
  cor.SKL <- corBrownian(phy = t.SKL, form = ~spp.SKL)
  
  #create the correlation matrix 
  #this step is performed since the tree isn't ultrametric
  
  mat.SKL <- vcv(t.SKL, corr = T)
  
  tip.heights.SKL <- diag(mat.SKL)
  
  #fit the PGLS model and save the result in moa
  
  moa <- gls(log.y ~ log.x, data = dodo, correlation = cor.SKL, weights = varFixed(~tip.heights.SKL))
  
  PE.SKL[i,1] <- 100*abs(10^(SKL_neck[i,2]) - 10^(predict(moa, data.frame(log.x = SKL_neck[i,1]))))/(10^(predict(moa, data.frame(log.x = SKL_neck[i,1]))))
  
  #fit the OLS model and save the result in moa
  
  moa <- lm(log.y ~ log.x, data = dodo)
  
  PE.SKL[i,2] <- 100*abs(10^(SKL_neck[i,2]) - 10^(predict(moa, data.frame(log.x = SKL_neck[i,1]))))/(10^(predict(moa, data.frame(log.x = SKL_neck[i,1]))))
  
  #fit the nonlinear model and save the result in moa
  
  moa <- drm(log.y ~ log.x, fct = LL.4(), data = dodo)
  
  PE.SKL[i,3] <- 100*abs(10^(SKL_neck[i,2]) - 10^(predict(moa, data.frame(log.x = SKL_neck[i,1]))))/(10^(predict(moa, data.frame(log.x = SKL_neck[i,1]))))
}

#check the mean values

setNames(c(mean(PE.SKL[,1]), mean(PE.SKL[,2]), mean(PE.SKL[,3])), c("PGLS", "OLS", "nonlinear"))


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

plot(t.rib, no.margin = T)

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
  
  #prune the MCCT
  
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

tail <- read.csv("trunk_tail.csv", row.names = 1, header = T)

#log10-transform the data

tail[,1] <- log10(tail[,1])
tail[,2] <- log10(tail[,2])

#model fitting-----------------------------------------------------------------#

#prune the tree first, since it contains much more taxa than the dataset

t.tail <- keep.tip(tree, row.names(tail))

#plot to check

plot(t.tail, no.margin = T)

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

#fit the PGLS model

tail.pgls <- gls(Tail ~ Trunk, data = tail, correlation = cor.tail, weights = varFixed(~tip.heights.tail))

#fit the OLS model

tail.ols <- lm(Tail ~ Trunk, data = tail)

#compare the AIC values of the two models

tail.aic <- setNames(c(AIC(tail.pgls), AIC(tail.ols)), c("PGLS", "OLS"))

tail.aic

#compute the AICc values
#(note that the two models have the same number of parameters[=2] and samples[=19])

tail.aicc <- tail.aic + (2*2*3/16)

#view the results

summary(tail.pgls)
summary(tail.ols)

#calculate the mean per cent prediction error of each model--------------------#

PE.tail <- matrix(ncol = 2, nrow = nrow(tail))
row.names(PE.tail) <- row.names(tail)
colnames(PE.tail) <- c("PGLS", "OLS")

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
  
  #fit the PGLS model and save the result in moa
  
  moa <- gls(Tail ~ Trunk, data = dodo, correlation = cor.tail, weights = varFixed(~tip.heights.tail))
  
  PE.tail[i,1] <- 100*abs(10^(tail[i,2]) - 10^(predict(moa, data.frame(Trunk = rib[i,1]))))/(10^(predict(moa, data.frame(Trunk = rib[i,1]))))
  
  #fit the OLS model and save the result in moa
  
  moa <- lm(Tail ~ Trunk, data = dodo)
  
  PE.tail[i,2] <- 100*abs(10^(tail[i,2]) - 10^(predict(moa, data.frame(Trunk = rib[i,1]))))/(10^(predict(moa, data.frame(Trunk = rib[i,1]))))
}

setNames(c(mean(PE.tail[,1]), mean(PE.tail[,2])), c("PGLS", "OLS"))

#plot the skull-neck regression model------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#manually define the log-logstic function 

eq = function(x){(1.2838)/(1+((x/1.5883)^10.3988))-1.4494}

SKL_plot <- 
  ggplot(data = SKL_neck, aes(x = log.x, y = log.y))+
  geom_point(color = "black", size = 2.5)+
  labs(x = "lg(cervical number)", y = "lg(skull/(skull+neck))")+
  stat_function(fun= eq, color = "#F8766D", lwd = 1)+
  geom_smooth(method = lm, se = F, color = "#619CFF", linetype = "dashed")+
  annotate("text", x = 1.1, y = 0.1, label = "OLS", color = "#619CFF", size = 5)+
  annotate("text", x = 1.1, y = -0.35, label = "log-logistic", color = "#F8766D", size = 5)+
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 15),
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.background = 
          element_rect(fill = "white",
                       color = "black", 
                       linewidth = 1), 
          panel.grid = element_blank())+
  ggtitle("skull-neck regression")

SKL_plot

#change the working pathway before saving

setwd()

ggsave(filename = "Skull-neck scatter.pdf", SKL_plot, width = 9, height = 6)
