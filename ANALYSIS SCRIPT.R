
# Clear workspace ---------------------------------------------------------

rm(list=ls())


# Load R libraries --------------------------------------------------------

library(eegkit)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(gdata)
library(psych)
library(polycor)
library(gtools)
library(calibrate)
library(stringr)
library(nlme)
library(irr)
library(gdata)
library(nparcomp)
library(lawstat)
library(car)
library(lcmm)
library(reshape2)
library(readxl)
library(mgcv)
library(mgcViz)
library(itsadug)
library(mgm)
library(multcomp)
library(emmeans)
library(bootnet)
library(NetworkComparisonTest)
library(igraph)
library(network)
library(qgraph)
library(networktools)
library(intergraph)
library(GGally)
library(grid)
library(gridExtra)
library(geomtextpath)
library(patchwork)
library(scales)
library(eegUtils)

# Load datasets -----------------------------------------------------------

setwd("~/CUSSB/CHIARA EEG/FINAL")
load("~/CUSSB/CHIARA EEG/FINAL/ALPHA_data.RData")

variabili=read_xlsx("esempio dati4rev.xlsx",sheet =1)
n.elec = length(Electrodes)
n.obs = nrow(ALPHA_data)

# Preprocessing -----------------------------------------------------------

Electrodes <- sapply(sapply(lapply(names(variabili), function(x) strsplit(x, "\\...")), function(y) y[1]), function(z) z[1])[2:26]
Variables_to_Model = Electrodes
Variables_to_Model_lag1 =paste(Variables_to_Model, ".lag1.",sep = "")
Variables_to_Model_lag2 =paste(Variables_to_Model, ".lag2.",sep = "")
data_to_modelALPHA = data.frame(do.call(rbind,lapply(ALPHA_MEAN_Complete,function(x) x$sec0.5)))
data_to_modelALPHA = merge(ALPHA_MEAN_BAS,data_to_modelALPHA,by="ID")
data_to_modelALPHA[is.na(data_to_modelALPHA)] = NA
datA = data.matrix(data_to_modelALPHA[,2:301])
datA[datA==0] <- 0.0001
data_to_modelALPHA[,2:301] = log(datA)
data_to_modelALPHA[,sapply(data_to_modelALPHA,class)=="list"] <- sapply(data_to_modelALPHA[,sapply(data_to_modelALPHA,class)=="list"],as.numeric)
indFirst = which(data_to_modelALPHA$Video=="FirstT")
data_to_modelALPHA = subset(data_to_modelALPHA,Video=="FirstT")
ALPHA_data = data_to_modelALPHA

save(ALPHA_data,file = "ALPHA_data.RData")

# Fitting the first level model to the alpha band data --------------------

set.seed(123)
{
  datA =ALPHA_data
  datA$Previous = factor(datA$Previous)
  datA$VideoO1 = factor(datA$VideoO1)
  datA$VideoVal = factor(datA$VideoVal)
  datA$Gender = factor(datA$Gender)
  datA$ID =  factor(datA$ID,levels = unique(datA$ID))
  rownames(datA) = NULL
  NOMIND_MODELA = vector("list",25)
  for (i in 1:25){
    cat("\n")
    cat(paste("Fitting variable :",Electrodes[i]))
    cat("\n")
    ARY1 =  paste(Electrodes[i],".lag1.",sep = "")
    BAS =  paste(Electrodes[i],"_BAS",sep = "")
    modelAi=as.formula(paste(Electrodes[i],"~VideoVal+",ARY1,"+",BAS," + s(Time) + 
    TAS.TOT + STAI.STATO.TOT +STAI.TRATTO.TOT+Age+Gender+s(ID, bs='re')",sep = ""))
    gamA1 = try(gam(modelAi, na.action = na.omit,  data=datA))
    NOMIND_MODELA[[i]] <- gamA1
  }
  FITTED_GAMS_MEAN_0_5_LOGB = NOMIND_MODELA
}

#save(FITTED_GAMS_MEAN_0_5_LOGB,file = "FITTED_GAMS_MEAN_0_5_LOGB.RData")


# Get the residuals from the first level model ----------------------------

Alpha_residuals = matrix(NA, nrow = n.obs, ncol=n.elec,dimnames = list(NULL,Electrodes))
Alpha_band = FITTED_GAMS_MEAN_0_5_LOGB
indomA = lapply(1:n.elec,function(i){
  ARY1 =  paste(Electrodes[i],".lag1.",sep = "")
  BAS =  paste(Electrodes[i],"_BAS",sep = "")
  as.numeric(rownames(gam(as.formula(paste(Electrodes[i],"~VideoVal+",ARY1,"+",BAS," + s(Time) + 
    TAS.TOT + STAI.STATO.TOT +STAI.TRATTO.TOT+Age+Gender+s(ID, bs='re')",sep = "")), na.action = na.omit,  data=datA, fit = FALSE)$mf))
})
aRES=data.matrix(sapply(Alpha_band, function(x) resid_gam(x,incl_na = TRUE)))
dimnames(aRES)[[2]] = Electrodes
for (i in 1:n.elec) Alpha_residuals[indomA[[i]],i] = aRES[,i]
Alpha_residualsB=Alpha_residuals[,c(4,5,7,6,2,3,1,8:16,21,20,19,18,17,22,24,25,23)]
Alpha_residualsBVID = split(data.frame(Alpha_residualsB), list(datA$VideoVal))



# Network Inference and Comparison ----------------------------------------

#### Set the seed for reproducibility
set.seed(123)

#### Compare campaigns 1-2
COMPARE_12_ALPHAP_MEANb = NCT(na.omit(Alpha_residualsBVID[[1]]),
                              na.omit(Alpha_residualsBVID[[2]]), it=5000,
                              test.edges = TRUE, test.centrality = TRUE,centrality = c("closeness", "betweenness", "strength", "expectedInfluence"),
                              p.adjust.methods = "bonferroni")


#### Compare campaigns 1-3

COMPARE_13_ALPHAP_MEANb = NCT(na.omit(Alpha_residualsBVID[[1]]),
                              na.omit(Alpha_residualsBVID[[3]]), it=5000,
                              test.edges = TRUE, test.centrality = TRUE,centrality = c("closeness", "betweenness", "strength", "expectedInfluence"),
                              p.adjust.methods = "bonferroni")

#### Compare campaigns 1-4
COMPARE_14_ALPHAP_MEANb = NCT(na.omit(Alpha_residualsBVID[[1]]),
                              na.omit(Alpha_residualsBVID[[4]]), it=5000,
                              test.edges = TRUE, test.centrality = TRUE,centrality = c("closeness", "betweenness", "strength", "expectedInfluence"),
                              p.adjust.methods = "bonferroni")


#### Compare campaigns 2-3
COMPARE_23_ALPHAP_MEANb = NCT(na.omit(Alpha_residualsBVID[[2]]),
                              na.omit(Alpha_residualsBVID[[3]]), it=5000,
                              test.edges = TRUE, test.centrality = TRUE,centrality = c("closeness", "betweenness", "strength", "expectedInfluence"),
                              p.adjust.methods = "bonferroni")


#### Compare campaigns 2-4
COMPARE_24_ALPHAP_MEANb = NCT(na.omit(Alpha_residualsBVID[[2]]),
                              na.omit(Alpha_residualsBVID[[4]]), it=5000,
                              test.edges = TRUE, test.centrality = TRUE,centrality = c("closeness", "betweenness", "strength", "expectedInfluence"),
                              p.adjust.methods = "bonferroni")


#### Compare campaigns 3-4
COMPARE_34_ALPHAP_MEANb = NCT(na.omit(Alpha_residualsBVID[[3]]),
                              na.omit(Alpha_residualsBVID[[4]]), it=5000,
                              test.edges = TRUE, test.centrality = TRUE,centrality = c("closeness", "betweenness", "strength", "expectedInfluence"),
                              p.adjust.methods = "bonferroni")


# save(COMPARE_12_ALPHAP_MEANb, file = "COMPARE_12_ALPHAP_MEANb.Rdata")
# save(COMPARE_13_ALPHAP_MEANb, file = "COMPARE_13_ALPHAP_MEANb.Rdata")
# save(COMPARE_14_ALPHAP_MEANb, file = "COMPARE_14_ALPHAP_MEANb.Rdata")
# save(COMPARE_23_ALPHAP_MEANb, file = "COMPARE_23_ALPHAP_MEANb.Rdata")
# save(COMPARE_24_ALPHAP_MEANb, file = "COMPARE_24_ALPHAP_MEANb.Rdata")
# save(COMPARE_34_ALPHAP_MEANb, file = "COMPARE_34_ALPHAP_MEANb.Rdata")




# Results -----------------------------------------------------------------



## Get the adjacency matrices of the networks
Alpha1 = COMPARE_12_ALPHAP_MEANb$nw1
Alpha2 = COMPARE_12_ALPHAP_MEANb$nw2
Alpha3 = COMPARE_13_ALPHAP_MEANb$nw2
Alpha4 = COMPARE_14_ALPHAP_MEANb$nw2
electrodes = names(Alpha_residualsBVID[[1]])

dimnames(Alpha1) = dimnames(Alpha2) = dimnames(Alpha3) = dimnames(Alpha4) = list(electrodes,electrodes)

## Create igraph objects
ALPHAVideo1net_ig <- graph_from_adjacency_matrix(Alpha1,
                                                 mode = "undirected",
                                                 weighted = TRUE,
                                                 diag = FALSE)

ALPHAVideo2net_ig <- graph_from_adjacency_matrix(Alpha2,
                                                 mode = "undirected",
                                                 weighted = TRUE,
                                                 diag = FALSE)
ALPHAVideo3net_ig <- graph_from_adjacency_matrix(Alpha3,
                                                 mode = "undirected",
                                                 weighted = TRUE,
                                                 diag = FALSE)
ALPHAVideo4net_ig <- graph_from_adjacency_matrix(Alpha4,
                                                 mode = "undirected",
                                                 weighted = TRUE,
                                                 diag = FALSE)

### Calculate centralities
cent_matA1 =matrix(NA, ncol = 7, nrow = 25, 
                   dimnames = list(NULL,c("Degree","Strength","Exp.Inf2",
                                          "Closeness","Betweenness","Video","Band")))
cent_matA1[,1] = scale(centr_degree(ALPHAVideo1net_ig)$res,center = T,scale = T)
cent_matA1[,2] = scale(strength(ALPHAVideo1net_ig),center = T,scale = T)
cent_matA1[,3] = scale(expectedInf(ALPHAVideo1net_ig)$step2,center = T,scale = T)
cent_matA1[,4] = scale(centralization.closeness(ALPHAVideo1net_ig)$res,center = T,scale = T)
cent_matA1[,5] = scale(centralization.betweenness(ALPHAVideo1net_ig)$res,center = T,scale = T)
cent_matA1[,6] = rep(1,25)
cent_matA1[,7] = rep(1,25)
cent_matA1 <-  data.frame(cent_matA1)
cent_matA1$Node = factor(V(ALPHAVideo1net_ig)$name)
centA1 = melt(cent_matA1[,-c(6,7)], id.vars = "Node")
centA1$Video= "Video 1"
centA1s = subset(centA1, variable %in% c("Degree", "Closeness", "Betweenness") )


cent_matA2 =matrix(NA, ncol = 7, nrow = 25, 
                   dimnames = list(NULL,c("Degree","Strength","Exp.Inf2",
                                          "Closeness","Betweenness","Video","Band")))
cent_matA2[,1] = scale(centr_degree(ALPHAVideo2net_ig)$res,center = T,scale = T)
cent_matA2[,2] = scale(strength(ALPHAVideo2net_ig),center = T,scale = T)
cent_matA2[,3] = scale(expectedInf(ALPHAVideo2net_ig)$step2,center = T,scale = T)
cent_matA2[,4] = scale(centralization.closeness(ALPHAVideo2net_ig)$res,center = T,scale = T)
cent_matA2[,5] = scale(centralization.betweenness(ALPHAVideo2net_ig)$res,center = T,scale = T)
cent_matA2[,6] = rep(2,25)
cent_matA2[,7] = rep(1,25)
cent_matA2 <-  data.frame(cent_matA2)
cent_matA2$Node = factor(V(ALPHAVideo2net_ig)$name)
centA2 = melt(cent_matA2[,-c(6,7)], id.vars = "Node")
centA2$Video= "Video 2"
centA2s = subset(centA2, variable %in% c("Degree", "Closeness", "Betweenness") )

cent_matA3 =matrix(NA, ncol = 7, nrow = 25, 
                   dimnames = list(NULL,c("Degree","Strength","Exp.Inf2",
                                          "Closeness","Betweenness","Video","Band")))

cent_matA3[,1] = scale(centr_degree(ALPHAVideo3net_ig)$res,center = T,scale = T)
cent_matA3[,2] = scale(strength(ALPHAVideo3net_ig),center = T,scale = T)
cent_matA3[,3] = scale(expectedInf(ALPHAVideo3net_ig)$step2,center = T,scale = T)
cent_matA3[,4] = scale(centralization.closeness(ALPHAVideo3net_ig)$res,center = T,scale = T)
cent_matA3[,5] = scale(centralization.betweenness(ALPHAVideo3net_ig)$res,center = T,scale = T)
cent_matA3[,6] = rep(3,25)
cent_matA3[,7] = rep(1,25)
cent_matA3 <-  data.frame(cent_matA3)
cent_matA3$Node = factor(V(ALPHAVideo3net_ig)$name)
centA3 = melt(cent_matA3[,-c(6,7)], id.vars = "Node")
centA3$Video= "Video 3"
centA3s = subset(centA3, variable %in% c("Degree", "Closeness", "Betweenness") )


cent_matA4 =matrix(NA, ncol = 7, nrow = 25, 
                   dimnames = list(NULL,c("Degree","Strength","Exp.Inf2",
                                          "Closeness","Betweenness","Video","Band")))
cent_matA4[,1] = scale(centr_degree(ALPHAVideo4net_ig)$res,center = T,scale = T)
cent_matA4[,2] = scale(strength(ALPHAVideo4net_ig),center = T,scale = T)
cent_matA4[,3] = scale(expectedInf(ALPHAVideo4net_ig)$step2,center = T,scale = T)
cent_matA4[,4] = scale(centralization.closeness(ALPHAVideo4net_ig)$res,center = T,scale = T)
cent_matA4[,5] = scale(centralization.betweenness(ALPHAVideo4net_ig)$res,center = T,scale = T)
cent_matA4[,6] = rep(4,25)
cent_matA4[,7] = rep(1,25)
cent_matA4 <-  data.frame(cent_matA4)
cent_matA4$Node = factor(V(ALPHAVideo4net_ig)$name)
centA4 = melt(cent_matA4[,-c(6,7)], id.vars = "Node")
centA4$Video= "Video 4"
centA4s = subset(centA4, variable %in% c("Degree", "Closeness", "Betweenness") )

centall = data.frame(rbind(centA1,centA2,centA3,centA4))
centall$Group = (centall$Group)
centalls = data.frame(rbind(centA1s,centA2s,centA3s,centA4s))

## Convert igraph objects to network class objects
undirnet_ALPHAVideo1net<- asNetwork(ALPHAVideo1net_ig)
undirnet_ALPHAVideo2net<- asNetwork(ALPHAVideo2net_ig)
undirnet_ALPHAVideo3net<- asNetwork(ALPHAVideo3net_ig)
undirnet_ALPHAVideo4net<- asNetwork(ALPHAVideo4net_ig)


# Plot the networks -------------------------------------------------------



# Define node colors 
plotcols = c("lightsteelblue3", "lightsteelblue3", "lightsteelblue3", "lightsteelblue3", "lightsteelblue3", 
             "lightsteelblue3", "lightsteelblue3", 
             "steelblue3", "steelblue3","steelblue3", "steelblue3","steelblue3", 
             "steelblue3","steelblue3", "steelblue3","steelblue3",
             "grey","lightgreen","lightgreen","lightgreen",
             "grey","grey","gold2",
             "gold2","grey")

# load the node coordinates
load("provaeeg.RData")


# Fix the node coordinates
undirnet_ALPHAVideo1net %v% "x" = undirnet_ALPHAVideo2net %v% "x" = undirnet_ALPHAVideo3net %v% "x" = undirnet_ALPHAVideo4net %v% "x" = dd[, 1]
undirnet_ALPHAVideo1net %v% "y" = undirnet_ALPHAVideo2net %v% "y" = undirnet_ALPHAVideo3net %v% "y" = undirnet_ALPHAVideo4net %v% "y" = dd[, 2]
undirnet_ALPHAVideo1net %v% "nodecol" <- undirnet_ALPHAVideo2net %v% "nodecol" <- undirnet_ALPHAVideo3net %v% "nodecol" <- undirnet_ALPHAVideo4net %v% "nodecol" <- plotcols
undirnet_ALPHAVideo1net %e% "aweight" <- abs(undirnet_ALPHAVideo1net %e% "weight")
undirnet_ALPHAVideo2net %e% "aweight" <- abs(undirnet_ALPHAVideo2net %e% "weight")
undirnet_ALPHAVideo3net %e% "aweight" <- abs(undirnet_ALPHAVideo3net %e% "weight")
undirnet_ALPHAVideo4net %e% "aweight" <- abs(undirnet_ALPHAVideo4net %e% "weight")


# Plot the network of each group
plotundir_ALPHAVideo1net <- ggnet2(undirnet_ALPHAVideo1net, 
                                   mode = c("x","y"),
                                   color = "nodecol",
                                   edge.color = c("color", "black"),
                                   size="degree", label.size = 5,max_size = 20,
                                   label=TRUE)+ guides(size = "none") + 
  annotate(geom="text", x=0.1, y=0, label=paste("density:",round(edge_density(ALPHAVideo1net_ig),3),sep = " "), color="black")
plotundir_ALPHAVideo1net 


plotundir_ALPHAVideo2net <- ggnet2(undirnet_ALPHAVideo2net, 
                                   mode = c("x","y"),
                                   color = "nodecol",
                                   edge.color = c("color", "black"),
                                   size="degree", label.size = 5,max_size = 20,
                                   label=TRUE)+ guides(size = "none") + 
  annotate(geom="text", x=0.1, y=0, label=paste("density:",round(edge_density(ALPHAVideo2net_ig),3),sep = " "), color="black")
plotundir_ALPHAVideo2net 

plotundir_ALPHAVideo3net <- ggnet2(undirnet_ALPHAVideo3net, 
                                   mode = c("x","y"),
                                   color = "nodecol",
                                   edge.color = c("color", "black"),
                                   size="degree", label.size = 5,max_size = 20,
                                   label=TRUE)+ guides(size = "none") + 
  annotate(geom="text", x=0.1, y=0, label=paste("density:",round(edge_density(ALPHAVideo3net_ig),3),sep = " "), color="black")
plotundir_ALPHAVideo3net


plotundir_ALPHAVideo4net <- ggnet2(undirnet_ALPHAVideo4net, 
                                   mode = c("x","y"),
                                   color = "nodecol",
                                   edge.color = c("color", "black"),
                                   size="degree", label.size = 5,max_size = 20,
                                   label=TRUE)+ guides(size = "none") + 
  annotate(geom="text", x=0.1, y=0, label=paste("density:",round(edge_density(ALPHAVideo4net_ig),3),sep = " "), color="black")
plotundir_ALPHAVideo4net 



# Plot the network centralities -------------------------------------------



source("altcoord_curvedpolar.R", echo=FALSE)
cp1 =coord_polar(theta = "x", start = 0,  clip = "off")
cp1$is_free = function() TRUE
cp2 = altcoord_curvedpolar(theta = "x", start = 0,  clip = "off")
cp2$is_free = function() TRUE

bfakedat1<- centA1s[centA1s$Node == "AF3", ]
bfakedat1$Node <- ""
centA1N = data.frame(rbind(centA1s,bfakedat1))


membershipCol = rep(NA, n.elec+1)
splitnames = strsplit(levels(centA1N$Node),"(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])",perl=TRUE)
splitnames[[(n.elec+1)]] = c("1","1")
intex = as.numeric(sapply(splitnames,function(x) x[2]))/2
membershipCol[is.na(intex)] = "#0072B2"
membershipCol[intex-floor(intex)==0] = "#009E73"
membershipCol[is.na(membershipCol)] = "#CC79A7"



plotggA1 = ggplot(centA1N, aes(x=(Node),y=value,group=variable))+ylim(c(-5,5)) + 
  geom_line(aes(color=variable),linewidth=1.2)+ coord_polar()+
  scale_color_colorblind()+ 
  ylab("Centrality (standardized)") + xlab("Electrodes")+ geom_hline(yintercept = 0,col="#D55E00", linewidth=1.1) + labs(col="Centrality:") + 
  theme_bw()+ theme(legend.key.height = unit(1.5, 'cm'),
                    legend.key.size = unit(1.5, 'cm'),
                    #legend.key.width = unit(2, 'cm'),
                    legend.text = element_text(size = 15,face = "bold"),
                    legend.title = element_text(size = 18,face = "bold",color = "aquamarine4"),
                    plot.title = element_text(face="bold",hjust = 0.5))+
  ggtitle("Centrality indices for each electrode \nin Video 1 EEG network")
plotggA1 = plotggA1 + scale_x_discrete(expand = c(0,0))+cp2
plotggA1 = plotggA1 + theme(axis.text.x = element_text(colour = membershipCol))
plotggA1$theme$axis.text.x$colour = plotggA1$theme$axis.text.x$colour[c(2:26,1)]
plotggA1 = plotggA1 + theme(axis.text.y  = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.x =  element_text(size=13,face = "bold"))
plotggA1

bfakedat2<- centA2s[centA2s$Node == "AF3", ]
bfakedat2$Node <- ""
centA2N = data.frame(rbind(centA2s,bfakedat2))

membershipCol = rep(NA, n.elec+1)
splitnames = strsplit(levels(centA2N$Node),"(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])",perl=TRUE)
splitnames[[(n.elec+1)]] = c("1","1")
intex = as.numeric(sapply(splitnames,function(x) x[2]))/2
membershipCol[is.na(intex)] = "#0072B2"
membershipCol[intex-floor(intex)==0] = "#009E73"
membershipCol[is.na(membershipCol)] = "#CC79A7"


plotggA2 = ggplot(centA2N, aes(x=(Node),y=value,group=variable))+ylim(c(-5,5)) + 
  geom_line(aes(color=variable),linewidth=1.2)+ coord_polar()+
  scale_color_colorblind()+ 
  ylab("Centrality (standardized)") + xlab("Electrodes")+ geom_hline(yintercept = 0,col="#D55E00", linewidth=1.1) + labs(col="Centrality:") + 
  theme_bw()+ theme(legend.key.height = unit(1.5, 'cm'),
                    legend.key.size = unit(1.5, 'cm'),
                    #legend.key.width = unit(2, 'cm'),
                    legend.text = element_text(size = 15,face = "bold"),
                    legend.title = element_text(size = 18,face = "bold",color = "aquamarine4"),
                    plot.title = element_text(face="bold",hjust = 0.5))+
  ggtitle("Centrality indices for each electrode \nin Video 1 EEG network")
plotggA2 = plotggA2 + scale_x_discrete(expand = c(0,0))+cp2
plotggA2 = plotggA2 + theme(axis.text.x = element_text(colour = membershipCol))
plotggA2$theme$axis.text.x$colour = plotggA2$theme$axis.text.x$colour[c(2:26,1)]
plotggA2 = plotggA2 + theme(axis.text.y  = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.x = element_blank(),
                            axis.text.x =  element_text(size=13,face = "bold"))
plotggA2


bfakedat3<- centA3s[centA3s$Node == "AF3", ]
bfakedat3$Node <- ""
centA3N = data.frame(rbind(centA3s,bfakedat3))

membershipCol = rep(NA, n.elec+1)
splitnames = strsplit(levels(centA3N$Node),"(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])",perl=TRUE)
splitnames[[(n.elec+1)]] = c("1","1")
intex = as.numeric(sapply(splitnames,function(x) x[2]))/2
membershipCol[is.na(intex)] = "#0072B2"
membershipCol[intex-floor(intex)==0] = "#009E73"
membershipCol[is.na(membershipCol)] = "#CC79A7"


plotggA3 = ggplot(centA3N, aes(x=(Node),y=value,group=variable))+ylim(c(-5,5)) + 
  geom_line(aes(color=variable),linewidth=1.2)+ coord_polar()+
  scale_color_colorblind()+ 
  ylab("Centrality (standardized)") + xlab("Electrodes")+ geom_hline(yintercept = 0,col="#D55E00", linewidth=1.1) + labs(col="Centrality:") + 
  theme_bw()+ theme(legend.key.height = unit(1.5, 'cm'),
                    legend.key.size = unit(1.5, 'cm'),
                    #legend.key.width = unit(2, 'cm'),
                    legend.text = element_text(size = 15,face = "bold"),
                    legend.title = element_text(size = 18,face = "bold",color = "aquamarine4"),
                    plot.title = element_text(face="bold",hjust = 0.5))+
  ggtitle("Centrality indices for each electrode \nin Video 1 EEG network")
plotggA3 = plotggA3 + scale_x_discrete(expand = c(0,0))+cp2
plotggA3 = plotggA3 + theme(axis.text.x = element_text(colour = membershipCol))
plotggA3$theme$axis.text.x$colour = plotggA3$theme$axis.text.x$colour[c(2:26,1)]
plotggA3 = plotggA3 + theme(axis.text.y  = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.x = element_blank(),
                            axis.text.x =  element_text(size=13,face = "bold"))

plotggA3

bfakedat4<- centA4s[centA4s$Node == "AF3", ]
bfakedat4$Node <- ""
centA4N = data.frame(rbind(centA4s,bfakedat4))

membershipCol = rep(NA, n.elec+1)
splitnames = strsplit(levels(centA4N$Node),"(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])",perl=TRUE)
splitnames[[(n.elec+1)]] = c("1","1")
intex = as.numeric(sapply(splitnames,function(x) x[2]))/2
membershipCol[is.na(intex)] = "#0072B2"
membershipCol[intex-floor(intex)==0] = "#009E73"
membershipCol[is.na(membershipCol)] = "#CC79A7"


plotggA4 = ggplot(centA4N, aes(x=(Node),y=value,group=variable))+ylim(c(-5,5)) + 
  geom_line(aes(color=variable),linewidth=1.2)+ coord_polar()+
  scale_color_colorblind()+ 
  ylab("Centrality (standardized)") + xlab("Electrodes")+ geom_hline(yintercept = 0,col="#D55E00", linewidth=1.1) + labs(col="Centrality:") + 
  theme_bw()+ theme(legend.key.height = unit(1.5, 'cm'),
                    legend.key.size = unit(1.5, 'cm'),
                    #legend.key.width = unit(2, 'cm'),
                    legend.text = element_text(size = 15,face = "bold"),
                    legend.title = element_text(size = 18,face = "bold",color = "aquamarine4"),
                    plot.title = element_text(face="bold",hjust = 0.5))+
  ggtitle("Centrality indices for each electrode \nin Video 1 EEG network")
plotggA4 = plotggA4 + scale_x_discrete(expand = c(0,0))+cp2
plotggA4 = plotggA4 + theme(axis.text.x = element_text(colour = membershipCol))
plotggA4$theme$axis.text.x$colour = plotggA4$theme$axis.text.x$colour[c(2:26,1)]
plotggA4 = plotggA4 + theme(axis.text.y  = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.x = element_blank(),
                            axis.text.x =  element_text(size=13,face = "bold"))

plotggA4



# Combine individual plots (networks and centralities)  into a single plot --------


## Combine the network plots
plotundir_ALPHAVideo1net = plotundir_ALPHAVideo1net+ggtitle("Video 1 EEG network")+ theme_linedraw()+theme(legend.key.height = unit(1.5, 'cm'),
                                                                                           legend.key.size = unit(1.5, 'cm'),
                                                                                           plot.margin = unit(c(0,0,0,0), "cm"),
                                                                                           axis.title = element_blank(),
                                                                                           plot.title = element_text(hjust = 0.5,face = "bold",size = 14))
plotundir_ALPHAVideo2net = plotundir_ALPHAVideo2net+ggtitle("Video 2 EEG network")+ theme_linedraw()+theme(legend.key.height = unit(1.5, 'cm'),
                                                                                           legend.key.size = unit(1.5, 'cm'),
                                                                                           plot.margin = unit(c(0,0,0,0), "cm"),
                                                                                           axis.title = element_blank(),
                                                                                           plot.title = element_text(hjust = 0.5,face = "bold",size = 14))
plotundir_ALPHAVideo3net = plotundir_ALPHAVideo3net+ggtitle("Video 3 EEG network")+ theme_linedraw()+theme(legend.key.height = unit(1.5, 'cm'),
                                                                                           legend.key.size = unit(1.5, 'cm'),
                                                                                           plot.margin = unit(c(0,0,0,0), "cm"),
                                                                                           axis.title = element_blank(),
                                                                                           
                                                                                           plot.title = element_text(hjust = 0.5,face = "bold",size = 14))
plotundir_ALPHAVideo4net = plotundir_ALPHAVideo4net+ggtitle("Video 4 EEG network")+ theme_linedraw()+theme(legend.key.height = unit(1.5, 'cm'),
                                                                                           legend.key.size = unit(1.5, 'cm'),
                                                                                           plot.margin = unit(c(0,0,0,0), "cm"),
                                                                                           axis.title = element_blank(),
                                                                                           
                                                                                           plot.title = element_text(hjust = 0.5,face = "bold",size = 14))
combnets  = (plotundir_ALPHAVideo1net + plotundir_ALPHAVideo2net)/(plotundir_ALPHAVideo3net + plotundir_ALPHAVideo4net) & theme(legend.position = "bottom")
combnets 


## Combine the plots for the centrality indices
plotggA1 = plotggA1 + ggtitle("Centralities for the electrodes\n in Video 1 EEG network")+ theme(axis.title.y = element_blank(),legend.position = "bottom")
plotggA2 = plotggA2 + ggtitle("Centralities for the electrodes\n in Video 2 EEG network")+ theme(axis.title.y = element_blank(),legend.position = "bottom")
plotggA3 = plotggA3 + ggtitle("Centralities for the electrodes\n in Video 3 EEG network")+ theme(axis.title.y = element_blank(),legend.position = "bottom")
plotggA4 = plotggA4 + ggtitle("Centralities for the electrodes\n in Video 4 EEG network")+ theme(axis.title.y = element_blank(),legend.position = "bottom")

combcent  = (plotggA1 + plotggA2)/(plotggA3 + plotggA4) & theme(legend.position = "bottom")
combcent + plot_layout(guides = "collect")


# Check central electrodes ------------------------------------------------

video1central = split(subset(centA1s, value>0), list(subset(centA1s, value>0)$variable))[c(1,4,5)]
video2central = split(subset(centA2s, value>0), list(subset(centA2s, value>0)$variable))[c(1,4,5)]
video3central = split(subset(centA3s, value>0), list(subset(centA3s, value>0)$variable))[c(1,4,5)]
video4central = split(subset(centA4s, value>0), list(subset(centA4s, value>0)$variable))[c(1,4,5)]


# Create topoplots from centrality indices --------------------------------


Degreevid1 = video1central[[1]]
rownames(Degreevid1) = NULL
Degreevid2 = video2central[[1]]
rownames(Degreevid2) = NULL
Degreevid3 = video3central[[1]]
rownames(Degreevid3) = NULL
Degreevid4 = video4central[[1]]
rownames(Degreevid4) = NULL



Closenvid1 = video1central[[2]]
rownames(Closenvid1) = NULL
Closenvid2 = video2central[[2]]
rownames(Closenvid2) = NULL
Closenvid3 = video3central[[2]]
rownames(Closenvid3) = NULL
Closenvid4 = video4central[[2]]
rownames(Closenvid4) = NULL


Betweenvid1 = video1central[[3]]
rownames(Betweenvid1) = NULL
Betweenvid2 = video2central[[3]]
rownames(Betweenvid2) = NULL
Betweenvid3 = video3central[[3]]
rownames(Betweenvid3) = NULL
Betweenvid4 = video4central[[3]]
rownames(Betweenvid4) = NULL

video1centraltop = split(centA1s, list(centA1s$variable))[c(1,4,5)]
video2centraltop = split(centA2s, list(centA2s$variable))[c(1,4,5)]
video3centraltop = split(centA3s, list(centA3s$variable))[c(1,4,5)]
video4centraltop = split(centA4s, list(centA4s$variable))[c(1,4,5)]


datatopo1 = dd
names(datatopo1) = c("x","y")
datatopo1$electrode = V(ALPHAVideo1net_ig)$name
datatopo1$amplitude = video1centraltop$Betweenness$value

datatopo2 = dd
names(datatopo2) = c("x","y")
datatopo2$electrode = V(ALPHAVideo1net_ig)$name
datatopo2$amplitude = video2centraltop$Betweenness$value

datatopo3 = dd
names(datatopo3) = c("x","y")
datatopo3$electrode = V(ALPHAVideo1net_ig)$name
datatopo3$amplitude = video3centraltop$Betweenness$value

datatopo4 = dd
names(datatopo4) = c("x","y")
datatopo4$electrode = V(ALPHAVideo1net_ig)$name
datatopo4$amplitude = video4centraltop$Betweenness$value

datatopo11 = dd
names(datatopo11) = c("x","y")
datatopo11$electrode = V(ALPHAVideo1net_ig)$name
datatopo11$amplitude = video1centraltop$Degree$value

datatopo22 = dd
names(datatopo22) = c("x","y")
datatopo22$electrode = V(ALPHAVideo1net_ig)$name
datatopo22$amplitude = video2centraltop$Degree$value

datatopo33 = dd
names(datatopo33) = c("x","y")
datatopo33$electrode = V(ALPHAVideo1net_ig)$name
datatopo33$amplitude = video3centraltop$Degree$value

datatopo44 = dd
names(datatopo44) = c("x","y")
datatopo44$electrode = V(ALPHAVideo1net_ig)$name
datatopo44$amplitude = video4centraltop$Degree$value

electrodeLocs <- read_delim("https://raw.githubusercontent.com/craddm/ExploringERPs/master/biosemi70elecs.loc",
                            "\t",
                            escape_double = FALSE,
                            col_names = c("chanNo","theta","radius","electrode"),
                            trim_ws = TRUE)


datatopo1$electrode[1] = "Fpz"
datatopo1$electrode[2] = "Fp2"
datatopo1$electrode[6] = "Fp1"
datatopo1$electrode[17] = "T8"
datatopo1$electrode[21] = "T7"
datatopo1$electrode[22] = "TP7"
datatopo1$electrode[25] = "TP8"


datatopo1$electrode[which(is.na(match(datatopo1$electrode, electrodeLocs$electrode)))]
electrodeLocs = electrodeLocs[match(datatopo1$electrode, electrodeLocs$electrode),]
electrodeLocs$radianTheta <- pi/180*electrodeLocs$theta
electrodeLocs <- electrodeLocs %>%
  mutate(x = .$radius*sin(.$radianTheta),
         y = .$radius*cos(.$radianTheta))


cartesian <- ggplot(electrodeLocs,
                    aes(x, y, label = electrode))+
  geom_text()+
  theme_bw()+
  coord_equal()


theme_topo <- function(base_size = 12)
{
  theme_bw(base_size = base_size) %+replace%
    theme(
      rect             = element_blank(),
      line             = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    )
}
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100) {
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
headShape <- circleFun(c(0, 0), round(max(electrodeLocs$x)), npoints = 100) # 0
nose <- data.frame(x = c(-0.075,0,.075),y=c(.495,.575,.495))

allData = electrodeLocs
allData$amplitude = datatopo1$amplitude
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
singleTimepoint <- allData
gridRes <- 100 # Specify the number of points for each grid dimension i.e. the resolution/smoothness of the interpolation
tmpTopo <- with(singleTimepoint,akima::interp(x = x, y = y, z = amplitude,
                                              xo = seq(min(x)*2,
                                                       max(x)*2,
                                                       length = gridRes),
                                              yo = seq(min(y)*2,
                                                       max(y)*2,
                                                       length = gridRes),linear = FALSE,extrap = TRUE)) 
interpTopo <- data.frame(x = tmpTopo$x, tmpTopo$z)
names(interpTopo)[1:length(tmpTopo$y)+1] <- tmpTopo$y
interpTopo <- gather(interpTopo,key = y,value = amplitude,-x,convert = TRUE)
interpTopo$incircle <- sqrt(interpTopo$x^2 + interpTopo$y^2) < .7 # mark grid elements that are outside of the plotting circle
interpTopo <- interpTopo[interpTopo$incircle,] #remove the elements outside the circle
maskRing <- circleFun(diameter = 1.42) #create a circle round the outside of the plotting area to mask the jagged edges of the interpolation

akimaPlot <- ggplot(interpTopo,aes(x = x, y = y, fill = amplitude)) +
  geom_raster() +
  theme_topo()+
  scale_fill_gradientn(colours = turbo(20),
                       limits = c(-5,5),
                       guide = "colourbar",
                       oob = squish) + 
  # stat_contour(aes(z = amplitude),
  #              colour = "black",
  #              binwidth = 0.5) +
  geom_path(data = maskRing,
            aes(x, y, z = NULL, fill =NULL),
            colour = "white",
            size = 6)+
  geom_point(data = singleTimepoint,
             aes(x, y),
             size = 1)+
  geom_path(data = headShape,
            aes(x, y, z = NULL, fill = NULL),
            size = 1.5)+
  geom_path(data = nose,
            aes(x, y, z = NULL, fill = NULL),
            size = 1.5)+
  coord_equal()

akimaPlot = akimaPlot + labs(fill="Degree")
#save(akimaPlot11,file = "akimaPlot.RData")
akimaPlot




akimaPlot1 = akimaPlot1+theme(legend.title = element_text(size=15,face = "bold"),legend.text  = element_text(size=10,face = "bold"))
akimaPlot2 = akimaPlot2+theme(legend.title = element_text(size=15,face = "bold"),legend.text  = element_text(size=10,face = "bold"))
akimaPlot3 = akimaPlot3+theme(legend.title = element_text(size=15,face = "bold"),legend.text  = element_text(size=10,face = "bold"))
akimaPlot4 = akimaPlot4+theme(legend.title = element_text(size=15,face = "bold"),legend.text  = element_text(size=10,face = "bold"))
combtopoB  = (akimaPlot1 + akimaPlot2)/(akimaPlot3 + akimaPlot4) & theme(legend.position = "bottom")
combtopoB + plot_layout(guides = "collect")




akimaPlot11 = akimaPlot11+theme(legend.title = element_text(size=15,face = "bold"),legend.text  = element_text(size=10,face = "bold"))
akimaPlot22 = akimaPlot22+theme(legend.title = element_text(size=15,face = "bold"),legend.text  = element_text(size=10,face = "bold"))
akimaPlot33 = akimaPlot33+theme(legend.title = element_text(size=15,face = "bold"),legend.text  = element_text(size=10,face = "bold"))
akimaPlot44 = akimaPlot44+theme(legend.title = element_text(size=15,face = "bold"),legend.text  = element_text(size=10,face = "bold"))
combtopoD  = (akimaPlot11 + akimaPlot22)/(akimaPlot33 + akimaPlot44) & theme(legend.position = "bottom")
combtopoD + plot_layout(guides = "collect")





p12n = SpyrosplotCCC.NCT(COMPARE_12_ALPHAP_MEANb,"network")+
  ggtitle("Network structure comparison \nVideo 1 - Video 2")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

p13n = SpyrosplotCCC.NCT(COMPARE_13_ALPHAP_MEANb,"network")+
  ggtitle("Network structure comparison \nVideo 1 - Video 3")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))


p14n = SpyrosplotCCC.NCT(COMPARE_14_ALPHAP_MEANb,"network")+
  ggtitle("Network structure comparison \nVideo 1 - Video 4")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

p23n = SpyrosplotCCC.NCT(COMPARE_23_ALPHAP_MEANb,"network")+
  ggtitle("Network structure comparison \nVideo 2 - Video 3")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

p24n = SpyrosplotCCC.NCT(COMPARE_24_ALPHAP_MEANb,"network")+
  ggtitle("Network structure comparison \nVideo 2 - Video 4")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

p34n = SpyrosplotCCC.NCT(COMPARE_34_ALPHAP_MEANb,"network")+
  ggtitle("Network structure comparison \nVideo 3 - Video 4")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

combtests  = (p12n + p13n + p14n)/(p23n + p24n + p34n) 
combtests 




p12c = SpyrosplotCCC.NCT(COMPARE_12_ALPHAP_MEANb,"strength")+
  ggtitle("Network strength comparison \nVideo 1 - Video 2")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

p13c = SpyrosplotCCC.NCT(COMPARE_13_ALPHAP_MEANb,"strength")+
  ggtitle("Network strength comparison \nVideo 1 - Video 3")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))


p14c = SpyrosplotCCC.NCT(COMPARE_14_ALPHAP_MEANb,"strength")+
  ggtitle("Network strength comparison \nVideo 1 - Video 4")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

p23c = SpyrosplotCCC.NCT(COMPARE_23_ALPHAP_MEANb,"strength")+
  ggtitle("Network strength comparison \nVideo 2 - Video 3")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

p24c = SpyrosplotCCC.NCT(COMPARE_24_ALPHAP_MEANb,"strength")+
  ggtitle("Network strength comparison \nVideo 2 - Video 4")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

p34c = SpyrosplotCCC.NCT(COMPARE_34_ALPHAP_MEANb,"strength")+
  ggtitle("Network strength comparison \nVideo 3 - Video 4")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

combtestst  = (p12c + p13c + p14c)/(p23c + p24c + p34c) 
combtestst
