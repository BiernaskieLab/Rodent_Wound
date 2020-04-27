library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(stats)
library(plotly)
library(tidyr)

########################################################################################################################
########## Plotting Interactome Profiles for Small Wound Day 8 cells: ##################################################
########################################################################################################################
# CSV provided under sarthaksinha123/Rodent_Wound/SmallWoundDay8_Interactome_Analysis
stage1.1 <- read.csv("Final_SW8.csv")
stage1.1_long <- gather(data = stage1.1, key = Interactions, value = Significant_Means, -c(1))
head(stage1.1_long)
ordered_y = stage1.1$interacting_pair
#View(stage1.1$interacting_pair)
stage1.1_heatmap <- ggplot(data = stage1.1_long, 
                           mapping = aes(x = Interactions, 
                                         y = factor(interacting_pair, level = ordered_y), 
                                         fill = Significant_Means)) +
  geom_tile() +
  xlab(label = "Cell Populations")+
  ylab(label = "Ligand-Receptor Pairs")+
  scale_fill_gradient(low = "white", high = "red", guide="colorbar", limits=c(0,3.16))

stage1.1_heatmap + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jpeg(file = "Epi_Mes_CPD_final_SW8.jpeg", width = 13, height = 20, units = "cm", res = 500)
stage1.1_heatmap + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


########################################################################################################################
########## Plotting Interactome Profiles for Large Wound Day 14 cells: #################################################
########################################################################################################################
# CSV provided under sarthaksinha123/Rodent_Wound/LargeWoundDay14_Interactome_Analysis/
stage1.1 <- read.csv("Final.csv")
stage1.1_long <- gather(data = stage1.1, key = Interactions, value = Significant_Means, -c(1))
head(stage1.1_long)
ordered_y = stage1.1$interacting_pair

#View(stage1.1$interacting_pair)
stage1.1_heatmap <- ggplot(data = stage1.1_long, 
                           mapping = aes(x = Interactions, 
                                         y = factor(interacting_pair, level = ordered_y), 
                                         fill = Significant_Means)) +
  geom_tile() +
  xlab(label = "Cell Populations")+
  ylab(label = "Ligand-Receptor Pairs")+
  scale_fill_gradient(low = "white", high = "red", guide="colorbar", limits=c(0,2.4))

stage1.1_heatmap + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jpeg(file = "Epi_Mes_CPD_final.jpeg", width = 13, height = 20, units = "cm", res = 500)
stage1.1_heatmap + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


########################################################################################################################
########## Plotting Interactome Profiles for Large Wound Center Dermal Condensate and Epithelial Placode cells: ########
########################################################################################################################
# CSV provided under sarthaksinha123/Rodent_Wound/DermalCondensate_EpithelialPlacode_Interactome_Analysis/
stage1.1 <- read.csv("Plac_Con.csv")
stage1.1_long <- gather(data = stage1.1, key = Interactions, value = Significant_Means, -c(1))
head(stage1.1_long)
ordered_y = stage1.1$interacting_pair
#View(stage1.1$interacting_pair)
stage1.1_heatmap <- ggplot(data = stage1.1_long, 
                           mapping = aes(x = Interactions, 
                                         y = factor(interacting_pair, level = ordered_y), 
                                         fill = Significant_Means)) +
  geom_tile() +
  xlab(label = "Cell Populations")+
  ylab(label = "Ligand-Receptor Pairs")+
  scale_fill_gradient(low = "white", high = "red", guide="colorbar", limits=c(0,2))

stage1.1_heatmap + theme(axis.text.x = element_text(angle = 90, hjust = 1))

jpeg(file = "Epi_Mes_CPD_final_Plac_Cond.jpeg", width = 13, height = 20, units = "cm", res = 500)
stage1.1_heatmap + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


