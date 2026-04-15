# Re-analysis of alpha diversity [13 March]

## Purpose:
1. To re-analyze the alpha diversity of CD and healthy individuals in the dataset by using filtered phyloseq object. 

## Materials:
1. R Studio 
2. filtered phyloseq object
## Methods:  
#Load all relevant packages
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(picante)

#Load filtered phloseq object 
load("Ryan_filt_new.RData")

#generating rarefaction curve, CANNOT run the alpha rarefaction curve
rarefaction_curve <- rarecurve(t(as.data.frame(otu_table(ryan_filt_new))), cex=0.1)
print(rarefaction_curve)

#Rarefaction 
ryan_rare_new <- rarefy_even_depth(ryan_filt_new, rngseed = 1, sample.size = 26522)

#Save the file
save(ryan_rare_new, file="ryan_rare_new.RData")

#Calculate alpha diversity and merge with the metadata
alphadiv <- estimate_richness(ryan_rare_new)
samp_dat <- sample_data(ryan_rare_new)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

#Create a formatted alpha diversity plot (Shannon)
alpha_plot_shannon <- ggplot(samp_dat_wdiv, aes(x = Condition,
                                        y = Shannon)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(aes(fill=Condition)) +
  scale_fill_manual(values = c(
    "Healthy" = "#7570B3",
    "Crohn's Disease" = "#D55E00"
  )) +
  labs(x = "Condition", y = "Shannon Diversity Measure") + 
  theme_classic(base_size =13) + 
  theme(legend.position = "none")
alpha_plot_shannon

ggsave(filename = "plot_Shannon.png"
       , alpha_plot_shannon
       , height=4, width=6)

#Create a formatted alpha diversity plot (Observed Features)
alpha_plot_observed <- ggplot(samp_dat_wdiv, aes(x = Condition,
                                                y = Observed)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(aes(fill=Condition)) +
  scale_fill_manual(values = c(
    "Healthy" = "#7570B3",
    "Crohn's Disease" = "#D55E00"
  )) +
  labs(x = "Condition", y = "Observed Features Diversity Measure") + 
  theme_classic(base_size =13) + 
  theme(legend.position = "none")
alpha_plot_observed

ggsave(filename = "plot_Observed.png"
       , alpha_plot_observed
       , height=4, width=6)

#Faith's Phylogenetic Diversity
phylo_dist <- pd(t(otu_table(ryan_rare_new)), phy_tree(ryan_rare_new),
                 include.root=F) 

sample_data(ryan_rare_new)$PD <- phylo_dist$PD

alpha_plot_pd <- ggplot(sample_data(ryan_rare_new), aes(x = Condition,
                                                 y = PD)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(aes(fill=Condition)) +
  scale_fill_manual(values = c(
    "Healthy" = "#7570B3",
    "Crohn's Disease" = "#D55E00"
  )) +
  labs(x = "Condition", y = "Faith's Phylogenetic Diversity Measure") + 
  theme_classic(base_size = 13) + 
  theme(legend.position = "none")

alpha_plot_pd

ggsave(filename = "plot_pd.png"
       , alpha_plot_pd
       , height=4, width=6)
