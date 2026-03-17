#### Re-do the diversity analysis by using filtered phyloseq object ####

# Load all relevant packages
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(picante)

# Load filtered phloseq object 
load("Ryan_filt_new.RData")

#### Rarefaction and Generating Rarefaction Curve####
# generating rarefaction curve, CANNOT run the alpha rarefaction curve
#rarefaction_curve <- rarecurve(t(as.data.frame(otu_table(ryan_filt_new))), cex=0.1)
#print(rarefaction_curve)

# Rarefaction 
ryan_rare_new <- rarefy_even_depth(ryan_filt_new, rngseed = 1, sample.size = 26522)

# Save the file
save(ryan_rare_new, file="ryan_rare_new.RData")

#### Alpha Diversity ####

#Alpha diversity analysis 
plot_richness(ryan_rare_new) 

# Observed features and Shannon 
plot_richness(ryan_rare_new, measures = c("Shannon","Observed")) 

gg_richness <- plot_richness(ryan_rare_new, x = "Condition", measures = c("Shannon","Observed")) +
  xlab("Subject ID") +
  geom_boxplot()
gg_richness

ggsave(filename = "plot_richness_new.png"
       , gg_richness
       , height=4, width=6)

# Faith's Phylogenetic Diversity 
phylo_dist <- pd(t(otu_table(ryan_rare_new)), phy_tree(ryan_rare_new),
                 include.root=F) 

sample_data(ryan_rare_new)$PD <- phylo_dist$PD

plot.pd <- ggplot(sample_data(ryan_rare_new), aes(Condition, PD)) + 
  geom_boxplot() +
  xlab("Subject ID") +
  ylab("Phylogenetic Diversity")

