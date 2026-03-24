# Significance Analysis for Alpha Diversity Metrics

# Load the packages 
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(picante)
library(ggsignif)
library(ggpubr)

# Load the relevant items in the RData file 
load("Alpha_formatted.RData")

# Calculate alpha diversity and merge with the metadata
alphadiv <- estimate_richness(ryan_rare_new)
samp_dat <- sample_data(ryan_rare_new)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

samp_dat_wdiv$Condition <- factor(
  samp_dat_wdiv$Condition,
  levels = c("Healthy", "Crohn's Disease")
)
#### Shannon Figure####  
# Wilcoxon Rank Sum Test 
wilcox.test(Shannon ~ Condition, data=samp_dat_wdiv, exact = FALSE)

# Figure generation 
alpha_plot_shannon <- ggplot(samp_dat_wdiv, aes(x = Condition,
                                                y = Shannon)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(aes(fill=Condition)) +
  scale_fill_manual(values = c(
    "Healthy" = "#7570B3",
    "Crohn's Disease" = "#D55E00"
  )) +
  labs(x = NULL, y = "Shannon Diversity Measure") + 
  geom_signif(comparisons = list(c("Healthy", "Crohn's Disease")), 
             test = "wilcox.test",
             annotations = c("p=0.03"), 
             y_position = 4.6) +
  theme_classic(base_size =13) + 
  theme(legend.position = "none")
alpha_plot_shannon

# Save the figure 
ggsave(filename = "plot_Shannon_sig.png"
       , alpha_plot_shannon
       , height=4, width=6)

#### Observed Features Figure####
# Wilcoxon Rank Sum Test 
wilcox.test(Observed ~ Condition, data=samp_dat_wdiv, exact = FALSE)

# Figure generation 
alpha_plot_observed <- ggplot(samp_dat_wdiv, aes(x = Condition,
                                                 y = Observed)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(aes(fill=Condition)) +
  scale_fill_manual(values = c(
    "Healthy" = "#7570B3",
    "Crohn's Disease" = "#D55E00"
  )) +
  labs(x = NULL, y = "Observed Features") + 
  geom_signif(comparisons = list(c("Healthy", "Crohn's Disease")), 
              test = "wilcox.test",
              annotations = c("p=0.0005"),, 
              y_position = 300) +
  theme_classic(base_size =13) + 
  theme(legend.position = "none")
alpha_plot_observed

# Export the figure 
ggsave(filename = "plot_Observed_sig.png"
       , alpha_plot_observed
       , height=4, width=6)

#### PD Figure ####
# faith PD calculation 
phylo_dist <- pd(t(otu_table(ryan_rare_new)), phy_tree(ryan_rare_new),
                 include.root=F) 

samp_dat_wdiv$PD <- phylo_dist$PD

# Wilcoxon Rank Sum Test 
wilcox.test(PD ~ Condition, data=samp_dat_wdiv, exact = FALSE)

# Figure generation 
alpha_plot_pd <- ggplot(samp_dat_wdiv, aes(x = Condition,
                                                        y = PD)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(aes(fill=Condition)) +
  scale_fill_manual(values = c(
    "Healthy" = "#7570B3",
    "Crohn's Disease" = "#D55E00"
  )) +
  labs(x = NULL, y = "Faith's Phylogenetic Diversity Measure") + 
  geom_signif(comparisons = list(c("Healthy", "Crohn's Disease")), 
              test = "wilcox.test",
              annotations = c("p=0.0007"),, 
              y_position = 16) +
  theme_classic(base_size = 13) + 
  theme(legend.position = "none")
alpha_plot_pd
