#### MICB 475 Project 2 ####
# Team 6 Ryan's dataset 

# Load all relevant packages
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(picante)

#### Load data####
# Load metadata file 
metafp <- "ryan_export/ryan_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

# Load feature-table file 
otufp <- "ryan_export/table_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

#load taxonomy file 
taxfp <- "ryan_export/taxonomy_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

#load phylogenetic tree file 
phylotreefp <- "ryan_export/rooted_tree_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Metadata filtering ####
# Filter out UC, non-inflamed CD from the metadata
meta_filt <- meta %>%
  filter(Condition == "Healthy" | 
           (Condition == "Crohn's Disease" & 
           Histological.status == "Inflamed tissue "))

#### Formatting relevant files####

# Format OUT table
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

# Format metadata 
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'sample-id'
SAMP <- sample_data(samp_df)
class(SAMP)

# Format filtered metadata
samp_df_filt <- as.data.frame(meta_filt[,-1])
rownames(samp_df_filt)<- meta_filt$'sample-id'
SAMP_filt <- sample_data(samp_df_filt)
class(SAMP_filt)

# Format taxonomy 
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)
class(TAX)


#### Create phyloseq object ####
# original phyloseq object
ryan_ps <- phyloseq(OTU, SAMP, TAX, phylotree)

# Filtered phyloseq object 
ryan_ps_filt <- phyloseq(OTU, SAMP_filt, TAX, phylotree)

#### Filtering ####
# Filter out the mitochondrial and chloroplast in the sample 
ryan_filt <- subset_taxa(ryan_ps,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
ryan_filt_new <- subset_taxa(ryan_ps_filt,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

#### Rarefaction and Generating Rarefaction Curve####
# generating rarefaction curve
rarefaction_curve <- rarecurve(t(as.data.frame(otu_table(ryan_filt))), cex=0.1)
print(rarefaction_curve)

# Rarefaction 
ryan_rare <- rarefy_even_depth(ryan_filt, rngseed = 1, sample.size = 26522)

#### Saving ####
save(ryan_filt, file="ryan_filt.RData")
save(ryan_rare, file="ryan_rare.RData")

#### Alpha Diversity ####
# Load data 
load("ryan_rare.RData")
load("ryan_filt.RData")

#Alpha diversity analysis 
plot_richness(ryan_rare) 

# Observed features and Shannon 
plot_richness(ryan_rare, measures = c("Shannon","Observed")) 

gg_richness <- plot_richness(ryan_rare, x = "Condition", measures = c("Shannon","Observed")) +
  xlab("Subject ID") +
  geom_boxplot()
gg_richness

ggsave(filename = "plot_richness.png"
       , gg_richness
       , height=4, width=6)

# Faith's Phylogenetic Diversity 
phylo_dist <- pd(t(otu_table(ryan_rare)), phy_tree(ryan_rare),
                 include.root=F) 

sample_data(ryan_rare)$PD <- phylo_dist$PD

plot.pd <- ggplot(sample_data(ryan_rare), aes(Condition, PD)) + 
  geom_boxplot() +
  xlab("Subject ID") +
  ylab("Phylogenetic Diversity")

