# Exporting, phyloseq creation and alpha diversity analysis [23 February]
## Purpose
To convert table.qza, phylogeentic tree, metadata, and taxonomy.qza to human readable files and export from QIIME2. 
To Create Phyloseq object and recreate rarefaction, alpha diversity in R. 
## Materials 
1. QIIME2 Group Server
  IP: root@10.19.139.189
  Password: Biome391
2. R Studio

## Method
### QIIME2 Portion 
1. Connect to server and go to the project working directory
2. Create an export folder by mkdir
3.convert and export table.qza into human readable file 
  qiime tools export \
  --input-path /data/project_2/Version_2/forwardryan-table.qza \
  --output-path table_export
4. Convert and export taxonomy.qza to human readable file
  qiime tools export \
  --input-path /data/project_2/Version_2/forwardryan-taxonomy.qza \
  --output-path taxonomy_export
5. Convert and export rooted tree
  qiime tools export \
  --input-path /data/project_2/Version_2/forwardryan-rooted-tree.qza \
  --output-path rooted_tree__export
6. Go to the export table directory and convert biom file to txt file.
  biom convert -i feature-table.biom --to-tsv -o feature-table.txt
7. Download the entire folder from server to local computer
  scp -r root@10.19.139.189:/data/project_2/Version_2/ryan_export .
8. Export metadata file
  scp root@10.19.139.189:/datasets/project_2/human_ibd/ryan_metadata.tsv .
### R Portion
1. Load all the relevant packages
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(picante)
2. Load data

2.1 Load metadata file 
metafp <- "ryan_export/ryan_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

2.2 Load feature-table file 
otufp <- "ryan_export/table_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

2.3 load taxonomy file 
taxfp <- "ryan_export/taxonomy_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

2.4 load phylogenetic tree file 
phylotreefp <- "ryan_export/rooted_tree_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

3. Formatting relevant files

3.1 Format OUT table
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

3.2 Format metadata 
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'sample-id'
SAMP <- sample_data(samp_df)
class(SAMP)

3.3 Format taxonomy 
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)
class(TAX)

4. Create phyloseq object
ryan_ps <- phyloseq(OTU, SAMP, TAX, phylotree)


5. Filtering
5.1 Filter out the mitochondrial and chloroplast in the sample 
ryan_filt <- subset_taxa(ryan_ps,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

6. Rarefaction and Generating Rarefaction Curve
6.1 generating rarefaction curve
rarefaction_curve <- rarecurve(t(as.data.frame(otu_table(ryan_filt))), cex=0.1)
print(rarefaction_curve)

6.2 Rarefaction 
ryan_rare <- rarefy_even_depth(ryan_filt, rngseed = 1, sample.size = 26522)

7. Saving
save(ryan_filt, file="ryan_filt.RData")
save(ryan_rare, file="ryan_rare.RData")

8. Load data for Alpha diversity analysis 
load("ryan_rare.RData")
load("ryan_filt.RData")

9. Alpha diversity analysis 
plot_richness(ryan_rare) 

10. Observed features and Shannon 
plot_richness(ryan_rare, measures = c("Shannon","Observed")) 

gg_richness <- plot_richness(ryan_rare, x = "Condition", measures = c("Shannon","Observed")) +
  xlab("Subject ID") +
  geom_boxplot()
gg_richness

ggsave(filename = "plot_richness.png"
       , gg_richness
       , height=4, width=6)

11. Faith's Phylogenetic Diversity 
phylo_dist <- pd(t(otu_table(ryan_rare)), phy_tree(ryan_rare),
                 include.root=F) 

sample_data(ryan_rare)$PD <- phylo_dist$PD

plot.pd <- ggplot(sample_data(ryan_rare), aes(Condition, PD)) + 
  geom_boxplot() +
  xlab("Subject ID") +
  ylab("Phylogenetic Diversity")


