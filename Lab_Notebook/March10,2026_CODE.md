# Re-create phyloseq object [10 March]
## Purpose
Conduct metadata filtering and recreate phyloseq object based on the filtered files

## Materials
- R Studio

## Method 
### Load all relevant packages
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(picante)

### Load metadata file 
metafp <- "ryan_export/ryan_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

### Load feature-table file 
otufp <- "ryan_export/table_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

#load taxonomy file 
taxfp <- "ryan_export/taxonomy_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

#load phylogenetic tree file 
phylotreefp <- "ryan_export/rooted_tree_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

### Filter out UC, non-inflamed CD from the metadata
meta_filt <- meta %>%
  filter(Condition == "Healthy" | 
           (Condition == "Crohn's Disease" & 
           Histological.status == "Inflamed tissue "))
### Format OUT table
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

### Format metadata 
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'sample-id'
SAMP <- sample_data(samp_df)
class(SAMP)

### Format filtered metadata
samp_df_filt <- as.data.frame(meta_filt[,-1])
rownames(samp_df_filt)<- meta_filt$'sample-id'
SAMP_filt <- sample_data(samp_df_filt)
class(SAMP_filt)

### Format taxonomy 
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)
class(TAX)

### original phyloseq object
ryan_ps <- phyloseq(OTU, SAMP, TAX, phylotree)

### Filtered phyloseq object 
ryan_ps_filt <- phyloseq(OTU, SAMP_filt, TAX, phylotree)

### Filter out the mitochondrial and chloroplast in the sample 
ryan_filt <- subset_taxa(ryan_ps,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
ryan_filt_new <- subset_taxa(ryan_ps_filt,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

