# Load all relevant packages
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

####Load data####
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
ryan_ps <- phyloseq(OTU, SAMP, TAX, phylotree)

#### Filtering ####
# Filter out the mitochondrial and chloroplast in the sample 
ryan_filt <- subset_taxa(ryan_ps,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

#### Rarefaction and Generating Rarefaction Curve####
# generating rarefaction curve
rarefaction_curve <- rarecurve(t(as.data.frame(otu_table(ryan_filt))), cex=0.1)
print(rarefaction_curve)

# Rarefaction 
ryan_rare <- rarefy_even_depth(ryan_filt, rngseed = 1, sample.size = 26522)

# Save the filtered and rarefied graphs 
save(ryan_filt, file="ryan_filt.RData")
save(ryan_rare, file="ryan_rare.RData")

#### Alpha Diversity ####
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggpubr)

# Load data 
load("ryan_rare.RData")
load("ryan_filt.RData")


# Find variables for condition and format the comparison list
unique(get_variable(ryan_rare, "Condition"))

alpha_comparisons <- list( c("Healthy", "Ulcerative Colitis"), 
                        c("Ulcerative Colitis", "Crohn's Disease"), 
                        c("Healthy", "Crohn's Disease") )

# Run Shannon and Observed features and generate alpha diversity plot 
gg_richness_alpha_with_stats <- plot_richness(ryan_rare, x = "Condition", 
                                              measures = c("Shannon", "Observed")) +
  geom_boxplot() +
  xlab("Condition") +
  stat_compare_means(method = "kruskal.test", label.y = 0) +    # Add Global Kruskal-Wallis p-value
  stat_compare_means(comparisons = alpha_comparisons,   #  Add Pairwise Wilcoxon p-values with brackets
                     method = "wilcox.test", 
                     p.adjust.method = "BH")

gg_richness_alpha_with_stats

# Save Shannon and Observed features plot
ggsave("plot_richness_with_stats.png",
       gg_richness_alpha_with_stats,
       height = 6,
       width = 8)

# Run Faith’s Phylogenetic Diversity
phylo_dist <- pd(t(otu_table(ryan_rare)),
                 phy_tree(ryan_rare),
                 include.root = FALSE)

sample_data(ryan_rare)$PD <- phylo_dist$PD

pd_df <- data.frame(
  PD = phylo_dist$PD,
  Condition = sample_data(ryan_rare)$Condition)

# Generate plot for Phylogenetic Diversity
plot_pd <- ggplot(pd_df, aes(x = Condition, y = PD)) +
  geom_boxplot() +
  xlab("Condition") +
  ylab("Faith's Phylogenetic Diversity") +
  stat_compare_means(method = "kruskal.test", label.y = 0) + 
  stat_compare_means(comparisons = alpha_comparisons, 
                     method = "wilcox.test", 
                     p.adjust.method = "BH")
plot_pd

# Save Phylogenetic Diversity plot
ggsave("plot_pd.png",
       plot_pd,
       height = 6,
       width = 8)

#### Beta diversity #####
# Run beta-diversity using bray-curtis as the metric
bc_dm <- distance(ryan_rare, method="bray")
pcoa_bc <- ordinate(ryan_rare, method="PCoA", distance=bc_dm)

# Run PERMANOVA for bray-curtis
perm_test <- adonis2(bc_dm ~ Condition, data = as(sample_data(ryan_rare), "data.frame"))
p_val <- perm_test$`Pr(>F)`[1]
f_val  <- perm_test$F[1]
stat_label <- paste0("PERMANOVA: F = ", round(f_val, 3), ", Pr(>F) = ", p_val)

# Generate PcoA plot for bray-curtis
gg_pcoa_bc_elipse <- plot_ordination(ryan_rare, pcoa_bc, color = "Condition") +
  stat_ellipse(aes(color = Condition), type = "t", level = 0.95) +
  labs(col = "Condition") +
  annotate("text", x = Inf, y = -Inf, label = stat_label, 
           hjust = 1.1, vjust = -0.5, size = 4, fontface = "bold")

gg_pcoa_bc_elipse

# Save PcoA plot for bray-curtis
ggsave("plot_pcoa_bc_eclipse.png",
       gg_pcoa_bc_elipse,
       height = 4,
       width = 6)


# Run beta-diversity using jaccard as the metric
j_dm <- distance(ryan_rare, method = "jaccard")
pcoa_j <- ordinate(ryan_rare, method="PCoA", distance=j_dm)

# Run PERMANOVA for jaccard
perm_test_j <- adonis2(j_dm ~ Condition, data = as(sample_data(ryan_rare), "data.frame"))
p_val_j <- perm_test$`Pr(>F)`[1]
f_val_j  <- perm_test$F[1]
stat_label_j <- paste0("PERMANOVA: F = ", round(f_val, 3), ", Pr(>F) = ", p_val)

# Generate PcoA plot for jaccard
gg_pcoa_j <- plot_ordination(ryan_rare, pcoa_j, color = "Condition") +
  stat_ellipse(aes(color = Condition), type = "t", level = 0.95) +
  labs(col = "Condition") +
  annotate("text", x = Inf, y = -Inf, label = stat_label, 
           hjust = 1.1, vjust = -0.5, size = 4, fontface = "bold")

gg_pcoa_j

# Save PcoA plot for jaccard
ggsave("plot_gg_pcoa_j.png",
       gg_pcoa_bc_elipse,
       height = 4,
       width = 6)

# Run beta-diversity using unweighted unifrac as the metric
uu_dm <- distance(ryan_rare, method = "uunifrac")
pcoa_uu <- ordinate(ryan_rare, method="PCoA", distance=uu_dm)

# Run PERMANOVA for unweighted unifrac
perm_test_uu <- adonis2(uu_dm ~ Condition, data = as(sample_data(ryan_rare), "data.frame"))
p_val_uu <- perm_test$`Pr(>F)`[1]
f_val_uu <- perm_test$F[1]
stat_label_uu <- paste0("PERMANOVA: F = ", round(f_val, 3), ", Pr(>F) = ", p_val)

# Generate PcoA plot for unweighted unifrac
gg_pcoa_uu <- plot_ordination(ryan_rare, pcoa_uu, color = "Condition") +
  stat_ellipse(aes(color = Condition), type = "t", level = 0.95) +
  labs(col = "Condition") +
  annotate("text", x = Inf, y = -Inf, label = stat_label, 
           hjust = 1.1, vjust = -0.5, size = 4, fontface = "bold")

gg_pcoa_uu

# Save PcoA plot for unweighted unifrac
ggsave("plot_gg_pcoa_uu.png",
       gg_pcoa_bc_elipse,
       height = 4,
       width = 6)
       
# Run beta-diversity using weighted unifrac as the metric
wu_dm <- distance(ryan_rare, method = "wunifrac")
pcoa_wu <- ordinate(ryan_rare, method="PCoA", distance=wu_dm)

# Run PERMANOVA for weighted unifrac
perm_test_wu <- adonis2(wu_dm ~ Condition, data = as(sample_data(ryan_rare), "data.frame"))
p_val_wu <- perm_test$`Pr(>F)`[1]
f_val_wu <- perm_test$F[1]
stat_label_wu <- paste0("PERMANOVA: F = ", round(f_val, 3), ", Pr(>F) = ", p_val)

# Generate PcoA plot for weighted unifrac
gg_pcoa_wu <- plot_ordination(ryan_rare, pcoa_wu, color = "Condition") +
  stat_ellipse(aes(color = Condition), type = "t", level = 0.95) +
  labs(col = "Condition") +
  annotate("text", x = Inf, y = -Inf, label = stat_label, 
           hjust = 1.1, vjust = -0.5, size = 4, fontface = "bold")

gg_pcoa_uu

# Save PcoA plot for weighted unifrac
ggsave("plot_gg_pcoa_wu.png",
       gg_pcoa_bc_elipse,
       height = 4,
       width = 6)

#### Taxonomy bar plots ####
# Plot bar plot of taxonomy. Group it by phylum
plot_bar(ryan_rare, fill="Phylum") 

# Convert to relative abundance.
ryan_RA <- transform_sample_counts(ryan_rare, function(x) x/sum(x))

# To remove black bars, "glom" by phylum first. We don't want to remove NAs
ryan_phylum <- tax_glom(ryan_RA, taxrank = "Phylum", NArm=FALSE)

# Create a bar plot based on phylum
gg_taxa <- plot_bar(ryan_phylum, fill="Phylum") + 
  facet_wrap(.~Condition, scales = "free_x")
gg_taxa

# Save the taxonomy bar plot
ggsave("plot_taxonomy.png"
       , gg_taxa
       , height=8, width =12)
