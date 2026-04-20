####Phyloseq creation####
# Load all relevant packages
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(picante)


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
# generating rarefaction curve, CANNOT run the alpha rarefaction curve
# rarefaction_curve <- rarecurve(t(as.data.frame(otu_table(ryan_filt_new))), cex=0.1)
# print(rarefaction_curve)

# Rarefaction 
ryan_rare_new <- rarefy_even_depth(ryan_filt_new, rngseed = 1, sample.size = 26522)

#### Saving ####
save(ryan_filt, file="ryan_filt.RData")
save(ryan_rare, file="ryan_rare.RData")

#### Alpha diversity analysis####
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

#### Beta diversity analysis####
# Load all relevant packages
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#### Load data####
# Load phloseq object
load("ryan_rare_new.RData")

#### Beta diversity #####
# Run beta-diversity using bray-curtis as the metric
bc_dm <- distance(ryan_rare_new, method="bray")
pcoa_bc <- ordinate(ryan_rare_new, method="PCoA", distance=bc_dm)

#PERMANOVA for bray-curtis
meta <- data.frame(sample_data(ryan_rare_new))
perm_test <- adonis2(
  bc_dm ~ Condition,
  data = meta,
  by = "terms")
p_val <- perm_test$`Pr(>F)`[1]
f_val  <- perm_test$F[1] 
stat_label <- paste0("PERMANOVA: F = ", round(f_val, 3), ", Pr(>F) = ", p_val)

#Generate PcoA plot for bray-curtis
# Set factor levels FIRST (controls legend order)
sample_data(ryan_rare_new)$Condition <- factor(
  sample_data(ryan_rare_new)$Condition,
  levels = c("Healthy", "Crohn's Disease")
)

# Extract variance explained
var_exp <- pcoa_bc$values$Relative_eig * 100

# Plot
gg_pcoa_bc <- plot_ordination(ryan_rare_new, pcoa_bc, color = "Condition") +
  geom_point(size = 2) +
  stat_ellipse(aes(color = Condition), type = "t", level = 0.95) +
  
  # Custom colors
  scale_color_manual(values = c(
    "Healthy" = "#7570B3",
    "Crohn's Disease" = "#D55E00"
  )) +
  
  # Axis labels with variance explained
  labs(
    x = paste0("PCo1 [", sprintf("%.1f", var_exp[1]), "%]"),
    y = paste0("PCo2 [", sprintf("%.1f", var_exp[2]), "%]"),
    color = "Condition"
  ) +
  
  # Theme
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12)
  ) +
  
  # Annotation
  annotate("text", x = Inf, y = -Inf, label = stat_label, 
           hjust = 1.1, vjust = -0.35, size = 5, fontface = "bold")

gg_pcoa_bc

#Save PcoA plot for bray-curtis
ggsave("plot_pcoa_bc.png",
       gg_pcoa_bc,
       height = 6,
       width = 10)


# Run beta-diversity using jaccard as the metric
j_dm <- distance(ryan_rare_new, method = "jaccard")
pcoa_j <- ordinate(ryan_rare_new, method="PCoA", distance=j_dm)

#PERMANOVA for jaccard
perm_test_j <- adonis2(
  j_dm ~ Condition,
  data = meta,
  by = "terms")
p_val_j <- perm_test_j$`Pr(>F)`[1]
f_val_j  <- perm_test_j$F[1]
stat_label_j <- paste0("PERMANOVA: F = ", round(f_val_j, 3), ", Pr(>F) = ", p_val_j)

# Set factor levels FIRST (controls legend order)
sample_data(ryan_rare_new)$Condition <- factor(
  sample_data(ryan_rare_new)$Condition,
  levels = c("Healthy", "Crohn's Disease")
)

# Extract variance explained
var_exp_j <- pcoa_j$values$Relative_eig * 100

# Plot
gg_pcoa_j <- plot_ordination(ryan_rare_new, pcoa_j, color = "Condition") +
  geom_point(size = 2) +
  stat_ellipse(aes(color = Condition), type = "t", level = 0.95) +
  
  # Custom colors
  scale_color_manual(values = c(
    "Healthy" = "#7570B3",
    "Crohn's Disease" = "#D55E00"
  )) +
  
  # Axis labels with variance explained
  labs(
    x = paste0("PCo1 [", sprintf("%.1f", var_exp_j[1]), "%]"),
    y = paste0("PCo2 [", sprintf("%.1f", var_exp_j[2]), "%]"),
    color = "Condition"
  ) +
  
  # Theme
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12)
  ) +
  
  # Annotation
  annotate("text", x = Inf, y = -Inf, label = stat_label_j, 
           hjust = 1.1, vjust = -0.35, size = 5, fontface = "bold")

gg_pcoa_j

#Save PcoA plot for jaccard
ggsave("plot_gg_pcoa_j.png",
       gg_pcoa_j,
       height = 6,
       width = 10)


# Run beta-diversity using unweighted unifrac as the metric
uu_dm <- distance(ryan_rare_new, method = "uunifrac")
pcoa_uu <- ordinate(ryan_rare_new, method="PCoA", distance=uu_dm)

#PERMANOVA for unweighted unifrac
perm_test_uu <- adonis2(
  uu_dm ~ Condition,
  data = meta,
  by = "terms")
p_val_uu <- perm_test_uu$`Pr(>F)`[1]
f_val_uu <- perm_test_uu$F[1]
stat_label_uu <- paste0("PERMANOVA: F = ", round(f_val_uu, 3), ", Pr(>F) = ", p_val_uu)


# Set factor levels FIRST (controls legend order)
sample_data(ryan_rare_new)$Condition <- factor(
  sample_data(ryan_rare_new)$Condition,
  levels = c("Healthy", "Crohn's Disease")
)

# Extract variance explained
var_exp_uu<- pcoa_uu$values$Relative_eig * 100

# Plot
gg_pcoa_uu<- plot_ordination(ryan_rare_new, pcoa_uu, color = "Condition") +
  geom_point(size = 2) +
  stat_ellipse(aes(color = Condition), type = "t", level = 0.95) +
  
  # Custom colors
  scale_color_manual(values = c(
    "Healthy" = "#7570B3",
    "Crohn's Disease" = "#D55E00"
  )) +
  
  # Axis labels with variance explained
  labs(
    x = paste0("PCo1 [", sprintf("%.1f", var_exp_uu[1]), "%]"),
    y = paste0("PCo2 [", sprintf("%.1f", var_exp_uu[2]), "%]"),
    color = "Condition"
  ) +
  
  # Theme
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12)
  ) +
  
  # Annotation
  annotate("text", x = Inf, y = -Inf, label = stat_label_uu, 
           hjust = 1.1, vjust = -0.35, size = 5, fontface = "bold")


gg_pcoa_uu

#Save PcoA plot for unweighted unifrac
ggsave("plot_gg_pcoa_uu.png",
       gg_pcoa_uu,
       height = 6,
       width = 10)

# Run beta-diversity using weighted unifrac as the metric
wu_dm <- distance(ryan_rare_new, method = "wunifrac")
pcoa_wu <- ordinate(ryan_rare_new, method="PCoA", distance=wu_dm)

#PERMANOVA for weighted unifrac
perm_test_wu <- adonis2(
  wu_dm ~ Condition,
  data = meta,
  by = "terms")
p_val_wu <- perm_test_wu$`Pr(>F)`[1]
f_val_wu <- perm_test_wu$F[1]
stat_label_wu <- paste0("PERMANOVA: F = ", round(f_val_wu, 3), ", Pr(>F) = ", p_val_wu)

# Set factor levels FIRST (controls legend order)
sample_data(ryan_rare_new)$Condition <- factor(
  sample_data(ryan_rare_new)$Condition,
  levels = c("Healthy", "Crohn's Disease")
)

# Extract variance explained
var_exp_wu<- pcoa_wu$values$Relative_eig * 100

# Plot
gg_pcoa_wu<- plot_ordination(ryan_rare_new, pcoa_wu, color = "Condition") +
  geom_point(size = 2) +
  stat_ellipse(aes(color = Condition), type = "t", level = 0.95) +
  
  # Custom colors
  scale_color_manual(values = c(
    "Healthy" = "#7570B3",
    "Crohn's Disease" = "#D55E00"
  )) +
  
  # Axis labels with variance explained
  labs(
    x = paste0("PCo1 [", sprintf("%.1f", var_exp_wu[1]), "%]"),
    y = paste0("PCo2 [", sprintf("%.1f", var_exp_wu[2]), "%]"),
    color = "Condition"
  ) +
  
  # Theme
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12)
  ) +
  
  # Annotation
  annotate("text", x = Inf, y = -Inf, label = stat_label_wu, 
           hjust = 1.1, vjust = -0.35, size = 5, fontface = "bold")

gg_pcoa_wu

#Save PcoA plot for weighted unifrac
ggsave("plot_gg_pcoa_wu.png",
       gg_pcoa_wu,
       height = 6,
       width = 10)

#### Taxonomy bar plots ####

# Plot bar plot of taxonomy. Group it by phylum
plot_bar(ryan_rare_new, fill="Phylum") 

# Convert to relative abundance.
ryan_RA <- transform_sample_counts(ryan_rare_new, function(x) x/sum(x))

# To remove black bars, "glom" by phylum first. We don't want to remove NAs
ryan_phylum <- tax_glom(ryan_RA, taxrank = "Phylum", NArm=FALSE)

# Create a bar plot based on phylum
df <- psmelt(ryan_phylum)

gg_taxa <- ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(
    stat = "identity",
    position = "stack",
    color = NA
  ) +
  facet_wrap(~Condition, scales = "free_x") +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0%", "25%", "50%", "75%", "100%")
  ) +
  scale_fill_viridis_d(option = "D") +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "white", colour = NA)
  ) +
  labs(y = "Relative Abundance(%)")
gg_taxa

ggsave("plot_taxonomy.png"
       , gg_taxa
       , height=8, width =12)
#### ISA####
library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
# Load the nonrarefied object
load("ryan_filt.RData")

#### Indicator Species/Taxa Analysis ####
# glom to Genus
# Group data based on particular taxanomic rank (genus)
ryan_genus <- tax_glom(ryan_filt, "Genus", NArm = FALSE)
# Convert counts to relative abundance
ryan_genus_RA <- transform_sample_counts(ryan_genus, fun=function(x) x/sum(x))

#ISA
# Calls th OTU table and transpose it
# Calculated that IVs for all ASVs
isa_ryan <- multipatt(t(otu_table(ryan_genus_RA)), cluster = sample_data(ryan_genus_RA)$`Condition`)
# Anything less than 0.05 will come up as indicator species
summary(isa_ryan)
# Create a taxanomy table, but ASVs are not the row names (separte column)
taxtable <- tax_table(ryan_filt) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# consider that your table is only going to be resolved up to the genus level, be wary of anything beyond the glomed taxa level
# Combine the taxanomic info to the list of indicator species
isa_table_all <- isa_ryan$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05)
View(isa_table_all)

write.csv(isa_table_all, "ISA_significant_taxa_all.csv", row.names = FALSE)

#### DeSeq####
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(Biostrings)

load("Ryan_filt_new.RData")
sequences <- readDNAStringSet("dna-sequences.fasta")
## NOTE: If you get a zeros error, then you need to add '1' count to all reads
# Convert taxonomy table to matrix
tax_mat <- as(tax_table(ryan_filt_new), "matrix")

tax_mat["f9283780d796ffcc12a283e56b806087", "Genus"] <- "g__Ruminococcus"
tax_mat["e96331d4191c5a32f8d1347fc3d295e2", "Genus"] <- "g__Harryflintia"
tax_mat["97351bb5e19f0e3d19bab552b7055653", "Genus"] <- "g__Faecalimonas"
tax_mat["f949c3d6392f7d23a7cbdb9e8f7180cd", "Genus"] <- "g__Coprococcus"
tax_mat["8c2297f07d90de0a6b60701c8ea376cb", "Genus"] <- "g__Anaerostipes"
tax_mat["8c5c1e33bdd2f244342f2acd6e5db57e", "Genus"] <- "g__Lachnoclostridium"
tax_mat[tax_mat[, "Genus"] == "g__UCG-005", "Genus"] <- "g__Faecousia"
tax_mat[tax_mat[, "Genus"] == "g__GCA-900066575", "Genus"] <- "g__Laedolimicola"

tax_table(ryan_filt_new) <- tax_table(tax_mat)
ryan_plus1 <- transform_sample_counts(ryan_filt_new, function(x) x+1)
ryan_deseq <- phyloseq_to_deseq2(ryan_plus1, ~`Condition`)
DESEQ_ryan <- DESeq(ryan_deseq)

# results for each condition against healthy
res_CD <- results(DESEQ_ryan, tidy=TRUE, 
                  #Healthy is the reference group
                  contrast = c("Condition","Crohn's Disease","Healthy"))
view(res_CD)

# Look at results 

# CD results 
## Volcano plot: effect size VS significance 
ggplot(res_CD) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
vol_plot_CD <- res_CD %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>3) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

vol_plot_CD
# Save file as png
ggsave(filename="vol_plot_CD.png",vol_plot_CD)

# To get table of results
sigASVs_CD <- res_CD %>% 
  filter(padj<0.01 & abs(log2FoldChange)>3) %>%
  dplyr::rename(ASV=row)
View(sigASVs_CD)

# Get only asv names
sigASVs_vec_CD <- sigASVs_CD %>%
  pull(ASV)

# Prune phyloseq file - CD (for bar plot) 
ryan_DESeq_CD <- prune_taxa(sigASVs_vec_CD,ryan_filt_new)
sigASVs_CD <- tax_table(ryan_DESeq_CD) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_CD) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

tax_table(ryan_DESeq_CD) %>% 
  as.data.frame() %>% 
  pull(Genus)

# Get ASVs where Genus is NA (only the ones that appear in DESeq)in order to BLAST them
na_asvs <- tax_table(ryan_DESeq_CD) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  filter(is.na(Genus) | Genus == "") %>%
  pull(ASV)

na_asvs

# Get ASVs where Genus is named in a different way in order to BLAST them
odd_asvs <- tax_table(ryan_DESeq_CD) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  filter(is.na(Genus) | Genus == "g__UCG-005", Genus == "g__GCA-900066575") %>%
  pull(ASV)

# Get their sequences from fasta file
na_seqs <- sequences[na_asvs]

# Export as fasta file for BLAST
writeXStringSet(na_seqs, filepath = "na_taxa_DESeq_sequences.fasta")


ggplot(sigASVs_CD) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

sigASVs_CD <- sigASVs_CD %>%
  mutate(fill_colour = case_when(
    grepl("Turicibacter", Genus) ~ "#D55E00",
    grepl("Laedolimicola", Genus) ~ "#D55E00",
    grepl("Faecousia", Genus) ~ "#D55E00",
    grepl("Fournierella", Genus) ~ "#D55E00",
    grepl("Howardella", Genus) ~ "#D55E00",
    grepl("Akkermansia", Genus) ~ "#D55E00",
    TRUE ~ "grey50"
  ))
DESeq_bar_plot <- ggplot(sigASVs_CD) +
  geom_bar(aes(x = Genus, y = log2FoldChange, fill = fill_colour), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

DESeq_bar_plot 
ggsave(filename="DESeq_bar_plot .png",DESeq_bar_plot )

#### Random Forest ####
