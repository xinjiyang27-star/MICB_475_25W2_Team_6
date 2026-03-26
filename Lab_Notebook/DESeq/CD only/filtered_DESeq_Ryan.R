# BiocManager::install("DeSeq2")
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(Biostrings)

load("Ryan_filt_new.RData")

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
  mutate(Genus = gsub("g__", "", Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

tax_table(ryan_DESeq_CD) %>% 
  as.data.frame() %>% 
  pull(Genus)

### the following code is for BLAST purposes- Don't run when doing DESeq runs again###

sequences <- readDNAStringSet("dna-sequences.fasta")
# Get ASVs where Genus is NA (only the ones that appear in DESeq)in order to BLAST them
na_asvs <- tax_table(ryan_DESeq_CD) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  filter(is.na(Genus) | Genus == "") %>%
  pull(ASV)

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

### End of BLAST-related code ###

ggplot(sigASVs_CD) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=0.5, vjust=0.5))

sigASVs_CD <- sigASVs_CD %>%
  mutate(fill_colour = case_when(
    grepl("Turicibacter", Genus) ~ "#ee55cc",
    grepl("Laedolimicola", Genus) ~ "#ee55cc",
    grepl("Faecousia", Genus) ~ "#ee55cc",
    grepl("Fournierella", Genus) ~ "#ee55cc",
    grepl("Howardella", Genus) ~ "#ee55cc",
    grepl("Akkermansia", Genus) ~ "#ee55cc",
    TRUE ~ "grey50"
  ))
DESeq_bar_plot <- ggplot(sigASVs_CD) +
  geom_bar(aes(x = Genus, y = log2FoldChange, fill = fill_colour), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  scale_fill_identity() +
  theme_classic(base_size=10) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
DESeq_bar_plot 
ggsave(filename="DESeq_bar_plot .png",DESeq_bar_plot )
