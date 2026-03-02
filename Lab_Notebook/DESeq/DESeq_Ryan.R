# BiocManager::install("DeSeq2")
library(tidyverse)
library(phyloseq)
library(DESeq2)

load("ryan_filt.RData")
## NOTE: If you get a zeros error, then you need to add '1' count to all reads
ryan_plus1 <- transform_sample_counts(ryan_filt, function(x) x+1)
ryan_deseq <- phyloseq_to_deseq2(ryan_plus1, ~`Condition`)
DESEQ_ryan <- DESeq(ryan_deseq)

# results for each condition against healthy
res_UC <- results(DESEQ_ryan, tidy=TRUE, 
               #Healthy is the reference group
               contrast = c("Condition","Ulcerative Colitis","Healthy"))
view(res_UC)

res_CD <- results(DESEQ_ryan, tidy=TRUE, 
                  #Healthy is the reference group
                  contrast = c("Condition","Crohn's Disease","Healthy"))
view(res_CD)

# Look at results 

# UC results 
## Volcano plot: effect size VS significance 
ggplot(res_UC) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
vol_plot_UC <- res_UC %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
# Save file as png
ggsave(filename="vol_plot_UC.png",vol_plot_UC)

# To get table of results
sigASVs_UC <- res_UC %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_UC)

# Get only asv names
sigASVs_vec_UC <- sigASVs_UC %>%
  pull(ASV)

# CD results 
## Volcano plot: effect size VS significance 
ggplot(res_CD) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
vol_plot_CD <- res_CD %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
# Save file as png
ggsave(filename="vol_plot_CD.png",vol_plot_CD)

# To get table of results
sigASVs_CD <- res_CD %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_CD)

# Get only asv names
sigASVs_vec_CD <- sigASVs_CD %>%
  pull(ASV)

# Prune phyloseq file - UC(for bar plot)
ryan_DESeq_UC <- prune_taxa(sigASVs_vec_UC,ryan_filt)
sigASVs_UC <- tax_table(ryan_DESeq_UC) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_UC) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(sigASVs_UC) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# Prune phyloseq file - CD (for bar plot) 
ryan_DESeq_CD <- prune_taxa(sigASVs_vec_CD,ryan_filt)
sigASVs_CD <- tax_table(ryan_DESeq_CD) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_CD) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(sigASVs_CD) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
