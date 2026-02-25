# Beta-diversity
#Run beta-diversity using bray-curtis as the metric

bc_dm <- distance(ryan_rare, method="bray")

pcoa_bc <- ordinate(ryan_rare, method="PCoA", distance=bc_dm)

gg_pcoa_bc <- plot_ordination(ryan_rare, pcoa_bc, color = "Condition") +
  labs(col = "Condition")
gg_pcoa_bc

#Run beta-diversity using jaccard as the metric

j_dm <- distance(ryan_rare, method="jaccard")

pcoa_j <- ordinate(ryan_rare, method="PCoA", distance=j_dm)

gg_pcoa_j <- plot_ordination(ryan_rare, pcoa_j, color = "Condition") +
  labs(col = "Condition")
gg_pcoa_j

#Run beta-diversity using unweighted unifrac as the metric

uu_dm <- distance(ryan_rare, method="uunifrac")

pcoa_uu <- ordinate(ryan_rare, method="PCoA", distance=uu_dm)

gg_pcoa_uu <- plot_ordination(ryan_rare, pcoa_uu, color = "Condition") +
  labs(col = "Condition")
gg_pcoa_uu

#Run beta-diversity using weighted unifrac as the metric

wu_dm <- distance(ryan_rare, method="wunifrac")

pcoa_wu <- ordinate(ryan_rare, method="PCoA", distance=wu_dm)

gg_pcoa_wu <- plot_ordination(ryan_rare, pcoa_wu, color = "Condition") +
  labs(col = "Condition")
gg_pcoa_wu

ggsave("plot_pcoa_bc.png"
       , gg_pcoa_bc
       , height=4, width=6)

# Taxonomy bar plots

#Plot bar plot of taxonomy. Group it by phylum

plot_bar(ryan_rare, fill="Phylum") 

#Convert to relative abundance.

ryan_RA <- transform_sample_counts(ryan_rare, function(x) x/sum(x))

#To remove black bars, "glom" by phylum first. We don't want to remove NAs

ryan_phylum <- tax_glom(ryan_RA, taxrank = "Phylum", NArm=FALSE)

#Create a bar plot based on phylum

gg_taxa <- plot_bar(ryan_phylum, fill="Phylum") + 
  facet_wrap(.~Condition, scales = "free_x")
gg_taxa

ggsave("plot_taxonomy.png"
       , gg_taxa
       , height=8, width =12)
