# March 2 [16:00 - 16:50] 
## Meeting Agenda 
1. Discuss alpha and beta diversity results, taxonomy bar plot, core microbiome, ISA, and DESeq results
2. Discuss Next step on building the random forest model

## Q and A
1. What are the requirements for the final figures in the paper?

## Experimental Results Summary
**1. Rarefaction Curve**
- [Rarefaction Curve](/Lab_Notebook/Alpha_and_Beta_Diversity/Rarefaction_Curve.png)
  
**2. Alpha diversity**
- [plot richness](/Lab_Notebook/Alpha_and_Beta_Diversity/plot_richness_with_stats.png)
  - For Observed: healthy vs CD, healthy vs UC, and UC vs CD is significant.
  - For Shannon: Only healthy vs CD shows significant differences. 
- [plot pd](/Lab_Notebook/Alpha_and_Beta_Diversity/plot_pd.png)
  - For phylogenetic distances: healthy vs CD, healthy vs UC, and UC vs CD is significant.
  
**3. Beta diversity**
- [bray curtis](/Lab_Notebook/Alpha_and_Beta_Diversity/plot_pcoa_bc.png)
    - For bray curtis: healthy vs UC vs CD is significant with PC at 17.7% 
- [jaccard](/Lab_Notebook/Alpha_and_Beta_Diversity/plot_gg_pcoa_j.png)
    - For jaccard: healthy vs UC vs CD is significant with PC at 11.1%. 
- [unweighted unifrac](/Lab_Notebook/Alpha_and_Beta_Diversity/plot_gg_pcoa_uu.png)
    - For unweighted unifrac: healthy vs UC vs CD is significant with PC at 9%. 
- [weighted unifrac](/Lab_Notebook/Alpha_and_Beta_Diversity/plot_gg_pcoa_wu.png)
    - For weighted unifrac: healthy vs UC vs CD is significant with PC at 21.9%. 
- But, we don't see distinct clusters for any groups.

**5. Taxonomy bar plot**
- [Taxonomy bar plot](/Lab_Notebook/Alpha_and_Beta_Diversity/plot_taxonomy.png)

**6. Core microbiome**
- [Core microbiome](/Lab_Notebook/Core_microbiome/Ryan_Venn_Condition.png)

**7. ISA**
- [ISA table](/Lab_Notebook/Indicator_Species_Analysis/ISA_significant_taxa.csv)

**8. DESeq**
- [DESeq sig ASVs UC vs Healthy](/Lab_Notebook/DESeq/DESeq_sigASVs_UC.png)
- [DESeq sig ASVs CD vs Healthy](/Lab_Notebook/DESeq/DESeq_sigASVs_CD.png)


## Meeting Minutes
1. Proposal revisions
- Hypothesis - Can we build a better model with biopsy samples rather than fecal samples
- Proposed approach - each dataset for each location, with the contingency that it has enough sample size. 
2. Alpha & beta diversity
- Alpha: Observed
  - Decrease for both UC vs CD, sig diff between all conditions (varying amount of species)
- Alpha: Shannon
  - Not sig Healthy vs UC and CD v sUC
  - Abundance changes are not diff → DESeq might not be great
- Alpha: Faith pd
  - All sig
  - Difference in taxa is what is diff, focus on ISA and Core Microbiome, rather than abundance.
- The beta diversity isn’t a very strong analysis (the % graphed for each analysis is very little) → Can’t trust the analysis. Be transparent → say you ran the metrics, but it only represents ~3-13%
- Beta: Jaccard
  - Not much diff
  - Permanova's new versions need to fix the p-value. 
- Beta: Weighted 
  - More spread than unweighted → takes into account Shannon, so obviously why not sig (same with Bay-curtis)

- Don’t combine CD and UC, as there are differences in the alpha metrics
- Possible issue with UC is that it looks very similar to healthy → may result in a poor model
- In manuscript → Diversity metrics are closer to healthy for UC so we don’t want to focus on that, hence why we dropped it and so only focus on CD 
3. Next steps:
- Keep all conditions
- Core microbiome
  - For core micro: make a Venn diagram, but get a list of microbes that are associated with each condition
- DESeq
  - Reconcile for how many are overlapping → only want ones that are unique
- ISA
  - Ones that are shared, remove them → only want one 1 across the 3 conditions
  - Stat indicates close to 1 is better → set as 0.4 and higher
  - Mention in the paper that we had to be lenient with the cut-off, as we didn’t get too many indicators. 
  - Keep any shared between CD and UC

- Keep everything we have rn and run the analyses again, but only for CD. (Possibly because UC is localised, and so certain biopsy locations would only have healthy or UC taxa)
- Reconcile taxa with all analyses and just pick those that agree all throughout → less emphasis on DESeq. Match the core and ISA first, then check how many match with DESeq
- 10-20 taxa is a good amount to feed the model, and maybe 5 metadata columns
- Pick metadata categories we want for the model → for example, biopsy location, histology status, medications, condition, gender (?) → can always remove these if it isn’t driving anything
- If we keep everything → more data = more noise → the more you see a decrease in accuracy

## TO-DO lists
1. Decide which taxa and metadata categories to use for building machine learning model
2. Watch random forest model and try with running the code



