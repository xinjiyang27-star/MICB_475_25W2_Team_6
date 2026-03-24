# Taxonomic analysis [14 February]
qiime feature-classifier extract-reads \
  --i-sequences /datasets/classifiers/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG \
  --p-r-primer GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC \
  --p-trunc-len 220 \
  --o-reads forward-ryan-ref-seqs-trimmed.qza

qiime feature-classifier extract-reads \
  --i-sequences /datasets/classifiers/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG \
  --p-r-primer GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC \
  --p-trunc-len 180 \
  --o-reads reverse-ryan-ref-seqs-trimmed.qza

qiime feature-table merge-seqs \
  --i-data forward-ryan-ref-seqs-trimmed.qza \
  --i-data reverse-ryan-ref-seqs-trimmed.qza \
  --o-merged-data ryan-ref-seqs.qza

# Train classifier with your new ref-seq file
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ryan-ref-seqs.qza \
  --i-reference-taxonomy /datasets/classifiers/silva_ref_files/silva-138-99-tax.qza \
  --o-classifier classifier.qza

# Use the trained classifier to assign taxonomy to your reads (rep-seqs.qza)
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads ryan-rep-seqs.qza \
  --o-classification ryan-taxonomy.qza

# Visualization 
qiime metadata tabulate \
  --m-input-file ryan-taxonomy.qza \
  --o-visualization ryan-taxonomy.qzv
  
# Taxonomy barplots
qiime taxa barplot \
  --i-table ryan-table.qza \
  --i-taxonomy ryan-taxonomy.qza \
  --m-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv \
  --o-visualization ryan-taxa-bar-plots.qzv

# Local 
scp root@10.19.139.189:/data/project_2/ryan-taxonomy.qzv .
scp root@10.19.139.189:/data/project_2/ryan-taxa-bar-plots.qzv .
