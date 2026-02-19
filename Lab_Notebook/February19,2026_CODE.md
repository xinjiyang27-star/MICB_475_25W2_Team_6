# Assign Taxonomy, Filtering, and Alpha Rarefaction

#Use trained classifier to create taxonomy file

qiime feature-classifier classify-sklearn \
--i-classifier ryan-classifier.qza \
--i-reads ryan-rep-seqs.qza \
--o-classification ryan-taxonomy.qza

#convert to a visualizable file

qiime metadata tabulate \
  --m-input-file ryan-taxonomy.qza \
  --o-visualization ryan-taxonomy.qzv

#create taxa barplot

qiime taxa barplot \
  --i-table ryan-table.qza \
  --i-taxonomy ryan-taxonomy.qza \
  --m-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv \
  --o-visualization ryan-taxa-bar-plots.qzv

#copy files to local

scp root@10.19.139.189:/data/project_2/Version_1/ryan-taxonomy.qzv .

scp root@10.19.139.189:/data/project_2/Version_1/ryan-taxa-bar-plots.qzv .

## Filtering
#Filtering out mitochondria and Chloroplasts

qiime taxa filter-table \
  --i-table ryan-table.qza \
  --i-taxonomy ryan-taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

#convert to a visualizable file

qiime feature-table summarize \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --o-visualization table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv

#Copy to local

scp root@10.19.139.189:/data/project_2/Version_1/table-no-mitochondria-no-chloroplast.qzv .

## Alpha Rarefaction
#Generate a tree for phylogenetic diversity analyses

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ryan-rep-seqs.qza \
  --o-alignment aligned-ryan-rep-seqs.qza \
  --o-masked-alignment masked-aligned-ryan-rep-seqs.qza \
  --o-tree ryan-unrooted-tree.qza \
  --o-rooted-tree ryan-rooted-tree.qza 

#Generate alpha rarefaction curve

qiime diversity alpha-rarefaction \
  --i-table ryan-table.qza \
  --i-phylogeny ryan-rooted-tree.qza \
  --p-max-depth 33000 \
  --m-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv \
  --o-visualization ryan-alpha-rarefaction.qzv

#Copy to local

scp root@10.19.139.189:/data/project_2/Version_1/ryan-alpha-rarefaction.qzv .

