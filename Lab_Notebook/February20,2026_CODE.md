# Forward Reads: Assign Taxonomy, Filtering, and Alpha Rarefaction [20 February]
## Purpose
To analyze the taxonomy, filter the ASVs and rarefy the samples by using the forward reads only

## Materials
MICB475 group server

## Methods 
1. Use trained classifier to create taxonomy file

qiime feature-classifier classify-sklearn \
--i-classifier forwardryan-classifier.qza \
--i-reads forwardryan-rep-seqs.qza \
--o-classification forwardryan-taxonomy.qza

2. convert to a visualizable file

qiime metadata tabulate \
  --m-input-file forwardryan-taxonomy.qza \
  --o-visualization forwardryan-taxonomy.qzv

3. create taxa barplot

qiime taxa barplot \
  --i-table forwardryan-table.qza \
  --i-taxonomy forwardryan-taxonomy.qza \
  --m-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv \
  --o-visualization forwardryan-taxa-bar-plots.qzv

4. copy files to local

scp root@10.19.139.189:/data/project_2/Version_2/forwardryan-taxonomy.qzv .

scp root@10.19.139.189:/data/project_2/Version_2/forwardryan-taxa-bar-plots.qzv .

5. Filtering
Filtering out mitochondria and Chloroplasts

qiime taxa filter-table \
  --i-table forwardryan-table.qza \
  --i-taxonomy forwardryan-taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table forwardtable-no-mitochondria-no-chloroplast.qza

6. convert to visualizable file

qiime feature-table summarize \
  --i-table forwardtable-no-mitochondria-no-chloroplast.qza \
  --o-visualization forwardtable-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv

7. Copy to local

scp root@10.19.139.189:/data/project_2/Version_2/forwardtable-no-mitochondria-no-chloroplast.qzv .

8. Alpha Rarefaction
Generate a tree for phylogenetic diversity analyses

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences forwardryan-rep-seqs.qza \
  --o-alignment aligned-forwardryan-rep-seqs.qza \
  --o-masked-alignment masked-aligned-forwardryan-rep-seqs.qza \
  --o-tree forwardryan-unrooted-tree.qza \
  --o-rooted-tree forwardryan-rooted-tree.qza 

9. Alpha rarefaction
qiime diversity alpha-rarefaction \
  --i-table forwardryan-table.qza \
  --i-phylogeny forwardryan-rooted-tree.qza \
  --p-max-depth 120000 \
  --m-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv \
  --o-visualization forward120000ryan-alpha-rarefaction.qzv

10. Copy to local

scp root@10.19.139.189:/data/project_2/Version_2/forward120000ryan-alpha-rarefaction.qzv .

11. Copy to local

scp root@10.19.139.189:/data/project_2/Version_1/ryan-alpha-rarefaction.qzv .
