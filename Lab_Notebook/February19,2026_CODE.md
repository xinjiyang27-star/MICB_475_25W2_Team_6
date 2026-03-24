# Paired Reads: Assign Taxonomy, Filtering, and Alpha Rarefaction [19 February]
## Purpose
To analyze the taxonomy, filter the sequence and do the alpha rarefaction 

## Material 
MICB475 group server 

## Method
### Paired end data analysis 
1. Use trained classifier to create taxonomy file

qiime feature-classifier classify-sklearn \
--i-classifier ryan-classifier.qza \
--i-reads ryan-rep-seqs.qza \
--o-classification ryan-taxonomy.qza

2. convert to a visualizable file

qiime metadata tabulate \
  --m-input-file ryan-taxonomy.qza \
  --o-visualization ryan-taxonomy.qzv

3. create taxa barplot

qiime taxa barplot \
  --i-table ryan-table.qza \
  --i-taxonomy ryan-taxonomy.qza \
  --m-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv \
  --o-visualization ryan-taxa-bar-plots.qzv

4. copy files to local

scp root@10.19.139.189:/data/project_2/Version_1/ryan-taxonomy.qzv .

scp root@10.19.139.189:/data/project_2/Version_1/ryan-taxa-bar-plots.qzv .

5. Filtering
#Filtering out mitochondria and Chloroplasts

qiime taxa filter-table \
  --i-table ryan-table.qza \
  --i-taxonomy ryan-taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

6. convert to a visualizable file

qiime feature-table summarize \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --o-visualization table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv

7. Copy to local

scp root@10.19.139.189:/data/project_2/Version_1/table-no-mitochondria-no-chloroplast.qzv .

8. Alpha Rarefaction
#Generate a tree for phylogenetic diversity analyses

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ryan-rep-seqs.qza \
  --o-alignment aligned-ryan-rep-seqs.qza \
  --o-masked-alignment masked-aligned-ryan-rep-seqs.qza \
  --o-tree ryan-unrooted-tree.qza \
  --o-rooted-tree ryan-rooted-tree.qza 

9. Generate alpha rarefaction curve

qiime diversity alpha-rarefaction \
  --i-table ryan-table.qza \
  --i-phylogeny ryan-rooted-tree.qza \
  --p-max-depth 33000 \
  --m-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv \
  --o-visualization ryan-alpha-rarefaction.qzv

10. Copy to local

scp root@10.19.139.189:/data/project_2/Version_1/ryan-alpha-rarefaction.qzv .

### Forward Sequences Only
#Code was rerun using only forward sequences in order to increase the quality of reads. Poor overlap of reads was causing read quality to suffer
1. Bring manifest.tsv file to local

scp root@10.19.139.189:/datasets/project_2/human_ibd/ryan_manifest.tsv .

2. Send edited manifest.tsv file back to server

scp ryan_manifest_forward.tsv root@10.19.139.189:/data/project_2

3. Import Data
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path ./data/project_2/Version_2/ryan_manifest_forward.tsv \
  --output-path ./data/project_2/Version_2/forwardryan_demux_seqs.qza

4. Trim
qiime cutadapt trim-single \
  --i-demultiplexed-sequences forwardryan_demux_seqs.qza \
  --p-front TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG \
  --p-match-read-wildcards \
  --p-match-adapter-wildcards \
  --p-discard-untrimmed \
  --o-trimmed-sequences forwardryan_demux-trimmed.qza \
  --verbose > cutadapt-log-2.txt

5. Visualize
qiime demux summarize \
  --i-data forwardryan_demux-trimmed.qza \
  --o-visualization forwardryan_demux-trimmed.qzv

6. Send to Local 
scp root@10.19.139.189:/data/project_2/Version_2/forwardryan_demux-trimmed.qzv .


7. Determine ASVs with DADA2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs forwardryan_demux-trimmed.qza \
  --p-trim-left 0 \
  --p-trunc-len 220 \
  --o-representative-sequences forwardryan-rep-seqs.qza \
  --o-table forwardryan-table.qza \
  --o-denoising-stats forwardryan-stats.qza


qiime feature-table summarize \
  --i-table forwardryan-table.qza \
  --o-visualization forwardryan-table.qzv \
  --m-sample-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data forwardryan-rep-seqs.qza \
  --o-visualization forwardryan-rep-seqs.qzv


8. Send to Local 
scp root@10.19.139.189:/data/project_2/Version_2/forwardryan-table.qzv .
scp root@10.19.139.189:/data/project_2/Version_2/forwardryan-rep-seqs.qzv .


9. Taxonomic analysis
qiime feature-classifier extract-reads \
  --i-sequences /datasets/classifiers/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-trunc-len 220 \
  --o-reads forward-ryan-ref-seqs-trimmed.qza

10. Train classifier with your new ref-seq file
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /data/project_2/Version_1/forward-ryan-ref-seqs-trimmed.qza \
  --i-reference-taxonomy /datasets/classifiers/silva_ref_files/silva-138-99-tax.qza \
  --o-classifier forwardryan-classifier.qza

