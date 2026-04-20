# Demultiplexing
#import the data from project 2 datasets
   qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path ./datasets/project_2/human_ibd/ryan_manifest.tsv \
  --output-path ./data/project_2/ryan_demux_seqs.qza

#Trim the adaptor within the reads (codes are provided by Hans)
   qiime cutadapt trim-paired \
 --i-demultiplexed-sequences ryan_demux_seqs.qza \
 --p-front-f TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG \
 --p-front-r GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC \
 --p-match-read-wildcards \
 --p-match-adapter-wildcards \
 --p-discard-untrimmed \
 --o-trimmed-sequences ryan_demux_seqs-trimmed.qza \
 --verbose > cutadapt-log-2.txt

#Create visualization of demultiplexed samples
   qiime demux summarize \
  --i-data ryan_demux_seqs-trimmed.qza \
  --o-visualization ryan_demux_seqs-trimmed.qzv
# Denosing
#Denoising and clustering 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ryan_demux_seqs-trimmed.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 180 \
  --o-representative-sequences ryan-rep-seqs.qza \
  --o-table ryan-table.qza \
  --o-denoising-stats ryan-stats.qza
#Ensure overlap is: forward length after truncation + reverse length after truncation - amplicon length ≥ 20–30 bp overlap

#Visualize Files
qiime feature-table summarize \
  --i-table ryan-table.qza \
  --o-visualization ryan-table.qzv \
  --m-sample-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv
   
qiime feature-table tabulate-seqs \
  --i-data ryan-rep-seqs.qza \
  --o-visualization ryan-rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file ryan-stats.qza \
  --o-visualization ryan-stats.qzv
# Taxonomy analysis 
#Extract the reads
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

#Train classifier with your new ref-seq file
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ryan-ref-seqs.qza \
  --i-reference-taxonomy /datasets/classifiers/silva_ref_files/silva-138-99-tax.qza \
  --o-classifier classifier.qza

#Use the trained classifier to assign taxonomy to your reads (rep-seqs.qza)
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads ryan-rep-seqs.qza \
  --o-classification ryan-taxonomy.qza

#Visualization 
qiime metadata tabulate \
  --m-input-file ryan-taxonomy.qza \
  --o-visualization ryan-taxonomy.qzv
  
#Taxonomy barplots
qiime taxa barplot \
  --i-table ryan-table.qza \
  --i-taxonomy ryan-taxonomy.qza \
  --m-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv \
  --o-visualization ryan-taxa-bar-plots.qzv

#Use trained classifier to create taxonomy file
qiime feature-classifier classify-sklearn \
--i-classifier forwardryan-classifier.qza \
--i-reads forwardryan-rep-seqs.qza \
--o-classification forwardryan-taxonomy.qza

#convert to a visualizable file
qiime metadata tabulate \
  --m-input-file forwardryan-taxonomy.qza \
  --o-visualization forwardryan-taxonomy.qzv

#create taxa barplot
qiime taxa barplot \
  --i-table forwardryan-table.qza \
  --i-taxonomy forwardryan-taxonomy.qza \
  --m-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv \
  --o-visualization forwardryan-taxa-bar-plots.qzv

# Filtering 
#Filtering out mitochondria and Chloroplasts

qiime taxa filter-table \
  --i-table forwardryan-table.qza \
  --i-taxonomy forwardryan-taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table forwardtable-no-mitochondria-no-chloroplast.qza

#convert to visualizable file

qiime feature-table summarize \
  --i-table forwardtable-no-mitochondria-no-chloroplast.qza \
  --o-visualization forwardtable-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv

# Rarefaction in QIIME2 
#Generate a tree for phylogenetic diversity analyses

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences forwardryan-rep-seqs.qza \
  --o-alignment aligned-forwardryan-rep-seqs.qza \
  --o-masked-alignment masked-aligned-forwardryan-rep-seqs.qza \
  --o-tree forwardryan-unrooted-tree.qza \
  --o-rooted-tree forwardryan-rooted-tree.qza 

#Alpha rarefaction
qiime diversity alpha-rarefaction \
  --i-table forwardryan-table.qza \
  --i-phylogeny forwardryan-rooted-tree.qza \
  --p-max-depth 120000 \
  --m-metadata-file /datasets/project_2/human_ibd/ryan_metadata.tsv \
  --o-visualization forward120000ryan-alpha-rarefaction.qzv
  
# Converting files to human-readable format and exporting 
#Create an export folder by mkdir
#convert and export table.qza into human readable file 
  qiime tools export \
  --input-path /data/project_2/Version_2/forwardryan-table.qza \
  --output-path table_export
#Convert and export taxonomy.qza to human readable file
  qiime tools export \
  --input-path /data/project_2/Version_2/forwardryan-taxonomy.qza \
  --output-path taxonomy_export
#Convert and export rooted tree
  qiime tools export \
  --input-path /data/project_2/Version_2/forwardryan-rooted-tree.qza \
  --output-path rooted_tree__export
#Go to the export table directory and convert biom file to txt file.
  biom convert -i feature-table.biom --to-tsv -o feature-table.txt
#Download the entire folder from server to local computer
  scp -r root@10.19.139.189:/data/project_2/Version_2/ryan_export .
#Export metadata file
  scp root@10.19.139.189:/datasets/project_2/human_ibd/ryan_metadata.tsv .
