qiime cutadapt trim-paired \
  --i-demultiplexed-sequences ryan_demux_seqs.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --p-adapter-f CTGTCTCTTATACACATCT \
  --p-adapter-r CTGTCTCTTATACACATCT \
  --p-overlap 12 \
  --p-match-read-wildcards \
  --p-match-adapter-wildcards \
  --o-trimmed-sequences ryan_demux_seqs-trimmed.qza \
  --verbose > cutadapt-log-2.txt

  qiime demux summarize \
  --i-data ryan_demux_seqs-trimmed.qza \
  --o-visualization ryan_demux_seqs-trimmed.qzv

  qiime vsearch merge-pairs \
  --i-demultiplexed-seqs ryan_demux_seqs-trimmed.qza \
  --p-minovlen 12 \
  --p-maxdiffs 10 \
  --p-threads 4 \
  --o-merged-sequences ryan_joined.qza \
  --o-unmerged-sequences ryan_unmerged.qza

  qiime quality-filter q-score \
  --i-demux ryan_joined.qza \
  --o-filtered-sequences ryan_joined_qfiltered.qza \
  --o-filter-stats ryan_joined_qfilter_stats.qza \
  --verbose

  qiime demux summarize --i-data ryan_joined.qza --o-visualization ryan_joined.qzv

  qiime deblur denoise-16S \
  --i-demultiplexed-seqs ryan_joined_qfiltered.qza \
  --p-trim-length 400 \
  --p-sample-stats \
  --p-jobs-to-start 4 \
  --o-representative-sequences ryan_rep_seqs.qza \
  --o-table ryan_table.qza \
  --o-stats ryan_deblur_stats.qza \
  --verbose

qiime deblur visualize-stats \
  --i-deblur-stats ryan_deblur_stats.qza \
  --o-visualization ryan_deblur_stats.qzv

qiime feature-table summarize \
  --i-table ryan_table.qza \
  --o-visualization ryan_table.qzv

qiime feature-table tabulate-seqs \
  --i-data ryan_rep_seqs.qza \
  --o-visualization ryan_rep_seqs.qzv

scp root@10.34.36.149:/data/ryan_deblur_stats.qzv ~/Desktop
scp root@10.34.36.149:/data/ryan_table.qzv ~/Desktop
scp root@10.34.36.149:/data/ryan_rep_seqs.qzv ~/Desktop


qiime feature-classifier classify-sklearn \
  --i-classifier /datasets/classifiers/silva-138-99-nb-classifier.qza \
  --i-reads ryan_rep_seqs.qza \
  --p-n-jobs 4 \
  --o-classification ryan_taxonomy.qza \
  --verbose

qiime metadata tabulate \
  --m-input-file ryan_taxonomy.qza \
  --o-visualization ryan_taxonomy.qzv

scp root@10.34.36.149:/data/ryan_taxonomy.qzv ~/Desktop

qiime taxa filter-table \
  --i-table ryan_table.qza \
  --i-taxonomy ryan_taxonomy.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table ryan_table_filt.qza

qiime taxa collapse \
  --i-table ryan_table_filt.qza \
  --i-taxonomy ryan_taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table ryan_genus_table.qza

qiime feature-table summarize \
  --i-table ryan_genus_table.qza \
  --o-visualization ryan_genus_table.qzv

qiime tools export --input-path ryan_genus_table.qza --output-path genus_export
biom convert -i genus_export/feature-table.biom \
  -o genus_export/ryan_genus_feature-table.tsv --to-tsv

qiime tools export --input-path ryan_taxonomy.qza --output-path taxonomy_export

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ryan_rep_seqs.qza \
  --o-alignment aligned-ryan-rep-seqs.qza \
  --o-masked-alignment masked-aligned-ryan-rep-seqs.qza \
  --o-tree ryan-unrooted-tree.qza \
  --o-rooted-tree ryan-rooted-tree.qza

mkdir -p ryan_export

qiime tools export --input-path ryan_table.qza --output-path ryan_export/table_export
biom convert -i ryan_export/table_export/feature-table.biom \
  -o ryan_export/table_export/feature-table.txt --to-tsv

qiime tools export --input-path ryan_taxonomy.qza --output-path ryan_export/taxonomy_export

qiime tools export --input-path ryan-rooted-tree.qza --output-path ryan_export/rooted_tree_export

cp /datasets/project_2/human_ibd/ryan_metadata.tsv ryan_export/ryan_metadata.tsv



scp root@10.34.36.149:/data/genus_export/ryan_genus_feature-table.tsv ~/Desktop
scp root@10.34.36.149:/datasets/project_2/human_ibd/ryan_metadata.tsv ~/Desktop
scp root@10.34.36.149:/data/ryan_genus_table.qzv ~/Desktop


