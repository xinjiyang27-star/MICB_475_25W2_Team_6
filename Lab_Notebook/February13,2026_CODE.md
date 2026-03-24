# Denoising and Clustering [13 February]
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

# Visualize Files
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
