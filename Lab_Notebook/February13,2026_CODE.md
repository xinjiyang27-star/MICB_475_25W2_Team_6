# Determine ASVs with DADA2
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
