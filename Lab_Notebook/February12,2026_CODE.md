Login to Server:

ssh root@10.19.139.189

Password: Biome391

#Open new screen
screen -S Ryan

#Import Data
qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path ./datasets/project_2/human_ibd/ryan_manifest.tsv \
  --output-path ./data/project_2/ryan_demux_seqs.qza

#Trim
qiime cutadapt trim-paired \
 --i-demultiplexed-sequences ryan_demux_seqs.qza \
 --p-front-f TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG \
 --p-front-r GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC \
 --p-match-read-wildcards \
 --p-match-adapter-wildcards \
 --p-discard-untrimmed \
 --o-trimmed-sequences ryan_demux_seqs-trimmed.qza \
 --verbose > cutadapt-log-2.txt

#Visualize
# Create visualization of demultiplexed samples
qiime demux summarize \
  --i-data ryan_demux_seqs-trimmed.qza \
  --o-visualization ryan_demux_seqs-trimmed.qzv

#Local 
scp root@10.19.139.189:/data/project_2/ryan_demux_seqs-trimmed.qzv .


# Determine ASVs with DADA2
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ryan_demux_seqs-trimmed.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 120 \
  --p-trunc-len-r 120 \
  --o-representative-sequences ryan-rep-seqs.qza \
  --o-table ryan-table.qza \
  --o-denoising-stats ryan-stats.qza

#Ensure overlap is: forward length after truncation + reverse length after truncation - amplicon length ≥ 20–30 bp overlap

