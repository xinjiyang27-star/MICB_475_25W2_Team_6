# Sequencing data process in QIIME2 12 February 
## Purpose: 
To import the data from Ryan's paper and demultiplex the data for downstream analysis. 

## Materials:
MICB475 Server
  - IP address: root@10.19.139.189
  - password: Biome391

## Method:
1. Login to server by provided IP address and password.
2. Open a separate screen for background job running named "Ryan"

   screen -S Ryan 
4. import the data from project 2 datasets

   qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path ./datasets/project_2/human_ibd/ryan_manifest.tsv \
  --output-path ./data/project_2/ryan_demux_seqs.qza

5. Trim the adaptor within the reads (codes are provided by Hans)

   qiime cutadapt trim-paired \
 --i-demultiplexed-sequences ryan_demux_seqs.qza \
 --p-front-f TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG \
 --p-front-r GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC \
 --p-match-read-wildcards \
 --p-match-adapter-wildcards \
 --p-discard-untrimmed \
 --o-trimmed-sequences ryan_demux_seqs-trimmed.qza \
 --verbose > cutadapt-log-2.txt

6. Create visualization of demultiplexed samples

   qiime demux summarize \
  --i-data ryan_demux_seqs-trimmed.qza \
  --o-visualization ryan_demux_seqs-trimmed.qzv

8. Download the qzv file to Local

   scp root@10.19.139.189:/data/project_2/ryan_demux_seqs-trimmed.qzv .

## Results:
[ryan_demux_seqs-trimmed.qzv](/Lan_Notebook/Experimental_Results/)

