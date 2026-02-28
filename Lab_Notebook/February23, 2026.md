# February 23 2026
## Purpose
To convert table.qza, phylogeentic tree, metadata, and taxonomy.qza to human readable files and export from QIIME2. 
To Create Phyloseq object and recreate rarefaction, alpha diversity in R. 
## Materials 
1. QIIME2 Group Server
  IP: root@10.19.139.189
  Password: Biome391
2. R Studio

## Method
1. Connect to server and go to the project working directory
2. Create an export folder by mkdir
3.convert and export table.qza into human readable file 
  qiime tools export \
  --input-path /data/project_2/Version_2/forwardryan-table.qza \
  --output-path table_export
4. Convert and export taxonomy.qza to human readable file
  qiime tools export \
  --input-path /data/project_2/Version_2/forwardryan-taxonomy.qza \
  --output-path taxonomy_export
5. Convert and export rooted tree
  qiime tools export \
  --input-path /data/project_2/Version_2/forwardryan-rooted-tree.qza \
  --output-path rooted_tree__export
6. Go to the export table directory and convert biom file to txt file.
  biom convert -i feature-table.biom --to-tsv -o feature-table.txt
7. Download the entire folder from server to local computer
  scp -r root@10.19.139.189:/data/project_2/Version_2/ryan_export .
8. Export metadata file
  scp root@10.19.139.189:/datasets/project_2/human_ibd/ryan_metadata.tsv .
