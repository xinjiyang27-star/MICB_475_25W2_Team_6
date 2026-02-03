# February 2 [16:00-17:00]
## Completed
1. Research topics were narrowed down to infant and chemotherapy
3. Dataset searching from published literature and Evelyn found 3 pontential datasets
   - 2 datasets for infant (one for infant antibiotics and one for pediatric UC)
   - 1 dataset for chemotherapy
5. Brainstorm of research questions based on the metadata and research papers
   
## Meeting Agenda 
1. Finalize the research topic and dataset
   - concern: infant and chemotherapy metadata don't contain too much information, how can we understand their data more?
      - Based on the papers, there should be more metadata such as ages, antibiotic combinations, etc. since the authors specified more "metadata like" information in the material and method section🤔
3. Discussion of research question ideas
   - pediatric UC: if combine this dataset with an adult IBD dataset: compare female and male pediatric UC and/or female/male pediatric UC with female/male adult UC.
   - Research questions about
        1) antibiotic exposure of infant: window specificity 
        2) UC: sex difference in UC
        3) combining 1) and 2): signature comparison (NOT causation) ➡️ taxa depletion by antibiotic
        4) combining 2) and Ryan's dataset (IBD in human) on Canvas ➡️ age-specific signature
        5) chemotherapy
        6) Differences in microbiota in samples collected from stool vs liquid biopsy
5. Finalize the potential research questions
6. Q&A about concepts, technicals, etc. (if there is)

## Meeting Notes 
1. infant antibiotics dataset:
   - maybe able to use the timeline/chart to annotate the metadata
   - issues: doesn't specify infant ID and which antibiotics the infant takes
2. pediatric UC dataset: best one with more information
   - small dataset: 60 samples
   - issue: age from 7-21 years
   - can combine with an adult IBD dataset: better to use fecal sample datasets
        - Halfvarson paper: doesn't contain age, but can assume adult
             - might not need age if we want to compare infant and adult
        - Human_IBD dataset: 20 samples and just the disease status
   - can bin the age
        - do a longitudinal study 
3. Combine dataset
   - combine at the beginning: import both datasets are together
        - can only do that if they sequence the same region and both either single or paired sequence
        - both datasets(pediatric UC and Halfvarson paper) sequenced V4
        - will do the combination in R
   - process the dataset and then generate two tables -> merge them afterwards
4. machine learning
   - train an pediatric UC dataset -> apply it to a new dataset -> predict the age of the dataset
   - 30% trainning and 70% testing
   - 1) build 2 models -> compare 2 models -> how comparable are these 2 models
           - advantage: there is significant between UC and healthy
           - select which microbes should be included -> choose those that are very associated with disease conditions
     3) build 1 model -> predict for another dataset
5. proposal
   - can change the proposal if we end up having one very strong model 
  
## Ideas
1) longitudianl study + funtional analysis -> which pathway is upregulated or downregulated
   - process pediatric UC only using QIIME2 -> longtitudinal alpha diversity using R -> bin ages (7-10,11-15,>16) -> indicator taxa, core bicrobiome, deseq in R 
        - if longtitudinal alpha diveristy doesn't show significance, not meaningful to do indicator taxa, core bicrobiome, deseq
2) compare pediatric and adult UC by combining 2 datasets
   - combine data in QIIME2 -> pediatric healthy and UC vs adualt healthy and UC in R -> diversity metrics, taxoomic bar in R 
3) machine learning
   - 2 poeple: pediatric UC in QIIME2 -> desep or indicator taxa in R -> random forest model (build the model) -> % accuracy
   - 3 people: Halfvarson in QIIME2 -> desep or indicator taxa in R -> random forest model (build the model) -> % accuracy
   - combine pediatric UC and Halfvarson model -> cross compare them
   - objective: Buid machine learning models to predict UC patients and healthy for pediatric and adult.
   
## Q&A
❓Q1 How do we combine datasets? What do we need to check to combine them?

❓Q2 Will the datasets and metadata files which we found be uploaded to the server?

❓Q3 If the research question that we came up with are unfortunately addressed or half-addressed by a paper from somewhere (Not the original dataset-containing paper nor UJIMI), how can we do? Do we need to switch to another question directly or add some conditions in the research question (since we would not use the same dataset, so can we add some limiting conditions based on the dataset and make the research question "novel")

## To-Do List 
1. go over the module for machine learning
2. read over the machine learning dataset from last term (opioid one)
