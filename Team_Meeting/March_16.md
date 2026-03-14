# March 16 [16:00 - ] 
## Meeting Agenda 
1. Discuss Random Forest Model and ISA results

## Q and A
1. For Random forest:
   - Why all the modes are overfitted?
   - Why the models are not consistent between runs with the same laptop and different laptops?
   - More information on mtry and min.node.size? How to set the most appropriate parameters?
   - Inconsistent between models(metadata only and metadata+taxa):
        - For unfiltered model, Age when sampled is MORE importance than biopsy location, which is consistent between metadata only and metadata+taxa.
        - For filtered model, Age when sampled is LESS important than biopsy location for metadata+taxa, but otherwise for metadata only.
        - We think it's driven by this. Which code is better:
           1) df_pivot = filtered_df %>% \
                    select(Sample,Condition,Genus,Abundance) %>% \
                    # Turn each Genus into its own column \
                    pivot_wider(names_from = Genus, values_from = Abundance) 
           2) predictors <- df_noNA %>% \
                 select(-Condition) %>% \
                 #Need to change categorical variables to factors \
                 mutate( \
                    Biopsy_location = as.factor(Biopsy_location), \
                    Smoking.status = as.factor(Smoking.status), \
                    Gender = as.factor(Gender), \
                    Age_when_sampled = as.numeric(Age_when_sampled))
   - For model with metadata only, should we load metadata or phloseq object?
   - Other other suggestions?
  

## Experimental Results Summary
- Background and Terminology:
  - non-filtered metadata: Keeping all samples, including UC, CD, and healthy; No filteration conducted
  - filtered metadata: Filtering out all UC and non-inflamed CD, only keep inflamed CD and all healthy
- Attempt 1-3 are trying/ playing around of RF with our data

### Attempt 4: Build the RF by using ISA taxa only
   - metadata: non-filtered 
   - predictors: the 7 taxa identified by ISA
     ![ROC](/images/Attempt4_ROC.png)
     ![Importance bar plot](/images/Attempt4_Importance.png)

 ### Attempt 5: Build the RF by using non-filtered metadata only
   - metadata: non-filtered 
   - predictors: Age, biopsy location, gender, smoking condition
     ![ROC](/images/Attempt5_ROC.png)
     ![Importance bar plot](/images/Attempt5_Importance.png)

### Attempt 6: Build the RF by combining nonfiltered metadata and non-filtered ISA taxa 
   - metadata: non-filtered 
   - predictors: 7 ISA taxa, Age, biopsy location, gender, smoking condition
    ![ROC](/images/Attempt6_ROC.jpeg)
    ![ROC](/images/Attempt6_Importance.jpeg)

### Attempt 7: Build the RF by using filtered ISA taxa only 
   - metadata: filtered 
   - predictors: 15 taxa identified by ISA
    ![ROC](/images/Attempt7_ROC.png)
    ![ROC](/images/Attempt7_Importance.png)

### Attempt 8: Build the RF by using filtered metadata only 
   - metadata: filtered 
   - predictors: Age, biopsy location, gender, smoking condition
    ![ROC](/images/Attempt8_ROC.png)
    ![ROC](/images/Attempt8_Importance.png)

### Attempt 9: Build the RF by combining filtered metadata and filtered ISA taxa 
   - metadata: filtered 
   - predictors: 15 taxa identified by ISA, age, biopsy location, gender, smoking condition
    ![ROC](/images/Attempt9_ROC.png)
    ![ROC](/images/Attempt9_Importance.png)
