# February 9 [16:00-]
## Completed
1. Finalized the dataset selection which will be used for research project
   - Halfvarson paper
2. Finalized research ideas: build a machine learning model to distinguish UC, CD, and healthy for diagnostic purposes
4. Draft/Skeleton of the research ideas and proposal 

## Meeting Agenda 
1. ❗(Please give Jiyang ~5 minutes to talk about the issues and his current ideas, thanks!) Halfvarson paper discussion: I found something really important about this paper, which may affect our research topic, let's discuss it in the meeting -- Jiyang 
2. Go through the proposal skeleton (Introduction&Background; Hypothesis; Experimental aims) and ask Evelyn for suggestions
3. Q and A:
   - For introduction, how much details should we include about UC and CD microbiota? Do we need to include how species are different or more general trend?
   - In terms of literature review: Do we have a least numbers of papers that need to be included in the literature review (I looked through the example proposal, which seems like they didn't include many papers for that section)?
   - In terms of the current knoweldge gap: my current idea 1) most ML is trained on metagenomic sequencing, not too many papers talks about training based on 16S sequencing. 2) Identify any gap from the original Halfvarson paper
   - ❗Halfvarson paper also did machine learning (random forest) to predict the subtype of IBD... We need to discuss a novel research idea/question
   - Would normalizing sequencing depths and number of samples be enough for aim 1
   - How specific do we have to be in the hypothesis? Should we mention known taxa that are elevated/reduced between the diseases?
This is the skeleton of Introduction&Background section(This is only a draft, which some contents may be repetitive, and will be revised later)
![Alt text](/images/Intro-1.png)
![Alt text](/images/Intro-2.png)
## Meeting Notes 
1. Halfvarson paper
  - They used machine learning by RF to do the prediction classifier
2. Ryan
   - No Machine learning
   - we can go with Ryan paper to use the dataset to do comparison with Halfvarson paper's model --> which one is better for prediction 
   - Need to control for location of biopsy
   - Do many small models for each location --> to see which one is the best one
        - control biopsy location
        - Put gender and medication in the model --> to see medication or gender is a good predictor
        - Can include all in the model
    - finding the key ASV to train the model

3. Contact Hans --> Ryan's dataset
   - Hans will know --> to cut the barcode
   - We have to remove the barcodes before processing. 


## To-do List
