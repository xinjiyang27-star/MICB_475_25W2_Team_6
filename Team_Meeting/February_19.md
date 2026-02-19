# February 19 [12:00-12:45]
## Questions to discuss 
1. quality score cut off:
   1) stay with what we have so far so far --> if we stay with it, what about the rarefaction?
<img width="1007" height="538" alt="Screen Shot 2026-02-18 at 10 04 23 PM" src="https://github.com/user-attachments/assets/3e65d716-89eb-4870-b970-219ff6d631eb" />
<img width="508" height="313" alt="Screen Shot 2026-02-18 at 10 04 58 PM" src="https://github.com/user-attachments/assets/e4d11504-1024-4216-b0c7-7ad294299b67" />
<img width="510" height="316" alt="Screen Shot 2026-02-18 at 10 06 08 PM" src="https://github.com/user-attachments/assets/eaba60f4-d6de-44b5-8636-c8e4a16a1805" />
<img width="1142" height="673" alt="Screen Shot 2026-02-18 at 10 09 50 PM" src="https://github.com/user-attachments/assets/e8d55da6-ba9f-4920-a67a-2726bdc2f1f7" />

   2) lower to 20
<img width="997" height="492" alt="Screen Shot 2026-02-18 at 10 06 42 PM" src="https://github.com/user-attachments/assets/425a3a9f-c874-43aa-8c5a-844f3bd83c24" />
   <img width="512" height="337" alt="Screen Shot 2026-02-18 at 10 06 59 PM" src="https://github.com/user-attachments/assets/c33f5c62-8dc0-4494-b83f-1060b5453fca" />
<img width="506" height="335" alt="Screen Shot 2026-02-18 at 10 07 22 PM" src="https://github.com/user-attachments/assets/2e620e3f-7e3f-4a99-8fdc-6f028219f54d" />
<img width="1141" height="723" alt="Screen Shot 2026-02-18 at 10 10 21 PM" src="https://github.com/user-attachments/assets/d016af41-1d4e-419d-bb33-7e01d3bcd9c5" />

   3) lower to 25, paper says
   4) do forward only, discard the reverse reads   
1. purpose of running alpha diversity and beta diversity and what statistical test should we use for the alpha.
   - KW test? permanova? Tukey post-hoc test
2. Explanation of running alpha and b diversity --> why we are doing this, doing that
3. Do we choose one of the A or B diversity metrics or running all of them; if we choose one of them, how can we choose? 
4. questions embedded in the chart
5. DeSeq: After we get the values from DeSeq, is the statitical test just for determining the p value cut off
   - Deseq will tell us the p value
   - should we run the stat test after we get the p values?
6. ISA: Results/Output is the chart?
7. We run all of the CoreMicrobiome, ISA, DeSeq --> determine the taxa used to train the model? but If each of them gives us a different result, what can we do?
8. What statistical test do CoreMicrobiome, ISA, DeSeq use? Since we don't choose the statistical test, do we need to find the default statistical test for them?
9. First line of chart --> Aim1 and Aim 2 sample using --> use all? use some? ...
10. research sub-aim: locations of the bipsy for CD, UC, Healthy classification? Main aim is using biopsy for ML IBD subtype classification. 
## Meeting Notes
**Answers to Questions:**
1. Proceed with just forward reads. Keep proposal as it is, as though paired reads, but in reality, run it as forward reads only.
  - adjust manifest, 3 columns rn, so what we have to do is delete reverse pathway, change forwards pathway to absolute. So it will only read forward reads; they are better quality anyway.
  - 180 or 220bp is fine for the cut-off
  - Benefit of paired, it's more beneficial IF THEY PAIR
  - Base rarefraction value based on alpha rarefraction, base it on a condition. Play with sampling depth to maximize where we don't see a decrease in frequency.
  - Alpha and Beta diversity are not super important. Important to still run the analysis to see how different the conditions are, best bet is to rerun as single-end reads, and see if we have the sampling depth issue. If there is an issue, rewrite data overview part.
  - Last thing we need is the rarefaction value: "if we try to rarefy, we lose a lot of samples, so we will re-run as single reads and do the rarefaction action, hence we cannot run the diversity metrics."
    
1'. Run all alpha diversity and Beta. If alpha ALWAYS KW test, if Beta ALWAYS Permanova. We want to see which metrics has significant differences, informs which downstream analysis is more important in picking taxa for the model.
  - We need to run the Post-hoc test too for the KW test
5. Stat test is embedded into the DESEQ, ISA, etc. When listing, we just say "Internal test determined by the actual analysis". Keep everything at 0.05 (p-value). Contingent on taxa - if we have a lot of indicator taxa then we need to be more stringent with the p-value. All the analyses highlight which taxa are different in the conditions - should be in experimental aims.
6. ISA output will be a table
9. Control for location, keep the samples for a specific location. For the proposal, say "depends on model, we may control for location or treatment, but we will apply it to the model first to determine if it is necessary" - for the first 3 analyses - will depend on the random forest model.
- For random forest, say "all samples, contingent on what we found as confounding variables."
10. Don't talk about comparing the locations of sampling in the proposal.
  
_Other notes:_
- Demux file looks good, no issues there. After trimming, however, reads went down A LOT (700ish from 19Million)
- Alignment of pair reads may just not be good because the barcodes may disrupt the alignment. 
- Don't have sub-aims.
- Start running single-end reads even if we don't put it in the proposal (oo munts did it alr)
- After we trim, and get the table, email Evelyn to see we get good reads. Not sure if 220 is stringent enough
- RF is a statistical test,
- Do mention what the cut off for the p-adj for volcano plot 


