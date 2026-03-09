# March 9 [16:00 - ] 
## Meeting Agenda 
1. Discuss candidate taxa to use for building the Random Forest model
2. Discuss metadata categories to use for building the Random Forest model
   - metadata categories: biopsy location, histology status, medications, condition,   gender 
4. Discuss RF parameters and resulting AUC values 

## Q and A
1. Is the random forest model working acceptably, or is it being overfit? If so, why and how can this be resolved?
2. Are there certain parameters we should control for, remove, or actually include in the model?

## Experimental Results Summary
1. Core microbiome
<img width="633" height="881" alt="image" src="https://github.com/user-attachments/assets/a892dad2-70e3-4d03-b6e6-27ae6387ad33" />

2. Indicator Taxa
<img width="964" height="1222" alt="image" src="https://github.com/user-attachments/assets/ba5b9ade-9bf2-4017-8d6d-3db3a8e19928" />

3. DESeq
<img width="418" height="556" alt="Screenshot 2026-03-08 at 10 00 27 PM" src="https://github.com/user-attachments/assets/abfffa92-5542-4c0a-89b2-4ef67aca5239" />

4. Random Forest
- First Attempt:
   - Retained:Gender,Condition,Genus,Abundance,Biopsy_location,Histological.status,Smoking.status. Other columns are either some kind of identifier or relatively not useful.
   - Altered the condition to just have crohns disease and healthy. Set the k crossfold as 10 as recommended for large datasets like ours. Started off with default mtry:5,7,10. Changed the      minnodules to be proportional to our sample size with half as max as recommended. 
   - Result: Suspiciously good with a training and test AUC of above 0.9. Overfitted?
   - Histological status most important variable, followed by genus Anaerostipes
   - Selected taxa based on highest abundance, not analyses.
   - Biopsy location actually came behind certain genera in importance.
     ![Importance Attempt 1](https://github.com/user-attachments/assets/1105deba-5054-4603-b418-61fd2066e426)
     ![Correlation Attempt 1](https://github.com/user-attachments/assets/fa7bd947-11b3-486e-b8a8-39b0152941e3)
     ![AUC Attempt 1](https://github.com/user-attachments/assets/78a9c131-3ac4-434b-b498-98e33d8d0eee)
- Second Attempt:
   - Removed histological status, no other changes
     ![Correlation Attempt 2](https://github.com/user-attachments/assets/4b3b3e15-6403-493b-b838-3504630f44e3)
     ![AUC Attempt 2](https://github.com/user-attachments/assets/7f2e9664-2d4e-49ea-8695-2927c723f58c)
     ![Importance Attempt 2](https://github.com/user-attachments/assets/68a5f649-7663-4f7d-b7a3-201f6595e08f)
- Third Attempt:
  - Altered hyperparameter tuning.
    ![Correlation Attempt 3](https://github.com/user-attachments/assets/d22937e9-92e1-4b0c-836a-a8f35721340a)
    ![AUC Attempt 3](https://github.com/user-attachments/assets/aa39a62a-857f-4cc6-a0c5-fc3f30381339)
    ![Importance Attempt 3](https://github.com/user-attachments/assets/504168ed-29da-48f7-98c2-bbad44e2171d)
## Meeting Minutes












## TO-DO lists
