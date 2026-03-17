#Load Libraries
library(randomForest)
library(caret)
library(ranger)
library(pROC)
library(boot)
library(ggplot2)
library(phyloseq)
library(tidyverse)

# Load the dataset and aggregate reads to the Genus level
filtered_metadata <- read.delim("meta_filt.tsv", sep = "\t")

# Selecting the metadata columns + outcomes
df <- filtered_metadata %>%
  select(Condition, Age_when_sampled, Biopsy_location, Smoking.status, Gender) %>%
  # Apostrophes causes errors in R
  mutate(Condition = str_replace_all(Condition, "Crohn's Disease", "Crohns_Disease")) %>%
  # Filter to exclude UC
  filter(Condition %in% c("Crohns_Disease", "Healthy")) %>%
  # Have to change Age_when_sampled to factor before the others
  # as it introduces new NAs if done after the actual NA removal step
  # this is because it has a string that cannot be converted to NA, the string
  # says "uncertain"
  mutate(Age_when_sampled = as.numeric(Age_when_sampled))

# Remove rows with NA values
df_noNA <- df %>% na.omit()

# Set Predictors
predictors <- df_noNA %>%
  select(-Condition) %>%
  # Need to change categorical variables to factors
  mutate(
    Biopsy_location = as.factor(Biopsy_location),
    Smoking.status = as.factor(Smoking.status),
    Gender = as.factor(Gender),
    Age_when_sampled = as.numeric(Age_when_sampled))

# Set Outcomes
outcome <- df_noNA %>% pull(Condition) %>%
  factor(levels = c("Healthy", "Crohns_Disease"))

# Randomly subsets the rows into k equal bins.
k = 10
set.seed(060)
folds = createFolds(outcome, k = k, list = TRUE)

# Each of these folds will take a turn being the test dataset.
str(folds)
#Hyperparameter Tuning
# mtry: number of variables that will be used per forest. 
#       High = overfitting, low = uninformative

# splitrule: affects how decision trees are calculated.
#            Use gini or extratrees for boolean outcomes (ex. subject)
#            Use variance for continuous outcomes (ex. age)

# min.node.size: Related to tree complexity. Larger = simpler tree.
# Often best as a proportional fraction of your sample size.

# These are generic values. Depending on your dataset, you may need to adjust the numeric ranges up or down.
tune_grid = expand.grid(mtry = c(1,2,3,4), 
                        # mtry has max 4 since we have 4 predictors will test it
                        splitrule = c("gini","extratrees"),
                        min.node.size = c(21,42,83))



#Run RF
source('randomforest_functions.R')
filt_clinical_model <- run_rf(
  X = predictors,
  y = outcome,
  fold_list = folds,
  hyper = tune_grid,
  
  rngseed = 060)

names(filt_clinical_model)

#Interpretation

roc_test = roc(filt_clinical_model$test_labels$true_labels,
               filt_clinical_model$test_labels$predicted_probabilities)

roc_train = roc(filt_clinical_model$train_labels$true_labels,
                filt_clinical_model$train_labels$predicted_probabilities)

ggplot() +
  # Training data: this is a type of control
  geom_line(aes(x = 1 - roc_train$specificities, 
                y = roc_train$sensitivities), 
            color = "red",size=1) +
  # Test data: tells us the strength of the prediction
  geom_line(aes(x = 1 - roc_test$specificities,
                y = roc_test$sensitivities), 
            color = "black",size=1) +
  geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dashed",size=1) +
  labs(x = "False Positive Rate", y = "True Positive Rate", 
       title = "Filtered Clinical Model") +
  annotate("text", x = 0.7, y = 0.2, 
           label = sprintf("Train (red): %.2f (%.2f-%.2f)\nTest (black): %.2f (%.2f-%.2f)",
                           auc(roc_train), filt_clinical_model$auc_train_ci[1], filt_clinical_model$auc_train_ci[2],
                           auc(roc_test), filt_clinical_model$auc_test_ci[1], filt_clinical_model$auc_test_ci[2]), 
           size = 6) +
  theme_minimal(base_size=18)

filt_clinical_model$importance

filt_clinical_model$importance %>% 
  # Data are automatically arranged by decreasing importance - turn it into a factor.
  # Otherwise the features will show up alphabetically in the plot.
  mutate(Feature = factor(.$Feature,levels = .$Feature)) %>% 
  ggplot(aes(Feature,MeanDecreaseGini,fill=MeanDecreaseGini)) +
  geom_col() +
  theme_classic(base_size=18) +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1),
        plot.margin = margin(10, 0, 10, 10)) +
  ylab('Importance (Gini)') + xlab(NULL) 

