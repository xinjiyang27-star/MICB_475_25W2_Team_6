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


roc_plot <- ggplot() +
  # Training ROC
  geom_line(aes(x = 1 - roc_train$specificities, 
                y = roc_train$sensitivities,
                color = "Training"),
            size = 1.2) +
  
  # Test ROC
  geom_line(aes(x = 1 - roc_test$specificities,
                y = roc_test$sensitivities,
                color = "Test"),
            size = 1.2) +
  
  
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", 
              color = "grey50", 
              size = 0.8) +
  
  
  labs(
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Dataset"
  ) +
  
  
  annotate("text", x = 0.65, y = 0.15, 
           label = sprintf("Train: %.2f (%.2f–%.2f)\nTest: %.2f (%.2f–%.2f)",
                           auc(roc_train), filt_clinical_model$auc_train_ci[1], filt_clinical_model$auc_train_ci[2],
                           auc(roc_test), filt_clinical_model$auc_test_ci[1], filt_clinical_model$auc_test_ci[2]), 
           size = 4.5, hjust = 0) +
  
  
  scale_color_manual(values = c(
    "Training" = "#22B4EE",  
    "Test" = "#ee55cc"       
  )) +
  
  
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

roc_plot


ggsave("ROC_Plot_8.png", roc_plot, width = 8, height = 6, dpi = 300)


importance_plot <- filt_clinical_model$importance %>% 
  mutate(Feature = factor(Feature, levels = Feature)) %>% 
  
  ggplot(aes(x = Feature, y = MeanDecreaseGini, fill = MeanDecreaseGini)) +
  geom_col(width = 0.75, color = "black", size = 0.2) +
  
  scale_fill_gradient(
    low = "#E0F5FD", 
    high = "#22B4EE",
    name = "Gini Importance"
  ) +
  
  labs(
    y = "Importance (Gini)",
    x = NULL
  ) +
  
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    
    axis.line = element_line(size = 0.6),
    axis.ticks = element_line(size = 0.6)
  )
importance_plot
ggsave("Importance_Plot_8.png", importance_plot, width = 8, height = 6, dpi = 300)
