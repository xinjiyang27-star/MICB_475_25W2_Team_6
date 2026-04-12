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
load("ryan_filt_new.RData")
ps <- ryan_filt_new

# Convert taxonomy table to matrix and assign genus to your ASV
tax_mat <- as(tax_table(ps), "matrix")

tax_mat["f9283780d796ffcc12a283e56b806087", "Genus"] <- "g__Ruminococcus"
tax_mat["e96331d4191c5a32f8d1347fc3d295e2", "Genus"] <- "g__Harryflintia"
tax_mat["97351bb5e19f0e3d19bab552b7055653", "Genus"] <- "g__Faecalimonas"
tax_mat["f949c3d6392f7d23a7cbdb9e8f7180cd", "Genus"] <- "g__Coprococcus"
tax_mat["8c2297f07d90de0a6b60701c8ea376cb", "Genus"] <- "g__Anaerostipes"
tax_mat["8c5c1e33bdd2f244342f2acd6e5db57e", "Genus"] <- "g__Lachnoclostridium"
tax_mat[tax_mat[, "Genus"] == "g__UCG-005", "Genus"] <- "g__Faecousia"
tax_mat[tax_mat[, "Genus"] == "g__GCA-900066575", "Genus"] <- "g__Laedolimicola"
tax_mat[tax_mat[, "Genus"] == "g__GCA-900066575 ", "Genus"] <- "g__Laedolimicola"
# Reassign the updated taxonomy table to the phyloseq object
tax_table(ps) <- tax_table(tax_mat)

tax_table(ps) %>%
  as.data.frame() %>%
  filter(grepl("Laedolimicola", Genus))

# Check that the taxonomy table is valid
tax_table(ps)["8c5c1e33bdd2f244342f2acd6e5db57e", ]

# Aggregate reads to Genus level
ps <- tax_glom(ps,"Genus")

# Load the DESeq taxa results
DESeq_results <- read.csv("DESeq_CD_ryanfilt_new.csv")

# Normalize microbiome data and select only the significant taxa
# CLR transformation
ps_clr = ps %>% microbiome::transform('clr') 
# Filter taxa
ps_filt = ps_filt = prune_taxa(DESeq_results$ASV, ps_clr)

# Melt the dataset
df = psmelt(ps_filt) %>% 
  # dashes will cause errors
  mutate(Condition = str_replace_all(Condition,"Crohn's Disease","Crohns_Disease"))

# Add Z transformation to each Genus individually
# (mean of zero, standard deviation of 1)
df_transformed = df %>% 
  group_by(Genus) %>% 
  mutate(Abundance = scale(Abundance)) %>% 
  ungroup()

# Filter For Conditions
target_conditions <- c("Crohns_Disease", "Healthy")
filtered_df <- df_transformed %>% 
  filter(Condition %in% target_conditions)

# Final table should ONLY contain the outcome and explanatory variables, each as their own column. For now we'll also include the sample id.
df_pivot = filtered_df %>% 
  select(Sample,Condition,Genus,Abundance) %>% 
  # Turn each Genus into its own column
  pivot_wider(names_from = Genus, values_from = Abundance)

# Remove rows with NA values in the metadata
df_noNA = df_pivot %>% na.omit()

# Remove the sample ID column - otherwise the code will try to use it as an explanatory variable (just like the microbial genera).
# We do this after pivoting (try doing it before pivoting, see what happens!)
df_final = df_noNA %>% select(-Sample)

#Set Predictors
predictors = df_final %>% select(-Condition)

# We will transform subject into a factor.
# The first level of the factor is automatically used as the reference group.
# We will put subject_1 first.
# The results will describe subject 2 relative to the reference (subject_1).
outcome = df_final %>% pull(Condition) %>% 
  factor(levels = c("Crohns_Disease","Healthy"))

# Randomly subsets the rows into k equal bins.
k = 10
set.seed(066)
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
tune_grid = expand.grid(mtry = c(5,7,10), 
                        splitrule = c("gini","extratrees"),
                        min.node.size = c(21,42,83))

#Run RF
source('randomforest_functions.R')
pd_model = run_rf(X = predictors, y = outcome, 
                  fold_list = folds,
                  hyper = tune_grid, 
                  rngseed = 066)
names(pd_model)

#Interpretation

roc_test = roc(pd_model$test_labels$true_labels,
               pd_model$test_labels$predicted_probabilities)
roc_train = roc(pd_model$train_labels$true_labels,
                pd_model$train_labels$predicted_probabilities)
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
                           auc(roc_train), pd_model$auc_train_ci[1], pd_model$auc_train_ci[2],
                           auc(roc_test), pd_model$auc_test_ci[1], pd_model$auc_test_ci[2]), 
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


ggsave("ROC_Plot_DESeq_taxa_only.png", roc_plot, width = 8, height = 6, dpi = 300)


importance_plot <- pd_model$importance %>% 
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
ggsave("Importance_Plot_DESeq_taxa_only.png", importance_plot, width = 8, height = 6, dpi = 300)

#Check Correlation

temp = ps_clr %>% 
  subset_taxa(Genus=='g__Laedolimicola') %>% 
  psmelt()

box_plot_RF <- temp %>% 
  ggplot(aes(Condition,Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height=2, width = 0.2) +
  theme_classic(base_size=18)

box_plot_RF
ggsave("Laedolimicola_abundance.png", box_plot_RF, width = 8, height = 6, dpi = 300)
pd_model$importance
