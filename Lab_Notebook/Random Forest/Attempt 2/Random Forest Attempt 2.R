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
ps = readRDS('ryan_phyloseq_CDvsHealthy.rds') %>% 
  tax_glom('Genus')

# Calculate average abundance~~~~~~~~~~~~
avg_abundance = taxa_sums(ps)/sum(taxa_sums(ps)) 
# Sort high to low
avg_abundance = sort(avg_abundance, decreasing = T)
# Take the top 10
top_10 = avg_abundance[1:10]
# Extract taxa names 
top_10 = names(top_10)

# Normalize microbiome data and select only the significant taxa
# CLR transformation
ps_clr = ps %>% microbiome::transform('clr') 
# Filter taxa
ps_filt = prune_taxa(top_10,ps_clr)
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
  select(Sample,Gender,Condition,Genus,Abundance,Biopsy_location,Smoking.status) %>% 
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
set.seed(061)
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
                  rngseed = 061)
names(pd_model)

#Interpretation

roc_test = roc(pd_model$test_labels$true_labels,
               pd_model$test_labels$predicted_probabilities)
roc_train = roc(pd_model$train_labels$true_labels,
                pd_model$train_labels$predicted_probabilities)
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
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  annotate("text", x = 0.7, y = 0.2, 
           label = sprintf("Train (red): %.2f (%.2f-%.2f)\nTest (black): %.2f (%.2f-%.2f)",
                           auc(roc_train), pd_model$auc_train_ci[1], pd_model$auc_train_ci[2],
                           auc(roc_test), pd_model$auc_test_ci[1], pd_model$auc_test_ci[2]), 
           size = 6) +
  theme_minimal(base_size=18)

pd_model$importance

pd_model$importance %>% 
  # Data are automatically arranged by decreasing importance - turn it into a factor.
  # Otherwise the features will show up alphabetically in the plot.
  mutate(Feature = factor(.$Feature,levels = .$Feature)) %>% 
  ggplot(aes(Feature,MeanDecreaseGini,fill=MeanDecreaseGini)) +
  geom_col() +
  theme_classic(base_size=18) +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1)) +
  ylab('Importance (Gini)') + xlab(NULL)

#Check Coorelation

temp = ps_clr %>% 
  subset_taxa(Genus=='g__Anaerostipes') %>% 
  psmelt()

temp %>% ggplot(aes(Condition,Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height=0, width = 0.2) +
  theme_classic(base_size=18)
