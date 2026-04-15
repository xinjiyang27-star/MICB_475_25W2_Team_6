# Building Random Forest Models [25 March]

## Purpose:
To build the Random Forest models based on: 
1. micribial taxa only
2. metadata only
3. microbial taxa and metadata combination

## Materials:
1. R Studio 
2. filtered phyloseq object
## Methods:  
### Microbial taxa only 
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
tax_mat["949d11292146e12904547a6db61ce024", "Genus"] <- "g__Laedolimicola"
tax_mat["a716b99232fa3a34ee30ae758b9216c1", "Genus"] <- "g__Beduinella"
tax_mat["0e6b48ffce8c52ae15a77ccd71ed5f31", "Genus"] <- "g__Pseudoruminococcus"
tax_mat["42a8749495523917f1ff19114d6f374b", "Genus"] <- "g__Bifidobacterium"
tax_mat["30753be97d9d33c8c68f245f6c0466e1", "Genus"] <- "g__Brotomerdimonas"
tax_mat["007f2738f882920803aab2320cbb613b", "Genus"] <- "g__Faecousia"

# Reassign the updated taxonomy table to the phyloseq object
tax_table(ps) <- tax_table(tax_mat)

# Check that the taxonomy table is valid
tax_table(ps)["949d11292146e12904547a6db61ce024", ]

# Aggregate reads to Genus level
ps <- tax_glom(ps,"Genus")

# Load the indicator taxa results
ISA_results<- read.csv("ISA_ryanfilt_new.csv")

# Filter the data
filtered_data <- ISA_results %>%
  filter(stat > 0.5 & (s.Crohn.s.Disease == 1 | s.Healthy == 1))

# Save the filtered file
write.csv(filtered_data, "filtered_data.csv", row.names = FALSE)

# Modify the Genus column
filtered_data_with_genus <- filtered_data %>%
  mutate(Genus = case_when(
    ASV == "949d11292146e12904547a6db61ce024" ~ "g__Laedolimicola",
    ASV == "a716b99232fa3a34ee30ae758b9216c1" ~ "g__Beduinella",
    ASV == "0e6b48ffce8c52ae15a77ccd71ed5f31" ~ "g__Pseudoruminococcus",
    ASV == "42a8749495523917f1ff19114d6f374b" ~ "g__Bifidobacterium",
    ASV == "30753be97d9d33c8c68f245f6c0466e1" ~ "g__Brotomerdimonas",
    ASV == "007f2738f882920803aab2320cbb613b" ~ "g__Faecousia",
    TRUE ~ Genus
  ))

# Save the revised file
write.csv(filtered_data_with_genus, "filtered_data_with_genus.csv", row.names = FALSE)

# Normalize microbiome data and select only the significant taxa
# CLR transformation
ps_clr = ps %>% microbiome::transform('clr') 
# Filter taxa
ps_filt = prune_taxa(filtered_data_with_genus$ASV,ps_clr)

# Make sure ps_filt uses the curated Genus names
tax_table(ps_filt)[, "Genus"] <- filtered_data_with_genus$Genus[
  match(rownames(tax_table(ps_filt)), filtered_data_with_genus$ASV)
]
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


ggsave("ROC_Plot_7.png", roc_plot, width = 8, height = 6, dpi = 300)


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
ggsave("Importance_Plot_7.png", importance_plot, width = 8, height = 6, dpi = 300)

pd_model$importance
#Check Correlation

temp1 = ps_clr %>% 
  subset_taxa(Genus=='g__Laedolimicola') %>% 
  psmelt()
temp1$Condition <- factor(temp1$Condition, 
                          levels = c("Healthy", "Crohn's Disease"))
box_plot_Laedolimicola <- temp1 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Laedolimicola_abundance.png", box_plot_Laedolimicola, width = 8, height = 6, dpi = 300)

temp2 = ps_clr %>% 
  subset_taxa(Genus=='g__Beduinella') %>% 
  psmelt()
temp2$Condition <- factor(temp2$Condition, 
                          levels = c("Healthy", "Crohn's Disease"))
box_plot_Beduinella <- temp2 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Beduinella_abundance.png", box_plot_Beduinella, width = 8, height = 6, dpi = 300)

temp3 = ps_clr %>% 
  subset_taxa(Genus=='g__Akkermansia') %>% 
  psmelt()
temp3$Condition <- factor(temp3$Condition, 
                          levels = c("Healthy", "Crohn's Disease"))
box_plot_Akkermansia <- temp3 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Akkermansia_abundance.png", box_plot_Akkermansia, width = 8, height = 6, dpi = 300)

temp4 = ps_clr %>% 
  subset_taxa(Genus=='g__Brotomerdimonas') %>% 
  psmelt()
temp4$Condition <- factor(temp4$Condition, 
                          levels = c("Healthy", "Crohn's Disease"))
box_plot_Brotomerdimonas <- temp4 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Brotomerdimonas_abundance.png", box_plot_Brotomerdimonas, width = 8, height = 6, dpi = 300)

temp5 = ps_clr %>% 
  subset_taxa(Genus=='g__Desulfovibrio') %>% 
  psmelt()
temp5$Condition <- factor(temp5$Condition, 
                          levels = c("Healthy", "Crohn's Disease"))
box_plot_Desulfovibrio <- temp5 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Desulfovibrio_abundance.png", box_plot_Desulfovibrio, width = 8, height = 6, dpi = 300)

temp6 = ps_clr %>% 
  subset_taxa(Genus=='g__Veillonella') %>% 
  psmelt()
temp6$Condition <- factor(temp6$Condition, 
                          levels = c("Healthy", "Crohn's Disease"))
box_plot_Veillonella <- temp6 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Veillonella_abundance.png", box_plot_Veillonella, width = 8, height = 6, dpi = 300)

temp7 = ps_clr %>% 
  subset_taxa(Genus=='g__Pseudoruminococcus') %>% 
  psmelt()
temp7$Condition <- factor(temp7$Condition, 
                          levels = c("Healthy", "Crohn's Disease"))
box_plot_Pseudoruminococcus <- temp7 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Pseudoruminococcus_abundance.png", box_plot_Pseudoruminococcus, width = 8, height = 6, dpi = 300)

temp8 = ps_clr %>% 
  subset_taxa(Genus=='g__Hydrogenoanaerobacterium') %>% 
  psmelt()
temp8$Condition <- factor(temp8$Condition, 
                          levels = c("Healthy", "Crohn's Disease"))
box_plot_Hydrogenoanaerobacterium <- temp8 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Hydrogenoanaerobacterium_abundance.png", box_plot_Hydrogenoanaerobacterium, width = 8, height = 6, dpi = 300)

temp9 = ps_clr %>% 
  subset_taxa(Genus=='g__Fournierella') %>% 
  psmelt()
temp9$Condition <- factor(temp9$Condition, 
                          levels = c("Healthy", "Crohn's Disease"))
box_plot_Fournierella <- temp9 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Fournierella_abundance.png", box_plot_Fournierella, width = 8, height = 6, dpi = 300)

temp10 = ps_clr %>% 
  subset_taxa(Genus=='g__Howardella') %>% 
  psmelt()
temp10$Condition <- factor(temp10$Condition, 
                          levels = c("Healthy", "Crohn's Disease"))
box_plot_Howardella <- temp10 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Howardella_abundance.png", box_plot_Howardella, width = 8, height = 6, dpi = 300)

temp11 = ps_clr %>% 
  subset_taxa(Genus=='g__Faecousia') %>% 
  psmelt()
temp11$Condition <- factor(temp11$Condition, 
                           levels = c("Healthy", "Crohn's Disease"))
box_plot_Faecousia <- temp11 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Faecousia_abundance.png", box_plot_Faecousia, width = 8, height = 6, dpi = 300)

temp12 = ps_clr %>% 
  subset_taxa(Genus=='g__Bifidobacterium') %>% 
  psmelt()
temp12$Condition <- factor(temp12$Condition, 
                           levels = c("Healthy", "Crohn's Disease"))
box_plot_Bifidobacterium <- temp12 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Bifidobacterium_abundance.png", box_plot_Bifidobacterium, width = 8, height = 6, dpi = 300)

temp13 = ps_clr %>% 
  subset_taxa(Genus=='g__Turicibacter') %>% 
  psmelt()
temp13$Condition <- factor(temp13$Condition, 
                           levels = c("Healthy", "Crohn's Disease"))
box_plot_Turicibacter <- temp13 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Turicibacter_abundance.png", box_plot_Turicibacter, width = 8, height = 6, dpi = 300)

temp14 = ps_clr %>% 
  subset_taxa(Genus=='g__Senegalimassilia') %>% 
  psmelt()
temp14$Condition <- factor(temp14$Condition, 
                           levels = c("Healthy", "Crohn's Disease"))
box_plot_Senegalimassilia <- temp14 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Senegalimassilia_abundance.png", box_plot_Senegalimassilia, width = 8, height = 6, dpi = 300)

temp15 = ps_clr %>% 
  subset_taxa(Genus=='g__Slackia') %>% 
  psmelt()
temp15$Condition <- factor(temp15$Condition, 
                           levels = c("Healthy", "Crohn's Disease"))
box_plot_Slackia <- temp15 %>% 
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Crohn's Disease" = "#D55E00",
    "Healthy" = "#7570B3"
  )) +
  labs(x = NULL, y = "Abundance", fill = NULL) +
  theme_classic(base_size = 18)
ggsave("Slackia_abundance.png", box_plot_Slackia, width = 8, height = 6, dpi = 300)

### Metadata only 
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
filt_clinical_model$importance


### Combined microbial taxa and metadata 
#### Loading packages and data####
library(caret)
library(ranger)
library(pROC)
library(boot)
library(ggplot2)
library(phyloseq)
library(tidyverse)
library(randomForest)

#### Formatting/Preparing files for RF####
# Load your saved phyloseq object
load ("ryan_filt_new.RData")
ps <- ryan_ps_filt
# Convert taxonomy table to matrix and assign genus to your ASV
tax_mat <- as(tax_table(ps), "matrix")
tax_mat["949d11292146e12904547a6db61ce024", "Genus"] <- "g__Laedolimicola"
tax_mat["a716b99232fa3a34ee30ae758b9216c1", "Genus"] <- "g__Beduinella"
tax_mat["0e6b48ffce8c52ae15a77ccd71ed5f31", "Genus"] <- "g__Pseudoruminococcus"
tax_mat["42a8749495523917f1ff19114d6f374b", "Genus"] <- "g__Bifidobacterium"
tax_mat["30753be97d9d33c8c68f245f6c0466e1", "Genus"] <- "g__Brotomerdimonas"
tax_mat["007f2738f882920803aab2320cbb613b", "Genus"] <- "g__Faecousia"

# Reassign the updated taxonomy table to the phyloseq object
tax_table(ps) <- tax_table(tax_mat)

# Check that the taxonomy table is valid
tax_table(ps)["949d11292146e12904547a6db61ce024", ]

# Aggregate reads to Genus level
ps <- tax_glom(ps,"Genus")

# Load the indicator taxa results
ISA_results<- read.csv("ISA_ryanfilt_new.csv")

# Filter the data
filtered_data <- ISA_results %>%
  filter(stat > 0.5 & (s.Crohn.s.Disease == 1 | s.Healthy == 1))

# Save the filtered file
write.csv(filtered_data, "filtered_data.csv", row.names = FALSE)

# Modify the Genus column
filtered_data_with_genus <- filtered_data %>%
  mutate(Genus = case_when(
    ASV == "949d11292146e12904547a6db61ce024" ~ "g__Laedolimicola",
    ASV == "a716b99232fa3a34ee30ae758b9216c1" ~ "g__Beduinella",
    ASV == "0e6b48ffce8c52ae15a77ccd71ed5f31" ~ "g__Pseudoruminococcus",
    ASV == "42a8749495523917f1ff19114d6f374b" ~ "g__Bifidobacterium",
    ASV == "30753be97d9d33c8c68f245f6c0466e1" ~ "g__Brotomerdimonas",
    ASV == "007f2738f882920803aab2320cbb613b" ~ "g__Faecousia",
    TRUE ~ Genus
  ))



# Save the revised file
write.csv(filtered_data_with_genus, "filtered_data_with_genus.csv", row.names = FALSE)

# Normalize microbiome data and select only the significant taxa
# CLR transformation
ps_clr = ps %>% microbiome::transform('clr') 
# Filter taxa
ps_filt = prune_taxa(filtered_data_with_genus$ASV,ps_clr)

# Make sure ps_filt uses the curated Genus names
tax_table(ps_filt)[, "Genus"] <- filtered_data_with_genus$Genus[
  match(rownames(tax_table(ps_filt)), filtered_data_with_genus$ASV)
]
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
  select(Sample,Condition,Genus,Abundance,Gender,Biopsy_location,Smoking.status,Age_when_sampled) %>% 
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
set.seed(065)
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

####Run RF####
source('randomforest_functions.R')
pd_model = run_rf(X = predictors, y = outcome, 
                  fold_list = folds,
                  hyper = tune_grid, 
                  rngseed = 065)
names(pd_model)

####Interpretation####

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


ggsave("ROC_Plot_9_Final.png", roc_plot, width = 8, height = 6, dpi = 300)


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
ggsave("Importance_Plot_9_Final.png", importance_plot, width = 8, height = 6, dpi = 300)
pd_model$importance
#Check Correlation

temp = ps_clr %>% 
  subset_taxa(Genus=='g__Laedolimicola') %>% 
  psmelt()

box_plot_RF <- temp %>% 
  ggplot(aes(Condition,Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height=2, width = 0.2) +
  theme_classic(base_size=18)

ggsave("Laedolimicola_abundance.png", box_plot_RF, width = 8, height = 6, dpi = 300)
