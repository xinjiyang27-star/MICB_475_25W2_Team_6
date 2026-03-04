library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
# Load the nonrarefied object
load("ryan_filt.RData")

#### Indicator Species/Taxa Analysis ####
# glom to Genus
# Group data based on particular taxanomic rank (genus)
ryan_genus <- tax_glom(ryan_filt, "Genus", NArm = FALSE)
# Convert counts to relative abundance
ryan_genus_RA <- transform_sample_counts(ryan_genus, fun=function(x) x/sum(x))

#ISA
# Calls th OTU table and transpose it
# Calculated that IVs for all ASVs
isa_ryan <- multipatt(t(otu_table(ryan_genus_RA)), cluster = sample_data(ryan_genus_RA)$`Condition`)
# Anything less than 0.05 will come up as indicator species
summary(isa_ryan)
# Create a taxanomy table, but ASVs are not the row names (separte column)
taxtable <- tax_table(ryan_filt) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# consider that your table is only going to be resolved up to the genus level, be wary of anything beyond the glomed taxa level
# Combine the taxanomic info to the list of indicator species
isa_table_all <- isa_ryan$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05)
View(isa_table_all)

# Filter out ASVs that are shared between UC+healthy and CD+healthy
isa_table_unique_UCCD <- isa_table_all %>%
  filter(
    !(
      (s.Healthy == 1 & `s.Ulcerative Colitis` == 1 & `s.Crohn's Disease` == 0) |  # Healthy+UC only
        (s.Healthy == 1 & `s.Crohn's Disease` == 1 & `s.Ulcerative Colitis` == 0)    # Healthy+CD only
    )
  ) 
View(isa_table_unique_UCCD)

# Create a new column that labels each ASV based on which condition(s) it belongs to.
isa_unique_uccd_labeled <- isa_table_unique_UCCD %>%
  mutate(group_label = case_when(
    s.Healthy == 0 & `s.Ulcerative Colitis` == 1 & `s.Crohn's Disease` == 0 ~ "UC_only",
    s.Healthy == 0 & `s.Ulcerative Colitis` == 0 & `s.Crohn's Disease` == 1 ~ "CD_only",
    s.Healthy == 0 & `s.Ulcerative Colitis` == 1 & `s.Crohn's Disease` == 1 ~ "UC+CD_shared",
    s.Healthy == 1 & `s.Ulcerative Colitis` == 1 & `s.Crohn's Disease` == 1 ~ "Healthy+UC+CD",
    s.Healthy == 1 & `s.Ulcerative Colitis` == 0 & `s.Crohn's Disease` == 0 ~ "Healthy_only",
    TRUE ~ "Other"
  ))

# Extracts the ASV IDs into separate lists based on the group category
asv_uc_only  <- isa_unique_uccd_labeled %>% filter(group_label == "UC_only") %>% pull(ASV)
asv_cd_only  <- isa_unique_uccd_labeled %>% filter(group_label == "CD_only") %>% pull(ASV)
asv_uc_cd    <- isa_unique_uccd_labeled %>% filter(group_label == "UC+CD_shared") %>% pull(ASV)
asv_all3     <- isa_unique_uccd_labeled %>% filter(group_label == "Healthy+UC+CD") %>% pull(ASV)

asv_uc_only
asv_cd_only
asv_uc_cd
asv_all3

# Save the table
write.csv(isa_table_unique_UCCD, "ISA_significant_taxa_unique_uccd.csv", row.names = FALSE)

# Create a table only contains ASVs related to CD and Healthy
isa_CD_Healthy <- isa_table_all %>%
  filter(`s.Ulcerative Colitis` == 0)
View(isa_CD_Healthy)

# Create a new column that labels each ASV based on which condition(s) it belongs to.
isa_CD_Healthy_labeled <- isa_CD_Healthy %>%
  mutate(group_label = case_when(
    s.Healthy == 1 & `s.Crohn's Disease` == 0 ~ "Healthy_only",
    s.Healthy == 0 & `s.Crohn's Disease` == 1 ~ "CD_only",
    s.Healthy == 1 & `s.Crohn's Disease` == 1 ~ "CD+Healthy_shared"
  ))

# Extracts the ASV IDs into separate lists based on the group category
asv_cd_only <- isa_CD_Healthy_labeled %>% filter(group_label == "CD_only") %>% pull(ASV)
asv_healthy_only <- isa_CD_Healthy_labeled %>% filter(group_label == "Healthy_only") %>% pull(ASV)
asv_cd_healthy <- isa_CD_Healthy_labeled %>% filter(group_label == "CD+Healthy_shared") %>% pull(ASV)

asv_cd_only
asv_healthy_only
asv_cd_healthy

# Save the table
write.csv(isa_CD_Healthy, "ISA_significant_taxa_cd_healthy.csv", row.names = FALSE)