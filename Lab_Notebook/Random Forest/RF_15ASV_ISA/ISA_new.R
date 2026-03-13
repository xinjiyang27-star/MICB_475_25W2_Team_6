library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
# Load the nonrarefied object
load("ryan_filt_new.RData")

#### Indicator Species/Taxa Analysis ####
# glom to Genus
# Group data based on particular taxanomic rank (genus)
ryan_genus <- tax_glom(ryan_filt_new, "Genus", NArm = FALSE)
# Convert counts to relative abundance
ryan_genus_RA <- transform_sample_counts(ryan_genus, fun=function(x) x/sum(x))

#ISA
# Calls th OTU table and transpose it
# Calculated that IVs for all ASVs
isa_ryan <- multipatt(t(otu_table(ryan_genus_RA)), cluster = sample_data(ryan_genus_RA)$`Condition`)
# Anything less than 0.05 will come up as indicator species
summary(isa_ryan)
# Create a taxanomy table, but ASVs are not the row names (separte column)
taxtable <- tax_table(ryan_filt_new) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# consider that your table is only going to be resolved up to the genus level, be wary of anything beyond the glomed taxa level
# Combine the taxanomic info to the list of indicator species
isa_table <- isa_ryan$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05)
View(isa_table)

isa_cd <- isa_table %>%
  filter(`s.Crohn's Disease` == 1)

# Save the tables
write.csv(isa_table, "ISA_ryanfilt_new.csv", row.names = FALSE)
write.csv(isa_cd, "ISA_onlyCD.csv", row.names = FALSE)