library(tidyverse)
library(phyloseq)
library(indicspecies)
library(dplyr)

#### Load data ####
# Load the nonrarefied object
load("Ryan_filt_new.RData")

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

isa_table_stat0.5 <- isa_table %>%
  filter(stat>0.5)

# Save the tables
write.csv(isa_table, "ISA_ryanfilt_new.csv", row.names = FALSE)
write.csv(isa_table_stat0.5, "ISA_15_taxa.csv", row.names = FALSE)

# Remove all g__, f__, p__, etc from the table
isa_table_stat0.5 <- isa_table_stat0.5 %>%
  mutate(across(
    c(Domain, Phylum, Class, Order, Family, Genus, Species),
    ~ stringr::str_remove(., "^[a-z]__")
  ))

# Replace the six ASVs's genus with formal genus name
isa_table_stat0.5 <- isa_table_stat0.5 %>%
  mutate(Genus = case_when(
    ASV == "949d11292146e12904547a6db61ce024" ~ "Laedolimicola",
    ASV == "a716b99232fa3a34ee30ae758b9216c1" ~ "Beduinella",
    ASV == "0e6b48ffce8c52ae15a77ccd71ed5f31" ~ "Pseudoruminococcus",
    ASV == "42a8749495523917f1ff19114d6f374b" ~ "Bifidobacterium",
    ASV == "30753be97d9d33c8c68f245f6c0466e1" ~ "Brotomerdimonas",
    ASV == "007f2738f882920803aab2320cbb613b" ~ "Faecousia",
    TRUE ~ Genus
  ))

isa_pub <- isa_table_stat0.5

# Rename indicator columns safely
isa_pub <- isa_pub %>%
  rename(
    CD = matches("Crohn"),
    Healthy = matches("Healthy")
  )

# Add group column
isa_pub <- isa_pub %>%
  mutate(Group = case_when(
    CD == 1 ~ "Crohn's Disease",
    Healthy == 1 ~ "Healthy",
    TRUE ~ NA_character_
  )) %>%
  select(Group, Phylum, Family, Genus, stat, p.value) %>%
  rename(
    `Indicator Value` = stat,
    `p-value` = p.value
  ) %>%
  arrange(Group, desc(`Indicator Value`))

# Split into two sections
cd_table <- isa_pub %>%
  filter(Group == "Crohn's Disease") %>%
  select(-Group)

healthy_table <- isa_pub %>%
  filter(Group == "Healthy") %>%
  select(-Group)

# Add group label only once
cd_table <- cd_table %>%
  mutate(Group = "") %>%
  relocate(Group) %>%
  as.data.frame()

healthy_table <- healthy_table %>%
  mutate(Group = "") %>%
  relocate(Group) %>%
  as.data.frame()

cd_table$Group[1] <- "Crohn's Disease"
healthy_table$Group[1] <- "Healthy"

# Combine into one large table
final_table <- bind_rows(cd_table, blank_row, healthy_table)

# Remove the two empty columns
final_table <- final_table %>%
  select(-Indicator.Value, -p.value)

# Remove blank row
final_table <- bind_rows(cd_table, healthy_table)

# Save the table
write.csv(final_table, "ISA_final_grouped_table.csv", row.names = FALSE)