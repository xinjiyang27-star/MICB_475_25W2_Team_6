# Indicator Species/Taxa Analysis [27 February]

## Purpose:
1. To run ISA to identify the indicator species for each condition

## Materials:
1. R-studio

## Methods:  
library(tidyverse)
library(phyloseq)
library(indicspecies)

#Load the nonrarefied object
load("ryan_filt.RData")

#### Indicator Species/Taxa Analysis ####
#glom to Genus

#Group data based on particular taxanomic rank (genus)

ryan_genus <- tax_glom(ryan_filt, "Genus", NArm = FALSE)

#Convert counts to relative abundance

ryan_genus_RA <- transform_sample_counts(ryan_genus, fun=function(x) x/sum(x))

#ISA

#Calls th OTU table and transpose it. Calculated that IVs for all ASVs

isa_ryan <- multipatt(t(otu_table(ryan_genus_RA)), cluster = sample_data(ryan_genus_RA)$`Condition`)

#Anything less than 0.05 will come up as indicator species

summary(isa_ryan)

#Create a taxanomy table, but ASVs are not the row names (separte column)

taxtable <- tax_table(ryan_filt) %>% as.data.frame() %>% rownames_to_column(var="ASV")

#Combine the taxanomic info to the list of indicator species

isa_table <- isa_ryan$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View()

#Save the table

write.csv(isa_table, "ISA_significant_taxa.csv", row.names = FALSE)
