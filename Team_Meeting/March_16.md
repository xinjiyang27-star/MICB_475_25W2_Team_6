# March 16 [16:00 - ] 
## Meeting Agenda 
1. Discuss Random Forest Model

## Q and A
1. For ISA:
   - Why are we not getting consistent ISA results between runs with the same labtop 
   - 
3. For Random forest:
   - Why all the modes are overfitted?
   - Why the models are not consistent between runs/laptops?
   - More information on mtry and min.node.size? How to set the most appropriate parameters?
   - Which code is better:
     1) df_pivot = filtered_df %>% 
                    select(Sample,Condition,Genus,Abundance) %>% 
                    # Turn each Genus into its own column
                    pivot_wider(names_from = Genus, values_from = Abundance)

     2) From Tara 
   - 
  
   

## Experimental Results Summary
