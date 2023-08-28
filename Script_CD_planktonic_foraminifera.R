# Script written by Raphael Morard, Christiane Hassenrück and Chiara Vanni
# Data tidying
library(tidyverse)
library(dplyr)
library(textshape)

library(igraph)
library(reshape)
library(CINNA)

library("dada2")
library("seqRFLP")

library("data")

# For PCA
library(factoextra)

# alternative ordination
require(umap)

# To produce maps
library(maps)
library(parallel)
library(hexbin)
library(scales)

library("janitor")

# Loading the data
Data_sanger                         <- readxl::read_xlsx(path = "Table S1.xlsx", sheet = 1,skip = 1, col_types="guess") %>% arrange(Internal_ID)
RFLP_data                           <- readxl::read_xlsx(path = "Table S2.xlsx", sheet = 1, col_types="guess") 

Geography_Curated                   <- readxl::read_xlsx(path = "Table S4.xlsx", sheet = 1, col_types="guess")

Molecular_Nomenclature              <- readxl::read_xlsx(path = "Table S5.xlsx", sheet = 1, col_types="guess")
Equivalence_Nomenclature            <- readxl::read_xlsx(path = "Table S5.xlsx", sheet = 3, col_types="guess")

ASV_Occurence                       <- readxl::read_xlsx(path = "Table S6.xlsx", sheet = 1, col_types="guess")
PF_ASV_seq                          <- readxl::read_xlsx(path = "Table S6.xlsx", sheet = 2, col_types="guess")
Molecular_Nomenclature_Novel_ASVs   <- readxl::read_xlsx(path = "Table S6.xlsx", sheet = 3, col_types="guess")

# importing dataset # The order of the row is sorted using the unique internal provided for each sequence to make sure that the table is processed in the same way at each iteration
# Exact same thing as Script 1 to rebuilt the same files
# Assess quality of each individual region to select Basetypes
# Subset the dataset for the coverage and sequence of each region

Data_sanger_subset_37F <- Data_sanger %>% 
  select(`Internal_ID`,`37F_coverage`,`37F_seq`) %>% 
  filter(`37F_coverage` == "Complete") %>% 
  group_by(`37F_seq`) %>%
  add_count() %>% 
  dplyr::rename("N_obs_37F" = "n") %>%
  ungroup(`37F_seq`) %>% 
  select(`Internal_ID`, N_obs_37F)

Data_sanger_subset_3741 <- Data_sanger %>% 
  select(`Internal_ID`,`37-41_coverage`,`37-41_seq`) %>% 
  filter(`37-41_coverage` == "Complete") %>% 
  group_by(`37-41_seq`) %>%
  add_count() %>% 
  dplyr::rename("N_obs_3741" = "n") %>%
  ungroup(`37-41_seq`) %>% 
  select(`Internal_ID`, N_obs_3741)

Data_sanger_subset_41F <- Data_sanger %>% 
  select(`Internal_ID`,`41F_coverage`,`41F_seq`) %>% 
  filter(`41F_coverage` == "Complete") %>% 
  group_by(`41F_seq`) %>%
  add_count() %>% 
  dplyr::rename("N_obs_41F" = "n") %>%
  ungroup(`41F_seq`) %>% 
  select(`Internal_ID`, N_obs_41F)

Data_sanger_subset_3943 <- Data_sanger %>% 
  select(`Internal_ID`,`39-43_coverage`,`39-43_seq`) %>% 
  filter(`39-43_coverage` == "Complete") %>% 
  group_by(`39-43_seq`) %>%
  add_count() %>% 
  dplyr::rename("N_obs_3943" = "n") %>%
  ungroup(`39-43_seq`) %>% 
  select(`Internal_ID`, N_obs_3943)

Data_sanger_subset_43E <- Data_sanger %>% 
  select(`Internal_ID`,`43E_coverage`,`43E_seq`) %>% 
  filter(`43E_coverage` == "Complete") %>% 
  group_by(`43E_seq`) %>%
  add_count() %>% 
  dplyr::rename("N_obs_43E" = "n") %>%
  ungroup(`43E_seq`) %>% 
  select(`Internal_ID`, N_obs_43E)

Data_sanger_subset_4445 <- Data_sanger %>% 
  select(`Internal_ID`,`44-45_coverage`,`44-45_seq`) %>% 
  filter(`44-45_coverage` == "Complete") %>% 
  group_by(`44-45_seq`) %>%
  add_count() %>% 
  dplyr::rename("N_obs_4445" = "n") %>%
  ungroup(`44-45_seq`) %>% 
  select(`Internal_ID`, N_obs_4445)

Data_sanger_subset_45E47F <- Data_sanger %>% 
  select(`Internal_ID`,`45E-47F_coverage`,`45E-47F_seq`) %>% 
  filter(`45E-47F_coverage` == "Complete") %>% 
  group_by(`45E-47F_seq`) %>%
  add_count() %>% 
  dplyr::rename("N_obs_45E47F" = "n") %>%
  ungroup(`45E-47F_seq`) %>% 
  select(`Internal_ID`, N_obs_45E47F)

Data_sanger_subset_4749 <- Data_sanger %>% 
  select(`Internal_ID`,`47-49_coverage`,`47-49_seq`) %>% 
  filter(`47-49_coverage` == "Complete") %>% 
  group_by(`47-49_seq`) %>%
  add_count() %>% 
  dplyr::rename("N_obs_4749" = "n") %>%
  ungroup(`47-49_seq`) %>% 
  select(`Internal_ID`, N_obs_4749)

Data_sanger_subset_49E <- Data_sanger %>% 
  select(`Internal_ID`,`49E_coverage`,`49E_seq`) %>% 
  filter(`49E_coverage` == "Complete") %>% 
  group_by(`49E_seq`) %>%
  add_count() %>% 
  dplyr::rename("N_obs_49E" = "n") %>%
  ungroup(`49E_seq`) %>% 
  select(`Internal_ID`, N_obs_49E)

# Assemble all table to obtain the observation of each region for each sequence

Data_sanger_quality <- full_join(
  full_join(
    full_join(full_join(Data_sanger_subset_37F,Data_sanger_subset_3741),
              full_join(Data_sanger_subset_41F,Data_sanger_subset_3943))
    ,
    full_join(full_join(Data_sanger_subset_43E,Data_sanger_subset_4445),
              full_join(Data_sanger_subset_45E47F,Data_sanger_subset_4749)))
  ,Data_sanger_subset_49E) %>% replace(is.na(.), 0)

# Count the number of regions observed at least 3 times 
counts <-  data.frame(`sum_2` =rowSums(Data_sanger_quality %>% select(-Internal_ID) > 2))

# Subset only sequences for which all regions have been observed thrice
Data_sanger_selection <- cbind(Data_sanger_quality,counts) %>% 
  filter(`sum_2` == 9) %>% 
  select(Internal_ID)

# Assemble the basetypes and subset relevant metadata for taxonomy
Selection_holotypes <- inner_join(Data_sanger_selection, Data_sanger) %>% 
  unite("Basetype", `37F_seq`:`49E_seq`, remove = FALSE,sep= "") %>% 
  select(Internal_ID,Voucher_code, Basetype, Clade_curated, Genus_curated, species_curated, Quality, Clade_delimitation)

# add to the data the representative sequences of the genotypes that did not go through the filtering
Selection_representative <-Data_sanger %>% 
  filter(Quality=="REPRESENTATIVE") %>% 
  unite("Basetype", `37F_seq`:`49E_seq`, remove = FALSE,sep= "") %>% 
  select(Internal_ID,Voucher_code, Basetype, Clade_curated, Genus_curated, species_curated, Quality, Clade_delimitation)

# Provide names for unique basetypes
Naming_Basetype <- Selection_holotypes %>% 
  select(Basetype) %>% 
  distinct() %>%
  mutate(rank = row_number()) %>% 
  as.data.frame() %>% 
  add_column(Dummy = "Basetype") %>% 
  unite("Basetype_ID", Dummy:rank)

# Provide names for representative sequences
Naming_representative <- Selection_representative %>% 
  select(Basetype) %>% 
  distinct() %>%
  mutate(rank = row_number()) %>% 
  as.data.frame() %>% 
  add_column(Dummy = "Representative") %>% 
  unite("Basetype_ID", Dummy:rank)

# Assign basetype ID to each select sequences and aggregates tables with BASETYPES and REPRESENTATIVE sequences
Occurence_basetypes_specimens <- left_join(Selection_holotypes, Naming_Basetype) %>% 
  rbind(left_join(Selection_representative, Naming_representative))

#Produce co occurence matrix of the basetypes between the specimens
Occurence_matrix_basetypes <- Occurence_basetypes_specimens %>% 
  select(`Voucher_code`, Basetype_ID) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(`Voucher_code`,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1)

# definition of basegroups
test_in <- Occurence_matrix_basetypes
colnames(test_in) <- paste0("S_", colnames(test_in))
test_in$bt <- rownames(test_in)
test_out <- reshape::melt(test_in, id.vars = "bt") %>% filter(value != 0)
test_graph <- graph_from_data_frame(test_out[, 1:2], directed = F)
plot(
  test_graph, 
  layout = layout_with_graphopt, 
  vertex.color = as.numeric(names(V(test_graph)) %in% rownames(Occurence_matrix_basetypes)) + 1,
  vertex.size = 2,
  vertex.label = NA, 
  edge.color = "grey",
  main = ""
)
test_split <- graph_extract_components(test_graph)
test_vnames <- lapply(test_split, function(x) names(V(x)))
test_bt <- lapply(test_vnames, function(x) grep("^S_", x, value = T, invert = T))
names(test_bt) <- paste0("basegroup_", 1:length(test_bt))
tmp <- data.frame(
  basegroup = rep(names(test_bt), sapply(test_bt, length)),
  basetype = unlist(test_bt)
) %>% column_to_rownames("basetype")
Occurence_basetypes_specimens$Basegroup_ID <- tmp[Occurence_basetypes_specimens$Basetype_ID, 1]

####### Subseting data to produce fasta file ######

Curated_sequences <- Occurence_basetypes_specimens %>% 
  distinct(Basetype_ID, .keep_all = TRUE) %>% 
  unite("Name_seq", Basetype_ID,Basegroup_ID,Genus_curated,species_curated, remove = FALSE,sep= "|")

        ################################################################
        #                                                              #
        # Beginning of script 3 with import of Molecular Nomenclature  #
        #                                                              #
        ################################################################


# Merge all table to have the molecular nomenclature together with the sequences
Key_assignement <- full_join(Curated_sequences, Molecular_Nomenclature %>% select(Basetype_ID, MOTU_lvl_1, MOTU_lvl_2, MOTU_lvl_3, MOTUs)) %>% 
  left_join(Data_sanger %>% select(Internal_ID,`37F_seq`))     %>% 
  left_join(Data_sanger %>% select(Internal_ID,`37-41_seq`))   %>%
  left_join(Data_sanger %>% select(Internal_ID,`41F_seq`))     %>% 
  left_join(Data_sanger %>% select(Internal_ID,`39-43_seq`))   %>% 
  left_join(Data_sanger %>% select(Internal_ID,`43E_seq`))     %>% 
  left_join(Data_sanger %>% select(Internal_ID,`44-45_seq`))   %>% 
  left_join(Data_sanger %>% select(Internal_ID,`45E-47F_seq`)) %>% 
  left_join(Data_sanger %>% select(Internal_ID,`47-49_seq`))   %>%
  left_join(Data_sanger %>% select(Internal_ID,`49E_seq`)) 

# Now the game is to find which sequence portion will be diagnostic for a given taxonomic unit and use it to re assign the rest of the database.
# We start from the finest level to the coarsest

# Now we produce a file with only necessary informations
Key_assignement_MOTUS_lvl3 <- Key_assignement %>% 
  unite("Taxo_lvl3", Clade_curated,Genus_curated,species_curated, MOTU_lvl_1, MOTU_lvl_2, MOTU_lvl_3, remove = FALSE,sep= "|") %>%
  select(Taxo_lvl3,  `37F_seq`:`49E_seq`)
  
# Identify individual sequence pattern that occur only in a single MOTU lvl3 and that is therefore diagnostic

List_37F_Diagnostic_MOTU3 <- Key_assignement_MOTUS_lvl3 %>% 
  select(Taxo_lvl3, `37F_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl3,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("37F_seq") %>% 
  select(`37F_seq`) %>% left_join(Key_assignement_MOTUS_lvl3 %>% select(Taxo_lvl3, `37F_seq`) %>% distinct(`37F_seq`, .keep_all =TRUE)) 

List_41F_Diagnostic_MOTU3 <- Key_assignement_MOTUS_lvl3 %>% 
  select(Taxo_lvl3, `41F_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl3,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("41F_seq") %>% 
  select(`41F_seq`) %>% left_join(Key_assignement_MOTUS_lvl3 %>% select(Taxo_lvl3, `41F_seq`) %>% distinct(`41F_seq`, .keep_all =TRUE)) 

List_43E_Diagnostic_MOTU3 <- Key_assignement_MOTUS_lvl3 %>% 
  select(Taxo_lvl3, `43E_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl3,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("43E_seq") %>% 
  select(`43E_seq`) %>% left_join(Key_assignement_MOTUS_lvl3 %>% select(Taxo_lvl3, `43E_seq`) %>% distinct(`43E_seq`, .keep_all =TRUE)) 

List_45E47F_Diagnostic_MOTU3 <- Key_assignement_MOTUS_lvl3 %>% 
  select(Taxo_lvl3, `45E-47F_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl3,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("45E-47F_seq") %>% 
  select(`45E-47F_seq`) %>% left_join(Key_assignement_MOTUS_lvl3 %>% select(Taxo_lvl3, `45E-47F_seq`) %>% distinct(`45E-47F_seq`, .keep_all =TRUE)) 

List_49E_Diagnostic_MOTU3 <- Key_assignement_MOTUS_lvl3 %>% 
  select(Taxo_lvl3, `49E_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl3,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("49E_seq") %>% 
  select(`49E_seq`) %>% left_join(Key_assignement_MOTUS_lvl3 %>% select(Taxo_lvl3, `49E_seq`) %>% distinct(`49E_seq`, .keep_all =TRUE))

# Re assignemtn of the database by matching diagnostic sequence (Either entire basetype or single region) incrementaly. Starting from basetype then each region, filter out the part of the database that is not assigned every time

##################### First Basetype level #################

DB_assigned_Basetype <-  Data_sanger %>%
  unite("Basetype", `37F_seq`:`49E_seq`, remove = FALSE,sep= "") %>%
  left_join(Key_assignement_MOTUS_lvl3 %>% unite("Basetype", `37F_seq`:`49E_seq`, remove = TRUE,sep= "") %>% 
              select(Basetype,Taxo_lvl3) %>% 
              add_column(Quality_assignement = "Basetype")) 

DB_assigned_Basetype$Quality_assignement <- DB_assigned_Basetype$Quality_assignement %>% replace(is.na(.), 0)

################### Assignement lvl 3 for each diagnostic sequence fragment ############

DB_assigned_lvl3_37F <- DB_assigned_Basetype %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl3, -Quality_assignement) %>% 
  left_join(List_37F_Diagnostic_MOTU3 %>% add_column(Quality_assignement = "Assignement_lvl3"))
DB_assigned_lvl3_37F$Quality_assignement <- DB_assigned_lvl3_37F$Quality_assignement %>% replace(is.na(.), 0)

DB_assigned_lvl3_41F <- DB_assigned_lvl3_37F %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl3, -Quality_assignement) %>% 
  left_join(List_41F_Diagnostic_MOTU3 %>% add_column(Quality_assignement = "Assignement_lvl3"))
DB_assigned_lvl3_41F$Quality_assignement <- DB_assigned_lvl3_41F$Quality_assignement %>% replace(is.na(.), 0)

DB_assigned_lvl3_43E <- DB_assigned_lvl3_41F %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl3, -Quality_assignement) %>% 
  left_join(List_43E_Diagnostic_MOTU3 %>% add_column(Quality_assignement = "Assignement_lvl3"))
DB_assigned_lvl3_43E$Quality_assignement <- DB_assigned_lvl3_43E$Quality_assignement %>% replace(is.na(.), 0)

DB_assigned_lvl3_45E47F <- DB_assigned_lvl3_43E %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl3, -Quality_assignement) %>% 
  left_join(List_45E47F_Diagnostic_MOTU3 %>% add_column(Quality_assignement = "Assignement_lvl3"))
DB_assigned_lvl3_45E47F$Quality_assignement <- DB_assigned_lvl3_45E47F$Quality_assignement %>% replace(is.na(.), 0)

DB_assigned_lvl3_49E <- DB_assigned_lvl3_45E47F %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl3, -Quality_assignement) %>% 
  left_join(List_49E_Diagnostic_MOTU3 %>% add_column(Quality_assignement = "Assignement_lvl3"))
DB_assigned_lvl3_49E$Quality_assignement <- DB_assigned_lvl3_49E$Quality_assignement %>% replace(is.na(.), 0)

###### Merge all data assigned to either the Basetype or MOTUs lvl-3 #########
DB_assigned_lvl3 <- DB_assigned_Basetype %>% filter(Quality_assignement != 0) %>% 
  rbind(DB_assigned_lvl3_37F) %>%    filter(Quality_assignement != 0) %>% 
  rbind(DB_assigned_lvl3_41F) %>%    filter(Quality_assignement != 0) %>% 
  rbind(DB_assigned_lvl3_43E) %>%    filter(Quality_assignement != 0) %>%
  rbind(DB_assigned_lvl3_45E47F) %>% filter(Quality_assignement != 0) %>% 
  rbind(DB_assigned_lvl3_49E) %>%    filter(Quality_assignement != 0)

# Same as before but for the level 2
# Now we produce a file with only necessary informations
Key_assignement_MOTUS_lvl2 <- Key_assignement %>% 
  unite("Taxo_lvl2", Clade_curated,Genus_curated,species_curated, MOTU_lvl_1, MOTU_lvl_2, remove = FALSE,sep= "|") %>%
  select(Taxo_lvl2,  `37F_seq`:`49E_seq`)

# Identify individual sequence pattern that occur only in a single MOTU lvl2 and that is therefore diagnostic

List_37F_Diagnostic_MOTU2 <- Key_assignement_MOTUS_lvl2 %>% 
  select(Taxo_lvl2, `37F_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl2,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("37F_seq") %>% 
  select(`37F_seq`) %>% left_join(Key_assignement_MOTUS_lvl2 %>% select(Taxo_lvl2, `37F_seq`) %>% distinct(`37F_seq`, .keep_all =TRUE)) 

List_41F_Diagnostic_MOTU2 <- Key_assignement_MOTUS_lvl2 %>% 
  select(Taxo_lvl2, `41F_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl2,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("41F_seq") %>% 
  select(`41F_seq`) %>% left_join(Key_assignement_MOTUS_lvl2 %>% select(Taxo_lvl2, `41F_seq`) %>% distinct(`41F_seq`, .keep_all =TRUE)) 

List_43E_Diagnostic_MOTU2 <- Key_assignement_MOTUS_lvl2 %>% 
  select(Taxo_lvl2, `43E_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl2,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("43E_seq") %>% 
  select(`43E_seq`) %>% left_join(Key_assignement_MOTUS_lvl2 %>% select(Taxo_lvl2, `43E_seq`) %>% distinct(`43E_seq`, .keep_all =TRUE)) 

List_45E47F_Diagnostic_MOTU2 <- Key_assignement_MOTUS_lvl2 %>% 
  select(Taxo_lvl2, `45E-47F_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl2,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("45E-47F_seq") %>% 
  select(`45E-47F_seq`) %>% left_join(Key_assignement_MOTUS_lvl2 %>% select(Taxo_lvl2, `45E-47F_seq`) %>% distinct(`45E-47F_seq`, .keep_all =TRUE)) 

List_49E_Diagnostic_MOTU2 <- Key_assignement_MOTUS_lvl2 %>% 
  select(Taxo_lvl2, `49E_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl2,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("49E_seq") %>% 
  select(`49E_seq`) %>% left_join(Key_assignement_MOTUS_lvl2 %>% select(Taxo_lvl2, `49E_seq`) %>% distinct(`49E_seq`, .keep_all =TRUE))

##################### First Basetype level #################

DB_assigned_Basetype <-  Data_sanger %>%
  unite("Basetype", `37F_seq`:`49E_seq`, remove = FALSE,sep= "") %>%
  left_join(Key_assignement_MOTUS_lvl3 %>% unite("Basetype", `37F_seq`:`49E_seq`, remove = TRUE,sep= "") %>% 
              select(Basetype,Taxo_lvl3) %>% 
              add_column(Quality_assignement = "Basetype")) 
DB_assigned_Basetype$Quality_assignement <- DB_assigned_Basetype$Quality_assignement %>% replace(is.na(.), 0)

################### Assignement lvl 2 for each diagnostic sequence fragment ############
# As a starting point take the last file used for the assignement at lvl 3 to extract unassigned sequences

DB_assigned_lvl2_37F <- DB_assigned_lvl3_49E %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl3, -Quality_assignement) %>% 
  left_join(List_37F_Diagnostic_MOTU2 %>% add_column(Quality_assignement = "Assignement_lvl2")) 
DB_assigned_lvl2_37F$Quality_assignement <- DB_assigned_lvl2_37F$Quality_assignement %>% replace(is.na(.), 0)

DB_assigned_lvl2_41F <- DB_assigned_lvl2_37F %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl2, -Quality_assignement) %>% 
  left_join(List_41F_Diagnostic_MOTU2 %>% add_column(Quality_assignement = "Assignement_lvl2"))
DB_assigned_lvl2_41F$Quality_assignement <- DB_assigned_lvl2_41F$Quality_assignement %>% replace(is.na(.), 0)

DB_assigned_lvl2_43E <- DB_assigned_lvl2_41F %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl2, -Quality_assignement) %>% 
  left_join(List_43E_Diagnostic_MOTU2 %>% add_column(Quality_assignement = "Assignement_lvl2"))
DB_assigned_lvl2_43E$Quality_assignement <- DB_assigned_lvl2_43E$Quality_assignement %>% replace(is.na(.), 0)

DB_assigned_lvl2_45E47F <- DB_assigned_lvl2_43E %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl2, -Quality_assignement) %>% 
  left_join(List_45E47F_Diagnostic_MOTU2 %>% add_column(Quality_assignement = "Assignement_lvl2"))
DB_assigned_lvl2_45E47F$Quality_assignement <- DB_assigned_lvl2_45E47F$Quality_assignement %>% replace(is.na(.), 0)

DB_assigned_lvl2_49E <- DB_assigned_lvl2_45E47F %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl2, -Quality_assignement) %>% 
  left_join(List_49E_Diagnostic_MOTU2 %>% add_column(Quality_assignement = "Assignement_lvl2"))
DB_assigned_lvl2_49E$Quality_assignement <- DB_assigned_lvl2_49E$Quality_assignement %>% replace(is.na(.), 0)

# Now merge all data at taxonomic level 2
DB_assigned_lvl2 <-  rbind(DB_assigned_lvl2_37F) %>%    filter(Quality_assignement != 0) %>% 
  rbind(DB_assigned_lvl2_41F) %>%    filter(Quality_assignement != 0) %>% 
  rbind(DB_assigned_lvl2_43E) %>%    filter(Quality_assignement != 0) %>%
  rbind(DB_assigned_lvl2_45E47F) %>% filter(Quality_assignement != 0) %>% 
  rbind(DB_assigned_lvl2_49E) %>%    filter(Quality_assignement != 0)
  
# Final round to assign all possible data to level 1
  
# Now we produce a file with only necessary informations
Key_assignement_MOTUS_lvl1 <- Key_assignement %>% 
  unite("Taxo_lvl1", Clade_curated,Genus_curated,species_curated, MOTU_lvl_1, remove = FALSE,sep= "|") %>%
  select(Taxo_lvl1,  `37F_seq`:`49E_seq`)

# Identify individual sequence pattern that occur only in a single MOTU lvl1 and that is therefore diagnostic

List_37F_Diagnostic_MOTU1 <- Key_assignement_MOTUS_lvl1 %>% 
  select(Taxo_lvl1, `37F_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl1,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("37F_seq") %>% 
  select(`37F_seq`) %>% left_join(Key_assignement_MOTUS_lvl1 %>% select(Taxo_lvl1, `37F_seq`) %>% distinct(`37F_seq`, .keep_all =TRUE)) 

List_41F_Diagnostic_MOTU1 <- Key_assignement_MOTUS_lvl1 %>% 
  select(Taxo_lvl1, `41F_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl1,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("41F_seq") %>% 
  select(`41F_seq`) %>% left_join(Key_assignement_MOTUS_lvl1 %>% select(Taxo_lvl1, `41F_seq`) %>% distinct(`41F_seq`, .keep_all =TRUE)) 

List_43E_Diagnostic_MOTU1 <- Key_assignement_MOTUS_lvl1 %>% 
  select(Taxo_lvl1, `43E_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl1,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("43E_seq") %>% 
  select(`43E_seq`) %>% left_join(Key_assignement_MOTUS_lvl1 %>% select(Taxo_lvl1, `43E_seq`) %>% distinct(`43E_seq`, .keep_all =TRUE)) 

List_45E47F_Diagnostic_MOTU1 <- Key_assignement_MOTUS_lvl1 %>% 
  select(Taxo_lvl1, `45E-47F_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl1,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("45E-47F_seq") %>% 
  select(`45E-47F_seq`) %>% left_join(Key_assignement_MOTUS_lvl1 %>% select(Taxo_lvl1, `45E-47F_seq`) %>% distinct(`45E-47F_seq`, .keep_all =TRUE)) 

List_49E_Diagnostic_MOTU1 <- Key_assignement_MOTUS_lvl1 %>% 
  select(Taxo_lvl1, `49E_seq`) %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Taxo_lvl1,pres,fill=0) %>% remove_rownames %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 1) %>% 
  select(-sum) %>% 
  rownames_to_column("49E_seq") %>% 
  select(`49E_seq`) %>% left_join(Key_assignement_MOTUS_lvl1 %>% select(Taxo_lvl1, `49E_seq`) %>% distinct(`49E_seq`, .keep_all =TRUE))

################### Assignement lvl 1 for each diagnostic sequence fragment ############
# As a starting point take the last file used for the assignement at lvl 2 to extract unassigned sequences

DB_assigned_lvl1_37F <- DB_assigned_lvl2_49E %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl2, -Quality_assignement) %>% 
  left_join(List_37F_Diagnostic_MOTU1 %>% add_column(Quality_assignement = "Assignement_lvl1")) 
DB_assigned_lvl1_37F$Quality_assignement <- DB_assigned_lvl1_37F$Quality_assignement %>% replace(is.na(.), 0)

DB_assigned_lvl1_41F <- DB_assigned_lvl1_37F %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl1, -Quality_assignement) %>% 
  left_join(List_41F_Diagnostic_MOTU1 %>% add_column(Quality_assignement = "Assignement_lvl1"))
DB_assigned_lvl1_41F$Quality_assignement <- DB_assigned_lvl1_41F$Quality_assignement %>% replace(is.na(.), 0)

DB_assigned_lvl1_43E <- DB_assigned_lvl1_41F %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl1, -Quality_assignement) %>% 
  left_join(List_43E_Diagnostic_MOTU1 %>% add_column(Quality_assignement = "Assignement_lvl1"))
DB_assigned_lvl1_43E$Quality_assignement <- DB_assigned_lvl1_43E$Quality_assignement %>% replace(is.na(.), 0)

DB_assigned_lvl1_45E47F <- DB_assigned_lvl1_43E %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl1, -Quality_assignement) %>% 
  left_join(List_45E47F_Diagnostic_MOTU1 %>% add_column(Quality_assignement = "Assignement_lvl1")) 
DB_assigned_lvl1_45E47F$Quality_assignement <- DB_assigned_lvl1_45E47F$Quality_assignement %>% replace(is.na(.), 0)

DB_assigned_lvl1_49E <- DB_assigned_lvl1_45E47F %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl1, -Quality_assignement) %>% 
  left_join(List_49E_Diagnostic_MOTU1 %>% add_column(Quality_assignement = "Assignement_lvl1"))
DB_assigned_lvl1_49E$Quality_assignement <- DB_assigned_lvl1_49E$Quality_assignement %>% replace(is.na(.), 0)
# Now merge all data at taxonomic level 1

DB_assigned_lvl1 <-  rbind(DB_assigned_lvl1_37F) %>%    filter(Quality_assignement != 0) %>% 
  rbind(DB_assigned_lvl1_41F) %>%    filter(Quality_assignement != 0) %>% 
  rbind(DB_assigned_lvl1_43E) %>%    filter(Quality_assignement != 0) %>%
  rbind(DB_assigned_lvl1_45E47F) %>% filter(Quality_assignement != 0) %>% 
  rbind(DB_assigned_lvl1_49E) %>%    filter(Quality_assignement != 0)

# remaining sequences that cannot be assigned using perfect matches

DB_No_assignement <- DB_assigned_lvl1_49E %>% filter(Quality_assignement == 0) %>% 
  select(-Taxo_lvl1, -Quality_assignement) %>%  add_column(Quality_assignement = "No_Assignement")

# Merging all datasets
# The aim is to have for each taxonomic level a taxonomic path (TP) when possible. For example, if a sequence is assigned only to level 1, then there will be only 1 TP available but
# First homogenize the 3 files with the different assignment and create the same amount of columns for each of them.

DB_assigned_lvl3_bis <- DB_assigned_lvl3 %>% 
  separate(col = Taxo_lvl3, into = c("Clade_curated_assignement","Genus_curated_assignement","species_curated_assignement", "MOTU_lvl_1", "MOTU_lvl_2", "MOTU_lvl_3"), sep = "\\|") %>% 
  unite("TP_lvl1", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, remove = FALSE,sep= "|") %>% 
  unite("TP_lvl2", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, MOTU_lvl_2, remove = FALSE,sep= "|") %>% 
  unite("TP_lvl3", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, MOTU_lvl_2, MOTU_lvl_3, remove = FALSE,sep= "|")

DB_assigned_lvl2_bis <- DB_assigned_lvl2 %>% 
  separate(col = Taxo_lvl2, into = c("Clade_curated_assignement","Genus_curated_assignement","species_curated_assignement", "MOTU_lvl_1", "MOTU_lvl_2"), sep = "\\|") %>% 
  unite("TP_lvl1", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, remove = FALSE,sep= "|") %>% 
  unite("TP_lvl2", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, MOTU_lvl_2, remove = FALSE,sep= "|") %>% 
  add_column(TP_lvl3 = "Not_available", MOTU_lvl_3="Not_available")

DB_assigned_lvl1_bis <- DB_assigned_lvl1 %>% 
  separate(col = Taxo_lvl1, into = c("Clade_curated_assignement","Genus_curated_assignement","species_curated_assignement", "MOTU_lvl_1"), sep = "\\|") %>% 
  unite("TP_lvl1", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, remove = FALSE,sep= "|") %>% 
  add_column(TP_lvl3 = "Not_available", TP_lvl2 = "Not_available", MOTU_lvl_3="Not_available", MOTU_lvl_2="Not_available")

# Now merge all 3 files. The sequences without assignment are left out for now.

DB_all_assignement <- rbind(DB_assigned_lvl3_bis,DB_assigned_lvl2_bis, DB_assigned_lvl1_bis)

# Prepare files for Manual check. Making sure that sequences belonging to the same specimens are not attributes to more than 1 TP (Violation of the taxonomic rules)
# start with lvl-1

Check_co_occurence_lvl1 <- DB_all_assignement %>% select(Voucher_code, TP_lvl1) %>% filter(TP_lvl1 != "Not_available") %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(TP_lvl1,pres,fill=0) %>% 
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum != 1)

Check_co_occurence_lvl2 <- DB_all_assignement %>% select(Voucher_code, TP_lvl2) %>% filter(TP_lvl2 != "Not_available") %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(TP_lvl2,pres,fill=0) %>% 
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum != 1)

# and lvl-3
Check_co_occurence_lvl3 <- DB_all_assignement %>% select(Voucher_code, TP_lvl3) %>% filter(TP_lvl3 != "Not_available") %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(TP_lvl3,pres,fill=0) %>% 
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum != 1)

# Examination revealed that most issues are related to misatribution to neighbour Basegroups and more rarely genotypes. Other problems such as attribution to several morphospecies should be solved with manual curation.
# This is the case of obliquiloculata and pachyderma and due to the fact that there is not enough clones relative to the IG variability and there is oversplit at this level of the taxonomy.
# Best way is to substract the specimens affected and use another attribution method

List_wrong_attribution <- rbind((Check_co_occurence_lvl3 %>% rownames_to_column("Voucher_code") %>% select(Voucher_code)),  
                                (Check_co_occurence_lvl2 %>% rownames_to_column("Voucher_code") %>% select(Voucher_code)), 
                                (Check_co_occurence_lvl1 %>% rownames_to_column("Voucher_code") %>% select(Voucher_code))) %>% distinct()
  
List_to_be_re_assigned <- inner_join(DB_all_assignement, List_wrong_attribution)
DB_correct_matches <- anti_join(DB_all_assignement, List_wrong_attribution)

##################### Reassigned part of the dataset that could not be assigned using perfect matches #################

# Subset from the DB the sequences that needs to be assigned using the DADA2 package
# Format the Assignemnt key to be used by DADA2

DB_for_assignement <- Key_assignement %>% 
  unite("Name_seq", Clade_curated,Genus_curated,species_curated, MOTU_lvl_1, remove = FALSE,sep= ";") %>% 
  select(Name_seq,Basetype)

# Gather sequenes that need to be assigned or reassigned

DB_to_be_assigned <- rbind(List_to_be_re_assigned %>% select(Internal_ID, Basetype), DB_No_assignement %>% select(Internal_ID, Basetype))
Query_to_be_assigned <- DB_to_be_assigned %>% select(Basetype)

# Code for a function to write FASTA file
# Piece of code shamlessely taken from: https://bootstrappers.umassmed.edu/guides/main/r_writeFasta.html
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"Name_seq"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"Basetype"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# write fasta for DADA2
writeFasta(DB_for_assignement, "DB_assignement.fasta")
# Assign the sequences not qualified as Basetype using a bootstrap 100
tax <- assignTaxonomy(Query_to_be_assigned$Basetype,
                      refFasta = "DB_assignement.fasta",
                      taxLevels = c("Clade_curated", "Genus_curated", "species_curated", "MOTU_lvl_1"),
                      minBoot = 100)

# Count the number of sequences not attributed to each taxonomic level
tax <- as.data.frame(tax)
sum(is.na(tax$MOTU_lvl_3))
sum(is.na(tax$MOTU_lvl_2))
sum(is.na(tax$MOTU_lvl_1))

DB_assigned_DADA2 <- cbind((DB_to_be_assigned %>%  
                              select(Internal_ID)), 
                           (tax %>% rownames_to_column("seq") %>% select(-seq)))

# Now qualify the level of attibution of each sequence

DB_assigned_DADA2$MOTU_lvl_1 <- DB_assigned_DADA2$MOTU_lvl_1 %>% replace(is.na(.), "Not_available")
DB_assigned_DADA2_lvl1 <- DB_assigned_DADA2  %>% filter(MOTU_lvl_1 != "Not_available") %>% add_column(Quality_assignement = "Assignement_lvl1")
DB_not_assigned_DADA2  <- DB_assigned_DADA2 %>% filter(MOTU_lvl_1 == "Not_available") %>% add_column(Quality_assignement = "No_Assignement")


DB_assigned_DADA2_merged <-  DB_assigned_DADA2_lvl1 %>% inner_join(Data_sanger %>% select(-Clade_curated, -Genus_curated, -species_curated)) %>% 
  dplyr::rename( "Clade_curated_assignement" = "Clade_curated",	"Genus_curated_assignement" = Genus_curated,	"species_curated_assignement" = species_curated)

################################################################################################################


DB_assignement_Final <- rbind(DB_correct_matches       %>% select(Internal_ID, Voucher_code,`Latitude [Decimal]`,`Longitude  [Decimal]`, Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement,MOTU_lvl_1,MOTU_lvl_2,MOTU_lvl_3,Quality_assignement),
                              DB_assigned_DADA2_merged %>% select(Internal_ID, Voucher_code,`Latitude [Decimal]`,`Longitude  [Decimal]`, Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement,MOTU_lvl_1,Quality_assignement) %>% add_column(MOTU_lvl_2 ="Not_Attributed",MOTU_lvl_3="Not_Attributed"))

# And add Taxonomic path for plotting

DB_assigned_TPlvl3 <- DB_assignement_Final %>% filter(Quality_assignement != "Assignement_lvl1") %>% filter(Quality_assignement != "Assignement_lvl2") %>% filter(Quality_assignement != "Not_Attributed") %>%filter(Quality_assignement != "No_Assignement") %>%
  unite("TP_MS", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, remove = FALSE,sep= "|") %>%
  unite("TP_lvl1", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, remove = FALSE,sep= "|") %>% 
  unite("TP_lvl2", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, MOTU_lvl_2, remove = FALSE,sep= "|") %>% 
  unite("TP_lvl3", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, MOTU_lvl_2, MOTU_lvl_3, remove = FALSE,sep= "|")

DB_assigned_TPlvl2 <- DB_assignement_Final %>% filter(Quality_assignement == "Assignement_lvl2") %>% 
  unite("TP_MS", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, remove = FALSE,sep= "|") %>%
  unite("TP_lvl1", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, remove = FALSE,sep= "|") %>% 
  unite("TP_lvl2", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, MOTU_lvl_2, remove = FALSE,sep= "|") %>% 
  add_column(TP_lvl3 = "Not_Attributed")

DB_assigned_TPlvl1 <- DB_assignement_Final %>% filter(Quality_assignement == "Assignement_lvl1") %>% 
  unite("TP_MS", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, remove = FALSE,sep= "|") %>%
  unite("TP_lvl1", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, remove = FALSE,sep= "|") %>% 
  add_column(TP_lvl2 = "Not_Attributed") %>% 
  add_column(TP_lvl3 = "Not_Attributed")

# And finally merge again all parts of the DB
DB_all_assignement <- rbind(DB_assigned_TPlvl3, DB_assigned_TPlvl2, DB_assigned_TPlvl1) %>% filter(TP_lvl1 != "Not_Attributed")


# Prepare files for Manual check. Making sure that sequences belonging to the same specimens are not attributed to more than 1 TP (Violation of the taxonomic rules)
# start with lvl-1
Check_co_occurence_lvl1 <- DB_all_assignement %>% select(Voucher_code, TP_lvl1) %>% filter(TP_lvl1 != "Not_Attributed") %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(TP_lvl1,pres,fill=0) %>% 
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum != 1)

# then lvl-2
# write.table(Check_co_occurence_lvl1, "checking_lvl1.txt")
Check_co_occurence_lvl2 <- DB_all_assignement %>% select(Voucher_code, TP_lvl2) %>% filter(TP_lvl2 != "Not_Attributed") %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(TP_lvl2,pres,fill=0) %>% 
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum != 1)

write.table(Check_co_occurence_lvl2, "checking_lvl2.txt")
# and lvl-3
Check_co_occurence_lvl3 <- DB_all_assignement %>% select(Voucher_code, TP_lvl3) %>% filter(TP_lvl3 != "Not_Attributed") %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(TP_lvl3,pres,fill=0) %>% 
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum != 1)

# Remove from the Sanger dataset the specimens that could not be attributed reliably by any method
# First list all specimens with obvious misatribution
List_wrong_attribution_final <- rbind(Check_co_occurence_lvl1 %>% rownames_to_column(var = "Voucher_code") %>% select(Voucher_code),
                                      Check_co_occurence_lvl2 %>% rownames_to_column(var = "Voucher_code") %>% select(Voucher_code),
                                      Check_co_occurence_lvl3 %>% rownames_to_column(var = "Voucher_code") %>% select(Voucher_code))

#Then remove them from the the final Sanger dataset
DB_all_assignement <- anti_join(DB_all_assignement, List_wrong_attribution_final)

Sumarize_Sanger_attribution <- rbind(DB_all_assignement %>%  filter(TP_lvl3 != "Not_Attributed") %>%  select(TP_lvl3) %>% count() %>% add_column(N_path = "MOTUs lvl-3"),
                                     DB_all_assignement %>%  filter(TP_lvl3 != "Not_Attributed") %>%  filter(TP_lvl2 != "Not_Attributed")  %>%  select(TP_lvl2)  %>% count() %>% add_column(N_path = "MOTUs lvl-2"),
                                     DB_all_assignement %>%  filter(TP_lvl3 != "Not_Attributed") %>%  filter(TP_lvl2 != "Not_Attributed")  %>%  filter(TP_lvl1 != "Not_Attributed") %>%  select(TP_lvl1) %>% count() %>% add_column(N_path = "MOTUs lvl-1"))


Sumarize_Sanger_attribution <- DB_all_assignement %>% count(Quality_assignement)

   ########################################################################
   #                  Import and Format Genotyping data                   #
   ########################################################################

# Here we translate the RFLP data into their taxonomic equivalent into the novel molecular Nomenclature
RFLP_data_Formated <- RFLP_data %>%
  select('ID-RFLP','Depth collection', 9:32) %>% 
  gather("Genotyping_taxonomy", n_observation,3:26) %>% 
  left_join(Equivalence_Nomenclature) %>% 
  unite("TP_MS", Clade_curated,Genus_curated,species_curated, remove = FALSE,sep= "|") %>%
  unite("TP_lvl1", Clade_curated,Genus_curated,species_curated, MOTU_lvl_1, remove = FALSE,sep= "|") %>% 
  unite("TP_lvl2", Clade_curated,Genus_curated,species_curated, MOTU_lvl_1, MOTU_lvl_2, remove = FALSE,sep= "|") %>% 
  add_column(TP_lvl3 = "Not_Attributed") %>% 
  filter(n_observation != 0)


    ##############################################################################
    #               Integrating V9 data from Metabarcoding dataset               #
    ##############################################################################

# Extract the number of observation for each PF ASV
N_Obs_ASV <-  ASV_Occurence %>% select(1:2084) %>%  gather("Samples", "reads", 2:2084) %>% filter(reads != "0") %>% select(ASV_ID) %>% count(ASV_ID) %>% dplyr::rename("N_obs_samples" = "n")

# Import ASV segmented and keep only ASVs that occured more than once in the entire dataset
ASV_data <- PF_ASV_seq %>% left_join(N_Obs_ASV) %>% filter(N_obs_samples >"1")
  
# Used same strategy as with the sanger to match the exact pattern first.
ASV_MOTU_lvl3 <- left_join(List_49E_Diagnostic_MOTU3,ASV_data %>% select(ASV_ID, `49E_seq`)) %>% filter(ASV_ID !="NA")    
ASV_MOTU_lvl2 <- left_join(List_49E_Diagnostic_MOTU2,ASV_data %>% select(ASV_ID, `49E_seq`)  %>% anti_join(ASV_MOTU_lvl3)) %>% filter(ASV_ID !="NA")
ASV_MOTU_lvl1 <- left_join(List_49E_Diagnostic_MOTU1,ASV_data %>% select(ASV_ID, `49E_seq`)  %>% anti_join(ASV_MOTU_lvl3)  %>% anti_join(ASV_MOTU_lvl2)) %>% filter(ASV_ID !="NA")

# Isolate ASVs that do not have a perfect match in the nomenclature
ASV_data_unknown <- ASV_data %>% anti_join(ASV_MOTU_lvl3) %>% anti_join(ASV_MOTU_lvl2)%>% anti_join(ASV_MOTU_lvl1)

# Here we prepare the file to compare the ASVs to MOTUs
# Recruit fragment of the BASETYPE covering the v9
Selection_holotypes_V9 <- Key_assignement %>% 
  unite("Taxo_lvl3",Basetype_ID, Clade_curated,Genus_curated,species_curated, MOTU_lvl_1, MOTU_lvl_2, MOTU_lvl_3, remove = FALSE,sep= "|") %>%
  add_column(additional_name = "BASETYPE") %>% 
  unite("Name_seq", additional_name,Taxo_lvl3, remove = FALSE,sep= "_") %>%
  select(Name_seq,  `47-49_seq`:`49E_seq`) %>% 
  unite("Sequence", `47-49_seq`:`49E_seq`, remove = TRUE,sep= "")

# Recruit the ASVs that need to be compared to the MOTUs 
ASV_unknown_sequences <- ASV_data_unknown %>% 
  add_column(additional_name = "ENVSEQ") %>%
  unite("Name_seq", additional_name,ASV_ID, remove = FALSE,sep= "_") %>%
  unite("Sequence", `47-49_seq`:`49E_seq`, remove = TRUE,sep= "") %>% 
  select(Name_seq, Sequence)


# Create a function to produce a fasta file from a 2 column dataframe
# Piece of code shamlessely taken from: https://bootstrappers.umassmed.edu/guides/main/r_writeFasta.html
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"Name_seq"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"Sequence"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# Concatanate infos to write fasta
# writeFasta(rbind(Selection_holotypes_V9,ASV_unknown_sequences), "Unknown_ASV_w_Basetype.fasta")

# See detailed methods for the explanation of how the ASVs not identical to the reference are analysed
# The resulting Taxonomy are provided in Table S4 and imported as Molecular_Nomenclature_Novel_ASVs object intot the working environment

# Now merge all ASVs taxonomy into a single dataset
# First group attriution of ASVs
ASV_assigned_TPlvl3 <- ASV_MOTU_lvl3 %>% select(-`49E_seq`) %>% 
  separate(col ="Taxo_lvl3", into =c("Clade_curated_assignement", "Genus_curated_assignement", "species_curated_assignement", "MOTU_lvl_1", "MOTU_lvl_2", "MOTU_lvl_3"), sep = "\\|") %>% 
  unite("TP_MS", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, remove = FALSE,sep= "|") %>%
  unite("TP_lvl1", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, remove = FALSE,sep= "|") %>% 
  unite("TP_lvl2", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, MOTU_lvl_2, remove = FALSE,sep= "|") %>% 
  unite("TP_lvl3", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, MOTU_lvl_2, MOTU_lvl_3, remove = FALSE,sep= "|")

ASV_assigned_TPlvl2 <- ASV_MOTU_lvl2 %>% select(-`49E_seq`) %>% 
  rbind(Molecular_Nomenclature_Novel_ASVs %>% filter(Level_attribution == 'MOTUlvl2') %>% select(ASV_ID,'Taxonomy attribution') %>% dplyr::rename("Taxo_lvl2" = "Taxonomy attribution")) %>% 
  separate(col ="Taxo_lvl2", into =c("Clade_curated_assignement", "Genus_curated_assignement", "species_curated_assignement", "MOTU_lvl_1", "MOTU_lvl_2", "MOTU_lvl_3"), sep = "\\|") %>% 
  unite("TP_MS", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, remove = FALSE,sep= "|") %>%
  unite("TP_lvl1", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, remove = FALSE,sep= "|") %>% 
  unite("TP_lvl2", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, MOTU_lvl_2, remove = FALSE,sep= "|") %>% 
  add_column(TP_lvl3 = "Not_Attributed")

ASV_assigned_TPlvl1 <- ASV_MOTU_lvl1 %>% select(-`49E_seq`) %>% 
  rbind(Molecular_Nomenclature_Novel_ASVs %>% filter(Level_attribution == 'MOTUlvl1') %>% select(ASV_ID,'Taxonomy attribution') %>% dplyr::rename("Taxo_lvl1" = "Taxonomy attribution")) %>% 
  separate(col ="Taxo_lvl1", into =c("Clade_curated_assignement", "Genus_curated_assignement", "species_curated_assignement", "MOTU_lvl_1", "MOTU_lvl_2", "MOTU_lvl_3"), sep = "\\|") %>% 
  unite("TP_MS", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, remove = FALSE,sep= "|") %>%
  unite("TP_lvl1", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, remove = FALSE,sep= "|") %>% 
  add_column(TP_lvl3 = "Not_Attributed")%>%  add_column(TP_lvl2 = "Not_Attributed")

ASV_assigned_MS <-  Molecular_Nomenclature_Novel_ASVs %>% filter(Level_attribution == 'MS') %>% select(ASV_ID,'Taxonomy attribution') %>% dplyr::rename("Taxo_MS" = "Taxonomy attribution") %>% 
  separate(col ="Taxo_MS", into =c("Clade_curated_assignement", "Genus_curated_assignement", "species_curated_assignement", "MOTU_lvl_1", "MOTU_lvl_2", "MOTU_lvl_3"), sep = "\\|") %>% 
  unite("TP_MS", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, remove = FALSE,sep= "|") %>%
  add_column(TP_lvl3 = "Not_Attributed") %>%  add_column(TP_lvl2 = "Not_Attributed") %>%  add_column(TP_lvl1 = "Not_Attributed") 

# Merge the Taxonomy of all ASV into a single file
ASV_all_Assignement <- rbind(ASV_assigned_TPlvl3, ASV_assigned_TPlvl2, ASV_assigned_TPlvl1, ASV_assigned_MS) %>% unique()

# First extract the Planktonic foraminifera ASVs to reduce the size of the table
PF_ASV_Occurence <- ASV_Occurence %>% inner_join(ASV_all_Assignement) %>% select(1:2084)
# Merge table with geographic coordinates
ASV_Biogeography <- PF_ASV_Occurence %>% 
  t() %>% 
  row_to_names(row_number = 1) %>% as.data.frame() %>% 
  rownames_to_column(var = "Original_Sample_ID")

ASV_Matrix <- ASV_Biogeography %>% select(-"Original_Sample_ID") %>% mutate_if(is.character, as.numeric)

ASV_Biogeography_final <- cbind((ASV_Biogeography %>% select(Original_Sample_ID)), ASV_Matrix) %>% 
                           group_by(Original_Sample_ID) %>%
                           summarise_each(funs(sum)) %>% 
                           drop_na()

# Now do a long list and merge the molecular taxonomy
ASV_Occurence_taxo <- ASV_Biogeography_final %>% gather("ASV_ID", "reads", 2:143) %>% 
  filter(reads != 0) %>% 
  full_join(ASV_all_Assignement %>% select(ASV_ID, TP_MS, TP_lvl1, TP_lvl2, TP_lvl3)) %>% 
  add_column(Dataset = "METABARCODING") %>% 
  dplyr::rename("n_observation" = "reads") %>% select(-ASV_ID)

# Now the Sanger dataset
Sanger_Occurence_taxo <- DB_all_assignement %>% 
  select(Internal_ID, TP_MS, TP_lvl1, TP_lvl2, TP_lvl3) %>%
  dplyr::rename("Original_Sample_ID" = "Internal_ID") %>%
  add_column(n_observation = 1) %>% 
  add_column(Dataset = "SANGER")

# And RFLP
RFLP_Occurence_taxo <- RFLP_data_Formated %>% 
  select('ID-RFLP', TP_MS, TP_lvl1, TP_lvl2, TP_lvl3, n_observation) %>%
  dplyr::rename("Original_Sample_ID"  = 'ID-RFLP') %>%
  add_column(Dataset = "GENOTYPING") %>% 
  filter(n_observation != 0 )

# Merge data
All_occurences_final <- rbind(RFLP_Occurence_taxo, 
                              Sanger_Occurence_taxo, 
                              ASV_Occurence_taxo) %>% distinct() %>%
                              dplyr::rename("Original Sample ID"  = "Original_Sample_ID")

# Export occurences to make them available without having to run the scipt
write.xlsx(
  All_occurences_final_geography %>%  
    as.data.frame(),
  file="../Table_S07bis.xlsx",
  sheetName = "Sheet1",
  col.names = TRUE,
  row.names = FALSE
)





# Merge with geographical metadata
All_occurences_final_geography <- All_occurences_final  %>% 
  select(-Dataset) %>% 
  left_join(Geography_Curated) %>%  
  dplyr::rename("latitude" = "Averaged_Latitude", "longitude" ="Averaged_Longitude") %>% distinct()



#################################################################################################
#                                                                                               #
#                                Plot MAPs For all MOTUs                                        #
#                                                                                               #
#################################################################################################


# Produce one dataset per taxonomic path 
Occurences_MS     <- All_occurences_final_geography %>% select(latitude, longitude, TP_MS)   %>% filter(TP_MS   != "Not_Attributed")   %>% unique()
Occurences_level1 <- All_occurences_final_geography %>% select(latitude, longitude, TP_lvl1) %>% filter(TP_lvl1 != "Not_Attributed") %>% unique()
Occurences_level2 <- All_occurences_final_geography %>% select(latitude, longitude, TP_lvl2) %>% filter(TP_lvl2 != "Not_Attributed") %>% unique()
Occurences_level3 <- All_occurences_final_geography %>% select(latitude, longitude, TP_lvl3) %>% filter(TP_lvl3 != "Not_Attributed") %>% unique()

# create list of tables based on the TP column
level1.list <- split(Occurences_level1, list(Occurences_level1$TP_lvl1), drop=T)
level2.list <- split(Occurences_level2, list(Occurences_level2$TP_lvl2), drop=T)
level3.list <- split(Occurences_level3, list(Occurences_level3$TP_lvl3), drop=T)
levelMS.list <- split(Occurences_MS, list(Occurences_MS$TP_MS), drop=T)
# Plot on a world map
library(maps)
library(parallel)
library(hexbin)
library(scales)
map_occurrences <- function(X,level){
  mapWorld <- borders("world", colour="#0a0c2a", fill="#0a0c2a") # create a layer of borders
  df <- as.tibble(X) %>% drop_na()
  df_coord <- df[, c("latitude", "longitude")] %>%
    mutate(latitude=as.numeric(latitude),
           longitude=as.numeric(longitude))
  
  fake_coord <- tibble(
    latitude = c(seq(-100, 100, 1), rep(200,401)),
    longitude = c(rep(400,201),seq(-200, 200, 1))
  )
  
  df_plot <- rbind(df_coord,fake_coord)
  
  name <- gsub("\\|","_",unique(df[,3]))
  
  zp1 <- ggplot(df_plot) + mapWorld
  zp1 <- zp1 + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.ticks = element_blank(), 
                     axis.text.x = element_text (size = 12, vjust = 0), 
                     axis.text.y = element_text (size = 12, hjust = 1.3)) + 
    coord_cartesian() + labs(y="",x="") 
  
  zp1 <- zp1 + 
    geom_hex(data = df_plot,
             aes(x = longitude, y = latitude,fill=..count..),
             alpha=0.8, col="black", bins=60) +
    stat_binhex(data = df_plot,
                aes(x = longitude, y = latitude,label=..count..), geom="text",
                colour="black", size=2, bins=60) +
    coord_cartesian(xlim=c(-180,180),ylim=c(-90,90))
  
  zp1 <- zp1 +
    scale_fill_distiller("# Occurrences",
                         palette = "YlOrRd",
                         direction = 1,
                         breaks = scales::pretty_breaks()) +
    theme_bw() + #coord_equal(ratio=1) + 
    theme(legend.position="bottom",
          legend.title = element_text(size=12,vjust = 0.8),
          legend.text = element_text(size=10),
          legend.spacing.x = unit(0.2,"cm")) +
    ggtitle(paste0(level,"_",name)) +
    xlab("") + ylab("")   
  ggsave(plot = zp1, paste0(level,"_",name,".jpg"), width = 8.8, height = 6.5,device = "jpg")
}

# The maps will be printed in the same folder as the scripts and the excel tables
# to change the folder add folder name in ggsave:
# ggsave(plot = zp1, paste0("new_folder/",level,"_",df$TP_lvl3,".pdf"), width = 8, height = 5,device = "pdf")
# To save in a different format then .pdf, for example .jpg:
# ggsave(plot = zp1, paste0(level,"_",df$TP_lvl1,".jpeg"), width = 8, height = 5,device = "jpeg")

maps_level1 <- lapply(level1.list, map_occurrences, level="Level1")
maps_level2 <- lapply(level2.list, map_occurrences, level="Level2")
maps_level3 <- lapply(level3.list, map_occurrences, level="Level3")
maps_levelMS <- lapply(levelMS.list, map_occurrences, level="LevelMS")







######################################### END of official script #################################################

# data for new figure
# Compile diversity vs number of observation
# At lvl-1

N_TP_lvl1 <- All_occurences_final_geography %>% 
  filter(Dataset != "Genotyping") %>% 
  select(TP_MS, TP_lvl1) %>% 
  filter(TP_lvl1 != "Not_Attributed") %>% 
  distinct() %>% 
  count(TP_MS) %>% dplyr::rename("n_TP"="n")%>% 
  add_column(Dataset = "TP_lvl1")

N_TP_lvl2 <- All_occurences_final_geography %>% 
  filter(Dataset != "Genotyping") %>% 
  select(TP_MS, TP_lvl2) %>%
  filter(TP_lvl2 != "Not_Attributed") %>% 
  distinct() %>% 
  count(TP_MS) %>% dplyr::rename("n_TP"="n")%>% 
  add_column(Dataset = "TP_lvl2")

N_TP_lvl3 <- All_occurences_final_geography %>% 
  filter(Dataset != "Genotyping") %>% 
  select(TP_MS, TP_lvl3) %>% 
  filter(TP_lvl3 != "Not_Attributed") %>% 
  distinct() %>% 
  count(TP_MS) %>% dplyr::rename("n_TP"="n")%>% 
  add_column(Dataset = "TP_lvl3")


N_Obs_lvl1 <- All_occurences_final_geography %>% 
  filter(Dataset != "Genotyping") %>% 
  select(TP_MS, TP_lvl1, Curated_Geography_ID) %>% 
  filter(TP_lvl1 != "Not_Attributed") %>%
  distinct() %>% 
  count(TP_MS) %>% dplyr::rename("n_Obs"="n")%>% 
  add_column(Dataset = "TP_lvl1")

N_Obs_lvl2 <- All_occurences_final_geography %>% 
  filter(Dataset != "Genotyping") %>% 
  select(TP_MS, TP_lvl2, Curated_Geography_ID) %>% 
  filter(TP_lvl2 != "Not_Attributed") %>%
  distinct() %>% 
  count(TP_MS) %>% dplyr::rename("n_Obs"="n") %>% 
  add_column(Dataset = "TP_lvl2")

N_Obs_lvl3 <- All_occurences_final_geography %>% 
  filter(Dataset != "Genotyping") %>% 
  select(TP_MS, TP_lvl3, Curated_Geography_ID) %>% 
  filter(TP_lvl3 != "Not_Attributed") %>%
  distinct() %>% 
  count(TP_MS) %>% dplyr::rename("n_Obs"="n") %>% 
  add_column(Dataset = "TP_lvl3")


N_seq_lvl1 <- All_occurences_final_geography %>% 
  filter(Dataset != "Genotyping") %>% 
  select(TP_MS, TP_lvl1, Curated_Geography_ID, n_observation) %>% 
  filter(TP_lvl1 != "Not_Attributed") %>%
  distinct() %>%
  group_by(TP_MS) %>% 
  summarise(N_seq = sum(n_observation)) %>% 
  ungroup() %>% distinct() %>%  
  add_column(Dataset = "TP_lvl1")

N_seq_lvl2 <- All_occurences_final_geography %>% 
  filter(Dataset != "Genotyping") %>% 
  select(TP_MS, TP_lvl2, Curated_Geography_ID, n_observation) %>% 
  filter(TP_lvl2 != "Not_Attributed") %>%
  distinct() %>%
  group_by(TP_MS) %>% 
  summarise(N_seq = sum(n_observation)) %>% 
  ungroup() %>% distinct() %>%  
  add_column(Dataset = "TP_lvl2")

N_seq_lvl3 <- All_occurences_final_geography %>% 
  filter(Dataset == "Sanger") %>% 
  select(TP_MS, TP_lvl3, Curated_Geography_ID, n_observation) %>% 
  filter(TP_lvl3 != "Not_Attributed") %>%
  distinct() %>%
  group_by(TP_MS) %>% 
  summarise(N_seq = sum(n_observation)) %>% 
  ungroup() %>% distinct() %>%  
  add_column(Dataset = "TP_lvl3")


  

Data_div_vs_sample <- rbind(N_TP_lvl1,
                            N_TP_lvl2,
                            N_TP_lvl3) %>% inner_join(
                      rbind(N_Obs_lvl1,
                            N_Obs_lvl2,
                            N_Obs_lvl3))%>% inner_join(
                      rbind(N_seq_lvl1,
                            N_seq_lvl2,
                            N_seq_lvl3))


library(ggpmisc)
library(ggpubr)

ggplot(Data_div_vs_sample %>% filter(n_Obs> 0), aes(y=n_Obs, x=n_TP))+
  geom_point()+
  facet_grid(.~Dataset, scales = "free")+
  geom_smooth(method='lm')+
  stat_regline_equation(label.x = 3, label.y = 2)+
  stat_cor(label.x = 3, label.y = 1)+
  scale_y_log10()




  ########################################################################
#                  Summarizing observation                             #
########################################################################

# Produce a file to have an overview of the amount of data available for each morphological species

# Number of specimens per Morphological species.
Summary_observation_specimens <- Data_sanger %>% 
  unite("TP_MS", Clade_curated,Genus_curated,species_curated, remove = FALSE,sep= "|") %>% 
  select(Voucher_code, TP_MS) %>%
  distinct() %>% 
  group_by(TP_MS) %>% 
  summarise(n_specimens_sequenced = n()) 

# Number of sequences per Morphological species.
Summary_observation_sequences <- Data_sanger %>% 
  unite("TP_MS", Clade_curated,Genus_curated,species_curated, remove = FALSE,sep= "|") %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_sequences = n())

# Number of Genotyping per Morphological species.
Summary_observation_Genotyping <- RFLP_data_Formated %>% 
  select(TP_MS, n_observation) %>%
  group_by(TP_MS) %>% 
  summarise(n_specimens_genotyped = sum(n_observation))

# Number of sequences having quality to be qualified as Basetype.
Summary_observation_basetypes <- Data_sanger %>% filter(Quality=="BASETYPE") %>% 
  unite("TP_MS", Clade_curated,Genus_curated,species_curated, remove = FALSE,sep= "|") %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_basetype = n())

# Number of sequences chosen to be representative of sequences of genotypes or MS presented in the literature
Summary_observation_representative <- Data_sanger %>% filter(Quality=="REPRESENTATIVE") %>% 
  unite("TP_MS", Clade_curated,Genus_curated,species_curated, remove = FALSE,sep= "|") %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_representative = n())

# Number of sequences attributed to the level 1  for Barcoding data
Summary_observation_TPlvl1_SC <- DB_all_assignement %>% select(TP_MS, TP_lvl1) %>%
  filter(TP_lvl1 != "Not_Attributed") %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_obs_TPlvl1 = n()) %>% 
  left_join(RFLP_data_Formated %>% select(TP_MS,n_observation) %>% group_by(TP_MS) %>% summarise(n_observation = sum(n_observation))) 

Summary_observation_TPlvl1[is.na(Summary_observation_TPlvl1)] <- 0
Summary_observation_TPlvl1 <- Summary_observation_TPlvl1 %>%  
  mutate(n_obs_TPlvl1_SC = n_obs_TPlvl1+n_obs) %>% 
  select(-n_obs)


# Number of sequences attributed to the level 2  for Barcoding data
Summary_observation_TPlvl2_SC <- DB_all_assignement %>% select(TP_MS, TP_lvl2) %>%
  filter(TP_lvl2 != "Not_Attributed") %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_obs_TPlvl2 = n()) %>% 
  left_join(RFLP_data_Formated %>% select(TP_MS,n_observation) %>% group_by(TP_MS) %>% summarise(n_observation = sum(n_observation)))

Summary_observation_TPlvl2[is.na(Summary_observation_TPlvl2)] <- 0
Summary_observation_TPlvl2 <- Summary_observation_TPlvl2 %>% 
  mutate(n_obs_TPlvl2_SC = n_obs_TPlvl2+n_obs) %>% 
  select(-n_obs)

# Number of sequences attributed to the level 3  for Barcoding data
Summary_observation_TPlvl3_SC <- DB_all_assignement %>% select(TP_MS, TP_lvl3) %>%
  filter(TP_lvl3 != "Not_Attributed") %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_obs_TPlvl3_SC = n())

# Number of unique TP at level 1  for Barcoding data
Summary_observation_N_TPlvl1_SC <- DB_all_assignement %>% select(TP_MS, TP_lvl1) %>%
  filter(TP_lvl1 != "Not_Attributed") %>% 
  distinct() %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_Taxonomic_Path_lvl1_SC = n())

# Number of unique TP at level 2 for Barcoding data
Summary_observation_N_TPlvl2_SC <- DB_all_assignement %>% select(TP_MS, TP_lvl2) %>%
  filter(TP_lvl2 != "Not_Attributed") %>% 
  distinct() %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_Taxonomic_Path_lvl2_SC = n())

# Number of unique TP at level 3 for Barcoding data
Summary_observation_N_TPlvl3_SC <- DB_all_assignement %>% select(TP_MS, TP_lvl3) %>%
  filter(TP_lvl3 != "Not_Attributed") %>% 
  distinct() %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_Taxonomic_Path_lvl3_SC = n())

# Number of ASVs
Summary_observation_N_ASVs <- ASV_all_Assignement %>% select(TP_MS, ASV_ID) %>%
  distinct() %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_ASVs = n())
Summary_observation_N_ASVs
# Number of reads

Summary_observation_N_reads <- PF_ASV_Occurence %>% column_to_rownames(1) %>% rowSums() %>% as.data.frame() %>% 
  dplyr::rename("n_reads"=".") %>% 
  rownames_to_column(var = "ASV_ID") %>% 
  left_join(ASV_all_Assignement %>% select(TP_MS, ASV_ID)) %>% select(-ASV_ID) %>% 
  group_by(TP_MS) %>%
  summarise_each(funs(sum)) 

# Total_Number of Taxonomic path
Summary_observation_N_TPlvl1_SC <- Key_assignement %>% 
  unite("TP_MS", Clade_curated,Genus_curated,species_curated, remove = FALSE,sep= "|") %>%
  unite("TP_lvl1", Clade_curated,Genus_curated,species_curated, MOTU_lvl_1, remove = FALSE,sep= "|") %>% 
  select(TP_MS, TP_lvl1) %>%
  filter(TP_lvl1 != "Not_Attributed") %>% 
  distinct() %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_Taxonomic_Path_lvl1_SC = n())

Summary_observation_N_TPlvl2_SC <- Key_assignement %>% 
  unite("TP_MS", Clade_curated,Genus_curated,species_curated, remove = FALSE,sep= "|") %>%
  unite("TP_lvl2", Clade_curated,Genus_curated,species_curated, MOTU_lvl_1,MOTU_lvl_2, remove = FALSE,sep= "|") %>% 
  select(TP_MS, TP_lvl2) %>%
  filter(TP_lvl2 != "Not_Attributed") %>% 
  distinct() %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_Taxonomic_Path_lvl2_SC = n())

Summary_observation_N_TPlvl3_SC <- Key_assignement %>% 
  unite("TP_MS", Clade_curated,Genus_curated,species_curated, remove = FALSE,sep= "|") %>%
  unite("TP_lvl3", Clade_curated,Genus_curated,species_curated, MOTU_lvl_1,MOTU_lvl_2,MOTU_lvl_3, remove = FALSE,sep= "|") %>% 
  select(TP_MS, TP_lvl3) %>%
  filter(TP_lvl3 != "Not_Attributed") %>% 
  distinct() %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_Taxonomic_Path_lvl3_SC = n())

Summary_observation_N_TPlvl1_ALL <- All_occurences_final %>% select(TP_MS, TP_lvl1) %>%
  filter(TP_lvl1 != "Not_Attributed") %>% 
  distinct() %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_Taxonomic_Path_lvl1_ALL = n())

Summary_observation_N_TPlvl2_ALL <- All_occurences_final %>% select(TP_MS, TP_lvl2) %>%
  filter(TP_lvl2 != "Not_Attributed") %>% 
  distinct() %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_Taxonomic_Path_lvl2_ALL = n())

Summary_observation_N_TPlvl3_ALL <- All_occurences_final %>% select(TP_MS, TP_lvl3) %>%
  filter(TP_lvl3 != "Not_Attributed") %>% 
  distinct() %>% 
  select(TP_MS) %>%
  group_by(TP_MS) %>% 
  summarise(n_Taxonomic_Path_lvl3_ALL = n())

# Merging all data

Summary_observation_all <- left_join(Summary_observation_specimens,
                                     Summary_observation_sequences) %>%
  left_join(Summary_observation_Genotyping)     %>%
  full_join(Summary_observation_N_reads)        %>%
  full_join(Summary_observation_N_ASVs)         %>%
  left_join(Summary_observation_basetypes)      %>% 
  left_join(Summary_observation_representative) %>% 
  left_join(Summary_observation_TPlvl1_SC)      %>% 
  left_join(Summary_observation_TPlvl2_SC)      %>% 
  left_join(Summary_observation_TPlvl3_SC)      %>% 
  left_join(Summary_observation_N_TPlvl1_SC)    %>% 
  left_join(Summary_observation_N_TPlvl2_SC)    %>% 
  left_join(Summary_observation_N_TPlvl3_SC)    %>% 
  left_join(Summary_observation_N_TPlvl1_ALL)   %>% 
  left_join(Summary_observation_N_TPlvl2_ALL)   %>% 
  left_join(Summary_observation_N_TPlvl3_ALL)  




View(Summary_observation_all)

write.table(Summary_observation_all, "Summary_observation_all.txt")


library(ggplot2)
library(hrbrthemes)
library(viridis)
library(scales)

Summary_obs_plot_Formated <- Summary_observation_all %>% 
  select(TP_MS, n_reads, n_sequences, n_Taxonomic_Path_lvl1_ALL, n_Taxonomic_Path_lvl2_ALL, n_Taxonomic_Path_lvl3_ALL) %>% 
  gather(key=MOTU_level, value = N_Taxonomic_Path, 4:6) %>% replace(is.na(.), 0)

ggplot(Summary_obs_plot_Formated, aes( x=n_sequences, y=n_reads, size=N_Taxonomic_Path,color=N_Taxonomic_Path))+
  geom_point(alpha =0.5) + 
  scale_size(range = c(1, 20), name="N_taxonomic_path")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1,200000)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1,200000)) +
  theme_bw()+ annotation_logticks() +coord_fixed(ratio = 1)+ facet_wrap(.~MOTU_level)

ggplot(Summary_observation_all, aes( x=n_ASVs, y=n_Taxonomic_Path_lvl3_ALL))+
  geom_point(alpha =0.5)

################################################################################
#                                                                              #
#                    Prepare elements for diversity figure                     #
#                                                                              #
################################################################################


# First rarefaction curves
library(iNEXT)
# We need to calculate one curve per taxonomic level because we do not have the same amount of observation at each station for all taxonomic level
# Extract information sequencially, we use the number of occrences at unique locality to calculate the rarefaction, and remove genotyping data
# Morphological species
data_rarefaction_MS <- All_occurences_final_geography %>% 
  filter(Dataset != "Genotyping") %>% 
  select(TP_MS, Curated_Geography_ID) %>% 
  unique() %>%
  select(TP_MS) %>% count(TP_MS) %>% select(-TP_MS)
# MOTUs lvl-1
data_rarefaction_lvl1 <- All_occurences_final_geography %>% 
  filter(Dataset != "Genotyping") %>% 
  select(TP_lvl1, Curated_Geography_ID) %>% 
  unique() %>%
  select(TP_lvl1) %>% count(TP_lvl1) %>% select(-TP_lvl1)
# MOTUs lvl-2
data_rarefaction_lvl2 <- All_occurences_final_geography %>% 
  filter(Dataset != "Genotyping") %>% 
  select(TP_lvl2, Curated_Geography_ID) %>% 
  unique() %>%
  select(TP_lvl2) %>% count(TP_lvl2) %>% select(-TP_lvl2)
# MOTUs lvl-3
data_rarefaction_lvl3 <- All_occurences_final_geography %>% 
  filter(Dataset != "Genotyping") %>% 
  select(TP_lvl3, Curated_Geography_ID) %>% 
  unique() %>%
  select(TP_lvl3) %>% count(TP_lvl3) %>% select(-TP_lvl3)
# ASVs
data_rarefaction_ASV <- ASV_Biogeography_3 %>% gather("ASV_ID", "reads", 2:142) %>% 
  filter(reads != 0) %>% 
  select(-reads) %>% 
  unique() %>% 
  select(ASV_ID) %>% count(ASV_ID) %>% select(-ASV_ID)
# And basetypes
data_rarefaction_Basetype <- Data_sanger %>%
  unite("Basetype", `37F_seq`:`49E_seq`, remove = FALSE,sep= "") %>%
  select(Internal_ID,Basetype)  %>% 
  right_join(Key_assignement %>% select(Basetype_ID, Basetype)) %>% 
  select(-Basetype) %>% 
  left_join(Geography_Curated %>% select('Original Sample ID',Curated_Geography_ID) %>% dplyr::rename("Internal_ID" ="Original Sample ID")) %>%
  select(-Internal_ID) %>% 
  unique() %>% 
  select(Basetype_ID) %>% count(Basetype_ID) %>% select(-Basetype_ID)

#Now calculate the rarefaction curve for each level
ggiNEXT(iNEXT(sapply(data_rarefaction_MS, as.numeric),
              q=0, datatype="abundance", knots=100),
              type=1, se = TRUE, grey =TRUE)+
              ylim(c(0,400)) +
              xlim(c(0,5000))+
              theme(aspect.ratio=1)

ggiNEXT(iNEXT(sapply(data_rarefaction_lvl1, as.numeric),
              q=0, datatype="abundance", knots=100),
              type=1, se = TRUE, grey =TRUE)+
              ylim(c(0,400)) +
              xlim(c(0,5000))+
              theme(aspect.ratio=1)

ggiNEXT(iNEXT(sapply(data_rarefaction_lvl2, as.numeric),
              q=0, datatype="abundance", knots=100),
              type=1, se = TRUE, grey =TRUE)+
              ylim(c(0,400)) +
              xlim(c(0,5000))+
              theme(aspect.ratio=1)

ggiNEXT(iNEXT(sapply(data_rarefaction_lvl3, as.numeric),
              q=0, datatype="abundance", knots=100),
              type=1, se = TRUE, grey =TRUE)+
              ylim(c(0,400)) +
              xlim(c(0,5000))+
              theme(aspect.ratio=1)

ggiNEXT(iNEXT(sapply(data_rarefaction_ASV, as.numeric),
              q=0, datatype="abundance", knots=1000),
              type=1, se = TRUE, grey =TRUE)+
              ylim(c(0,400)) +
              xlim(c(0,5000))+
              theme(aspect.ratio=1)

ggiNEXT(iNEXT(sapply(data_rarefaction_Basetype, as.numeric),
              q=0, datatype="abundance", knots=1000),
              type=1, se = TRUE, grey =TRUE)+
              ylim(c(0,400)) +
              xlim(c(0,5000))+
              theme(aspect.ratio=1)

# Summarizing the number of taxonomic path revealed by each dataset

# For morphospecies
N_TP_MS_All   <- All_occurences_final %>% select(TP_MS)   %>% unique() %>% count()

N_TP_MS_common <- ASV_Occurence_taxo %>% select(TP_MS)   %>% unique() %>% 
  inner_join(Sanger_Occurence_taxo %>% select(TP_MS)   %>% unique()) %>% count()

N_TP_MS_unique_env <- ASV_Occurence_taxo %>% select(TP_MS)   %>% unique() %>% 
  anti_join(Sanger_Occurence_taxo %>% select(TP_MS)   %>% unique()) %>% count()

N_TP_MS_unique_SC <- Sanger_Occurence_taxo %>% select(TP_MS)   %>% unique() %>% 
  anti_join(ASV_Occurence_taxo %>% select(TP_MS)   %>% unique()) %>% count()

# For Level-1
N_TP_lvl1_All   <- All_occurences_final %>% select(TP_lvl1) %>% filter(TP_lvl1 != "Not_Attributed") %>% unique() %>% count()

N_TP_lvl1_common <- ASV_Occurence_taxo %>% select(TP_lvl1)  %>% filter(TP_lvl1 != "Not_Attributed") %>% unique() %>% 
  inner_join(Key_assignement_MOTUS_lvl1 %>% select(Taxo_lvl1) %>% dplyr::rename("TP_lvl1"="Taxo_lvl1") %>% unique()) %>% count()

N_TP_lvl1_unique_env <- ASV_Occurence_taxo %>% select(TP_lvl1)  %>% filter(TP_lvl1 != "Not_Attributed") %>% unique() %>% 
  anti_join(Key_assignement_MOTUS_lvl1 %>% select(Taxo_lvl1) %>% dplyr::rename("TP_lvl1"="Taxo_lvl1") %>% unique()) %>% count()

N_TP_lvl1_unique_SC <-Key_assignement_MOTUS_lvl1 %>% select(Taxo_lvl1) %>% dplyr::rename("TP_lvl1"="Taxo_lvl1") %>% unique() %>% 
  anti_join(ASV_Occurence_taxo %>% select(TP_lvl1)  %>% filter(TP_lvl1 != "Not_Attributed") %>% unique()) %>% count()

# For Level-2
N_TP_lvl2_All   <- All_occurences_final %>% select(TP_lvl2) %>% filter(TP_lvl2 != "Not_Attributed") %>% unique() %>% count()

N_TP_lvl2_common <- ASV_Occurence_taxo %>% select(TP_lvl2) %>% filter(TP_lvl2 != "Not_Attributed")   %>% unique() %>% 
  inner_join(Key_assignement_MOTUS_lvl2 %>% select(Taxo_lvl2) %>% dplyr::rename("TP_lvl2"="Taxo_lvl2") %>% unique()) %>% count()

N_TP_lvl2_unique_env <- ASV_Occurence_taxo %>% select(TP_lvl2) %>% filter(TP_lvl2 != "Not_Attributed")  %>% unique() %>% 
  anti_join(Key_assignement_MOTUS_lvl2 %>% select(Taxo_lvl2) %>% dplyr::rename("TP_lvl2"="Taxo_lvl2") %>% unique()) %>% count()

N_TP_lvl2_unique_SC <- Key_assignement_MOTUS_lvl2 %>% select(Taxo_lvl2) %>% dplyr::rename("TP_lvl2"="Taxo_lvl2") %>% unique() %>% 
  anti_join(ASV_Occurence_taxo %>% select(TP_lvl2) %>% filter(TP_lvl2 != "Not_Attributed")  %>% unique()) %>% count()

# For Level-3
N_TP_lvl3_All   <- Key_assignement_MOTUS_lvl3 %>% select(Taxo_lvl3) %>% dplyr::rename("TP_lvl3"="Taxo_lvl3") %>% unique() %>% count()

N_TP_lvl3_common <- ASV_Occurence_taxo %>% select(TP_lvl3) %>% filter(TP_lvl3 != "Not_Attributed")  %>% unique() %>% 
  inner_join(Key_assignement_MOTUS_lvl3 %>% select(Taxo_lvl3) %>% dplyr::rename("TP_lvl3"="Taxo_lvl3") %>% unique()) %>% count()

N_TP_lvl3_unique_env <- ASV_Occurence_taxo %>% select(TP_lvl3) %>% filter(TP_lvl3 != "Not_Attributed")    %>% unique() %>% 
  anti_join(Key_assignement_MOTUS_lvl3 %>% select(Taxo_lvl3) %>% dplyr::rename("TP_lvl3"="Taxo_lvl3") %>% unique()) %>% count()

N_TP_lvl3_unique_SC <- Key_assignement_MOTUS_lvl3 %>% select(Taxo_lvl3) %>% dplyr::rename("TP_lvl3"="Taxo_lvl3") %>% unique() %>% 
  anti_join(ASV_Occurence_taxo %>% select(TP_lvl3) %>% filter(TP_lvl3 != "Not_Attributed")    %>% unique()) %>% count()

# ASV
N_ASV      <- ASV_data %>% select(ASV) %>% count()
# Basetype
N_Basetype <- Key_assignement %>% select(Basetype_ID) %>% unique() %>% count()


#Create a matrix

Matrix_diversity <-
  rbind(
    c("Occurence" , "01_TP_MS"            ,  "02_TP_lvl1"           , "03_TP_lvl2"           , "04_TP_lvl3"           , "05_ASV" , "06_Basetype" ),
    c("03_SC_ENV"  , N_TP_MS_common     ,  N_TP_lvl1_common    , N_TP_lvl2_common    , N_TP_lvl3_common    , 0     , 0          ),
    c("02_Env_only", N_TP_MS_unique_env ,  N_TP_lvl1_unique_env, N_TP_lvl2_unique_env, N_TP_lvl3_unique_env, N_ASV , 0          ),
    c("01_SC_only" , N_TP_MS_unique_SC  ,  N_TP_lvl1_unique_SC , N_TP_lvl2_unique_SC,  N_TP_lvl3_unique_SC , 0     , N_Basetype )
  )%>% row_to_names(row_number = 1) %>% as.data.frame() %>% gather("Taxo_level","N_Taxa",2:7) %>% as.data.frame()


Matrix_diversity$N_Taxa <- as.numeric(Matrix_diversity$N_Taxa)
Matrix_diversity$Occurence <- as.character(Matrix_diversity$Occurence)

#Now Plot the data
ggplot(Matrix_diversity, aes(x=Taxo_level, y=N_Taxa, fill=Occurence))+
  geom_bar(stat="identity")+
  theme_bw()

# Now look at the rate of discovery depending of time

Occurence_year_MS <- All_occurences_final %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Year_of_Collection)) %>% 
  select(Year_of_Collection, TP_MS) %>% unique() %>% 
  mutate(pres=1) %>% 
  spread(TP_MS,pres,fill=0)

Occurence_year_lvl1 <- All_occurences_final %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Year_of_Collection)) %>% 
  select(Year_of_Collection, TP_lvl1) %>% unique() %>% 
  mutate(pres=1) %>% 
  spread(TP_lvl1,pres,fill=0)

Occurence_year_lvl2 <- All_occurences_final %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Year_of_Collection)) %>% 
  select(Year_of_Collection, TP_lvl2) %>% unique() %>% 
  mutate(pres=1) %>% 
  spread(TP_lvl2,pres,fill=0)

Occurence_year_lvl3 <- All_occurences_final %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Year_of_Collection)) %>% 
  select(Year_of_Collection, TP_lvl3) %>% unique() %>% 
  mutate(pres=1) %>% 
  spread(TP_lvl3,pres,fill=0)

Occurence_year_ASV <- ASV_Biogeography_3 %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Year_of_Collection)) %>%
  select(-Curated_Geography_ID) %>% 
  unique() %>% 
  relocate(Year_of_Collection) %>% 
  gather(ASV, n_obs,2:250) %>% 
  filter(n_obs != 0) %>% 
  select(-n_obs) %>%
  unique() %>% 
  mutate(pres=1) %>% 
  spread(ASV, pres,fill=0)
  
Occurence_year_Basetype <- Data_sanger %>%
  unite("Basetype", `37F_seq`:`49E_seq`, remove = FALSE,sep= "") %>%
  select(Internal_ID,Basetype)  %>% 
  right_join(Key_assignement %>% select(Basetype_ID, Basetype)) %>% 
  select(-Basetype) %>% 
  left_join(Geography_Curated %>% select('Original Sample ID',Curated_Geography_ID,Year_of_Collection) %>% dplyr::rename("Internal_ID" ="Original Sample ID")) %>% 
  select(-Curated_Geography_ID,-Internal_ID) %>%
  unique() %>% 
  mutate(pres=1) %>% 
  spread(Basetype_ID, pres,fill=0) 


# Export tables to do reformating
write.table(Occurence_year_MS  ,      "Occurence_year_MS.txt")
write.table(Occurence_year_lvl1,      "Occurence_year_lvl1.txt")
write.table(Occurence_year_lvl2,      "Occurence_year_lvl2.txt")
write.table(Occurence_year_lvl3,      "Occurence_year_lvl3.txt")
write.table(Occurence_year_ASV ,      "Occurence_year_ASV.txt")
write.table(Occurence_year_Basetype,  "Occurence_year_Basetype.txt")

#import formated data

Occurence_year_final <- readxl::read_xlsx(path = "Collection_per_year.xlsx", sheet = 1, col_types="guess")

ggplot(Occurence_year_final, aes(x=Year, y=N_obs, colour = Dataset))+
  geom_point()+
  geom_step()+
  theme_bw()

# Now look into the quality of attribution of the ASVs

Match_ASV_lvl3 <- ASv_MOTU_lvl3 %>% select(ASV) %>% add_column(Quality_assignement = "11_Match_lvl3")
Match_ASV_lvl2 <- ASv_MOTU_lvl2 %>% select(ASV) %>% add_column(Quality_assignement = "10_Match_lvl2")
Match_ASV_lvl1 <- ASv_MOTU_lvl1 %>% select(ASV) %>% add_column(Quality_assignement = "09_Match_lvl1")

attributed_ASV_lvl2 <- Molecular_Nomenclature_Novel_ASVs %>%  filter(`Putative new Genotype` == "NO") %>% filter(Level_attribution == "MOTUlvl2") %>% select(ASV_Name) %>% dplyr::rename("ASV" = "ASV_Name") %>% add_column(Quality_assignement = "08_Attributed_lvl2")
attributed_ASV_lvl1 <- Molecular_Nomenclature_Novel_ASVs %>%  filter(`Putative new Genotype` == "NO") %>% filter(Level_attribution == "MOTUlvl1") %>% select(ASV_Name) %>% dplyr::rename("ASV" = "ASV_Name") %>% add_column(Quality_assignement = "07_Attributed_lvl1")
attributed_ASV_MS   <- Molecular_Nomenclature_Novel_ASVs %>%  filter(`Putative new Genotype` == "NO") %>% filter(Level_attribution == "MS") %>% select(ASV_Name) %>% dplyr::rename("ASV" = "ASV_Name") %>% add_column(Quality_assignement = "06_Attributed_MS")

Unclear_ASV_Attribution <- Molecular_Nomenclature_Novel_ASVs %>%  filter(`Putative new Genotype` == "UNCLEAR") %>%  select(ASV_Name) %>% dplyr::rename("ASV" = "ASV_Name") %>% add_column(Quality_assignement = "05_Unclear")

New_ASV_lvl2 <- Molecular_Nomenclature_Novel_ASVs %>%  filter(`Level of novelty` == "MOTUlvl2")  %>% select(ASV_Name) %>% dplyr::rename("ASV" = "ASV_Name") %>% add_column(Quality_assignement = "04_New_lvl2") 
New_ASV_lvl1 <- Molecular_Nomenclature_Novel_ASVs %>%  filter(`Level of novelty` == "MOTUlvl1")  %>% select(ASV_Name) %>% dplyr::rename("ASV" = "ASV_Name") %>% add_column(Quality_assignement = "03_New_lvl1")
New_ASV_MS   <- Molecular_Nomenclature_Novel_ASVs %>%  filter(`Level of novelty` == "MS")  %>% select(ASV_Name) %>% dplyr::rename("ASV" = "ASV_Name") %>% add_column(Quality_assignement = "02_New_MS")

Not_attributed_ASV <- ASV_data %>% select(ASV) %>% anti_join(ASV_all_Assignement %>% select(ASV_ID) %>% dplyr::rename("ASV" = "ASV_ID")) %>% add_column(Quality_assignement = "01_Unresolved")

Single_occurence_ASV <- readxl::read_xlsx("../../02_Datasets/07_ASVs_PF/PF_ASVs_Segmented.xlsx", sheet = 1, col_types="guess") %>% 
  left_join(N_Obs_ASV %>% dplyr::rename("ASV" ="ASV_ID")) %>% anti_join(ASV_data) %>% 
  select(ASV)  %>% add_column(Quality_assignement = "00_single_occurence")


view(Level_attribution_ASV)

Level_attribution_ASV <- rbind(Match_ASV_lvl3,
                               Match_ASV_lvl2,
                               Match_ASV_lvl1,
                               attributed_ASV_lvl2,
                               attributed_ASV_lvl1,
                               attributed_ASV_MS,
                               Unclear_ASV_Attribution,
                               New_ASV_lvl2,
                               New_ASV_lvl1,
                               New_ASV_MS,
                               Not_attributed_ASV,
                               Single_occurence_ASV) %>% count(Quality_assignement)
  
ggplot(Level_attribution_ASV, aes(x=2, y=n, fill=Quality_assignement))+
  geom_bar(stat="identity")+
  xlim(0.5, 2.5) +
  coord_polar(theta = "y")+
  labs(x=NULL, y=NULL)+
  labs(fill="") +
  theme_bw()

# Compare number of ASV and basetype per morphospecies
view(comparison_ASV_Basetype)
comparison_ASV_Basetype <- Key_assignement %>% select(Clade_curated,Genus_curated,species_curated, Basetype_ID) %>% 
unite("TP_MS", `Clade_curated`:`species_curated`, remove = TRUE,sep= "|") %>% 
count(TP_MS) %>% dplyr::rename("n_Basetype" = "n") %>% 
  full_join(Summary_observation_N_ASVs) %>% replace(is.na(.), 0) %>% 
  separate(col = TP_MS, into = c("Clade","Genus","species"), sep = "\\|")
  
ggplot(comparison_ASV_Basetype, aes(x=n_Basetype,y=n_ASVs, colour = Clade))+
  geom_point()+
  geom_abline(intercept = 0,
              slope = 1)+
  scale_y_sqrt(limits =c(0,70))+
  scale_x_sqrt(limits =c(0,70))+
  theme_bw()+
  coord_fixed()

###################################################################################
#                                                                                 #
#                                   Ecology                                       #
#                                                                                 #
###################################################################################

# Try to find fun and creative ways to look at the diversity of all the cryptic species....
library("UpSetR")

#Plot the preference for temperature
ggplot(All_occurences_final_geography, aes(x=TP_lvl1,y=Temperature_Monthly))+
  geom_point(
    shape = 1,
    #size = 10,
    alpha = .2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_flip()

library("ggplot2")
library("ggproto")


ggplot(All_occurences_final_geography, aes(x=TP_lvl1,y=Temperature_Monthly))+
  geom_vline(data = All_occurences_final_geography, aes(xintercept=Temperature_Monthly), alpha=0.2)+
  facet_wrap(~TP_lvl1)
  
  geom_point(alpha = .1, shape = "|")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_flip()



ggplot(All_occurences_final_geography, aes(x=TP_lvl1,y=Temperature_Monthly))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_flip()




ggplot(All_occurences_final_geography, aes(x=TP_lvl1,y=Temperature_Monthly))+
  geom_jitter(alpha = 0.2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_flip()




g_interval <- 
  ggplot(All_occurences_final_geography, aes(TP_MS, Temperature_Monthly)) +
  scale_color_viridis_d(
    option = "mako", name = "Level:", direction = -1, 
    begin = .15, end = .9
  ) +
  guides(
    color = guide_legend(reverse = TRUE, title.position = "top")
  ) +
  theme(
    legend.position = c(.75, .95), legend.direction = "horizontal",
    legend.text = element_text(family = "Roboto Mono", size = 18),
    legend.title = element_text(face = "bold", size = 22, hjust = .5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )+
  coord_flip()


g_interval +
  ggdist::stat_interval(size = 2)


ggplot(All_occurences_final_geography, aes(TP_lvl2, Temperature_Monthly))+ ggdist::stat_gradientinterval(
  width = .3, color = "black"
)


library(tidyverse)     ## data wrangling + ggplot2
library(colorspace)    ## adjust colors
library(rcartocolor)   ## Carto palettes
library(ggforce)       ## sina plots
library(ggdist)        ## halfeye plots
library(ggridges)      ## ridgeline plots
library(ggbeeswarm)    ## beeswarm plots
library(gghalves)      ## off-set jitter
library(systemfonts)   ## custom fonts


g_interval +
  ggdist::stat_interval(
    .width = c(.25, .5, .95), 
    size = 1
  ) +
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    width = .7, fill = "grey85",
    interval_colour = NA, point_colour = "black",
    shape = 23, stroke = 1.5, point_size = 2, point_fill = "white",
    position = position_nudge(x = .03),
    aes(thickness = stat(f*n))
  ) +
  scale_color_viridis_d(
    option = "mako", name = "Level:", direction = -1, 
    begin = .15, end = .9,
    labels = function(x) paste0(as.numeric(x)*100, "%")
  )



biogeo_simplified <- All_occurences_final_geography %>% 
  select(TP_MS, Region_Curated) %>% 
  unique() %>% 
  mutate(pres = 1)

ggplot(biogeo_simplified, aes(y=TP_MS, x=Region_Curated, size = pres))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



biogeo_simplified <- All_occurences_final_geography %>% 
  select(TP_lvl1, Region_Curated) %>% 
  unique() %>% 
  mutate(pres = 1)

ggplot(biogeo_simplified, aes(y=TP_lvl1, x=Region_Curated, size = pres))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

biogeo_simplified <- All_occurences_final_geography %>% 
  select(TP_lvl2, Region_Curated) %>% 
  unique() %>% 
  mutate(pres = 1)

ggplot(biogeo_simplified, aes(y=TP_lvl1, x=Region_Curated, size = pres))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme_bw()


###################################################################
# Taxa proportion reads in size fractions

#N Reads per taxonomic paths
N_Reads_TP_Lvl_1 <- ASV_all_Assignement %>% select(TP_lvl1, ASV_ID) %>% 
  full_join(PF_ASV_Occurence) %>% select(-ASV_ID) %>%group_by(TP_lvl1) %>%
  summarise_each(funs(sum)) %>%
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.))  %>% 
  select(sum) %>% 
  rownames_to_column(var = "TP_lvl1")

ggplot(N_Reads_TP_Lvl_1, aes(y=TP_lvl1, x=sum)) + 
  geom_col()+
  scale_x_log10()
  
#Proportion Reads environments
Prop_reads_env <- ASV_all_Assignement %>% select(TP_lvl1, ASV_ID) %>% 
  full_join(PF_ASV_Occurence) %>% select(-ASV_ID) %>% group_by(TP_lvl1) %>%
  summarise_each(funs(sum)) %>%
  column_to_rownames(1) %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "Original Sample ID") %>%  
  left_join(Geography_Curated %>% select('Original Sample ID', Curated_Depth)) %>% 
  select(-'Original Sample ID') %>% 
  group_by(Curated_Depth) %>%
  summarise_each(funs(sum)) %>% 
  filter(Curated_Depth != "NA") %>% 
  column_to_rownames(1) %>% as.matrix() %>% 
  prop.table(margin = 2) %>% as.data.frame() %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "TP_lvl1") %>% 
  gather(key = "Depth", value = "prop", 2:4)
  
ggplot(Prop_reads_env, aes(y=TP_lvl1, x=prop, fill= fct_rev(Depth))) + 
  geom_col()

#Proportion Reads size fraction
Prop_reads_SizeFraction <- ASV_all_Assignement %>% select(TP_lvl1, ASV_ID) %>% 
  full_join(PF_ASV_Occurence) %>% select(-ASV_ID) %>% group_by(TP_lvl1) %>%
  summarise_each(funs(sum)) %>%
  column_to_rownames(1) %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "Original Sample ID") %>%  
  left_join(Geography_Curated %>% select('Original Sample ID', Curated_size_fraction)) %>% 
  select(-'Original Sample ID') %>% 
  group_by(Curated_size_fraction) %>%
  summarise_each(funs(sum)) %>% 
  filter(Curated_size_fraction != "NA") %>% 
  column_to_rownames(1) %>% as.matrix() %>% 
  prop.table(margin = 2) %>% as.data.frame() %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "TP_lvl1") %>% 
  gather(key = "SizeFraction", value = "prop", 2:7)

ggplot(Prop_reads_SizeFraction, aes(y=TP_lvl1, x=prop, fill= fct_rev(SizeFraction))) + 
  geom_col()

#Intensity sampling
N_stations_Sampled <- All_occurences_final %>% select(TP_lvl1,Curated_Geography_ID, Dataset) %>%  
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Dataset, pres, fill=0) %>% 
  select(-Curated_Geography_ID) %>% 
  group_by(TP_lvl1) %>%
  summarise_each(funs(sum)) %>% 
  gather(key = "Dataset", value = "N_Obs", 2:4)

ggplot(N_stations_Sampled, aes(y=TP_lvl1, x=N_Obs, fill= fct_rev(Dataset))) + 
  geom_col()

################# Now arrange the plot together
library("ggpubr")
library("scales")
library("viridis")

complete_list_TPlvl1 <- N_stations_Sampled %>% select(TP_lvl1) %>% distinct()


Plot_N_station <- ggplot(N_stations_Sampled %>% filter(TP_lvl1 != "Not_Attributed"),
                         aes(y=fct_rev(TP_lvl1), x=N_Obs, fill= fct_rev(Dataset))) + 
                         geom_col()+
                         theme(legend.position="none")

Plot_N_station2 <- ggplot(N_stations_Sampled %>% filter(TP_lvl1 != "Not_Attributed"),
                         aes(y=fct_rev(TP_lvl1), x=N_Obs, fill= fct_rev(Dataset))) + 
                         geom_col()+ 
                         theme(legend.position="none",
                         axis.text.y =element_blank(),
                         axis.title.y =element_blank())


Plot_N_reads <- ggplot(N_Reads_TP_Lvl_1 %>% right_join(complete_list_TPlvl1)%>% filter(TP_lvl1 != "Not_Attributed"),
                       aes(y=fct_rev(TP_lvl1), x=sum)) + 
                  geom_col()+
                  scale_x_sqrt()+ 
                  theme(legend.position="none",
                  axis.text.y =element_blank(),
                  axis.title.y =element_blank())

Plot_prop_env <- ggplot(Prop_reads_env %>% right_join(complete_list_TPlvl1)%>% filter(TP_lvl1 != "Not_Attributed"),
                        aes(y=fct_rev(TP_lvl1), x=prop, fill= fct_rev(Depth))) +  
                  geom_col()+ 
                  theme(legend.position="none",
                  axis.text.y =element_blank(),
                  axis.title.y =element_blank())+
                  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

Plot_prop_SizeFraction <- ggplot(Prop_reads_SizeFraction %>% right_join(complete_list_TPlvl1)%>% filter(TP_lvl1 != "Not_Attributed"), 
                          aes(y=fct_rev(TP_lvl1), x=prop, fill= fct_rev(SizeFraction))) + 
                          geom_col() + 
                          theme(legend.position="none",
                          axis.text.y =element_blank(),
                          axis.title.y =element_blank())+ 
                          scale_fill_brewer(palette="Dark2")

Plot_SST <- ggplot(All_occurences_final_geography %>% filter(TP_lvl1 != "Not_Attributed"),
                   aes(y=fct_rev(TP_lvl1),x=Temperature_Monthly, color = Temperature_Monthly))+
                   geom_point(alpha = 0.2)+ 
                   theme(legend.position="none",
                         axis.text.y =element_blank(),
                         axis.title.y =element_blank())+ 
  scale_color_viridis()


ggarrange(Plot_N_station,
          Plot_N_station2,
          Plot_N_reads,
          Plot_prop_env, 
          Plot_prop_SizeFraction,
          Plot_SST,
          nrow = 1)
Plot_prop_SizeFraction

ggplot(Prop_reads_SizeFraction %>% right_join(complete_list_TPlvl1)%>% filter(TP_lvl1 != "Not_Attributed"), 
       aes(y=fct_rev(TP_lvl1), x=prop, fill= fct_rev(SizeFraction))) + 
  geom_col() + 
  theme(
        axis.text.y =element_blank(),
        axis.title.y =element_blank())+ 
  scale_fill_brewer(palette="Dark2")

#### Checking lenght foram seq
view(ASV_data)

lenght_ASV <- ASV_data %>% mutate(ASV_Lenght = str_length(Full_sequence)) %>% 
                           select(ASV, ASV_Lenght) %>% 
                           dplyr::rename("ASV_ID" = "ASV") %>% 
                           left_join(ASV_all_Assignement)

ggplot(lenght_ASV %>% filter(Clade_curated_assignement != "NA"), aes(x=ASV_Lenght, fill=Clade_curated_assignement)) + 
  geom_histogram(alpha=0.5)

#################### Ecology ####################

################### Write a clean script for the environmental data analyses
# Subset Annual env variables
Env_data_selection_Annual <- All_occurences_final_geography %>% 
  select(Curated_Geography_ID, 	
         Salinity_Annual,	
         Depth_Annual, 
         Chlorophyll_Annual_Aquamodis, 
         SST_Annual_Aquamodis, 
         PAR_Annual_Aquamodis, 
         POC_Annual_Aquamodis, 
         PIC_Annual_Aquamodis, 
         EuphoticDepth_Annual_Aquamodis) %>% 
  filter(Salinity_Annual!= "NaN",	
         Depth_Annual!= "NaN",	
         Chlorophyll_Annual_Aquamodis != "NaN", 
         SST_Annual_Aquamodis!= "NaN", 
         PAR_Annual_Aquamodis!= "NaN", 
         POC_Annual_Aquamodis!= "NaN", 
         PIC_Annual_Aquamodis!= "NaN",
         EuphoticDepth_Annual_Aquamodis!= "NaN") %>% 
  distinct() %>% 
  group_by(Curated_Geography_ID) %>% 
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("Curated_Geography_ID")  
  
# Subset Monthly env variables
Env_data_selection_Monthly <- All_occurences_final_geography %>% 
  select(Curated_Geography_ID, 	
         Salinity_Monthly,	
         Depth_Monthly, 
         Chlorophyll_Monthly_Aquamodis, 
         SST_Monthly_Aquamodis, 
         PAR_Monthly_Aquamodis, 
         POC_Monthly_Aquamodis, 
         PIC_Monthly_Aquamodis, 
         EuphoticDepth_Monthly_Aquamodis) %>% 
  filter(Salinity_Monthly != "NaN",	
         Depth_Monthly != "NaN",	
         Chlorophyll_Monthly_Aquamodis != "NaN", 
         SST_Monthly_Aquamodis != "NaN", 
         PAR_Monthly_Aquamodis != "NaN", 
         POC_Monthly_Aquamodis != "NaN", 
         PIC_Monthly_Aquamodis != "NaN",
         EuphoticDepth_Monthly_Aquamodis!= "NaN") %>% 
  distinct() %>% 
  mutate_at(c('Salinity_Monthly', 'Depth_Monthly', 'Chlorophyll_Monthly_Aquamodis','SST_Monthly_Aquamodis','PAR_Monthly_Aquamodis','PAR_Monthly_Aquamodis','POC_Monthly_Aquamodis','PIC_Monthly_Aquamodis','EuphoticDepth_Monthly_Aquamodis'), as.numeric) %>% 
  group_by(Curated_Geography_ID) %>% 
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("Curated_Geography_ID")  




# Now calculate PCA
# Annual data
#First do the PCA
res.pca.annual <- prcomp(Env_data_selection_annual, scale = TRUE)
# Visualize the eingenvalues
fviz_eig(res.pca.annual)
# Visualize the contribution of the parameters to the principal components
fviz_pca_var(res.pca.annual,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_biplot(res.pca.annual, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)


# Monthly data
#First do the PCA
res.pca.monthly <- prcomp(Env_data_selection_Monthly, scale = TRUE)
# Visualize the eingenvalues
fviz_eig(res.pca.monthly)
# Visualize the contribution of the parameters to the principal components
fviz_pca_var(res.pca.monthly,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_biplot(res.pca.monthly, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)



# Remove Outliers in POC values that cause the PCA to be unbalanded
# Script taken from: https://www.geeksforgeeks.org/how-to-remove-outliers-from-multiple-columns-in-r-dataframe/
# create detect outlier function
detect_outlier <- function(x) {
  
  # calculate first quantile
  Quantile1 <- quantile(x, probs=.00)
  
  # calculate third quantile
  Quantile3 <- quantile(x, probs=.95)
  
  # calculate inter quartile range
  IQR = Quantile3-Quantile1
  
  # return true or false
  x > Quantile3 | x < Quantile1
}

# create remove outlier function
# create remove outlier function
remove_outlier <-  function(dataframe,
                            columns=names(dataframe)) {
  
  # for loop to traverse in columns vector
  for (col in columns) {
    
    # remove observation if it satisfies outlier function
    dataframe <- dataframe[!detect_outlier(dataframe[[col]]), ]
  }
  
  # return dataframe
  print("Remove outliers")
  print(dataframe)
}


# Now remove outliers
Env_data_selection_Annual_no_outliers  <- remove_outlier(Env_data_selection_Annual, c('POC_Annual_Aquamodis'))
Env_data_selection_Monthly_no_outliers <- remove_outlier(Env_data_selection_Monthly, c('POC_Monthly_Aquamodis'))



# Check that PCA results are now more balanced:
# For Annual data
res.pca.Annual.no.outliers <- prcomp(Env_data_selection_Annual_no_outliers, scale = TRUE)
fviz_eig(res.pca.Annual.no.outliers)

fviz_pca_var(res.pca.Annual.no.outliers,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_biplot(res.pca.Annual.no.outliers, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color
                label = "var"
)


# For Monthly data
res.pca.Monthly.no.outliers <- prcomp(Env_data_selection_Monthly_no_outliers, scale = TRUE)
fviz_eig(res.pca.Monthly.no.outliers)

fviz_pca_var(res.pca.Monthly.no.outliers,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_biplot(res.pca.Monthly.no.outliers, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color
                label = "var"
)







#Recover PCA corrdinates
# For Annual data
res.PCA.Annual <- get_pca_ind(res.pca.Annual.no.outliers)
PCA_Coordinates_Env_parameters_Annual <- res.PCA.Annual$coord %>%as.data.frame() %>% rownames_to_column(var = 'Curated_Geography_ID')

# For Monthly data
res.PCA.Monthly <- get_pca_ind(res.pca.Monthly.no.outliers)
PCA_Coordinates_Env_parameters_Monthly <- res.PCA.Monthly$coord %>%as.data.frame() %>% rownames_to_column(var = 'Curated_Geography_ID')

#

#Now merge with occurences of species
All_occurences_vs_PCA <- All_occurences_final %>% 
  left_join(PCA_Coordinates_Env_parameters_Annual %>% select(Curated_Geography_ID, Dim.1, Dim.2) %>% dplyr::rename("PCA_1_Annual" = "Dim.1", "PCA_2_Annual" = "Dim.2")) %>%
  left_join(PCA_Coordinates_Env_parameters_Monthly %>% select(Curated_Geography_ID, Dim.1, Dim.2) %>% dplyr::rename("PCA_1_Monthly" = "Dim.1", "PCA_2_Monthly" = "Dim.2")) %>%
  filter(TP_lvl2 != "Not_Attributed") %>% 
  separate(col = TP_lvl2, into = c("Clade","Genus","species", "MOTU_lvl_1", "MOTU_lvl_2"), sep = "\\|") %>% 
  unite("Morphotaxa", Genus,species, remove = FALSE,sep= " ") %>%
  unite("MOTU", MOTU_lvl_1,MOTU_lvl_2, remove = FALSE,sep= "")
            
#Now subset species with cryptic diversity
Data_format_ecology_CD_only <- rbind(All_occurences_vs_PCA %>% filter(Morphotaxa == "Candeina nitida"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Globigerinita glutinata"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Globigerinita minuta"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Globorotalia truncatulinoides"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Globorotalia inflata"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Neogloboquadrina incompta"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Neogloboquadrina pachyderma"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Pulleniatina obliquiloculata"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Globigerina bulloides"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Globigerinella siphonifera"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Globigerinoides ruber-albus"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Globigerinoides conglobatus"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Orbulina universa"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Hastigerina pelagica"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Turborotalita clarkei"),
                                     All_occurences_vs_PCA %>% filter(Morphotaxa == "Turborotalita quinqueloba"))
                                     
                                     
# And plot the results
ggplot(Data_format_ecology_CD_only,
       aes(x=PCA_1_Annual, y=PCA_2_Annual,col=MOTU))+
       geom_point(size = 2) +
       facet_wrap(.~Morphotaxa) +
       coord_equal()+scale_colour_viridis_d()+
      geom_convexhull(alpha=0)
  
ggplot(Data_format_ecology_CD_only,
       aes(x=PCA_1_Monthly, y=PCA_2_Monthly,col=MOTU))+
  geom_point(size = 3) +
  facet_wrap(.~Morphotaxa, nrow = 3) +
  coord_fixed()+scale_colour_viridis_d()+
  geom_convexhull(alpha=0)+
  theme_light()




write.table(Data_format_ecology_CD_only, "Data_for_the_magician.txt")

# Ask a competent friend to do the stats....
# Now calculate ANOVA
library(vegan)


lm(Data_format_ecology_CD_only~factor(PCA_1_Annual)*factor(PCA_2_Annual))



PERMANOVA.res <- adonis(PCA_1_Annual ~ PCA_2_Annual + MOTU, data = (Data_format_ecology_CD_only %>% filter(PCA_1_Annual != "NA") %>% filter(Morphotaxa == "Neogloboquadrina pachyderma")))

summary(PERMANOVA.res)


test <- unlist(PERMANOVA.res)



PERMANOVA.res2 <- adonis2(PCA_1_Annual ~ PCA_2_Annual + MOTU, data = Data_format_ecology_CD_only %>% filter(PCA_1_Annual != "NA"))


########################### Find a simple way to disply the biogeography

Biogeo_simplified_data <- All_occurences_final_geography %>% 
  select(TP_lvl2, Curated_Geography_ID, Region_Curated_3) %>% 
  distinct() %>% 
  select(-Curated_Geography_ID) %>% 
  group_by(TP_lvl2, Region_Curated_3) %>% 
  count() %>% 
  filter(TP_lvl2 != "Not_Attributed" ) %>%
  filter(Region_Curated_3 != "NA") %>% 
  spread(TP_lvl2, n, fill=0) %>% 
  gather(TP_lvl2, n_obs, 2:95) %>% 
  mutate(N_Obs_log = log10(n_obs+1))
  
ggplot(Biogeo_simplified_data, aes(x=Region_Curated_3, y=TP_lvl2, col=N_Obs_log))+
  geom_point(size= 1)+
  scale_colour_viridis_c()+
  theme_bw()

ggplot(Biogeo_simplified_data, aes(x=n))+
  geom_histogram(bins = 120)

Biogeo_T_simplified_data <- All_occurences_final_geography %>% 
  select(TP_lvl2, Curated_Geography_ID, SST_Monthly_Aquamodis) %>% 
  distinct() %>% 
  select(-Curated_Geography_ID) %>% 
  filter(TP_lvl2 != "Not_Attributed" ) %>%
  filter(SST_Monthly_Aquamodis != "NA") %>%
  filter(SST_Monthly_Aquamodis != "NaN") %>%
  mutate_at(c('SST_Monthly_Aquamodis'), as.numeric)


ggplot(Biogeo_T_simplified_data , aes(x=SST_Monthly_Aquamodis, y=TP_lvl2))+
  geom_point(size= 1)+
  theme_light()+
  xlim(-5, 35)



### Do a fig with FORCENS to show how species with CD are distributed in respect to latitude
FORCENS <- readxl::read_xlsx(path = "../../02_Datasets/10_FORCENS/ForCenS_reduced.xlsx", sheet = 1, col_types="guess") %>% select(1:48, unidentified)
FORCENS_taxo <- readxl::read_xlsx(path = "../../02_Datasets/10_FORCENS/ForCenS_reduced.xlsx", sheet = 2, col_types="guess")

library(tidyverse)

Forcens_curated <- FORCENS %>% 
  gather(Taxonomy, Proportions,8:49) %>% 
  filter(Proportions != "N/A") %>% 
  mutate_at(c('Proportions'), as.numeric) %>% 
  spread(Taxonomy, Proportions, fill=0)

Forcens_N_taxa <- Forcens_curated %>% 
  select(Sample_ID, 8:49) %>%
  gather(Taxonomy, Proportions,2:43) %>% 
  group_by(Sample_ID) %>% 
  summarise(N_Morphotaxa = sum(Proportions > 0, na.rm = TRUE))%>% 
  add_column(Cryptic_Diversity = "ALL")
  
Forcens_N_taxa_CD <- Forcens_curated %>% 
  select(Sample_ID, 8:49) %>%
  gather(Taxonomy, Proportions,2:43) %>%
  left_join(FORCENS_taxo %>% dplyr::rename("Taxonomy" = "Morphotaxa")) %>% 
  filter(Identified_cryptic_diversity == "YES") %>% 
  group_by(Sample_ID) %>% 
  summarise(N_Morphotaxa = sum(Proportions > 0, na.rm = TRUE)) %>% 
  add_column(Cryptic_Diversity = "YES")

Forcens_N_taxa_noCD <- Forcens_curated %>% 
  select(Sample_ID, 8:49) %>%
  gather(Taxonomy, Proportions, 2:43) %>%
  left_join(FORCENS_taxo %>% dplyr::rename("Taxonomy" = "Morphotaxa")) %>% 
  filter(Identified_cryptic_diversity == "NO") %>% 
  group_by(Sample_ID) %>% 
  summarise(N_Morphotaxa = sum(Proportions > 0, na.rm = TRUE)) %>% 
  add_column(Cryptic_Diversity = "NO")

Forcens_distrib_CD <- rbind(Forcens_N_taxa, Forcens_N_taxa_CD,Forcens_N_taxa_noCD) %>% left_join(FORCENS %>% select(Sample_ID, Latitude))

ggplot(Forcens_distrib_CD, aes(x=Latitude, y=N_Morphotaxa, colour = Cryptic_Diversity))+
  geom_point()+
  geom_smooth(method = "loess")+
  ylim(0, 30)+
  theme_bw()+
  scale_x_continuous(limits =c(-90,90), breaks=c(-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90))+
  scale_fill_viridis_c()


Forcens_Prop_CD <- Forcens_curated %>% 
  select(Sample_ID, 8:49) %>%
  gather(Morphotaxa, Proportions,2:43) %>% 
  left_join(FORCENS_taxo) %>%
  select(-Morphotaxa) %>% 
  group_by(Sample_ID, Identified_cryptic_diversity) %>%
  summarise_each(funs(sum)) %>% 
  filter(Identified_cryptic_diversity == "YES") %>% 
  select(-Identified_cryptic_diversity) 


Forcens_Prop_CD <- Forcens_curated %>% 
  select(Sample_ID, 8:49) %>%
  column_to_rownames("Sample_ID") %>%
  as.matrix() %>% 
  prop.table(margin = 1) %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample_ID") %>% 
  gather(Morphotaxa, Proportions,2:43) %>% 
  left_join(FORCENS %>% select(Sample_ID, Class_latitude)) %>% 
  select(-Sample_ID) %>% 
  group_by(Class_latitude, Morphotaxa) %>%
  summarise_each(funs(mean)) %>% 
  left_join(FORCENS_taxo) %>% 
  drop_na()


library(viridis)

ggplot(Forcens_Prop_CD, aes(x=Class_latitude, y=Proportions, fill = Identified_cryptic_diversity))+
  geom_col(color = "black")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_viridis_d(option="inferno")


ggplot(Forcens_Prop_CD, aes(x=Latitude, y=Proportions, colour = Identified_cryptic_diversity))+
  geom_point(alpha = 0.05)+
  geom_smooth(method = "loess")+
  #  xlim(-90,90)+
  ylim(0, 1)+
  scale_x_continuous(limits =c(-90,90), breaks=c(-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90))


ggplot(Forcens_Prop_CD, aes(x=Latitude, y=Proportions, colour = Identified_cryptic_diversity))+
  geom_point(alpha = 0.05)+
  geom_smooth(method = "loess")+
  #  xlim(-90,90)+
  ylim(0, 1)+
  scale_x_continuous(limits =c(-90,90), breaks=c(-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90))









 
  filter(Identified_cryptic_diversity == "YES") %>% 
  select(-Identified_cryptic_diversity) 







Forcens_comp <- Forcens_curated %>% 
  select(Sample_ID, Class_latitude) %>% 
  left_join(Forcens_N_taxa) %>% 
  left_join(Forcens_Prop_CD) %>% 
  gather(Dataset, data, 3:4)
  
ggplot(Forcens_comp, aes(x=Class_latitude,y=data))+
  geom_boxplot()+
  facet_grid(Dataset~., scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



library(tidyverse)


spread(TP_MS, pres, fill=0)







#######
# ASV_all_assignement -> Taxonomy of all AVS
# ASV_Occurence -> Matrix with occurence of all ASV in all samples
# Biogeography_Cordier -> All environmental samples and their associated metadata and data
# Geography_curated -> ALL env data

view(ASV_all_Assignement)
view(ASV_Occurence)
view(Geography_Curated)

# Prepare the file with the relevant variables to test
Metadata_Tara_Lite <- Geography_Curated %>% 
  select(`Original Sample ID`, Region_Curated_2, Curated_Geography_ID, Curated_size_fraction, Curated_Depth)

ASV_all_Assignement_modif <- ASV_all_Assignement %>% 
  select(-TP_lvl1, -TP_lvl2, -TP_lvl3) %>% 
  unite("TP_lvl1", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, remove = FALSE,sep= "|") %>%
  unite("TP_lvl2", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, MOTU_lvl_2, remove = FALSE,sep= "|")

################################################################################################

# Look at the biogeography specifically
# Create observation dataset at the Morphospecies level
Matrix_occurence_MS <- ASV_all_Assignement_modif %>% 
  select(TP_MS, ASV_ID) %>% 
  left_join(ASV_Occurence) %>%
  group_by(TP_MS) %>% 
  summarise(across(2:2084, sum)) %>% 
  column_to_rownames(1) %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'Original Sample ID') %>% 
  left_join(Metadata_Tara_Lite %>% select(`Original Sample ID`, Curated_Geography_ID)) %>% 
  select(-`Original Sample ID`) %>% 
  group_by(Curated_Geography_ID) %>% 
  summarise_all(sum) %>% 
  inner_join(Metadata_Tara_Lite %>% select(Curated_Geography_ID, Region_Curated_2) %>% distinct()) %>% 
  select(-Curated_Geography_ID) %>% 
  rowwise() %>%
  mutate(sum_cols = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
  filter(sum_cols > 0) %>% select(-sum_cols)

tmp_MS <- Matrix_occurence_MS %>% select(-Region_Curated_2) %>% 
  as.matrix() %>% 
  prop.table(margin = 1) %>% 
  as.data.frame()

#Now calculate NMDS at MS level
dissimilarity_matrix_MS <- vegdist(decostand(tmp_MS, "pa"), method = "jaccard", binary = T)
nmds_result_MS <- metaMDS(dissimilarity_matrix_MS, k = 2, trymax=100)
plot(nmds_result_MS$points, col = as.numeric(as.factor(Matrix_occurence_MS$Region_Curated_2)), pch = 16, asp = 1)

# umap ordination
umap.config <- umap.defaults
umap.config$input <- "dist"
umap_result_MS <- umap(as.matrix(dissimilarity_matrix_MS), config = umap.config)
plot(umap_result_MS$layout, col = as.numeric(as.factor(Matrix_occurence_MS$Region_Curated_2)), pch = 16, asp = 1)
# very prone to changes in parameters (without any clear guidelines on how to tune them)

# PERMANOVA
table(Matrix_occurence_MS$Region_Curated_2)
# potentially quite unbalanced, but leave as is for now
tmp_permanova <- adonis2(dissimilarity_matrix_MS ~ Matrix_occurence_MS$Region_Curated_2, sqrt.dist = T)
source("C:/Users/morard/Desktop/permanova_pairwise.R")
tmp_permanova_pw <- permanova_pairwise(dissimilarity_matrix_MS, Matrix_occurence_MS$Region_Curated_2, sqrt.dist = T, padj = "fdr")
tmp_permanova_pw


# PERMANOVA
table(Matrix_occurence_lvl1$Region_Curated_2)
# potentially quite unbalanced, but leave as is for now
tmp_permanova_2 <- adonis2(dissimilarity_matrix_lvl1 ~ Matrix_occurence_lvl1$Region_Curated_2, sqrt.dist = T)
tmp_permanova_pw_2 <- permanova_pairwise(dissimilarity_matrix_lvl1, Matrix_occurence_lvl1$Region_Curated_2, sqrt.dist = T, padj = "fdr")



tmp_permanova_pw_2



# PERMANOVA
table(Matrix_occurence_lvl2$Region_Curated_2)
# potentially quite unbalanced, but leave as is for now
tmp_permanova_3 <- adonis2(dissimilarity_matrix_lvl2 ~ Matrix_occurence_lvl1$Region_Curated_2, sqrt.dist = T)
tmp_permanova_pw_3 <- permanova_pairwise(dissimilarity_matrix_lvl2, Matrix_occurence_lvl2$Region_Curated_2, sqrt.dist = T, padj = "fdr")



tmp_permanova_pw_2




# betadispersion
tmp_betadisper <- betadisper(dissimilarity_matrix_MS, Matrix_occurence_MS$Region_Curated_2, sqrt.dist = T)
plot(tmp_betadisper)
boxplot(tmp_betadisper)

# ANOSIM
tmp_anosim <- anosim(dissimilarity_matrix_MS, Matrix_occurence_MS$Region_Curated_2)
tmp_anosim_pw <- anosim_pairwise(dissimilarity_matrix_MS, Matrix_occurence_MS$Region_Curated_2, padj = "fdr")

tmpumap <- umap(decostand(tmp_MS, "pa"))
str(tmpumap)
plot(tmpumap$layout, col = as.numeric(as.factor(Matrix_occurence_MS$Region_Curated_2)), pch = 16)




#Extract coordinates to put them in a dataframe
Coord_NMDS_MS <- nmds_result_MS$points %>% as.data.frame() %>% cbind(Matrix_occurence_MS %>% select(Region_Curated_2))

#Now calculate ANOVA
summary(aov(MDS1 ~ MDS2 * Region_Curated_2, data = Coord_NMDS_MS))

# Create observation dataset at the MOTU lvl 1
Matrix_occurence_lvl1 <- ASV_all_Assignement_modif %>% 
  select(TP_lvl1, ASV_ID) %>% 
  left_join(ASV_Occurence) %>%
  group_by(TP_lvl1) %>% 
  summarise(across(2:2084, sum)) %>% 
  column_to_rownames(1) %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'Original Sample ID') %>% 
  left_join(Metadata_Tara_Lite %>% select(`Original Sample ID`, Curated_Geography_ID)) %>% 
  select(-`Original Sample ID`) %>% 
  group_by(Curated_Geography_ID) %>% 
  summarise_all(sum) %>% 
  inner_join(Metadata_Tara_Lite %>% select(Curated_Geography_ID, Region_Curated_2) %>% distinct()) %>% 
  select(-Curated_Geography_ID) %>% 
  rowwise() %>%
  mutate(sum_cols = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
  filter(sum_cols > 0) %>% select(-sum_cols)

tmp_TP_lvl1 <- Matrix_occurence_lvl1 %>% select(-Region_Curated_2) %>% 
  as.matrix() %>% 
  prop.table(margin = 1) %>% 
  as.data.frame()

#Now calculate NMDS at MS level
dissimilarity_matrix_lvl1 <- vegdist(decostand(tmp_TP_lvl1, "pa"), method = "jaccard", binary = T)
nmds_result_lvl1 <- metaMDS(dissimilarity_matrix_lvl1, k = 2,trymax=100)


tmpumap <- umap(decostand(tmp_TP_lvl1, "pa"))
str(tmpumap)
plot(tmpumap$layout, col = as.numeric(as.factor(Matrix_occurence_lvl1$Region_Curated_2)), pch = 16, type = "n")




#Extract coordinates to put them in a dataframe
Coord_NMDS_lvl1 <- nmds_result_lvl1$points %>% as.data.frame() %>% cbind(Matrix_occurence_lvl1 %>% select(Region_Curated_2))

#Now calculate ANOVA
summary(aov(MDS1 ~ MDS2 * Region_Curated_2, data = Coord_NMDS_lvl1))

# Create observation dataset at the MOTU lvl 2
Matrix_occurence_lvl2 <- ASV_all_Assignement_modif %>% 
  select(TP_lvl2, ASV_ID) %>% 
  left_join(ASV_Occurence) %>%
  group_by(TP_lvl2) %>% 
  summarise(across(2:2084, sum)) %>% 
  column_to_rownames(1) %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'Original Sample ID') %>% 
  left_join(Metadata_Tara_Lite %>% select(`Original Sample ID`, Curated_Geography_ID)) %>% 
  select(-`Original Sample ID`) %>% 
  group_by(Curated_Geography_ID) %>% 
  summarise_all(sum) %>% 
  inner_join(Metadata_Tara_Lite %>% select(Curated_Geography_ID, Region_Curated_2) %>% distinct()) %>% 
  select(-Curated_Geography_ID) %>% 
  rowwise() %>%
  mutate(sum_cols = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
  filter(sum_cols > 100) %>% select(-sum_cols)

tmp_TP_lvl2 <- Matrix_occurence_lvl1 %>% select(-Region_Curated_2) %>% 
  as.matrix() %>% 
  prop.table(margin = 1) %>% 
  as.data.frame()


#Now calculate NMDS at MS level
dissimilarity_matrix_lvl2 <- vegdist(tmp_TP_lvl2, method = "jaccard")
nmds_result_lvl2 <- metaMDS(dissimilarity_matrix_lvl2, k = 2,trymax=100)

#Extract coordinates to put them in a dataframe
Coord_NMDS_lvl2 <- nmds_result_lvl2$points %>% as.data.frame() %>% cbind(Matrix_occurence_lvl2 %>% select(Region_Curated_2))

#Now calculate ANOVA
summary(aov(MDS1 ~ MDS2 * Region_Curated_2, data = Coord_NMDS_lvl2))











# Now Look at the size fraction
# Create observation dataset at the Morphospecies level
Matrix_occurence_MS_Size <- ASV_all_Assignement_modif %>% 
  select(TP_MS, ASV_ID) %>% 
  left_join(ASV_Occurence) %>%
  group_by(TP_MS) %>% 
  summarise(across(2:2084, sum)) %>% 
  column_to_rownames(1) %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'Original Sample ID') %>% 
  left_join(Metadata_Tara_Lite %>% select(`Original Sample ID`, Curated_Geography_ID, Curated_size_fraction)) %>% 
  select(-`Original Sample ID`) %>% 
  group_by(Curated_Geography_ID, Curated_size_fraction) %>% 
  summarise_all(sum) %>%
  ungroup() %>% 
  distinct() %>% 
  select(-Curated_Geography_ID) %>% 
  rowwise() %>%
  mutate(sum_cols = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
  filter(sum_cols > 100) %>% select(-sum_cols)

tmp_MS_size <- Matrix_occurence_MS_Size %>% select(-Curated_size_fraction) %>% 
  as.matrix() %>% 
  prop.table(margin = 1) %>% 
  as.data.frame()

#Now calculate NMDS at MS level
dissimilarity_matrix_MS_size <- vegdist(tmp_MS_size, method = "jaccard")
nmds_result_MS_size <- metaMDS(dissimilarity_matrix_MS_size, k = 2,trymax=100)

#Extract coordinates to put them in a dataframe
Coord_NMDS_MS_size <- nmds_result_MS_size$points %>% as.data.frame() %>% cbind(Matrix_occurence_MS_Size %>% select(Curated_size_fraction))

#Now calculate ANOVA
summary(aov(MDS1 ~ MDS2 * Curated_size_fraction, data = Coord_NMDS_MS_size))

# Create observation dataset at the MOTU lvl1
Matrix_occurence_lvl1_Size <- ASV_all_Assignement_modif %>% 
  select(TP_lvl1, ASV_ID) %>% 
  left_join(ASV_Occurence) %>%
  group_by(TP_lvl1) %>% 
  summarise(across(2:2084, sum)) %>% 
  column_to_rownames(1) %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'Original Sample ID') %>% 
  left_join(Metadata_Tara_Lite %>% select(`Original Sample ID`, Curated_Geography_ID, Curated_size_fraction)) %>% 
  select(-`Original Sample ID`) %>% 
  group_by(Curated_Geography_ID, Curated_size_fraction) %>% 
  summarise_all(sum) %>%
  ungroup() %>% 
  distinct() %>% 
  select(-Curated_Geography_ID) %>% 
  rowwise() %>%
  mutate(sum_cols = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
  filter(sum_cols > 0) %>% select(-sum_cols)

tmp_lvl1_size <- Matrix_occurence_lvl1_Size %>% select(-Curated_size_fraction) %>% 
  as.matrix() %>% 
  prop.table(margin = 1) %>% 
  as.data.frame()

#Now calculate NMDS at MS level
dissimilarity_matrix_lvl1_size <- vegdist(decostand(tmp_lvl1_size, "pa"), method = "jaccard", binary = T)
nmds_result_lvl1_size <- metaMDS(dissimilarity_matrix_lvl1_size, k = 2, trymax=20)

# PERMANOVA
table(Matrix_occurence_lvl1_Size$Curated_size_fraction)
# potentially quite unbalanced, but leave as is for now
tmp_permanova <- adonis2(dissimilarity_matrix_lvl1_size ~ Matrix_occurence_lvl1_Size$Curated_size_fraction, sqrt.dist = T)
tmp_permanova_pw <- permanova_pairwise(dissimilarity_matrix_lvl1_size, Matrix_occurence_lvl1_Size$Curated_size_fraction, sqrt.dist = T, padj = "fdr")

# betadispersion
tmp_betadisper <- betadisper(dissimilarity_matrix_lvl1_size, Matrix_occurence_lvl1_Size$Curated_size_fraction, sqrt.dist = T)
plot(tmp_betadisper)
boxplot(tmp_betadisper)
barplot(cumsum(tmp_betadisper$eig/sum(tmp_betadisper$eig))[1:20])

# ANOSIM
tmp_anosim <- anosim(dissimilarity_matrix_lvl1_size, Matrix_occurence_lvl1_Size$Curated_size_fraction)
tmp_anosim_pw <- anosim_pairwise(dissimilarity_matrix_lvl1_size, Matrix_occurence_lvl1_Size$Curated_size_fraction, padj = "fdr")


tmp_pa <- decostand(Matrix_occurence_lvl1_Size[, -1], "pa")
tmp_check <- tmp_pa[rowSums(tmp_pa) == 1, ]
tmp_check <- tmp_check[, colSums(tmp_check) > 0]
colSums(tmp_pa[rowSums(tmp_pa) > 1, colnames(tmp_check)])


View(t(tmp_pa[tmp_pa[, "Non-spinose|Neogloboquadrina|pachyderma|IV"] != 0, ]))


require(umap)

?umap

tmpumap <- umap(decostand(tmp_lvl1_size, "pa"))
str(tmpumap)
plot(tmpumap$layout, col = as.numeric(as.factor(Matrix_occurence_lvl1_Size$Curated_size_fraction)))








#Extract coordinates to put them in a dataframe
Coord_NMDS_lvl1_size <- nmds_result_lvl1_size$points %>% as.data.frame() %>% cbind(Matrix_occurence_lvl1_Size %>% select(Curated_size_fraction))

#Now calculate ANOVA
summary(aov(MDS1 ~ MDS2 * Curated_size_fraction, data = Coord_NMDS_lvl1_size))


# Create observation dataset at the MOTU lvl2
Matrix_occurence_lvl2_Size <- ASV_all_Assignement_modif %>% 
  select(TP_lvl2, ASV_ID) %>% 
  left_join(ASV_Occurence) %>%
  group_by(TP_lvl2) %>% 
  summarise(across(2:2084, sum)) %>% 
  column_to_rownames(1) %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'Original Sample ID') %>% 
  left_join(Metadata_Tara_Lite %>% select(`Original Sample ID`, Curated_Geography_ID, Curated_size_fraction)) %>% 
  select(-`Original Sample ID`) %>% 
  group_by(Curated_Geography_ID, Curated_size_fraction) %>% 
  summarise_all(sum) %>%
  ungroup() %>% 
  distinct() %>% 
  select(-Curated_Geography_ID) %>% 
  rowwise() %>%
  mutate(sum_cols = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
  filter(sum_cols > 100) %>% select(-sum_cols)

tmp_lvl2_size <- Matrix_occurence_lvl2_Size %>% select(-Curated_size_fraction) %>% 
  as.matrix() %>% 
  prop.table(margin = 1) %>% 
  as.data.frame()

#Now calculate NMDS at MS level
dissimilarity_matrix_lvl2_size <- vegdist(tmp_lvl2_size, method = "jaccard")
nmds_result_lvl2_size <- metaMDS(dissimilarity_matrix_lvl2_size, k = 2,trymax=100)

#Extract coordinates to put them in a dataframe
Coord_NMDS_lvl2_size <- nmds_result_lvl2_size$points %>% as.data.frame() %>% cbind(Matrix_occurence_lvl2_Size %>% select(Curated_size_fraction))

#Now calculate ANOVA
summary(aov(MDS1 ~ MDS2 * Curated_size_fraction, data = Coord_NMDS_lvl2_size))




# Now Look at the Depth
# Create observation dataset at the Morphospecies level
Matrix_occurence_MS_Depth <- ASV_all_Assignement_modif %>% 
  select(TP_MS, ASV_ID) %>% 
  left_join(ASV_Occurence) %>%
  group_by(TP_MS) %>% 
  summarise(across(2:2084, sum)) %>% 
  column_to_rownames(1) %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'Original Sample ID') %>% 
  left_join(Metadata_Tara_Lite %>% select(`Original Sample ID`, Curated_Geography_ID, Curated_Depth)) %>% 
  select(-`Original Sample ID`) %>% 
  group_by(Curated_Geography_ID, Curated_Depth) %>% 
  summarise_all(sum) %>%
  ungroup() %>% 
  distinct() %>% 
  select(-Curated_Geography_ID) %>% 
  rowwise() %>%
  mutate(sum_cols = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
  filter(sum_cols > 100) %>% select(-sum_cols)

tmp_MS_Depth <- Matrix_occurence_MS_Depth %>% select(-Curated_Depth) %>% 
  as.matrix() %>% 
  prop.table(margin = 1) %>% 
  as.data.frame()

#Now calculate NMDS at MS level
dissimilarity_matrix_MS_Depth <- vegdist(tmp_MS_Depth, method = "jaccard")
nmds_result_MS_Depth <- metaMDS(dissimilarity_matrix_MS_Depth, k = 2,trymax=100)

#Extract coordinates to put them in a dataframe
Coord_NMDS_MS_Depth <- nmds_result_MS_Depth$points %>% as.data.frame() %>% cbind(Matrix_occurence_MS_Depth %>% select(Curated_Depth))

#Now calculate ANOVA
summary(aov(MDS1 ~ MDS2 * Curated_Depth, data = Coord_NMDS_MS_Depth))



# Create observation dataset at the MOTU lvl1
Matrix_occurence_lvl1_Depth <- ASV_all_Assignement_modif %>% 
  select(TP_lvl1, ASV_ID) %>% 
  left_join(ASV_Occurence) %>%
  group_by(TP_lvl1) %>% 
  summarise(across(2:2084, sum)) %>% 
  column_to_rownames(1) %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'Original Sample ID') %>% 
  left_join(Metadata_Tara_Lite %>% select(`Original Sample ID`, Curated_Geography_ID, Curated_Depth)) %>% 
  select(-`Original Sample ID`) %>% 
  group_by(Curated_Geography_ID, Curated_Depth) %>% 
  summarise_all(sum) %>%
  ungroup() %>% 
  distinct() %>% 
  select(-Curated_Geography_ID) %>% 
  rowwise() %>%
  mutate(sum_cols = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
  filter(sum_cols > 100) %>% select(-sum_cols)

tmp_lvl1_Depth <- Matrix_occurence_lvl1_Depth %>% select(-Curated_Depth) %>% 
  as.matrix() %>% 
  prop.table(margin = 1) %>% 
  as.data.frame()

#Now calculate NMDS at MOTU lvl1
dissimilarity_matrix_lvl1_Depth <- vegdist(tmp_lvl1_Depth, method = "jaccard")
nmds_result_lvl1_Depth <- metaMDS(dissimilarity_matrix_lvl1_Depth, k = 2,trymax=100)

#Extract coordinates to put them in a dataframe
Coord_NMDS_lvl1_Depth <- nmds_result_lvl1_Depth$points %>% as.data.frame() %>% cbind(Matrix_occurence_lvl1_Depth %>% select(Curated_Depth))

#Now calculate ANOVA
summary(aov(MDS1 ~ MDS2 * Curated_Depth, data = Coord_NMDS_lvl1_Depth))




# Create observation dataset at the MOTU lvl2
Matrix_occurence_lvl2_Depth <- ASV_all_Assignement_modif %>% 
  select(TP_lvl2, ASV_ID) %>% 
  left_join(ASV_Occurence) %>%
  group_by(TP_lvl2) %>% 
  summarise(across(2:2084, sum)) %>% 
  column_to_rownames(1) %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'Original Sample ID') %>% 
  left_join(Metadata_Tara_Lite %>% select(`Original Sample ID`, Curated_Geography_ID, Curated_Depth)) %>% 
  select(-`Original Sample ID`) %>% 
  group_by(Curated_Geography_ID, Curated_Depth) %>% 
  summarise_all(sum) %>%
  ungroup() %>% 
  distinct() %>% 
  select(-Curated_Geography_ID) %>% 
  rowwise() %>%
  mutate(sum_cols = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
  filter(sum_cols > 10) %>% select(-sum_cols)

tmp_lvl2_Depth <- Matrix_occurence_lvl2_Depth %>% select(-Curated_Depth) %>% 
  as.matrix() %>% 
  prop.table(margin = 1) %>% 
  as.data.frame()

#Now calculate NMDS at MOTU lvl2
dissimilarity_matrix_lvl2_Depth <- vegdist(tmp_lvl2_Depth, method = "jaccard")
nmds_result_lvl2_Depth <- metaMDS(dissimilarity_matrix_lvl2_Depth, k = 2,trymax=100)

#Extract coordinates to put them in a dataframe
Coord_NMDS_lvl2_Depth <- nmds_result_lvl2_Depth$points %>% as.data.frame() %>% cbind(Matrix_occurence_lvl2_Depth %>% select(Curated_Depth))

#Now calculate ANOVA
summary(aov(MDS1 ~ MDS2 * Curated_Depth, data = Coord_NMDS_lvl2_Depth))



#Now plot the result of all NMDS

NMDS_all_results <- (Coord_NMDS_MS %>% add_column(Dataset = "Morphospecies") %>% add_column(Variable = "Geography")      %>% dplyr::rename("Data" = "Region_Curated_2")) %>% 
  rbind(Coord_NMDS_lvl1            %>% add_column(Dataset = "MOTU_lvl1")     %>% add_column(Variable = "Geography")      %>% dplyr::rename("Data" = "Region_Curated_2")) %>%
  rbind(Coord_NMDS_lvl2            %>% add_column(Dataset = "MOTU_lvl2")     %>% add_column(Variable = "Geography")      %>% dplyr::rename("Data" = "Region_Curated_2")) %>% 
  rbind(Coord_NMDS_MS_size         %>% add_column(Dataset = "Morphospecies") %>% add_column(Variable = "Size-fraction")  %>% dplyr::rename("Data" = "Curated_size_fraction")) %>% 
  rbind(Coord_NMDS_lvl1_size       %>% add_column(Dataset = "MOTU_lvl1")     %>% add_column(Variable = "Size-fraction")  %>% dplyr::rename("Data" = "Curated_size_fraction")) %>%
  rbind(Coord_NMDS_lvl2_size       %>% add_column(Dataset = "MOTU_lvl2")     %>% add_column(Variable = "Size-fraction")  %>% dplyr::rename("Data" = "Curated_size_fraction")) %>%  
  rbind(Coord_NMDS_MS_Depth        %>% add_column(Dataset = "Morphospecies") %>% add_column(Variable = "Depth")          %>% dplyr::rename("Data" = "Curated_Depth")) %>% 
  rbind(Coord_NMDS_lvl1_Depth      %>% add_column(Dataset = "MOTU_lvl1")     %>% add_column(Variable = "Depth")          %>% dplyr::rename("Data" = "Curated_Depth")) %>%
  rbind(Coord_NMDS_lvl2_Depth      %>% add_column(Dataset = "MOTU_lvl2")     %>% add_column(Variable = "Depth")          %>% dplyr::rename("Data" = "Curated_Depth"))



ggplot(NMDS_all_results, aes(x = MDS1, y = MDS2, col = Data))+
  geom_point(size = 3)+
  facet_grid(Variable~Dataset)+
  coord_fixed()+
  scale_colour_viridis_d()


################################## Do a figure to show LDG at each level of the taxonomy #######################

# All_occurences_final -> File with every occurences in the dataset.
# Geography_curated    -> File with classification per latitudinal bin.
library(reshape2)

Geography_Curated     <- readxl::read_xlsx("../../02_Datasets/09_Geographic_metadata/09_Curated_Geographic_Data_20230224.xlsx", sheet = 1, col_types="guess")



Latitudinal_bin <- Geography_Curated %>% select(Curated_Geography_ID, Class_latitude,Region_Curated_2)

view(All_occurences_final)

Occurences_MS_Geography <- All_occurences_final %>% 
  select(TP_MS, Curated_Geography_ID) %>% 
  left_join(Latitudinal_bin %>% select(-Region_Curated_2)) %>% 
  select(-Curated_Geography_ID) %>% 
  distinct() %>%
  add_column(pres=1) %>% 
  spread(TP_MS, pres, fill=0) %>%
  filter(Class_latitude != "NA") %>% 
  column_to_rownames(1) 

Occurences_MS_Geography$sum <- rowSums(Occurences_MS_Geography)

Occurences_lvl1_Geography <- All_occurences_final %>% 
  select(TP_lvl1, Curated_Geography_ID) %>% filter(TP_lvl1 != "Not_Attributed") %>% 
  left_join(Latitudinal_bin %>% select(-Region_Curated_2)) %>% 
  select(-Curated_Geography_ID) %>% 
  distinct() %>%
  add_column(pres=1) %>% 
  spread(TP_lvl1, pres, fill=0) %>%
  filter(Class_latitude != "NA") %>% 
  column_to_rownames(1) 

Occurences_lvl1_Geography$sum <- rowSums(Occurences_lvl1_Geography)

Occurences_lvl2_Geography <- All_occurences_final %>% 
  select(TP_lvl2, Curated_Geography_ID) %>% filter(TP_lvl2 != "Not_Attributed") %>% 
  left_join(Latitudinal_bin%>% select(-Region_Curated_2)) %>% 
  select(-Curated_Geography_ID) %>% 
  distinct() %>%
  add_column(pres=1) %>% 
  spread(TP_lvl2, pres, fill=0) %>%
  filter(Class_latitude != "NA") %>% 
  column_to_rownames(1) 

Occurences_lvl2_Geography$sum <- rowSums(Occurences_lvl2_Geography)

Occurences_lvl3_Geography <- All_occurences_final %>% 
  select(TP_lvl3, Curated_Geography_ID) %>% filter(TP_lvl3 != "Not_Attributed") %>% 
  left_join(Latitudinal_bin%>% select(-Region_Curated_2)) %>% 
  select(-Curated_Geography_ID) %>% 
  distinct() %>%
  add_column(pres=1) %>% 
  spread(TP_lvl3, pres, fill=0) %>%
  filter(Class_latitude != "NA") %>% 
  column_to_rownames(1) 

Occurences_lvl3_Geography$sum <- rowSums(Occurences_lvl3_Geography)


Occurences_MS_Geography_summary   <- Occurences_MS_Geography     %>% select(sum) %>% dplyr::rename("N_MS" = "sum") %>% rownames_to_column(var = "Latitude_class")
Occurences_lvl1_Geography_summary <- Occurences_lvl1_Geography %>% select(sum) %>% dplyr::rename("N_lvl1" = "sum") %>% rownames_to_column(var = "Latitude_class")
Occurences_lvl2_Geography_summary <- Occurences_lvl2_Geography %>% select(sum) %>% dplyr::rename("N_lvl2" = "sum") %>% rownames_to_column(var = "Latitude_class")
Occurences_lvl3_Geography_summary <- Occurences_lvl3_Geography %>% select(sum) %>% dplyr::rename("N_lvl3" = "sum") %>% rownames_to_column(var = "Latitude_class")


LDG_all_data <- Occurences_MS_Geography_summary %>% 
  left_join(Occurences_lvl1_Geography_summary)  %>% 
  left_join(Occurences_lvl2_Geography_summary)  %>% 
  left_join(Occurences_lvl3_Geography_summary)  %>% 
  gather(Taxonomy, N_taxa,2:5) %>% 
  add_column("Biogeography" = "00_World")


ggplot(LDG_all_data, aes(x =Latitude_class, y = N_taxa, fill = Taxonomy))+
  geom_col(position="dodge")



Latitudinal_bin <- Geography_Curated %>% select(Curated_Geography_ID, Class_latitude,Region_Curated_2)


Occurences_MS_Geography_2 <- All_occurences_final %>% 
  select(TP_MS, Curated_Geography_ID) %>% 
  left_join(Latitudinal_bin) %>% 
  select(-Curated_Geography_ID) %>%
  unite("Biogeo", Class_latitude, Region_Curated_2, remove = TRUE,sep= "|") %>%
  distinct() %>%
  add_column(pres=1) %>% 
  spread(TP_MS, pres, fill=0) %>%
  filter(Biogeo != "NA|NA") %>% 
  column_to_rownames(1) 

Occurences_MS_Geography_2$sum <- rowSums(Occurences_MS_Geography_2)


Occurences_lvl1_Geography_2 <- All_occurences_final %>% 
  select(TP_lvl1, Curated_Geography_ID) %>% 
  left_join(Latitudinal_bin) %>% 
  select(-Curated_Geography_ID) %>%
  unite("Biogeo", Class_latitude, Region_Curated_2, remove = TRUE,sep= "|") %>%
  distinct() %>%
  add_column(pres=1) %>% 
  spread(TP_lvl1, pres, fill=0) %>%
  filter(Biogeo != "NA|NA") %>% 
  column_to_rownames(1) 

Occurences_lvl1_Geography_2$sum <- rowSums(Occurences_lvl1_Geography_2)



Occurences_lvl2_Geography_2 <- All_occurences_final %>% 
  select(TP_lvl2, Curated_Geography_ID) %>% 
  left_join(Latitudinal_bin) %>% 
  select(-Curated_Geography_ID) %>%
  unite("Biogeo", Class_latitude, Region_Curated_2, remove = TRUE,sep= "|") %>%
  distinct() %>%
  add_column(pres=1) %>% 
  spread(TP_lvl2, pres, fill=0) %>%
  filter(Biogeo != "NA|NA") %>% 
  column_to_rownames(1) 

Occurences_lvl2_Geography_2$sum <- rowSums(Occurences_lvl2_Geography_2)




Occurences_lvl3_Geography_2 <- All_occurences_final %>% 
  select(TP_lvl2, Curated_Geography_ID) %>% 
  left_join(Latitudinal_bin) %>% 
  select(-Curated_Geography_ID) %>%
  unite("Biogeo", Class_latitude, Region_Curated_2, remove = TRUE,sep= "|") %>%
  distinct() %>%
  add_column(pres=1) %>% 
  spread(TP_lvl2, pres, fill=0) %>%
  filter(Biogeo != "NA|NA") %>% 
  column_to_rownames(1) 

Occurences_lvl3_Geography_2$sum <- rowSums(Occurences_lvl3_Geography_2)





Occurences_MS_Geography_summary_2   <- Occurences_MS_Geography_2   %>% select(sum) %>% dplyr::rename("N_MS" = "sum") %>% rownames_to_column(var = "Latitude_class")
Occurences_lvl1_Geography_summary_2 <- Occurences_lvl1_Geography_2 %>% select(sum) %>% dplyr::rename("N_lvl1" = "sum") %>% rownames_to_column(var = "Latitude_class")
Occurences_lvl2_Geography_summary_2 <- Occurences_lvl2_Geography_2 %>% select(sum) %>% dplyr::rename("N_lvl2" = "sum") %>% rownames_to_column(var = "Latitude_class")
Occurences_lvl3_Geography_summary_2 <- Occurences_lvl3_Geography_2 %>% select(sum) %>% dplyr::rename("N_lvl3" = "sum") %>% rownames_to_column(var = "Latitude_class")



LDG_all_data_2 <- Occurences_MS_Geography_summary_2 %>% left_join(Occurences_lvl1_Geography_summary_2) %>% left_join(Occurences_lvl2_Geography_summary_2) %>% left_join(Occurences_lvl3_Geography_summary_2) %>% gather(Taxonomy, N_taxa,2:5) %>% 
  separate(col = Latitude_class, into = c("Latitude_class","Biogeography"), sep = "\\|")




LDG_all_data_final <- rbind(LDG_all_data,LDG_all_data_2)



ggplot(LDG_all_data_final, aes(x = N_taxa, y = Latitude_class, fill = Taxonomy))+
  geom_col(position="dodge") + 
  facet_wrap(Biogeography ~., nrow=1)  +
  scale_y_discrete(limits=rev)



view(diversity_cryptic_latitude)

diversity_cryptic_latitude<- LDG_all_data_final %>% 
   spread(Taxonomy, N_taxa, fill=0) %>% 
   mutate(P_lvl1 = N_lvl1/N_MS) %>% 
   mutate(P_lvl2 = N_lvl2/N_MS) %>%
   mutate(P_lvl3 = N_lvl3/N_MS) %>%
  gather("Taxonomic_Level", Prop_Taxa_per_MS, 7:9)

view(diversity_cryptic_latitude)


diversity_cryptic_latitude_cumulative <- LDG_all_data_final %>% 
  spread(Taxonomy, N_taxa, fill=0) %>% 
  mutate('00_MS'   = N_MS) %>% 
  mutate('01_lvl1' = N_lvl1-N_MS) %>% 
  mutate('02_lvl2' = N_lvl2-N_lvl1) %>% 
  mutate('03_lvl3' = N_lvl3-N_lvl2) %>%
  gather("Taxonomic_Level", N_Taxa, 7:10) %>% 
  filter(N_Taxa >0)
  

graph1<-ggplot(diversity_cryptic_latitude_cumulative %>% filter(Biogeography == "00_World") %>% filter(Taxonomic_Level != "03_lvl3"), aes(x = N_Taxa, y = Latitude_class, fill = Taxonomic_Level))+
  geom_col(position = position_stack(reverse = TRUE),width = 0.5) + 
  scale_y_discrete(limits=rev)+
  scale_fill_viridis_d()+
  theme_minimal()


graph2<-ggplot(diversity_cryptic_latitude%>% filter(Biogeography == "00_World") %>% filter(Taxonomic_Level != "P_lvl3"), aes(x = Prop_Taxa_per_MS, y = Latitude_class, colour = Taxonomic_Level))+
  geom_point(size=6) +
  scale_y_discrete(limits=rev)+
  scale_colour_viridis_d()+
  theme_minimal()+
  xlim(1,2)



library(cowplot)

plot_grid(graph1,graph2)


################################ Try upset plot

 library(UpSetR)
# Prepare the data as Presence absence table 
# First format data as presence/absence data
 
Data_Upset_Geography_MS <- All_occurences_final %>% 
   select(TP_MS, Curated_Geography_ID) %>% filter(TP_MS != "Not_Attributed") %>% 
   left_join(Geography_Curated %>% select(Curated_Geography_ID, Region_Curated_4) %>% filter(Curated_Geography_ID != "NA")) %>% select(-Curated_Geography_ID) %>%
   distinct() %>% 
   mutate(pres=1) %>% 
   spread(Region_Curated_2, pres, fill=0) %>%
   column_to_rownames(1) %>% 
   select(-"<NA>")

upset(Data_Upset_Geography_MS,
      nsets = 6, 
      number.angles = 30, 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "N Morphospecies", 
      sets.x.label = "N Taxa per Oceanic basin", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 1), 
      order.by = c("freq"))

Data_Upset_Geography_lvl1 <- All_occurences_final %>% 
  select(TP_lvl1, Curated_Geography_ID) %>% filter(TP_lvl1 != "Not_Attributed") %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Region_Curated_2) %>% filter(Curated_Geography_ID != "NA")) %>% select(-Curated_Geography_ID) %>%
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Region_Curated_2, pres, fill=0) %>%
  column_to_rownames(1) %>% 
  select(-"<NA>")

upset(Data_Upset_Geography_lvl1,
      nsets = 6, 
      number.angles = 30, 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "N MOTUs lvl1", 
      sets.x.label = "N Taxa per Oceanic basin", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 1), 
      order.by = c("freq"))

Data_Upset_Geography_lvl2 <- All_occurences_final %>% 
  select(TP_lvl2, Curated_Geography_ID) %>% filter(TP_lvl2 != "Not_Attributed") %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Region_Curated_4) %>% filter(Curated_Geography_ID != "NA")) %>% select(-Curated_Geography_ID) %>%
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Region_Curated_4, pres, fill=0) %>%
  column_to_rownames("TP_lvl2") %>% 
  select(-"<NA>")

upset(Data_Upset_Geography_lvl2,
      nsets = 10, 
      number.angles = 15, 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "N MOTUs lvl2", 
      sets.x.label = "N Taxa per Oceanic basin", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 1), 
      order.by = c("freq"))

Data_Upset_Geography_lvl3 <- All_occurences_final %>% 
  select(TP_lvl3, Curated_Geography_ID) %>% filter(TP_lvl3 != "Not_Attributed") %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Region_Curated_2) %>% filter(Curated_Geography_ID != "NA")) %>% select(-Curated_Geography_ID) %>%
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Region_Curated_2, pres, fill=0) %>%
  column_to_rownames(1) %>% 
  select(-"<NA>")

upset(Data_Upset_Geography_lvl3,
      nsets = 6, 
      number.angles = 30, 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "N MOTUs lvl3", 
      sets.x.label = "N Taxa per Oceanic basin", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 1), 
      order.by = c("freq"))



library(ggpubr)

ggarrange(UpS_BioGeography_MS,
          UpS_BioGeography_lvl1,
          UpS_BioGeography_lvl2,
          UpS_BioGeography_lvl3)


UpS_BioGeography_MS



#Now Size fraction
Data_Upset_SF_MS <- All_occurences_final %>% 
  select(TP_MS, Curated_Geography_ID) %>% filter(TP_MS != "Not_Attributed") %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Curated_size_fraction) %>% filter(Curated_Geography_ID != "NA")) %>% select(-Curated_Geography_ID) %>%
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Curated_size_fraction, pres, fill=0) %>%
  column_to_rownames(1) %>% 
  select(-"<NA>")

upset(Data_Upset_SF_MS,
      nsets = 6, 
      number.angles = 30, 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "N taxa", 
      sets.x.label = "N Taxa per Oceanic basin", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 1), 
      order.by = c("freq"))

Data_Upset_SF_lvl1 <- All_occurences_final %>% 
  select(TP_lvl1, Curated_Geography_ID) %>% filter(TP_lvl1 != "Not_Attributed") %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Curated_size_fraction) %>% filter(Curated_Geography_ID != "NA")) %>% select(-Curated_Geography_ID) %>%
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Curated_size_fraction, pres, fill=0) %>%
  column_to_rownames(1) %>% 
  select(-"<NA>")

upset(Data_Upset_SF_lvl1,
      nsets = 6, 
      number.angles = 30, 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "N taxa", 
      sets.x.label = "N Taxa per Oceanic basin", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 1), 
      order.by = c("freq"))

Data_Upset_SF_lvl2 <- All_occurences_final %>% 
  select(TP_lvl2, Curated_Geography_ID) %>% filter(TP_lvl2 != "Not_Attributed") %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Curated_size_fraction) %>% filter(Curated_Geography_ID != "NA")) %>% select(-Curated_Geography_ID) %>%
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Curated_size_fraction, pres, fill=0) %>%
  column_to_rownames(1) %>% 
  select(-"<NA>")

upset(Data_Upset_SF_lvl2,
      nsets = 10, 
      number.angles = 30, 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "N taxa", 
      sets.x.label = "N Taxa per Oceanic basin", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 1), 
      order.by = c("freq"))

Data_Upset_SF_lvl3 <- All_occurences_final %>% 
  select(TP_lvl3, Curated_Geography_ID) %>% filter(TP_lvl3 != "Not_Attributed") %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Curated_size_fraction) %>% filter(Curated_Geography_ID != "NA")) %>% select(-Curated_Geography_ID) %>%
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Curated_size_fraction, pres, fill=0) %>%
  column_to_rownames(1) %>% 
  select(-"<NA>")

upset(Data_Upset_SF_lvl3,
      nsets = 6, 
      number.angles = 30, 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "N taxa", 
      sets.x.label = "N Taxa per Oceanic basin", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 1), 
      order.by = c("freq"))

# And now Depth
#Now Size fraction
Data_Upset_Depth_MS <- All_occurences_final %>% 
  select(TP_MS, Curated_Geography_ID) %>% filter(TP_MS != "Not_Attributed") %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Curated_Depth) %>% filter(Curated_Geography_ID != "NA")) %>% select(-Curated_Geography_ID) %>%
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Curated_Depth, pres, fill=0) %>%
  column_to_rownames(1) %>% 
  select(-"<NA>")

upset(Data_Upset_Depth_MS,
      nsets = 6, 
      number.angles = 30, 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "N taxa", 
      sets.x.label = "N Taxa per Oceanic basin", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 1), 
      order.by = c("freq"))

Data_Upset_Depth_lvl1 <- All_occurences_final %>% 
  select(TP_lvl1, Curated_Geography_ID) %>% filter(TP_lvl1 != "Not_Attributed") %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Curated_Depth) %>% filter(Curated_Geography_ID != "NA")) %>% select(-Curated_Geography_ID) %>%
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Curated_Depth, pres, fill=0) %>%
  column_to_rownames(1) %>% 
  select(-"<NA>")

upset(Data_Upset_Depth_lvl1,
      nsets = 6, 
      number.angles = 30, 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "N taxa", 
      sets.x.label = "N Taxa per Oceanic basin", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 1), 
      order.by = c("freq"))

Data_Upset_Depth_lvl2 <- All_occurences_final %>% 
  select(TP_lvl2, Curated_Geography_ID) %>% filter(TP_lvl2 != "Not_Attributed") %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Curated_Depth) %>% filter(Curated_Geography_ID != "NA")) %>% select(-Curated_Geography_ID) %>%
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Curated_Depth, pres, fill=0) %>%
  column_to_rownames(1) %>% 
  select(-"<NA>")

upset(Data_Upset_Depth_lvl2,
      nsets = 6, 
      number.angles = 30, 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "N taxa", 
      sets.x.label = "N Taxa per Oceanic basin", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 1), 
      order.by = c("freq"))

Data_Upset_Depth_lvl3 <- All_occurences_final %>% 
  select(TP_lvl3, Curated_Geography_ID) %>% filter(TP_lvl3 != "Not_Attributed") %>% 
  left_join(Geography_Curated %>% select(Curated_Geography_ID, Curated_Depth) %>% filter(Curated_Geography_ID != "NA")) %>% select(-Curated_Geography_ID) %>%
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(Curated_Depth, pres, fill=0) %>%
  column_to_rownames(1) %>% 
  select(-"<NA>")

upset(Data_Upset_Depth_lvl3,
      nsets = 6, 
      number.angles = 30, 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "N taxa", 
      sets.x.label = "N Taxa per Oceanic basin", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 1), 
      order.by = c("freq"))














library(UpSetR)
library(tidyverse)





















 
 
 
 
 
 Matrix_occurence_MS <- ASV_all_Assignement_modif %>% 
   select(TP_MS, ASV_ID) %>% 
   left_join(ASV_Occurence) %>%
   group_by(TP_MS) %>% 
   summarise(across(2:2084, sum)) %>% 
   column_to_rownames(1) %>% 
   t() %>%
   as.data.frame() %>% 
   rownames_to_column(var = 'Original Sample ID') %>% 
   left_join(Metadata_Tara_Lite %>% select(`Original Sample ID`, Curated_Geography_ID)) %>% 
   select(-`Original Sample ID`) %>% 
   group_by(Curated_Geography_ID) %>% 
   summarise_all(sum) %>% 
   inner_join(Metadata_Tara_Lite %>% select(Curated_Geography_ID, Region_Curated_2) %>% distinct()) %>% 
   select(-Curated_Geography_ID) %>% 
   rowwise() %>%
   mutate(sum_cols = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
   filter(sum_cols > 0) %>% select(-sum_cols)
 
 tmp_MS <- Matrix_occurence_MS %>% select(-Region_Curated_2) %>% 
   as.matrix() %>% 
   prop.table(margin = 1) %>% 
   as.data.frame()
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 Upset_Geography_lvl2 <- Matrix_occurence_lvl2 %>%
   group_by(Region_Curated_2) %>% 
   summarise_all(sum) %>%
   ungroup() %>% 
   column_to_rownames(loc = 1) %>% 
   t() %>% 
   as.data.frame()
   
 Upset_Geography_lvl2[Upset_Geography_lvl2>0] <- 1
 
 
 Upset_Geography_lvl1 <- Matrix_occurence_lvl1 %>%
   group_by(Region_Curated_2) %>% 
   summarise_all(sum) %>%
   ungroup() %>% 
   column_to_rownames(loc = 1) %>% 
   t() %>% 
   as.data.frame()
 
 Upset_Geography_lvl1[Upset_Geography_lvl1>0] <- 1
 
 Upset_Geography_MS <- Matrix_occurence_MS %>%
   group_by(Region_Curated_2) %>% 
   summarise_all(sum) %>%
   ungroup() %>% 
   column_to_rownames(loc = 1) %>% 
   t() %>% 
   as.data.frame()
 
 Upset_Geography_MS[Upset_Geography_MS>0] <- 1
 
 
 upset(Upset_Geography_MS,
       nsets = 5, 
       number.angles = 30, 
       point.size = 3.5, 
       line.size = 2, 
       mainbar.y.label = "N taxa", 
       sets.x.label = "N Taxa per Oceanic basin", 
       text.scale = c(1.3, 1.3, 1, 1, 2, 0.75), 
       order.by = c("freq"))
 
 
 
 upset(Upset_Geography_lvl1,
       nsets = 5, 
       number.angles = 30, 
       point.size = 3.5, 
       line.size = 2, 
       mainbar.y.label = "N taxa", 
       sets.x.label = "N Taxa per Oceanic basin", 
       text.scale = c(1.3, 1.3, 1, 1, 2, 0.75), 
       order.by = c("freq"))
 
 
 upset(Upset_Geography_lvl2,
       nsets = 5, 
       number.angles = 30, 
       point.size = 3.5, 
       line.size = 2, 
       mainbar.y.label = "N taxa", 
       sets.x.label = "N Taxa per Oceanic basin", 
       text.scale = c(1.3, 1.3, 1, 1, 2, 0.75), 
       order.by = c("freq"))
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 view(ASV_all_Assignement)
 view(ASV_Occurence)
 view(Geography_Curated)
 
 # Prepare the file with the relevant variables to test
 Metadata_Tara_Lite <- Geography_Curated %>% 
   select(`Original Sample ID`, Region_Curated_2, Curated_Geography_ID, Curated_size_fraction, Curated_Depth)
 
 ASV_all_Assignement_modif <- ASV_all_Assignement %>% 
   select(-TP_lvl1, -TP_lvl2, -TP_lvl3) %>% 
   unite("TP_lvl1", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, remove = FALSE,sep= "|") %>%
   unite("TP_lvl2", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, MOTU_lvl_2, remove = FALSE,sep= "|")
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 #######
 # ASV_all_assignement -> Taxonomy of all AVS
 # ASV_Occurence -> Matrix with occurence of all ASV in all samples
 # Biogeography_Cordier -> All environmental samples and their associated metadata and data
 # Geography_curated -> ALL env data
 
 view(ASV_all_Assignement)
 view(ASV_Occurence)
 view(Geography_Curated)
 
 # Prepare the file with the relevant variables to test
 Metadata_Tara_Lite <- Geography_Curated %>% 
   select(`Original Sample ID`, Region_Curated_2, Curated_Geography_ID, Curated_size_fraction, Curated_Depth)
 
 ASV_all_Assignement_modif <- ASV_all_Assignement %>% 
   select(-TP_lvl1, -TP_lvl2, -TP_lvl3) %>% 
   unite("TP_lvl1", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, remove = FALSE,sep= "|") %>%
   unite("TP_lvl2", Clade_curated_assignement,Genus_curated_assignement,species_curated_assignement, MOTU_lvl_1, MOTU_lvl_2, remove = FALSE,sep= "|")
 
 ################################################################################################
 
 # Look at the biogeography specifically
 # Create observation dataset at the Morphospecies level
 Matrix_occurence_MS <- ASV_all_Assignement_modif %>% 
   select(TP_MS, ASV_ID) %>% 
   left_join(ASV_Occurence) %>%
   group_by(TP_MS) %>% 
   summarise(across(2:2084, sum)) %>% 
   column_to_rownames(1) %>% 
   t() %>%
   as.data.frame() %>% 
   rownames_to_column(var = 'Original Sample ID') %>% 
   left_join(Metadata_Tara_Lite %>% select(`Original Sample ID`, Curated_Geography_ID)) %>% 
   select(-`Original Sample ID`) %>% 
   group_by(Curated_Geography_ID) %>% 
   summarise_all(sum) %>% 
   inner_join(Metadata_Tara_Lite %>% select(Curated_Geography_ID, Region_Curated_2) %>% distinct()) %>% 
   select(-Curated_Geography_ID) %>% 
   rowwise() %>%
   mutate(sum_cols = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
   filter(sum_cols > 0) %>% select(-sum_cols)
 
 tmp_MS <- Matrix_occurence_MS %>% select(-Region_Curated_2) %>% 
   as.matrix() %>% 
   prop.table(margin = 1) %>% 
   as.data.frame()
 
 
 ###################################################### ANOSIM Script
 
 #setwd("C:/Users/chassenrueck/Documents/Additional_projects/RaphaelMorard/Molecular_Taxonomy_Raphael/")
 # save.image("Env_pattern_motus.Rdata")
 #load("Env_pattern_motus.Rdata")
 
 dat <- read.table(
   "Data_for_the_magician.txt",
   sep = " ",
   quote = "\""
 )
 
require(tidyverse)
require(vegan)
require(reshape)
require(caret)
require(geometry)
 
 
install.packages("geometry")
 # remove NAs from data set (why are there any?)
 dat <- dat %>% filter(!is.na(Curated_Geography_ID), !is.na(PCA_1_Monthly), !is.na(PCA_2_Monthly))

 # format contextual data (PC values)
 env <- unique(dat[, c("Curated_Geography_ID", "PCA_1_Monthly", "PCA_2_Monthly")]) %>% 
   remove_rownames() %>% 
   column_to_rownames("Curated_Geography_ID")
 nrow(env)
 
 ### summary:
 # The below approaches are mutually exclusive in their justification of suitability or because of redundancy
 # Option 1: variation in motu composition explained by env --> straightforward, no duplication of env data to match MOTU occurrences, no need to remove motus, one value per morphospecies (R2)
 # Option 2: prediction of motu based on env --> biased by highly unbalanced motu observations, exclusion of rare motus required, one value per motu to quantify predictive success (i.e. env dependency)
 # Option 3: not run due to too many violations of assumptions
 # Option 4: analysis of similarity between motus based on their env preferences --> data duplication of env values may be an issue, this approach is somewhat backward to Option 1, exclusion of rare motus required (>=3), either one value for strength of separation in general or pairwise comparisons between motus
 # Option 5: percentage overlap in 2D of env preferences --> exclusion of rare motus required (>=3), pairwise comparisons between motus, crude approach based solely on convex hulls of point clouds in 2D
 
 
 ### Option 1: test if jaccard dissimilarity of cryptic specied per morphospecies is influenced by PC values
 # strategy: 
 #   cast presence/absence community matrix per morphospecies with curated ID as sample identifer
 #   run permanova with pa matrix as response
 
 # testing
 x <- dat %>% filter(Morphotaxa == "Candeina nitida")
 
 # full loop
 permanova_out <- dat %>% 
   group_by(Morphotaxa) %>% 
   group_map(
     function(x, ...) {
       pa_mat <- cast(x, "MOTU ~ Curated_Geography_ID", value = "n_observation", fun.aggregate = length) %>% 
         remove_rownames() %>% 
         column_to_rownames("MOTU") %>% 
         t() %>% 
         as.data.frame() %>% 
         filter(rownames(.) %in% rownames(env))
       # weirdly not all entries are always 1
       # that means that for the same morphospecies and motu there are more than 1 row per Curated_Geography_ID in the input table
       pa_mat <- decostand(pa_mat, "pa")
       out <- adonis2(pa_mat ~ ., data = env[rownames(pa_mat), ], method = "jaccard", sqrt.dist = T)
       # should I use sequential or marginal SS?
       # theoretically, PCs should not be colinear by definition, but removing samples compared to the full data set may introduce some again
       # if you argue that PC1 is the most important pattern and all others should just be viewed additional to that,
       # sequential SS should be ok
       return(out)
     }
   )
 names(permanova_out) <- unique(sort(dat$Morphotaxa))
 
 # format output
 permanova_out_df <- map_dfr(
   permanova_out,
   function(x) {
     x$R2[1:3]
   }
 ) %>% t()
 colnames(permanova_out_df) <- c("R2_PC1", "R2_PC2", "R2_resid")
 write.table(permanova_out_df, "env_effect_permanova_out_df.txt", quote = F, sep = "\t")

 
 ### Option 2: use PC1 and 2 values for each occurrence of MOTU per morphospecies as response, 
 #             to test grouping effect of MOTU factor
 # this analysis would be closer to the figure
 # potential issues: PC1+2 values are being duplicated if at one location more than 1 MOTU is found

 # checking nobs per MOTU
 dat %>% 
   group_by(Morphotaxa) %>% 
   group_map(
     function(x, ...) {
       table(x$MOTU)
     }
   )
 # the detection of MOTUs is too unbalanced
 # model outcome is not reliable

 # run case-by-case
 x <- dat %>% filter(Morphotaxa == "Turborotalita clarkei")
 rf_cv <- train( 
   x[, c("PCA_1_Monthly", "PCA_2_Monthly")], 
   y = x$MOTU, 
   method = "rf", 
   trControl = trainControl(method = 'LOOCV')
 )
 confusionMatrix(rf_cv$pred$pred, rf_cv$pred$obs, mode = "everything", positive = NULL)
 
 # remove rare occurrences (needs to be detected at least 4 times)
 dat_sub <- dat %>% 
   mutate(morpho_motu = paste(Morphotaxa, MOTU)) %>% 
   filter(morpho_motu %in% names(which(table(.$morpho_motu) > 10))) %>% 
   filter(Morphotaxa %in% names(which(sapply(by(.$MOTU, .$Morphotaxa, unique), length) > 1)))
 dat_sub %>% 
   group_by(Morphotaxa) %>% 
   group_map(
     function(x, ...) {
       table(x$MOTU)
     }
   )

 # run full loop
 rf_out <- dat_sub %>% 
   group_by(Morphotaxa) %>% 
   group_map(
     function(x, ...) {
       train( 
         x[, c("PCA_1_Monthly", "PCA_2_Monthly")], 
         y = x$MOTU, 
         method = "rf", 
         trControl = trainControl(method = 'LOOCV')
       )
     }
   )
 names(rf_out) <- unique(sort(dat_sub$Morphotaxa))
 cm_out <- map(
   rf_out,
   function(x) {
     if(length(levels(x$pred$obs)) == 2) {
       sapply(
         levels(x$pred$obs),
         function(y) {
           caret::confusionMatrix(x$pred$pred, x$pred$obs, mode = "everything", positive = y)$byClass
         }
       ) %>% t()
     } else {
       cm_tmp <- caret::confusionMatrix(x$pred$pred, x$pred$obs, mode = "everything", positive = NULL)$byClass
       rownames(cm_tmp) <- gsub("Class: ", "", rownames(cm_tmp))
       return(cm_tmp)
     }
   }
 )
 cm_out_df <- (
   data.frame(
     Morphotaxa = rep(names(cm_out), sapply(cm_out, nrow)),
     MOTU = unlist(sapply(cm_out, rownames)),
     Nobs = unlist(sapply(rf_out, function(x) table(x$pred$obs))),
     do.call("rbind", cm_out)
   )
 )
 rownames(cm_out_df) <- NULL
 cm_out_df$F1[is.na(cm_out_df$F1)] <- 0
 write.table(cm_out_df, "env_effect_rf_out_df.txt", quote = F, sep = "\t")
 
 
 ### Option3: MANOVA or similar
 # I am reluctant to use any method with the PCs as response variable
 # those values are already very derived and I am not sure, if this will not interfere with the analysis
 # e.g. if I used RDA with the PCs as response it would calculate another ordination on top of the existing one
 # also I would avoid methods that require normality...
 # I am not pursuing this further at this stage
 
 
 ### Option4: ANOSIM 
 # use PC scores as variables and motu occurrences as samples (observations)
 # this is a bit backwards compared to the PERMANOVA approach (option 1)
 source("C:/Users/chassenrueck/Documents/Repos/ARISA/anosimPosthoc.R")
 
 # testing
 x <- dat %>% 
   filter(Morphotaxa == "Globigerina bulloides") %>% 
   select(MOTU, PCA_1_Monthly, PCA_2_Monthly) %>% 
   filter(MOTU %in% names(which(table(.$MOTU) > 1)))
 
 # also remove motus with less than N (N >= 2) observations
 dat_sub <- dat %>% 
   mutate(morpho_motu = paste(Morphotaxa, MOTU)) %>% 
   filter(morpho_motu %in% names(which(table(.$morpho_motu) > 2))) %>% 
   filter(Morphotaxa %in% names(which(sapply(by(.$MOTU, .$Morphotaxa, unique), length) > 1)))
 
 anosim_global <- dat_sub %>% 
   group_by(Morphotaxa) %>% 
   group_map(
     function(x, ...) {
       anosim(x[, c("PCA_1_Monthly", "PCA_2_Monthly")], x$MOTU, distance = "euclidean")
     }
   )
 names(anosim_global) <- unique(sort(dat_sub$Morphotaxa))
 
 anosim_pw <- dat_sub %>% 
   group_by(Morphotaxa) %>% 
   group_map(
     function(x, ...) {
       ANOSIMposthoc(x[, c("PCA_1_Monthly", "PCA_2_Monthly")], as.factor(x$MOTU), distance = "euclidean")
     }
   )
 names(anosim_pw) <- unique(sort(dat_sub$Morphotaxa))
 
 
 ### Option5: percentage overlap in 2D space
 # this is maybe a very crude approach based purely on the data visualization
 # although, maybe simpler is better?
 # as in option 4 I had to remove singleton and doubleton occurrences

 x <- dat %>% 
   filter(Morphotaxa == "Neogloboquadrina pachyderma") %>% 
   filter(MOTU %in% names(which(table(.$MOTU) > 2)))
 pw_motu <- data.frame(t(combn(unique(sort(x$MOTU)), 2)))
 pw_motu$overlap <- map_dbl(
   1:nrow(pw_motu),
   function(y) {
     tmp <- intersectn(
       x[x$MOTU == pw_motu[y, 1], c("PCA_1_Monthly", "PCA_2_Monthly")],
       x[x$MOTU == pw_motu[y, 2], c("PCA_1_Monthly", "PCA_2_Monthly")]
     )
     if(is.null(tmp$ps)) {
       0
     } else {
       tmp$ch$area/min(c(tmp$ch1$area, tmp$ch2$area)) * 100
     }
   }
 )
 
 

 
 
 
 
 
 
 
 
 
 
 
 #####################################################################################
 
 Data_formating_for_Siccha <- All_occurences_final %>% 
   select(Curated_Geography_ID, n_observation, TP_lvl2) %>% 
   group_by(Curated_Geography_ID,TP_lvl2) %>%
   summarise_each(funs(sum)) %>% 
   filter(TP_lvl2 != "Not_Attributed") %>% 
   spread(TP_lvl2, n_observation, fill=0)
   
   
write.table(Data_formating_for_Siccha, "Occurence_MOTUs_lvl2.txt")
 
