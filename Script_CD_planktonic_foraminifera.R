# Script written by Raphael Morard, Christiane Hassenrück and Chiara Vanni
# Loading packages 
# Data tidying
library(tidyverse)
library(dplyr)
library(textshape)
library("janitor")

# Networks
library(igraph)
library(reshape)
library(CINNA)

# Bioinformatics
library("dada2")
library("seqRFLP")

# To produce maps
library(maps)
library(parallel)
library(hexbin)
library(scales)

# Loading the datasets
Data_sanger                         <- readxl::read_xlsx(path = "Table S1.xlsx", sheet = 1,skip = 1, col_types="guess") %>% arrange(Internal_ID)
RFLP_data                           <- readxl::read_xlsx(path = "Table S2.xlsx", sheet = 1, col_types="guess") 
Geography_Curated                   <- readxl::read_xlsx(path = "Table S4.xlsx", sheet = 1, col_types="guess")
Molecular_Nomenclature              <- readxl::read_xlsx(path = "Table S5.xlsx", sheet = 1, col_types="guess")
Equivalence_Nomenclature            <- readxl::read_xlsx(path = "Table S5.xlsx", sheet = 3, col_types="guess")
ASV_Occurence                       <- readxl::read_xlsx(path = "Table S6.xlsx", sheet = 1, col_types="guess")
PF_ASV_seq                          <- readxl::read_xlsx(path = "Table S6.xlsx", sheet = 2, col_types="guess")
Molecular_Nomenclature_Novel_ASVs   <- readxl::read_xlsx(path = "Table S6.xlsx", sheet = 3, col_types="guess")

# importing dataset 
# The order of the row is sorted using the unique internal provided for each sequence to make sure that the table is processed in the same way at each iteration
# Assess the quality of each individual region to select Basetypes
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
Check_co_occurence_lvl2 <- DB_all_assignement %>% select(Voucher_code, TP_lvl2) %>% filter(TP_lvl2 != "Not_Attributed") %>% 
  distinct() %>% 
  mutate(pres=1) %>% 
  spread(TP_lvl2,pres,fill=0) %>% 
  column_to_rownames(1) %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum != 1)

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
