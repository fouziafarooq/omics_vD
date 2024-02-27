library(BiocManager)
library(omePath) # chemical to pathways
library(readxl)
library(dplyr)

# need to ask:
# how to normalize data - recs?
# what to do when there are 2 HMDB numbers.  Do I keep the Chem_ID also b/c those are unique. 
# next steps? 
meta_data <- readxl::read_excel("~/box/1.GW/Career MODE/OMICS TANZANIA PILOT/vD Metabalomics/R/omics_vD/data/CHHB-01-21VW+ SERUM DATA TABLES.XLSX",
                                sheet = "Sample Meta Data")

chem_anno <- read_excel("~/box/1.GW/Career MODE/OMICS TANZANIA PILOT/vD Metabalomics/R/omics_vD/data/CHHB-01-21VW+ SERUM DATA TABLES.XLSX",
                        sheet = "Chemical Annotation")                                

peak_area <- read_excel("~/box/1.GW/Career MODE/OMICS TANZANIA PILOT/vD Metabalomics/R/omics_vD/data/CHHB-01-21VW+ SERUM DATA TABLES.XLSX",
                        sheet = "Peak Area Data")

sort(names(chem_anno))
chem_anno <- chem_anno %>%
  rowwise() %>%
  mutate(HMDB_id = strsplit(HMDB, ",")[[1]][1]) 

table(chem_anno$HMDB_id)

chem_ids <- names(peak_area)[-1] # -1 is don't take the first column. 
# chem_id_pairs <- combn(chem_ids, 2, simplify = FALSE) # don't need this. dont need pairs:-(((((((

group_membership <- meta_data %>% select(PARENT_SAMPLE_NAME, TREATMENT)
peak_full <- left_join(peak_area, group_membership)

results_df <- data.frame()

for(chem_id in chem_ids){
  x<- peak_full[, chem_id] %>% unlist()
  na_tab <- table(is.na(x), peak_full$TREATMENT) # b/c there are NAs in the data.
  if(na_tab[1,1]<2 | na_tab[1,2]<2){
    print(paste("Skipping b/c missing chem_id data in chem id:", chem_id))
    next
  }
  t_result <- t.test(x ~ peak_full$TREATMENT) # value of chem_id like 35 as a string. 
  results_row <- list(chem_id = chem_id, p.value = t_result$p.value)
  results_df <- rbind(results_df, results_row)
}

results_df <- left_join(results_df, 
                        chem_anno %>% select(chem_id = CHEM_ID, CHEMICAL_NAME) %>% 
                          mutate(chem_id = as.character(chem_id)))

