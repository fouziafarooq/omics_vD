library(BiocManager)
library(omePath) # chemical to pathways
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# need to ask:
# how to normalize data - recs?
# what to do when there are 2 HMDB numbers.  Do I keep the Chem_ID also b/c those are unique. 
# next steps? 
meta_data <- readxl::read_excel("~/box/1.GW/Career MODE/OMICS TANZANIA PILOT/vD Metabalomics/R/omics_vD/data/CHHB-01-21VW+ SERUM DATA TABLES.XLSX",
                                sheet = "Sample Meta Data")

chem_anno <- read_excel("~/box/1.GW/Career MODE/OMICS TANZANIA PILOT/vD Metabalomics/R/omics_vD/data/CHHB-01-21VW+ SERUM DATA TABLES.XLSX",
                        sheet = "Chemical Annotation")                                

peak_area <- read_excel("~/box/1.GW/Career MODE/OMICS TANZANIA PILOT/vD Metabalomics/R/omics_vD/data/CHHB-01-21VW+ SERUM DATA TABLES.XLSX",
                        sheet = "Batch-normalized Data")

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
  x<- peak_full[, chem_id] %>% unlist() # grab one column/chemical of data. 
  chem_df <- data.frame(feature_level = x, group = peak_full$TREATMENT)
  na_tab <- table(is.na(x), peak_full$TREATMENT) # b/c there are NAs in the data.
  if(na_tab[1,1]<2 | na_tab[1,2]<2){
    print(paste("Skipping b/c missing chem_id data in chem id:", chem_id))
    next
  }
  t_result <- t.test(x ~ peak_full$TREATMENT) # value of chem_id like 35 as a string. 
  mean_diff <- t_result$estimate[2] - t_result$estimate[1] # mean in group 1 minus mean in group 0. 
  # For each chemical compute the mean difference in group 1 vs. group 2. 
  mean_diff_ratio <- t_result$estimate[2] / t_result$estimate[1] 
  
  # minimum magnitude of chemicals
  min_mag <- chem_df %>% 
    group_by(group) %>% 
    summarize(min_level = min(feature_level, na.rm = TRUE)) %>% 
    pivot_wider(names_from = group, values_from = min_level) %>% 
    unlist()
  
  results_row <- list(chem_id = chem_id, 
                      p.value = t_result$p.value, 
                      mean.diff = mean_diff,
                      mean.diff.ratio = mean_diff_ratio,
                      min.level.G0 = min_mag["G0"], 
                      min.level.G1 = min_mag["G1"])
  
  results_df <- rbind(results_df, results_row)
}

results_df <- left_join(results_df, 
                        chem_anno %>% select(chem_id = CHEM_ID, CHEMICAL_NAME) %>% 
                          mutate(chem_id = as.character(chem_id)))

# ADJUSTED P-VALUE:
results_df$p.adj <- p.adjust(results_df$p.value, method = "BH")

# VUSUALIZATION
# Create a subset with checimcal that are the top hit based on p-value (note I am not used p-adjusted atm)
top_pval <- results_df %>% 
  filter(p.value <= 0.05) %>%
  select(chem_id, CHEMICAL_NAME, p.value)
str(top_pval$chem_id)

peak_full_top <- peak_full %>% 
  select(TREATMENT, top_pval$chem_id) 

peak_full_top_long <- peak_full_top %>% 
  pivot_longer(cols = -TREATMENT, 
               names_to = "feature",
               values_to = "level")

peak_full_top_long <- left_join(peak_full_top_long, 
                                chem_anno %>% select(chem_id = CHEM_ID, CHEMICAL_NAME) %>% 
                                  mutate(chem_id = as.character(chem_id)), 
                                by = c("feature" = "chem_id"))
peak_full_top_long %>% 
  ggplot() + 
  aes(y = CHEMICAL_NAME, 
      x = level,
      color = TREATMENT) + 
  geom_boxplot() + 
  xlim(0, 15)
