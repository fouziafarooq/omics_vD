library(BiocManager)
library(omePath) # chemical to pathways
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(tibble)
library(ggrepel)
library(circlize)
library(vegan)
library(ggnewscale)

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

# Circular plots:
# https://r-graph-gallery.com/297-circular-barplot-with-groups.html

names(chem_anno) <- tolower(names(chem_anno))
str(chem_anno$chem_id)

sort(names(chem_anno))
chem_anno <- chem_anno %>%
  rowwise() %>%
  mutate(hmdb_id = strsplit(hmdb, ",")[[1]][1]) 

table(chem_anno$hmdb_id)

chem_ids <- names(peak_area)[-1] # -1 is don't take the first column. 
# chem_id_pairs <- combn(chem_ids, 2, simplify = FALSE) # don't need this. dont need pairs:-(((((((

#########################################################
# Process peak area data (impute with median)
# peak_area <- peak_area %>% 
#   mutate_at(vars(-PARENT_SAMPLE_NAME), ~ifelse(is.na(.), median(., na.rm=TRUE), .))# . means this one we are looking at.  ~ b/c we are creating a lambda fn. Second . means it is the else value (meaning keep it) 

# Impute with 1/2 
peak_area <- peak_area %>% 
  mutate_at(vars(-PARENT_SAMPLE_NAME), ~ifelse(is.na(.), 0.5*min(., na.rm = TRUE), .))# . means this one we are looking at.  ~ b/c we are creating a lambda fn. Second . means it is the else value (meaning keep it) 

# Transform the imputed data - log transformation.
peak_area <- peak_area %>% 
  mutate_at(vars(-PARENT_SAMPLE_NAME), ~log(., base = 10)) # log transforming the data. 

#########################################################

group_membership <- meta_data %>% select(PARENT_SAMPLE_NAME, TREATMENT)
peak_full <- left_join(peak_area, group_membership)

results_df <- data.frame()

for(chem_id in chem_ids){
  x<- peak_full[, chem_id] %>% unlist() # grab one column/chemical of data. 
  #Does this for a single chemical/feature
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
  # Don't need this. 
  min_mag <- chem_df %>% 
    group_by(group) %>% 
    summarize(min_level = min(feature_level, na.rm = TRUE)) %>% 
    pivot_wider(names_from = group, values_from = min_level) %>% 
    unlist()
  
  # group means
 group_means <- chem_df %>% 
    group_by(group) %>% 
    summarize(mean_level = mean(feature_level, na.rm = TRUE)) %>% 
    pivot_wider(names_from = group, values_from = mean_level) %>% 
    unlist()
  
  results_row <- list(chem_id = chem_id, 
                      p.value = t_result$p.value, 
                      mean.diff = mean_diff,
                      mean.diff.ratio = mean_diff_ratio,
                      min.level.G0 = min_mag["G0"], 
                      min.level.G1 = min_mag["G1"],
                      mean.level.G0 = group_means["G0"],
                      mean.level.G1 = group_means["G1"],
                      SE = t_result$stderr)
  
  results_df <- rbind(results_df, results_row)
}

results_df <- left_join(results_df, 
                        chem_anno %>% select(chem_id, chemical_name) %>% 
                          mutate(chem_id = as.character(chem_id)))

# ADJUSTED P-VALUE:
results_df$p.adj <- p.adjust(results_df$p.value, method = "BH")

# VUSUALIZATION
# Create a subset with checimcal that are the top hit based on p-value (note I am not used p-adjusted atm)
top_pval <- results_df %>% 
  filter(p.value <= 0.05) %>%
  select(chem_id, chemical_name, p.value)
str(top_pval$chem_id)

peak_full_top <- peak_full %>% 
  select(TREATMENT, top_pval$chem_id) 

peak_full_top_long <- peak_full_top %>% 
  pivot_longer(cols = -TREATMENT, 
               names_to = "feature",
               values_to = "level")

peak_full_top_long <- left_join(peak_full_top_long, 
                                chem_anno %>% select(chem_id, chemical_name) %>% 
                                  mutate(chem_id = as.character(chem_id)), 
                                by = c("feature" = "chem_id"))
peak_full_top_long %>% 
  ggplot() + 
  aes(y = chemical_name, 
      x = level,
      color = TREATMENT) + 
  geom_boxplot() + 
  xlim(0, 5)

################
# VOLCANO PLOT:
################
volcano_df <- results_df %>% 
  select(mean.diff, p.value, chem_id, chemical_name) %>% 
  rowwise() %>% 
  mutate(neg_log_p_value = -log(p.value, base = 10)) %>% # doing negative log of p-value to get the values to be inverted. 
  mutate(chem_name_label = ifelse(p.value<0.05 & ((mean.diff < -0.5) | (mean.diff > 0.5)), 
                                   chemical_name, ""))

volcano_plot <- volcano_df %>% 
  ggplot() + 
  geom_text_repel(aes(x = mean.diff, y = neg_log_p_value, label = chem_name_label),
                  min.segment.length = 0) +
  geom_point(aes(x = mean.diff, y = neg_log_p_value), alpha = 0.5) + # will take the negative log of p-value
 #  geom_text(aes(x = mean.diff, y = neg_log_p_value, label = chem_id), nudge_y = 0.05) + 
  geom_hline(yintercept = -log(0.05, base=10), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log(0.1, base=10), color = "red", linetype = "dashed") +   # just adding another line. 
  xlim(-1.2, 1.2)
  
volcano_plot


######
# Vector of names for chemicals that are significant
chems_interest <- volcano_df %>% 
  filter(chem_name_label != "") %>%
  select(chem_name_label) %>% 
  pull()

# Let's make boxplots where chems are p<0.05 and mean difference between the two groups is < -0.5 or >0.5.
box_plot_df <- peak_full_top_long %>% 
  filter(chemical_name %in% chems_interest)

box_plot_df %>% 
  ggplot() + 
  aes(y = chemical_name, 
      x = level,
      color = TREATMENT) + 
  geom_boxplot() + 
  xlim(0, 5)

# Adding pathway info on to the box_plot_df:
temp.df <- chem_anno %>% 
  filter(chemical_name %in% chems_interest)
temp.df <- left_join(box_plot_df, temp.df, by = "chemical_name")

#########################
# EXPORT FILE
# write_xlsx(temp.df, path = '~/box/1.GW/Career MODE/OMICS TANZANIA PILOT/vD Metabalomics/R/omics_vD/data_out/pathways_20240501.xlsx')
#####################
########################
# PCA:
########################
peak_area_pca <- peak_area %>% column_to_rownames(var = "PARENT_SAMPLE_NAME") # NOTE: tibble also allows numbers as col names. 
# peak_area_pca <- peak_area_pca %>% replace(is.na(.), 0) #Ask if this was ok to do. She said use imputed. 

pca <- prcomp(peak_area_pca, scale. = FALSE) #TODO ask mentors! 
#TODO PCA with abundance levels with NAs - do i first center the data then scale it. Then replace with 0s? 
pca_summary <- summary(pca)
pca_importance <- pca_summary$importance # this is in a matrix.

group_membership_rownames <- group_membership %>% column_to_rownames(var = "PARENT_SAMPLE_NAME") 

pca_plot_df <- merge(group_membership_rownames, as.data.frame(pca$x), by = 0)

pca_plot_df %>% ggplot() + 
  geom_point(aes(x=PC2, y=PC1, color = TREATMENT)) + 
  theme_light()

pca_plot_df %>% ggplot() + 
  geom_point(aes(x=PC3, y=PC1, color = TREATMENT)) + 
  theme_light()

##########################################
# Another plot for Proportion of Variance
# Looking for an 'elbow' in this plot. 

importance_plot_df <- as.data.frame(t(pca_importance))
importance_plot_df$PC <- 1:nrow(importance_plot_df)

importance_plot_df %>% 
  head(10) %>% # plotting the first 10 PCs. 
  ggplot() + 
  geom_line(aes(x=PC, y=`Proportion of Variance`)) + 
  geom_point(aes(x=PC, y=`Proportion of Variance`)) + 
  scale_x_continuous(breaks = 1:10)

######################################
# MERGE ON THE PATHWAYS USING CHEMID #
######################################
chem_anno$chem_id <- as.character(chem_anno$chem_id)
results_df <- left_join(results_df, chem_anno, by = "chem_id")

################
# VOLCANO PLOT:
################
# Doing the volcano plot again b/c now i have the pathways merged on. 
volcano_df <- results_df %>% 
  select(mean.diff, p.value, chem_id, super_pathway, sub_pathway) %>% 
  mutate(neg_log_p_value = -log(p.value, base = 10)) # doing negative log of p-value to get the values to be inverted. 

volcano_plot <- volcano_df %>% 
  ggplot() + 
  geom_point(aes(x = mean.diff, y = neg_log_p_value), alpha = 0.5) + # will take the negative log of p-value
   # geom_text(aes(x = mean.diff, y = neg_log_p_value, label = chem_id), nudge_y = 0.05) + 
  geom_hline(yintercept = -log(0.05, base=10), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log(0.1, base=10), color = "red", linetype = "dashed")   # just adding another line. 
volcano_plot

#########################
# EXPORT FILE
# write_xlsx(results_df, path = '~/box/1.GW/Career MODE/OMICS TANZANIA PILOT/vD Metabalomics/R/omics_vD/data_out/omics_results_202405017.xlsx')
#########################

########################
# PERMANOVA TEST
# We are testing to see if a group of metabolites (in a subpathway) is different among the two groups - it's a multivariate anova. 
# we do this for each subpathway. 
#######################

# First get the dataset in the right shape. 
subpathway_pval <- data.frame()
subpathway_names <- unique(chem_anno$sub_pathway)
subpathway_names <- subpathway_names[!is.na(subpathway_names)]

for (subpathway in subpathway_names) {
  sp_chem_ids <- chem_anno %>% 
    filter(sub_pathway == subpathway) %>%
    pull(chem_id) # pulls out this variable and makes it a vector.
  
  # new data frame
  peak_sp <- peak_full %>%
    select(sp_chem_ids)
  
  treatment_sp <- peak_full %>%
    select(TREATMENT)
  
  if (length(sp_chem_ids)==1) {
    t.test_permanova <- t.test(unlist(peak_sp) ~ treatment_sp$TREATMENT)
    subpathway_row <- data.frame(sub_pathway = subpathway, pval = t.test_permanova$p.value)
    subpathway_pval <- rbind(subpathway_pval, subpathway_row)
  } else {
    sp_permanova <- adonis2(peak_sp ~ TREATMENT, data = treatment_sp)
    subpathway_row <- data.frame(sub_pathway = subpathway, pval = sp_permanova$`Pr(>F)`[1])
    subpathway_pval <- rbind(subpathway_pval, subpathway_row) # we will join this subpathway_pval dataframe with the subpathway_df that we have in the circular plots.
  }
}



########################
# CIRCULAR PLOTS
#######################
lipids_df <- results_df %>% 
  select(mean.diff, chemical_name.x, super_pathway, sub_pathway, p.value) %>% 
  filter(super_pathway == "Lipid") %>% 
#  filter(p.value<=0.05) %>% # head(20) %>%
  arrange(sub_pathway) %>% 
  mutate(i = row_number())

# Lipids
lipids_df$mean.diff <- unname(lipids_df$mean.diff) # removes the names

subpathway_df <- lipids_df %>% 
  group_by(sub_pathway) %>% 
  summarize(min_x = min(i), max_x = max(i)) %>%
  left_join(subpathway_pval) %>%
  mutate(significant = ifelse(pval <= 0.05, "p ≤ 0.05", "p > 0.05")) %>%
  mutate(subpathway_label = ifelse(significant=='p ≤ 0.05', sub_pathway, ""))

# calculate angle
xmax <- max(subpathway_df$max_x)


ggplot() +
  ylim(-1.5, 0.8) + 
  geom_col(data = lipids_df, # data coming from lipids_df
           aes(x = i,
               y = mean.diff*1.8,
               fill = mean.diff),# Identity forces the bar plot to not count on it's own. You have to specify what it will be. 
           color = 'black', linewidth = 0.05) + 
  coord_polar() + 
  scale_fill_distiller(palette = "Spectral") + # for continuous.
#  scale_fill_gradient2(low = "blue",
 #                      mid = "orange",
  #                     high = "red") + 
  theme_minimal() + 
  theme(legend.position = "none",
        axis.text = element_blank(), 
        axis.title = element_blank(), 
        panel.grid = element_blank()) +
  new_scale_fill() + # allows you to use a new color scale fill. 
  geom_rect(data=subpathway_df, # data coming from subpathway_df
               aes(xmin = min_x-0.1, ymin = -0.9, xmax = max_x+0.1, ymax = -0.8,
                   fill = significant)) +
  scale_fill_manual(values = c('p > 0.05' = 'lightblue', 'p ≤ 0.05' = 'red')) +
  #geom_vline(data = subpathway_df, aes(xintercept = min_x-0.5))
  geom_segment(data = subpathway_df, aes(x = min_x-0.5, xend=min_x-0.5, y = -0.9, yend = 0.5), linewidth = 0.1) + 
  geom_text(data = subpathway_df, aes(x = 0.5*(min_x + max_x),
                                      label = subpathway_label,
                                      angle = 90-0.5*(min_x + max_x)*360/xmax), 
            y = 0.55,
            hjust = 0,
            size = 2) + 
  geom_text(data = subpathway_df, x = 0, y = -1.45, label = "Lipids")
  
#############################################
# Trying to do a loop for all superpathways

results_df2 <- results_df %>% 
  drop_na(super_pathway) %>%
  filter(super_pathway!='Xenobiotics') %>%
  filter(super_pathway!='Partially Characterized Molecules')
superpathway_list <- unique(results_df2$super_pathway)

for(superpathway in superpathway_list) {
  superpathway_df <- results_df2 %>% 
    select(mean.diff, chemical_name.x, super_pathway, sub_pathway, p.value) %>% 
    filter(super_pathway == superpathway) %>% #superpathway is the looping variable that is in the for()
    #  filter(p.value<=0.05) %>% # head(20) %>%
    arrange(sub_pathway) %>% 
    mutate(i = row_number())
  
  superpathway_df$mean.diff <- unname(superpathway_df$mean.diff) # removes the names
  
  subpathway_df <- superpathway_df %>% 
    group_by(sub_pathway) %>% 
    summarize(min_x = min(i), max_x = max(i)) %>%
    left_join(subpathway_pval) %>%
    mutate(significant = ifelse(pval <= 0.05, "p ≤ 0.05", "p > 0.05")) %>%
    mutate(subpathway_label = ifelse(significant=='p ≤ 0.05', sub_pathway, sub_pathway),
           angle = 90-0.5*(min_x + max_x)*360/max(max_x))
  
  # calculate angle
  xmax <- max(subpathway_df$max_x)
  
  #######
  # ggplot
  #######
  
  p_plot <- ggplot() +
    ylim(-1.5, 1.0) + 
    geom_col(data = superpathway_df, # data coming from lipids_df
             aes(x = i,
                 y = mean.diff*1.8,
                 fill = mean.diff),# Identity forces the bar plot to not count on it's own. You have to specify what it will be. 
             color = 'black', linewidth = 0.05) + 
    coord_polar() + 
    scale_fill_distiller(palette = "Spectral") + # for continuous.
    #  scale_fill_gradient2(low = "blue",
    #                      mid = "orange",
    #                     high = "red") + 
    theme_minimal() + 
    theme(legend.position = "none",
          axis.text = element_blank(), 
          axis.title = element_blank(), 
          panel.grid = element_blank()) +
    new_scale_fill() + # allows you to use a new color scale fill. 
    geom_rect(data=subpathway_df, # data coming from subpathway_df
              aes(xmin = min_x-0.1, ymin = -0.9, xmax = max_x+0.1, ymax = -0.8,
                  fill = significant)) +
    scale_fill_manual(values = c('p > 0.05' = 'lightblue', 'p ≤ 0.05' = 'red')) +
    #geom_vline(data = subpathway_df, aes(xintercept = min_x-0.5))
    geom_segment(data = subpathway_df, aes(x = min_x-0.5, xend=min_x-0.5, y = -0.9, yend = 0.5), linewidth = 0.1) + 
    geom_text(data = subpathway_df, aes(x = 0.5*(min_x + max_x),
                                        label = subpathway_label,
                                        angle = ifelse(between(angle, -270, -90), angle+180, angle),
                                        hjust = ifelse(between(angle, -270, -90), 1, 0)),
              y = 0.55,
              size = 1.7) + 
    geom_text(data = subpathway_df, x = 0, y = -1.45, label = superpathway)
  
  ggsave(filename = paste0('plots/',superpathway, '.pdf'), p_plot)
  
}


