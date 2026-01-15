source(paste0(dirname(getwd()),'/map.r'))
source(paste0(dirname(getwd()),'/shortcuts.r'))

summary_fusion <- fread(paste0(TMP_DIR, "summary_fusion.csv"))

share_reported_fusions <- 
summary_fusion %>% 
 fi(reported) %>%
 ar(desc(mutated_samples)) %>% 
 fi(cohort != "Pan-Cancer")

fwrite(share_reported_fusions, paste0(SHARE_DIR, "supplement_share_reported_fusions.csv"))
