source(paste0(dirname(getwd()),'/map.r'))
source(paste0(dirname(getwd()),'/shortcuts.r'))

base <- fread(paste0(TMP_DIR, "base.csv"))

number_events <- 
base %>% 
 fi(cohort != "Pan-Cancer", grepl("FUSION", VariantType)) %>% 
 gb(cohort,sampleId) %>% 
 su(unique_fusions = n_distinct(name),
    unique_fusion_derived_peptides = n_distinct(Peptide)) 

reported_events <- 
base %>% 
 fi(cohort != "Pan-Cancer", grepl("FUSION", VariantType), reported) %>% 
 gb(cohort, sampleId) %>% 
 su(reported_unique_fusions = n_distinct(name), peptides_from_reported_unique_fusions = n_distinct(Peptide)) 

supplement_share_events <- 
number_events %>% 
 lj(reported_events, by = c("cohort", "sampleId")) %>% 
  mutate(across(everything(), ~replace_na(.x, 0)))

fwrite(supplement_share_events, paste0(SHARE_DIR, "supplement_events_per_sample.csv"))
