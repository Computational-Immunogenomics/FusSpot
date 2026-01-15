source(paste0(dirname(getwd()),'/map.r'))
source(paste0(dirname(getwd()),'/shortcuts.r'))

base <- fread(paste0(TMP_DIR, "base.csv"))

share <-  
base %>% 
 se(-cohort_group, -contains("Flank"), -contains("Type"), -cohort_samples) %>% 
 relocate(cohort, acronym) %>% 
 fi(cohort != "Pan-Cancer")

fwrite( share , paste0(SHARE_DIR, "supplement_big_table.csv"))
