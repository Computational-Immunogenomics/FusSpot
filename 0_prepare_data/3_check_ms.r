source(paste0(dirname(getwd()),'/map.r'))
source(paste0(dirname(getwd()),'/shortcuts.r'))

cohorts <- fread(paste0(TMP_DIR, "cohorts.csv"))
identified_ms <- fread(paste0(TMP_DIR, "identified_pep_ms.csv"))

ready <-
cohorts %>%
 se(sampleId, cohort, acronym) %>%
 ij(identified_ms, by = "sampleId") %>%
 fi(cohort != "Pan-Cancer")

validated <- 
expr((Gene == "EML4_ALK" & Peptide == "DVIINQVYR" & Allele %in% c("A6801", "A6601", "A3303", "A3401")) |
     (Gene == "MYB_NFIB" & Peptide == "NELKGQQSW" & Allele %in% c("B4403", "B4402")) | 
     (Gene == "FUS_DDIT3" & Peptide == "SPHLKADVL" & Allele %in% c("B0801", "B0702")))

ms_driven <- 
expr((Gene == "EML4_ALK" & Peptide == "DVIINQVYRR" & Allele %in% c("A6801")) |
     (Gene == "EML4_ALK" & Peptide == "VIINQVYR" & Allele %in% c("A6801", "A6601")) | 
     (Gene == "EML4_ALK" & Peptide == "VIINQVYRR" & Allele %in% c("A6801")))

share <- ready %>% 
 mu(label = case_when(
    !!validated ~ "Top Sample Score",
    !!ms_driven ~ "Mass Spectrometry Validated",
    TRUE ~ "C")) %>%
 fi(label != "C") %>%
 ar(sampleId, desc(Score)) %>%
 se(sampleId, label, cohort, acronym, VariantType, Gene, Allele, Peptide, label, Score, LikelihoodRank) %>%
 gb(label, cohort, Gene,  Peptide, Allele, VariantType,  Score, LikelihoodRank) %>%
 su(number_of_patients = n()) %>%
 ar(cohort, desc(label), desc(number_of_patients)) %>%
 mu(disease_entity =
    case_when(
     cohort == "Lung Adenocarcinoma" ~ "NSCLC",
     cohort == "Liposarcoma" ~ "MLS",       
     cohort == "Salivary Gland Carcinoma" ~ "ACC",
     TRUE ~ "Other")) %>%
 fi(disease_entity != "Other") %>%
 relocate(disease_entity) %>%
 ar(disease_entity, desc(label), desc(number_of_patients)) %>% ug() %>% rename(Fusion = Gene, NEO_Score = Score) %>% relocate(cohort)

fwrite(share, paste0(SHARE_DIR, "mass_spectrometry_validation_algorithm_scores.csv"))
