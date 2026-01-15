source(paste0(dirname(getwd()),'/map.r'))n
source(paste0(dirname(getwd()),'/shortcuts.r'))

linx_fusions <- fread(paste0(TMP_DIR, "fusions.csv")) 
neo_pep <- fread(paste0(TMP_DIR, "neo_pep.csv")) 
neo <- fread(paste0(TMP_DIR, "neo.csv"))
fragile_sites <- fread(paste0(dirname(getwd()), "/ref/fragile_sites.37.csv"))

cohorts <- fread(paste0(TMP_DIR, "cohorts.csv"))
cohort_summary <- fread(paste0(TMP_DIR, "cohorts_summary.csv"))
cohort_summary_group <- fread(paste0(TMP_DIR, "cohorts_summary_group.csv"))

#ggplot(neo_pep %>%
# gb(sampleId) %>%
# su(tot = n()) %>%
# ar(desc(tot)), aes(x = 1, y= tot)) + geom_boxplot()

neo_pep_ready <- 
neo_pep %>% 
 tm(sampleId, name = Gene, VariantType, Allele, Peptide, Score, Rank, Likelihood, LikelihoodRank, FlankUp, FlankDown) %>% 
 unique()

neo_ready <- 
neo %>% 
 tm(sampleId, name = GeneName, RnaFrags, RnaDepth, VariantInfo, 
    TranscriptsUp, TranscriptsDown, VariantCopyNumber,	CopyNumber,	SubclonalLikelihood) %>% 
 unique()
neo_rna <- neo %>% gb(sampleId) %>% su(tot = max(RnaDepth)) %>% ug() %>% tm(sampleId, hasRNA = (tot > 0))

neo_rna <- 
neo %>% 
 gb(sampleId) %>% 
 su(tot = max(RnaDepth)) %>%
 ug() %>% 
 tm(sampleId, hasRNA = (tot > 0))

linx_fusions_ready <- 
linx_fusions %>% 
 se(name, reported, reportedType, junctionCopyNumber, sampleId) %>% 
 rw() %>% mu(gene_up = strsplit(name, "_")[[1]][1], gene_down = strsplit(name, "_")[[1]][2]) %>% ug() %>% 
 mu(name = ifelse(gene_up == gene_down, gene_down, name), 
    fragile = (gene_down %in% fragile_sites$Gene  | gene_up %in% fragile_sites$Gene)) %>% 
 unique() %>% 
 se(-gene_up, -gene_down)

base <- 
neo_pep_ready %>% 
 lj(cohorts, by = c("sampleId")) %>%  
 lj(neo_ready, by = c("sampleId", "name")) %>% 
 lj(neo_rna, by = "sampleId") %>%
 lj(linx_fusions_ready, by = c("sampleId", "name")) %>% 
 unique() %>% 
 mu(Type = 
     case_when(
      VariantType %in% c("OUT_OF_FRAME_FUSION", "INFRAME_FUSION") ~ "Fusion",
      VariantType %in% c("MISSENSE", "INFRAME_DELETION", "INFRAME_INSERTION") ~ "Missense",
      VariantType %in% c("FRAMESHIFT", "STOP_LOST") ~ "Frameshift"), 
    Type2 = 
     case_when(
      VariantType %in% c("OUT_OF_FRAME_FUSION") ~ "Out-of-frame Fusion",
      VariantType %in% c("INFRAME_FUSION") ~ "Inframe Fusion",
      VariantType %in% c("MISSENSE") ~ "Missense",
      VariantType %in% c("INFRAME_DELETION", "INFRAME_INSERTION") ~ "Inframe Indel", 
      VariantType %in% c("FRAMESHIFT", "STOP_LOST") ~ "Frameshift"))

base_peptide <- 
base  %>% 
 gb(cohort, acronym, sampleId, Peptide) %>% 
 mu(rk = row_number(desc(Score)), alleles = paste0(unique(Allele), collapse = ",")) %>% 
 ug() %>% 
 fi(rk == 1) %>% se(-rk)

dim(base_peptide)

summary_peptide <- 
base_peptide %>% 
 gb(cohort, acronym, Peptide, reported) %>% 
 su(mutated_samples = n(), 
    fusion_partners = paste0(unique(name), collapse = ","), 
    n_alleles = n_distinct(Allele),
    top_Allele = paste0(unique(Allele), collapse = ","),
    all_Alleles = paste0(unique(alleles), collapse = ","),
    has_rna_exp = sum(RnaFrags > 0), 
    has_rna = sum(hasRNA > 0)) %>% 
 ug() %>% 
 rj(cohort_summary, by = "cohort") %>% 
 mu(pct_peptide_rna_exp = has_rna_exp/has_rna, pct_peptide = mutated_samples/cohort_samples) 

base_fusion <- 
base  %>% 
 gb(cohort, acronym, sampleId, name) %>% 
 mu(rk = row_number(desc(Score)), alleles = paste0(unique(Allele), collapse = ",")) %>% 
 ug() %>% 
 fi(rk == 1) %>% se(-rk) ### Only select Peptide with top score

summary_fusion <- 
base_fusion  %>% 
 gb(cohort, acronym, name, reported) %>% 
 su(mutated_samples = n(), 
    peptides = paste0(unique(Peptide), collapse = ","), 
    n_alleles = n_distinct(Allele),
    top_Allele = paste0(unique(Allele), collapse = ","),
    all_Alleles = paste0(unique(alleles), collapse = ","),
    has_rna_exp = sum(RnaFrags > 0), 
    has_rna = sum(hasRNA > 0)) %>% 
 ug() %>% 
 lj(cohort_summary, by = "cohort") %>% 
 mu(pct_peptide_rna_exp = has_rna_exp/has_rna, pct_peptide = mutated_samples/cohort_samples) 

fwrite( base , paste0(TMP_DIR, "base.csv"))
fwrite( base_peptide , paste0(TMP_DIR, "base_peptide.csv"))
fwrite( summary_peptide , paste0(TMP_DIR, "summary_peptide"))
fwrite( base_fusion , paste0(TMP_DIR, "base_fusion.csv"))
fwrite( summary_fusion , paste0(TMP_DIR, "summary_fusion.csv"))
