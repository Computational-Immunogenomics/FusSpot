source(paste0(dirname(getwd()),'/map.r'))
source(paste0(dirname(getwd()),'/shortcuts.r'))

### NEO_DIR <- INPUT DIRECTORY SPECIFIED IN map.r
### SOM_DIR <- INPUT DIRECTORY SPECIFIED IN map.r
### TMP <- OUTPUT DIRECTORY SPECIFIED IN map.r

get_fp <- function (i, type = "purity") {
  if (type == "neo_pep") { ac(paste0(NEO_DIR, i, "/", i, ".neo.peptide_scores.tsv.gz"))} 
  else if (type == "neo") { ac(paste0(NEO_DIR, i, "/", i, ".neo.neoepitope.tsv.gz"))} 
  else if (type == "purity") { ac(paste0(SOM_DIR, i, "/purple/", i, ".purple.purity.tsv"))} 
  else if (type == "fusion") { ac(paste0(SOM_DIR, i, "/linx/", i, ".linx.fusion.tsv"))}
}

reader <- function (i_file = "ab", sample = "cd") { fread(i_file) %>% mutate(sampleId = sample)}

patients <- list.files(NEO_DIR)

score_threshold <- 5
lr_threshold <- .001

tmp <- list()
system.time(
for( i in patients ){
  i_file <- get_fp(i, type = "neo_pep")
  if(file.exists(i_file)){    
    oo <- tryCatch({reader( i_file, sample = i ) %>% fi(Score > score_threshold, LikelihoodRank < lr_threshold)}, error = function(e) {print("No!"); NA})
    if(is.data.frame(oo)){ tmp[[i]] <- oo }}}
)
fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "neo_pep.csv"))

tmp <- list(); 
system.time(
for( i in patients ){
  i_file <- get_fp(i, type = "neo")
  if(file.exists(i_file)){    
      tmp[[i]] <- reader( i_file = i_file, sample = i)
  }})
fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "neo.csv"))

tmp <- list(); 
system.time(
for( i in patients ){
  i_file <- get_fp(i, type = "fusion")
  if(file.exists(i_file)){    
      tmp[[i]] <- reader( i_file = i_file, sample = i)
  }})
fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "fusions.csv"))

ms_peptides <- 
expr((Gene == "EML4_ALK" & Peptide %in% c("DVIINQVYR", "DVIINQVYRR", "VIINQVYR", "VIINQVYRR")) |
(Gene == "MYB_NFIB" & Peptide %in% c("NELKGQQSW")) |
(Gene == "FUS_DDIT3" & Peptide %in% c("SPHLKADVL")))

tmp <- list()
system.time(
for( i in patients ){
  i_file <- get_fp(i, type = "neo_pep")
  if(file.exists(i_file)){    
    oo <- tryCatch({reader( i_file, sample = i ) %>% mu(med_score = median(Score), med_Likelihood = median(Likelihood)) %>% fi(!!ms_peptides)}, error = function(e) {print("No!"); NA})
    if(is.data.frame(oo)){ tmp[[i]] <- oo }}}
)
fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "identified_pep_ms.csv"))

identified_ms <- fread(paste0(TMP_DIR, "identified_pep_ms.csv"))

identified_ms %>% 
 fi(Peptide == "DVIINQVYRR", Allele == "A6801") %>% 
 ar(desc(Score))
