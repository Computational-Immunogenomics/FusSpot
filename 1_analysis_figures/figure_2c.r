source(paste0(dirname(getwd()),'/map.r'))
source(paste0(dirname(getwd()),'/shortcuts.r'))
source(paste0(dirname(getwd()),'/settings.r'))

library(ggrepel)

base_fusion <- fread(paste0(TMP_DIR, "base_fusion.csv"))
base_peptide <- fread(paste0(TMP_DIR, "base_peptide.csv"))

peptide <- base_peptide %>% fi(!!cohort_filters)
fusion <- base_fusion %>% fi(!!cohort_filters)
peptide_rna <- base_peptide %>% fi(RnaFrags > 0, hasRNA, !!cohort_filters)

go_go <- function( df ) {
 df %>% 
  gb(cohort, acronym, Type, sampleId) %>% 
  su(ct = n_distinct(Peptide)) %>% 
  gb(cohort, acronym, Type) %>% 
  su(mn_events = median(ct)) %>%   
  ug() %>% 
  pivot_wider(names_from = Type, values_from = mn_events) 
}

fig_ready <- 
rbind(go_go(fusion) %>% mu(gp = "Events"), go_go(peptide) %>% mu(gp = "Peptides"), 
      go_go(peptide_rna) %>% mu(gp = "Peptides RNA")) %>% 
 mu( diff = Fusion - Missense, 
     gp2 = case_when((diff >= 15 & gp == "Peptides") | (diff >= 5 & gp != "Peptides") ~ "More Fusion",
                     (diff <= -15 & gp == "Peptides") | (diff <= -5 & gp != "Peptides") ~ "More Missense",  
                     (abs(diff) <= 15 & gp == "Peptides") | (abs(diff) <= 5 & gp != "Peptides") ~ "Similar")) 

theme_go <- 
base_theme + 
 theme(legend.text = element_text(size = 9), legend.position = c(.8, .8), 
       axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), 
       plot.title = element_text(size = 14), 
       axis.title.y = element_text(size = 13), 
       axis.title.x = element_text(size = 13) )

fig1c <- 
fig_ready %>% fi(gp == "Peptides") %>% 
 ggplot(aes( x = Fusion, y = Missense, color = gp2)) + 
 geom_point() + 
 scale_color_manual(values = color_map) + 
 geom_text_repel(aes(label = acronym), size = 3, show.legend = FALSE) +
 theme_go + 
 geom_abline(intercept = -15, slope = 1, alpha = .1) + 
 geom_abline(intercept = 15, slope = 1, alpha = .1) + 
 xlim(0,75) + ylim(0,75) +  
 labs(color = NULL, x = "Median Fusion Derived Peptides", y = "Median Missense Peptides", title = "Peptides Presented with High Likelihood")

options(repr.plot.height = 4.5, repr.plot.width = 4.5)
fig1c
ggsave(file = paste0(SHARE_DIR,"fig_1c.pdf"), width = 4.5, height = 4.5)
