source(paste0(dirname(getwd()),'/map.r'))
source(paste0(dirname(getwd()),'/shortcuts.r'))
source(paste0(dirname(getwd()),'/settings.r'))

base_peptide <- fread(paste0(TMP_DIR, "base_peptide.csv"))
cohort_summary_group <- fread(paste0(TMP_DIR, "cohorts_summary_group.csv"))

tmp <- 
base_peptide %>%
 fi(!!cohort_filters)  %>%  
 gb(cohort_group, acronym, sampleId) %>%
 su(Peptides = n(), `Peptides RNA Expressed`= sum(RnaFrags > 0), rna_avail = sum(hasRNA), `Fusion Events` = n_distinct(name)) %>% 
 ug() %>% 
 mu( `Expressed Peptides` = ifelse(rna_avail == 0, NA, `Peptides RNA Expressed`)) %>% 
 ij( cohort_summary_group, by = "cohort_group") %>% 
 se( cohort_group, acronym, Peptides , `Peptides RNA Expressed`, `Fusion Events`) %>% 
 ga(field, val, -cohort_group, -acronym) %>% 
 mu(field = factor(field, levels = c( "Peptides", "Peptides RNA Expressed"))) 

tmp_all <- tmp %>% fi(field == "Peptides") %>% mu( cohort_group = fct_reorder(cohort_group, val, .fun = median, .desc = TRUE), acronym = fct_reorder(acronym, val, .fun = median, .desc = TRUE))
tmp_rna <- tmp %>% fi(field == "Peptides RNA Expressed") %>% mu( cohort_group = factor(cohort_group, levels = levels(tmp_all$cohort_group)), acronym = factor(acronym, levels = levels(tmp_all$acronym)))

fig_ready <- rbind(tmp_all, tmp_rna)

theme_go <-  base_theme + theme( legend.position = c(.85, .85), axis.text.x = element_text(angle = 0, hjust = .5, vjust = 0, size = 12), axis.text.y = element_text(angle = 0, hjust = 0, vjust = 0, size = 15), axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16), plot.title = element_text(hjust = .5, size = 20) ) 

fig_2a <- 
fig_ready %>% 
 fi(field == "Peptides") %>% 
 ggplot( aes(acronym, y = val, fill = field )) +
 geom_boxplot(outlier.colour = "white", outlier.fill = "white", outlier.shape = 21, fill = color_map['Fusion']) + 
 geom_jitter(shape = 21,position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.75),alpha = .1, size = 2, aes(fill = color_map['Fusion'])) +
 theme_go + 
 labs(fill = NULL, x = "Cohorts", y = "High Likelihood Peptides Per Sample", 
      title = "Distribution of Number of Peptides with High Presentation Likelihood Generated from Gene Fusion Events") + 
 scale_y_log10() + scale_x_discrete(expand = expansion(mult = c(0.01, .01))) + 
 theme(legend.position = "none") + 
 scale_x_discrete(guide = guide_axis(n.dodge=2)) 

options(repr.plot.height = 6, repr.plot.width = 18)
fig_2a
ggsave(fig_1a, file = paste0(SHARE_DIR, "fig_2a.pdf"), width = 15, height = 6)

fig_1a_supplement <- 
fig_ready %>% 
 ggplot( aes(acronym, y = val, fill = "a", alpha = field )) +
 geom_boxplot(outlier.colour = "white", outlier.fill = "white", outlier.shape = 21, fill = color_map['Fusion']) + 
 geom_jitter( shape = 21,position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.75),alpha = 0.1, size = 1, show.legend = FALSE) +
 theme_go + 
 scale_alpha_manual( values = c(1, .2)) + 
 labs(fill = NULL, alpha = NULL, color = NULL, shape = NULL, x = "Cohorts", y = "High Likelihood Peptides Per Sample", 
      title = "Distribution of Number of Peptides with High Presentation Likelihood Generated from Gene Fusion Events") + 
 scale_y_log10() + 
 scale_x_discrete(guide = guide_axis(n.dodge=2)) 

options(repr.plot.height = 6, repr.plot.width = 15)
fig_1a_supplement
ggsave(fig_1a_supplement, file = "supplement_fig_1a.pdf", width = 15, height = 6)
