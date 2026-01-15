source(paste0(dirname(getwd()),'/map.r'))
source(paste0(dirname(getwd()),'/shortcuts.r'))
source(paste0(dirname(getwd()),'/settings.r'))

library(scales)

summary_fusion <- fread(paste0(TMP_DIR, "summary_fusion.csv"))

fig_ready <- 
summary_fusion %>% 
 fi(reported, !!cohort_filters, cohort_samples > 20, mutated_samples > 1) %>% 
 ar(desc(pct_peptide)) %>% 
 head(15) %>% 
 mu(go = paste0(name, "\n", cohort, "\n")) %>% 
 se(go, mutated_samples, has_rna_exp, has_rna, cohort_samples, n_alleles, top_Allele) %>% 
 mu(`No RNA Available` = mutated_samples - has_rna, `RNA No Evidence` = has_rna - has_rna_exp) %>% 
 rename(`RNA Expressed` = has_rna_exp) %>% 
 se(-mutated_samples, -has_rna) %>% 
 ga(gp, val, -go, -cohort_samples, -n_alleles, -top_Allele) %>% 
 mu(pct = val/cohort_samples) %>% 
 gb(go) %>% mu(tot_pct = sum(pct), tot_samples = sum(val)) %>% ug() %>% 
 mu(go = fct_reorder(go, desc(tot_pct)), 
    gp = factor(gp, levels = rev(c("RNA Expressed", "No RNA Available", "RNA No Evidence"))), 
    go2 = ifelse(n_alleles > 4, paste0(go, n_alleles, " HLA Alleles\n"), paste0(go, top_Allele, "\n")), 
    go2 = fct_reorder(go2, desc(tot_pct))) %>%
 mu(gp = factor(gp, levels = rev(c("RNA Expressed", "RNA No Evidence", "No RNA Available"))))

fig_1d <- 
fig_ready %>% 
 ggplot(aes(x = go2, y = pct, alpha = gp, fill = gp)) +
 geom_bar(stat = "identity", color = "black", width = .7) + 
 geom_text(aes(label = paste0(tot_samples, " total"), y = tot_pct), vjust = -0.5, size = 3, show.legend = FALSE) + 
 scale_alpha_manual(values = alpha_map_1d) + 
 scale_fill_manual(values = color_map_1d) + 
 guides(
    fill = guide_legend(override.aes = list(alpha = alpha_map_1d)),
    alpha = "none"  # hide the separate alpha legend
  ) +
 base_theme + 
 scale_y_continuous(labels = label_percent(), expand = expansion(mult = c(0, 0.1))) + 
 labs(alpha = NULL, fill = NULL, label = NULL, text = NULL, x = "Fusion Event within Cohort", y = "Fusion Cohort Prevalence", title = "Top Reported Gene Fusions with Predicted Presented Peptides") + 
 scale_x_discrete(guide = guide_axis(n.dodge=2)) + 
 theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))

options(repr.plot.width = 10.5, repr.plot.height = 5)
fig_1d
ggsave(fig_1d, file = paste0(SHARE_DIR, "fig_1d.pdf"), width = 10.5, height = 5)
