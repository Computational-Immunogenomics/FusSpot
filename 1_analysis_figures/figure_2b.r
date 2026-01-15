source(paste0(dirname(getwd()),'/map.r'))
source(paste0(dirname(getwd()),'/shortcuts.r'))
source(paste0(dirname(getwd()),'/settings.r'))

base <- fread(paste0(TMP_DIR, "base.csv"))

tmp <- 
base %>% 
 fi(cohort == "Pan-Cancer") %>% 
 gb(Type2, VariantInfo, RNA_available = hasRNA,  RNA_Expressed = (RnaFrags > 0)) %>% 
 su(ct = n_distinct(Peptide)) %>% 
 ug()

fig_ready <- 
rbind(
 tmp %>% fi(RNA_available) %>% gb(Type2, VariantInfo) %>% su(ct = sum(ct)) %>% ug() %>% mu(rna = "All Events"), 
 tmp %>% fi(RNA_available, RNA_Expressed) %>% gb(Type2, VariantInfo) %>% su(ct = sum(ct)) %>% ug() %>% mu(rna = "RNA Expressed")
) %>%
 gb(Type2, rna) %>% 
 su(mn_peptides = mean(ct), sd_peptides = sd(ct), tot = n(), 
    moe = sd_peptides/sqrt(tot), pe = 1.96*sd_peptides + moe, 
    moe_sd = sd_peptides/sqrt(tot*2), 
    sd_ci_low = sd_peptides - 1.96*moe_sd, 
    sd_ci_high = sd_peptides + 1.96*moe_sd,     
    ci_low = mn_peptides  - 1.96*moe, ci_high = mn_peptides  + 1.96*moe, 
    pe_low = mn_peptides - 1.96*pe, pe_high = mn_peptides + 1.96*pe) %>% 
 ug() %>% 
 mu( Type =  case_when(
     Type2 %in% c("Frameshift" ) ~ "Frameshift",
     Type2 %in% c("Inframe Fusion", "Out-of-frame Fusion") ~ "Fusion",
     Type2 %in% c("Inframe Indel", "Missense") ~ "Missense")) %>% 
 mu( Type = factor(Type, levels = c("Frameshift", "Fusion", "Missense")),
     Type2 = factor(Type2, levels = c("Out-of-frame Fusion", "Frameshift", "Inframe Fusion", "Missense", "Inframe Indel")))

theme_go <- 
base_theme + 
 theme(
    axis.text.x = element_text(angle = .5, hjust = .5, size = 12),
    axis.text.y = element_text(angle = .5, hjust = .5, size = 12), 
    legend.position = c(.8, .65), legend.text = element_text(size = 10), 
    legend.spacing.y = unit(0.08, 'cm') )

fig_2b <- 
ggplot(fig_ready, aes(x = Type2, y = mn_peptides, alpha = rna, color = Type)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.5), width = .3) +
  scale_alpha_manual(values = c(.5, 1)) +
  scale_color_manual(values = unlist(color_map)) + 
  theme_go + ylim(0,5) + 
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  labs(color = NULL, alpha = NULL, x = "Event Type", y = "Mean Presented Peptides per Event", title = "Mean Peptides Presented Per Event Type") 

options(repr.plot.height = 5, repr.plot.width = 6)
fig_2b
ggsave(file = paste0(SHARE_DIR, "fig_2b.pdf"), width = 6, height = 5)
