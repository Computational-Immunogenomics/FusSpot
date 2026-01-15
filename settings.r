### Theme for figures
base_theme <- 
theme_bw() + 
theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    plot.title = element_text(hjust = .5, size = 16), 
    axis.text.x = element_text(angle = .5, hjust = .5, size = 7),
    axis.text.y = element_text(angle = .5, hjust = .5, size = 11),
    axis.title.y = element_text(size = 14), 
    axis.title.x = element_text(size = 14), 
    legend.position = c(.8, .8), legend.text = element_text(size = 14)
)

### Color settings for figures
color_map <-
c("Fusion" = "#A55A8A",
  "Missense" = "#5AA575",
  "Frameshift" = "#888888",
  "More Missense" = "#5AA575",
  "More Fusion" = "#A55A8A",
  "Similar" = "#888888")

alpha_map_1d = c("No RNA Available" = .5, "RNA Expressed" = 1, "RNA No Evidence" = .5)
color_map_1d = c("No RNA Available" = "grey", "RNA Expressed" = '#A55A8A', "RNA No Evidence" = '#A55A8A')

### Expressions 
cohort_filters <-
expr(!grepl("CUPPA", cohort_group) &
     !grepl("Unknown", cohort_group) &
     !grepl("Other", cohort_group) &
     !grepl("Unspecified", cohort_group) &
     cohort_group != "Pan-Cancer")

typer1 <- 
expr( 
 case_when(
     Type2 %in% c("Frameshift" ) ~ "Frameshift",
     Type2 %in% c("Inframe Fusion", "Out-of-frame Fusion") ~ "Fusion",
     Type2 %in% c("Inframe Indel", "Missense") ~ "Missense")
)

ct_gps <-
expr(
 case_when(
    ct == 1 ~ "1",
    ct <= 5 ~ "2-5",
    ct <= 20 ~ "6-20",
    ct > 20 ~ "20+")
)

### Highlight the cohorts
highlight <- c("Pan-Cancer", "ER+/HER+ Breast Carcinoma", "Gastrointestinal Stromal Tumor", "Triple Negative Breast Carcinoma")