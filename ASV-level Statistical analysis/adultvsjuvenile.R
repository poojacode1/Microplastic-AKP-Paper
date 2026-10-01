
# ASV-LEVEL DEVELOPMENTAL-STAGE ANALYSIS
#
# Comparison:
#   23-day exposure – adult
#   23-day exposure – juvenile





library(phyloseq)
library(vegan)
library(tidyverse)




setwd(
  "C:/Users/u0175886/Documents/Microbiome_Microplastic/Exposure_Microplastic"
)

pseq <- readRDS("daphnia_phyloseq.rds")


cat("Exposure_Condition labels in original phyloseq object\n")


print(
  table(
    as.character(sample_data(pseq)$Exposure_Condition),
    useNA = "ifany"
  )
)


# KEEP ONLY ADULT AND JUVENILE SAMPLES AFTER 23 DAYS

keep_labels <- c(
  "Post_exposure",
  "Juvenile"
)

cond <- as.character(
  sample_data(pseq)$Exposure_Condition
)

ps <- prune_samples(
  sample_names(pseq)[cond %in% keep_labels],
  pseq
)


# TAXON FILTERING

ps <- filter_taxa(
  ps,
  function(x) {
    sum(x > 0) >= 3 && sum(x) >= 10
  },
  prune = TRUE
)



#  RELATIVE ABUNDANCE


ps <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)


# BRAY–CURTIS DISSIMILARITIES


d <- phyloseq::distance(
  ps,
  method = "bray"
)


# 6. PREPARE METADATA


md <- data.frame(
  sample_data(ps),
  stringsAsFactors = FALSE
)

# Ensure metadata order matches the Bray–Curtis distance object
md <- md[
  labels(d),
  ,
  drop = FALSE
]


# figures

md$Stage <- factor(
  as.character(md$Exposure_Condition),
  
  levels = c(
    "Post_exposure",
    "Juvenile"
  ),
  
  labels = c(
    "23-day exposure – adult",
    "23-day exposure – juvenile"
  )
)



cat("Group sizes used in developmental-stage analysis\n")

print(table(md$Stage))


# Safety checks
stopifnot(
  nlevels(md$Stage) == 2
)

stopifnot(
  !anyNA(md$Stage)
)

stopifnot(
  identical(
    rownames(md),
    labels(d)
  )
)



#
# Question:
# Does microbiome COMPOSITION differ between adults and juveniles
# after 23 days of exposure?

set.seed(1)

ad <- vegan::adonis2(
  d ~ Stage,
  data = md,
  permutations = 9999
)

cat("Adult vs juvenile after 23 days of exposure\n")

print(ad)



# 8. PERMDISP


bd <- vegan::betadisper(
  d,
  group = md$Stage,
  type = "centroid"
)


# Permutation test
set.seed(1)

pt <- permutest(
  bd,
  permutations = 9999
)



cat("PERMDISP\n")
cat("Adult vs juvenile after 23 days of exposure\n")


print(pt)


cat("\nMean distance to group centroid:\n")
print(
  round(
    bd$group.distances,
    4
  )
)


cat("\nMedian distance to group centroid:\n")
print(
  tapply(
    bd$distances,
    md$Stage,
    median
  )
)



# EXTRACT HEADLINE STATISTICS

permanova_R2 <- ad$R2[1]
permanova_F  <- ad$F[1]
permanova_p  <- ad$`Pr(>F)`[1]

permdisp_F <- pt$tab$F[1]
permdisp_p <- pt$tab$`Pr(>F)`[1]


adult_mean <- bd$group.distances[
  "23-day exposure – adult"
]

juvenile_mean <- bd$group.distances[
  "23-day exposure – juvenile"
]




cat(
  sprintf(
    "PERMANOVA: R2 = %.4f, F = %.3f, p = %.4f\n",
    permanova_R2,
    permanova_F,
    permanova_p
  )
)

cat(
  sprintf(
    "PERMDISP: F = %.3f, p = %.4f\n",
    permdisp_F,
    permdisp_p
  )
)

cat(
  sprintf(
    paste0(
      "Mean distance to centroid:\n",
      "  Adult    = %.4f\n",
      "  Juvenile = %.4f\n"
    ),
    adult_mean,
    juvenile_mean
  )
)


# PER-SAMPLE DISPERSION TABLE


disp_df <- data.frame(
  
  Sample = names(
    bd$distances
  ),
  
  Stage = md[
    names(bd$distances),
    "Stage"
  ],
  
  Distance_to_Centroid = as.numeric(
    bd$distances
  )
  
)


disp_df$Stage <- factor(
  disp_df$Stage,
  
  levels = c(
    "23-day exposure – adult",
    "23-day exposure – juvenile"
  )
)

# DISPERSION SUMMARY


disp_summary <- disp_df %>%
  
  group_by(Stage) %>%
  
  summarise(
    
    n = n(),
    
    mean_distance = mean(
      Distance_to_Centroid
    ),
    
    median_distance = median(
      Distance_to_Centroid
    ),
    
    sd_distance = sd(
      Distance_to_Centroid
    ),
    
    se_distance =
      sd_distance / sqrt(n),
    
    .groups = "drop"
  )


cat("Beta-dispersion summary\n")


print(disp_summary)


# Export for supplementary statistics
write.csv(
  disp_summary,
  "DevelopmentalStage_BetaDispersion_Summary.csv",
  row.names = FALSE
)

write.csv(
  disp_df,
  "DevelopmentalStage_BetaDispersion_PerSample.csv",
  row.names = FALSE
)

# STATISTICS TABLE FOR SUPPLEMENTARY FILE


stats_summary <- tibble(
  
  Analysis = c(
    "PERMANOVA",
    "PERMDISP"
  ),
  
  Comparison = c(
    "23-day exposure adult vs juvenile",
    "23-day exposure adult vs juvenile"
  ),
  
  Metric = c(
    "Community composition",
    "Beta dispersion"
  ),
  
  F = c(
    permanova_F,
    permdisp_F
  ),
  
  R2 = c(
    permanova_R2,
    NA_real_
  ),
  
  p_value = c(
    permanova_p,
    permdisp_p
  ),
  
  Permutations = c(
    9999,
    9999
  )
)


write.csv(
  stats_summary,
  "DevelopmentalStage_Statistical_Results.csv",
  row.names = FALSE
)


# COLOURS


stage_cols <- c(
  
  "23-day exposure – adult" =
    "#55A868",
  
  "23-day exposure – juvenile" =
    "#4C72B0"
  
)


# SIGNIFICANCE LABEL FOR PERMDISP FIGURE


sig_label <- case_when(
  
  permdisp_p < 0.001 ~ "***",
  
  permdisp_p < 0.01 ~ "**",
  
  permdisp_p < 0.05 ~ "*",
  
  TRUE ~ "ns"
)


# Position significance bracket automatically
y_max <- max(
  disp_df$Distance_to_Centroid,
  na.rm = TRUE
)

bracket_y <- y_max + 0.025

text_y <- y_max + 0.045


# BETA DISPERSION


p_disp <- ggplot(
  
  disp_df,
  
  aes(
    x = Stage,
    y = Distance_to_Centroid,
    fill = Stage
  )
  
) +
  
  geom_boxplot(
    outlier.shape = NA,
    width = 0.60,
    alpha = 0.80
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.65,
    size = 1.9
  ) +
  
  scale_fill_manual(
    values = stage_cols
  ) +
  

# Significance bracket

annotate(
  "segment",
  x = 1,
  xend = 2,
  y = bracket_y,
  yend = bracket_y
) +
  
  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = bracket_y - 0.012,
    yend = bracket_y
  ) +
  
  annotate(
    "segment",
    x = 2,
    xend = 2,
    y = bracket_y - 0.012,
    yend = bracket_y
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = text_y,
    label = sig_label,
    size = 5
  ) +
  
  labs(
    
    title =
      "Beta dispersion after 23 days of exposure",
    
    subtitle =
      sprintf(
        "PERMDISP: F = %.2f, p = %.4f",
        permdisp_F,
        permdisp_p
      ),
    
    x = NULL,
    
    y =
      "Distance to group centroid"
  ) +
  
  expand_limits(
    y = text_y + 0.02
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    
    legend.position = "none",
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    axis.text.x =
      element_text(
        size = 11
      )
  )


print(p_disp)


ggsave(
  
  "AdultVsJuvenile_BetaDispersion.png",
  
  p_disp,
  
  width = 7,
  
  height = 5.5,
  
  dpi = 300
  
)



#  PCoA

ord <- ordinate(
  ps,
  method = "PCoA",
  distance = d
)


# Percentage explained by first two axes
ve <- round(
  100 *
    ord$values$Relative_eig[1:2],
  1
)


# Extract sample coordinates
odf <- as.data.frame(
  ord$vectors[, 1:2]
)


colnames(odf) <- c(
  "PCoA1",
  "PCoA2"
)


odf$Sample <- rownames(
  odf
)


# Add developmental-stage information
odf <- left_join(
  
  odf,
  
  md %>%
    tibble::rownames_to_column(
      "Sample"
    ) %>%
    select(
      Sample,
      Stage
    ),
  
  by = "Sample"
  
)

# FIGURE B — PCoA


p_pcoa <- ggplot(
  
  odf,
  
  aes(
    x = PCoA1,
    y = PCoA2,
    color = Stage
  )
  
) +
  
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  
  stat_ellipse(
    aes(
      group = Stage
    ),
    type = "norm",
    linewidth = 0.8
  ) +
  
  scale_color_manual(
    values = stage_cols
  ) +
  
  labs(
    
    title =
      "Adult and juvenile microbiome composition after 23 days of exposure",
    
    subtitle =
      sprintf(
        "Bray–Curtis PCoA · PERMANOVA R² = %.3f, p = %.4f",
        permanova_R2,
        permanova_p
      ),
    
    x =
      paste0(
        "PCoA 1 (",
        ve[1],
        "%)"
      ),
    
    y =
      paste0(
        "PCoA 2 (",
        ve[2],
        "%)"
      ),
    
    color =
      "Developmental stage"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold"
      )
  )


print(p_pcoa)


ggsave(
  
  "AdultVsJuvenile_BrayCurtis_PCoA.png",
  
  p_pcoa,
  
  width = 7.5,
  
  height = 5.8,
  
  dpi = 300
  
)



#  BETADISPER ORDINATION/SPIDER PLOT



png(
  "AdultVsJuvenile_Betadisper_Ordination.png",
  width = 900,
  height = 700,
  res = 130
)

plot(
  bd,
  hull = FALSE,
  ellipse = TRUE,
  seg.col = "grey60",
  main =
    "Beta dispersion: adult vs juvenile after 23 days of exposure"
)

dev.off()



#FINAL MESSAGE



cat(
  paste0(
    "\nFiles written:\n",
    "  AdultVsJuvenile_BetaDispersion.png\n",
    "  AdultVsJuvenile_BrayCurtis_PCoA.png\n",
    "  AdultVsJuvenile_Betadisper_Ordination.png\n",
    "  DevelopmentalStage_BetaDispersion_Summary.csv\n",
    "  DevelopmentalStage_BetaDispersion_PerSample.csv\n",
    "  DevelopmentalStage_Statistical_Results.csv\n"
  )
)



# MICROPLASTIC EFFECT ON BETA DISPERSION WITHIN LIFE STAGE



run_mp_dispersion <- function(stage_name) {
  
  # samples belonging to this developmental stage
  ids <- rownames(md)[md$Stage == stage_name]
  
  # subset Bray-Curtis matrix
  d_stage <- as.dist(
    as.matrix(d)[ids, ids]
  )
  
  # metadata aligned to same samples
  md_stage <- md[ids, , drop = FALSE]
  
  # Convert raw plastic-content values to publication labels
  md_stage$Plastic_Group <- factor(
    ifelse(
      as.character(md_stage$Plastic_content) == "YES",
      "MP",
      "Control"
    ),
    levels = c("Control", "MP")
  )
  

  # PERMDISP
 
  
  bd_stage <- vegan::betadisper(
    d_stage,
    group = md_stage$Plastic_Group,
    type = "centroid"
  )
  
  set.seed(1)
  
  pt_stage <- permutest(
    bd_stage,
    permutations = 9999
  )
  

  # PERMANOVA
  
  
  set.seed(1)
  
  ad_stage <- vegan::adonis2(
    d_stage ~ Plastic_Group,
    data = md_stage,
    permutations = 9999
  )
  
  
  # Per-sample distances
 
  
  out_df <- data.frame(
    Sample = names(bd_stage$distances),
    Life_stage = stage_name,
    Plastic_Group =
      md_stage[names(bd_stage$distances), "Plastic_Group"],
    Distance_to_Centroid =
      as.numeric(bd_stage$distances)
  )
  

  # Statistics

  
  stats <- tibble(
    Life_stage = stage_name,
    
    PERMDISP_F =
      pt_stage$tab$F[1],
    
    PERMDISP_p =
      pt_stage$tab$`Pr(>F)`[1],
    
    PERMANOVA_R2 =
      ad_stage$R2[1],
    
    PERMANOVA_F =
      ad_stage$F[1],
    
    PERMANOVA_p =
      ad_stage$`Pr(>F)`[1],
    
    Control_mean =
      bd_stage$group.distances["Control"],
    
    MP_mean =
      bd_stage$group.distances["MP"]
  )
  
  list(
    dispersion = out_df,
    stats = stats
  )
}


# RUN ADULT AND JUVENILE ANALYSES


adult_mp <- run_mp_dispersion(
  "23-day exposure – adult"
)

juvenile_mp <- run_mp_dispersion(
  "23-day exposure – juvenile"
)


mp_disp_df <- bind_rows(
  adult_mp$dispersion,
  juvenile_mp$dispersion
)


mp_stats <- bind_rows(
  adult_mp$stats,
  juvenile_mp$stats
)



cat("CONTROL vs MP WITHIN DEVELOPMENTAL STAGE\n")


print(mp_stats)


write.csv(
  mp_stats,
  "DevelopmentalStage_MP_Control_Statistics.csv",
  row.names = FALSE
)




mp_disp_df$Life_stage_display <- factor(
  mp_disp_df$Life_stage,
  
  levels = c(
    "23-day exposure – adult",
    "23-day exposure – juvenile"
  ),
  
  labels = c(
    "Adult",
    "Juvenile"
  )
)




mp_stats <- mp_stats %>%
  mutate(
    
    Life_stage_display = factor(
      Life_stage,
      levels = c(
        "23-day exposure – adult",
        "23-day exposure – juvenile"
      ),
      labels = c(
        "Adult",
        "Juvenile"
      )
    ),
    
    significance = case_when(
      PERMDISP_p < 0.001 ~ "***",
      PERMDISP_p < 0.01  ~ "**",
      PERMDISP_p < 0.05  ~ "*",
      TRUE               ~ "ns"
    )
  )


# position annotation separately for each facet
ann_df <- mp_disp_df %>%
  group_by(Life_stage_display) %>%
  summarise(
    y = max(Distance_to_Centroid, na.rm = TRUE) + 0.04,
    .groups = "drop"
  ) %>%
  left_join(
    mp_stats %>%
      select(
        Life_stage_display,
        significance,
        PERMDISP_p
      ),
    by = "Life_stage_display"
  )



#FIGURE


plastic_cols <- c(
  "Control" = "#999999",
  "MP" = "#E67E22"
)


p_mp <- ggplot(
  mp_disp_df,
  aes(
    x = Plastic_Group,
    y = Distance_to_Centroid,
    fill = Plastic_Group
  )
) +
  
  geom_boxplot(
    outlier.shape = NA,
    width = 0.58,
    alpha = 0.8
  ) +
  
  geom_jitter(
    width = 0.12,
    alpha = 0.65,
    size = 1.8
  ) +
  
  facet_wrap(
    ~ Life_stage_display,
    nrow = 1
  ) +
  
  scale_fill_manual(
    values = plastic_cols
  ) +
  
  geom_text(
    data = ann_df,
    aes(
      x = 1.5,
      y = y,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 5
  ) +
  
  labs(
    title =
      "B. Microplastic effect on beta dispersion",
    
    subtitle =
      "Control vs MP within developmental stage",
    
    x = NULL,
    
    y =
      "Distance to group centroid"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    legend.position = "none",
    
    plot.title =
      element_text(face = "bold"),
    
    strip.text =
      element_text(face = "bold")
  )


print(p_mp)


ggsave(
  "Microplastic_BetaDispersion_Within_LifeStage.png",
  p_mp,
  width = 7,
  height = 5,
  dpi = 300
)



# ADULT vs JUVENILE DISPERSION


p_stage <- ggplot(
  disp_df,
  aes(
    x = Stage,
    y = Distance_to_Centroid,
    fill = Stage
  )
) +
  
  geom_boxplot(
    outlier.shape = NA,
    width = 0.60,
    alpha = 0.8
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.65,
    size = 1.9
  ) +
  
  scale_fill_manual(
    values = stage_cols
  ) +
  
  labs(
    title =
      "A. Beta dispersion between developmental stages",
    
    subtitle =
      sprintf(
        "PERMDISP: F = %.2f, p = %.4f",
        permdisp_F,
        permdisp_p
      ),
    
    x = NULL,
    
    y =
      "Distance to group centroid"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 10)
  )


# combinne images


library(patchwork)

fig_development <- p_stage + p_mp +
  
  plot_layout(
    widths = c(0.9, 1.3)
  ) +
  
  plot_annotation(
    title =
      "Developmental-stage differences in microbiome beta dispersion",
    
    theme = theme(
      plot.title =
        element_text(
          face = "bold",
          size = 17
        )
    )
  )


print(fig_development)


ggsave(
  "Figure_DevelopmentalStage_BetaDispersion.png",
  fig_development,
  width = 14,
  height = 6,
  dpi = 300
)




cat(
  "\nWrote:\n",
  "  Microplastic_BetaDispersion_Within_LifeStage.png\n",
  "  Figure_DevelopmentalStage_BetaDispersion.png\n",
  "  DevelopmentalStage_MP_Control_Statistics.csv\n"
)



# DEVELOPMENTAL-STAGE ANALYSIS OF MICROBIOME BETA DISPERSION



library(phyloseq)
library(vegan)
library(tidyverse)
library(patchwork)



setwd(
  "C:/Users/u0175886/Documents/Microbiome_Microplastic/Exposure_Microplastic"
)

pseq <- readRDS("daphnia_phyloseq.rds")




print(
  table(
    as.character(sample_data(pseq)$Exposure_Condition),
    useNA = "ifany"
  )
)



print(
  table(
    as.character(sample_data(pseq)$Plastic_content),
    useNA = "ifany"
  )
)



# KEEP ADULT AND JUVENILE SAMPLES AFTER 23 DAYS


# Raw labels in your metadata
keep_labels <- c(
  "Post_exposure",
  "Juvenile"
)

cond <- as.character(
  sample_data(pseq)$Exposure_Condition
)

keep <- cond %in% keep_labels

ps <- prune_samples(
  sample_names(pseq)[keep],
  pseq
)


# TAXON FILTERING


# Same ASV filtering rule used elsewhere:
# detected in >= 3 samples
# total abundance >= 10 reads

ps <- filter_taxa(
  ps,
  function(x) {
    sum(x > 0) >= 3 &&
      sum(x) >= 10
  },
  prune = TRUE
)


#  RELATIVE ABUNDANCE


ps <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)

# BRAY–CURTIS DISSIMILARITIES


d <- phyloseq::distance(
  ps,
  method = "bray"
)



# PREPARE METADATA


md <- data.frame(
  sample_data(ps),
  stringsAsFactors = FALSE
)

md <- md[
  labels(d),
  ,
  drop = FALSE
]


# Publication-friendly developmental-stage labels
md$Stage <- factor(
  as.character(md$Exposure_Condition),
  
  levels = c(
    "Post_exposure",
    "Juvenile"
  ),
  
  labels = c(
    "23-day exposure – adult",
    "23-day exposure – juvenile"
  )
)


print(
  table(md$Stage)
)


# Safety checks
stopifnot(
  !anyNA(md$Stage)
)

stopifnot(
  nlevels(md$Stage) == 2
)

stopifnot(
  identical(
    rownames(md),
    labels(d)
  )
)





# PERMANOVA


set.seed(1)

ad_stage <- vegan::adonis2(
  d ~ Stage,
  data = md,
  permutations = 9999
)

cat("\n============================================================\n")
cat("PERMANOVA — adult vs juvenile\n")
cat("============================================================\n")

print(ad_stage)


#  PERMDISP


bd_stage <- vegan::betadisper(
  d,
  group = md$Stage,
  type = "centroid"
)

set.seed(1)

pt_stage <- permutest(
  bd_stage,
  permutations = 9999
)



print(pt_stage)


cat("\nMean distance to group centroid:\n")

print(
  round(
    bd_stage$group.distances,
    4
  )
)


cat("\nMedian distance to group centroid:\n")

print(
  tapply(
    bd_stage$distances,
    md$Stage,
    median
  )
)



# STATISTICS



stage_permanova_R2 <- ad_stage$R2[1]
stage_permanova_F  <- ad_stage$F[1]
stage_permanova_p  <- ad_stage$`Pr(>F)`[1]

stage_permdisp_F <- pt_stage$tab$F[1]
stage_permdisp_p <- pt_stage$tab$`Pr(>F)`[1]


adult_mean <- bd_stage$group.distances[
  "23-day exposure – adult"
]

juvenile_mean <- bd_stage$group.distances[
  "23-day exposure – juvenile"
]





cat(
  sprintf(
    "PERMANOVA: R2 = %.4f, F = %.3f, p = %.4f\n",
    stage_permanova_R2,
    stage_permanova_F,
    stage_permanova_p
  )
)

cat(
  sprintf(
    "PERMDISP: F = %.3f, p = %.4f\n",
    stage_permdisp_F,
    stage_permdisp_p
  )
)

cat(
  sprintf(
    "Mean distance to centroid: Adult = %.4f, Juvenile = %.4f\n",
    adult_mean,
    juvenile_mean
  )
)



# PANEL A PER-SAMPLE DISPERSION TABLE


disp_df <- data.frame(
  
  Sample =
    names(
      bd_stage$distances
    ),
  
  Stage =
    md[
      names(bd_stage$distances),
      "Stage"
    ],
  
  Distance_to_Centroid =
    as.numeric(
      bd_stage$distances
    )
)


disp_df$Stage <- factor(
  
  disp_df$Stage,
  
  levels = c(
    "23-day exposure – adult",
    "23-day exposure – juvenile"
  )
)


# Short labels for plotting
disp_df$Stage_short <- factor(
  
  disp_df$Stage,
  
  levels = c(
    "23-day exposure – adult",
    "23-day exposure – juvenile"
  ),
  
  labels = c(
    "Adult",
    "Juvenile"
  )
)



# DISPERSION SUMMARY


disp_summary <- disp_df %>%
  
  group_by(Stage) %>%
  
  summarise(
    
    n = n(),
    
    mean_distance =
      mean(
        Distance_to_Centroid
      ),
    
    median_distance =
      median(
        Distance_to_Centroid
      ),
    
    sd_distance =
      sd(
        Distance_to_Centroid
      ),
    
    se_distance =
      sd_distance / sqrt(n),
    
    .groups = "drop"
  )



print(disp_summary)


write.csv(
  disp_summary,
  "DevelopmentalStage_BetaDispersion_Summary.csv",
  row.names = FALSE
)


write.csv(
  disp_df,
  "DevelopmentalStage_BetaDispersion_PerSample.csv",
  row.names = FALSE
)



# PANEL A STATISTICS TABLE


stage_stats <- tibble(
  
  Analysis = c(
    "PERMANOVA",
    "PERMDISP"
  ),
  
  Comparison = c(
    "Adult vs juvenile after 23-day exposure",
    "Adult vs juvenile after 23-day exposure"
  ),
  
  Metric = c(
    "Community composition",
    "Beta dispersion"
  ),
  
  F = c(
    stage_permanova_F,
    stage_permdisp_F
  ),
  
  R2 = c(
    stage_permanova_R2,
    NA_real_
  ),
  
  p_value = c(
    stage_permanova_p,
    stage_permdisp_p
  ),
  
  Permutations = c(
    9999,
    9999
  )
)


write.csv(
  stage_stats,
  "DevelopmentalStage_Statistical_Results.csv",
  row.names = FALSE
)




#  PANEL A SIGNIFICANCE LABEL



stage_sig <- case_when(
  
  stage_permdisp_p < 0.001 ~ "***",
  
  stage_permdisp_p < 0.01 ~ "**",
  
  stage_permdisp_p < 0.05 ~ "*",
  
  TRUE ~ "ns"
)


ymax_stage <- max(
  disp_df$Distance_to_Centroid,
  na.rm = TRUE
)

bracket_y_stage <- ymax_stage + 0.025
text_y_stage    <- ymax_stage + 0.045




# PANEL A FIGURE



stage_cols_short <- c(
  
  "Adult" =
    "#55A868",
  
  "Juvenile" =
    "#4C72B0"
)


p_stage <- ggplot(
  
  disp_df,
  
  aes(
    x = Stage_short,
    y = Distance_to_Centroid,
    fill = Stage_short
  )
  
) +
  
  geom_boxplot(
    outlier.shape = NA,
    width = 0.60,
    alpha = 0.80
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.65,
    size = 1.9
  ) +
  
  scale_fill_manual(
    values = stage_cols_short
  ) +
  
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = bracket_y_stage,
    yend = bracket_y_stage
  ) +
  
  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = bracket_y_stage - 0.012,
    yend = bracket_y_stage
  ) +
  
  annotate(
    "segment",
    x = 2,
    xend = 2,
    y = bracket_y_stage - 0.012,
    yend = bracket_y_stage
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = text_y_stage,
    label = stage_sig,
    size = 5.5
  ) +
  
  labs(
    
    title =
      "A. Beta dispersion between developmental stages",
    
    subtitle =
      sprintf(
        "PERMDISP: F = %.2f, p = %.4f",
        stage_permdisp_F,
        stage_permdisp_p
      ),
    
    x = NULL,
    
    y =
      "Distance to group centroid"
  ) +
  
  expand_limits(
    y = text_y_stage + 0.02
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    
    legend.position =
      "none",
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    axis.text.x =
      element_text(
        size = 11
      )
  )



# PANEL B ANALYSIS FUNCTION
#     CONTROL vs MP WITHIN EACH DEVELOPMENTAL STAGE



run_mp_dispersion <- function(stage_name) {
  
  
 
 
  # Samples in this life stage
  
  
  
  ids <- rownames(md)[
    md$Stage == stage_name
  ]
  
 
 
  # Stage-specific Bray–Curtis matrix


  
  d_stage <- as.dist(
    as.matrix(d)[
      ids,
      ids
    ]
  )
  
  
  
  
  # Stage-specific metadata
  
  
  md_stage <- md[
    ids,
    ,
    drop = FALSE
  ]
  
  
 
 
  # Keep only known Control / MP samples
  
  
  plastic_raw <- as.character(
    md_stage$Plastic_content
  )
  
  known <- plastic_raw %in% c(
    "YES",
    "NO"
  )
  
  md_stage <- md_stage[
    known,
    ,
    drop = FALSE
  ]
  
  
  ids2 <- rownames(
    md_stage
  )
  
  
  d_stage <- as.dist(
    as.matrix(d_stage)[
      ids2,
      ids2
    ]
  )
  
  
 
 
  md_stage$Plastic_Group <- factor(
    
    ifelse(
      as.character(
        md_stage$Plastic_content
      ) == "YES",
      
      "MP",
      
      "Control"
    ),
    
    levels = c(
      "Control",
      "MP"
    )
  )
  
  
  
  
  cat(
    "Stage: ",
    stage_name,
    "\n",
    sep = ""
  )
  
  cat(
    "Control / MP sample sizes:\n"
  )
  
  print(
    table(
      md_stage$Plastic_Group
    )
  )
  
  # PERMDISP
  
  
  
  bd_mp <- vegan::betadisper(
    
    d_stage,
    
    group =
      md_stage$Plastic_Group,
    
    type =
      "centroid"
  )
  
  
  set.seed(1)
  
  pt_mp <- permutest(
    
    bd_mp,
    
    permutations =
      9999
  )
  
  
 
 
  # PERMANOVA
  
  
  
  set.seed(1)
  
  ad_mp <- vegan::adonis2(
    
    d_stage ~ Plastic_Group,
    
    data =
      md_stage,
    
    permutations =
      9999
  )
  
  
  cat(
    "\nPERMDISP:\n"
  )
  
  print(
    pt_mp
  )
  
  
  cat(
    "\nPERMANOVA:\n"
  )
  
  print(
    ad_mp
  )
  
  
  cat(
    "\nMean distance to centroid:\n"
  )
  
  print(
    round(
      bd_mp$group.distances,
      4
    )
  )
  
  
  
  
  # Per-sample table
  
  out_df <- data.frame(
    
    Sample =
      names(
        bd_mp$distances
      ),
    
    Life_stage =
      stage_name,
    
    Plastic_Group =
      md_stage[
        names(
          bd_mp$distances
        ),
        "Plastic_Group"
      ],
    
    Distance_to_Centroid =
      as.numeric(
        bd_mp$distances
      )
  )
  
  
  
  # Statistics table
  
  
  stats <- tibble(
    
    Life_stage =
      stage_name,
    
    PERMDISP_F =
      pt_mp$tab$F[1],
    
    PERMDISP_p =
      pt_mp$tab$`Pr(>F)`[1],
    
    PERMANOVA_R2 =
      ad_mp$R2[1],
    
    PERMANOVA_F =
      ad_mp$F[1],
    
    PERMANOVA_p =
      ad_mp$`Pr(>F)`[1],
    
    Control_mean =
      unname(
        bd_mp$group.distances[
          "Control"
        ]
      ),
    
    MP_mean =
      unname(
        bd_mp$group.distances[
          "MP"
        ]
      )
  )
  
  
  return(
    list(
      dispersion =
        out_df,
      
      stats =
        stats
    )
  )
}


# 15. RUN CONTROL vs MP ANALYSES

adult_mp <- run_mp_dispersion(
  "23-day exposure – adult"
)

juvenile_mp <- run_mp_dispersion(
  "23-day exposure – juvenile"
)


mp_disp_df <- bind_rows(
  adult_mp$dispersion,
  juvenile_mp$dispersion
)


mp_stats <- bind_rows(
  adult_mp$stats,
  juvenile_mp$stats
)




print(
  mp_stats
)


write.csv(
  mp_stats,
  "DevelopmentalStage_MP_Control_Statistics.csv",
  row.names = FALSE
)


write.csv(
  mp_disp_df,
  "DevelopmentalStage_MP_Control_PerSample.csv",
  row.names = FALSE
)


#  CLEAN LABELS


mp_disp_df$Life_stage_display <- factor(
  
  mp_disp_df$Life_stage,
  
  levels = c(
    "23-day exposure – adult",
    "23-day exposure – juvenile"
  ),
  
  labels = c(
    "Adult",
    "Juvenile"
  )
)


mp_stats <- mp_stats %>%
  
  mutate(
    
    Life_stage_display = factor(
      
      Life_stage,
      
      levels = c(
        "23-day exposure – adult",
        "23-day exposure – juvenile"
      ),
      
      labels = c(
        "Adult",
        "Juvenile"
      )
    ),
    
    significance = case_when(
      
      PERMDISP_p < 0.001 ~ "***",
      
      PERMDISP_p < 0.01 ~ "**",
      
      PERMDISP_p < 0.05 ~ "*",
      
      TRUE ~ "ns"
    )
  )


#PANEL B SIGNIFICANCE POSITIONS

ann_df <- mp_disp_df %>%
  
  group_by(
    Life_stage_display
  ) %>%
  
  summarise(
    
    y =
      max(
        Distance_to_Centroid,
        na.rm = TRUE
      ) + 0.045,
    
    .groups =
      "drop"
  ) %>%
  
  left_join(
    
    mp_stats %>%
      
      select(
        Life_stage_display,
        significance,
        PERMDISP_p
      ),
    
    by =
      "Life_stage_display"
  )



#B FIGURE


plastic_cols <- c(
  
  "Control" =
    "#999999",
  
  "MP" =
    "#E67E22"
)


p_mp <- ggplot(
  
  mp_disp_df,
  
  aes(
    x = Plastic_Group,
    y = Distance_to_Centroid,
    fill = Plastic_Group
  )
  
) +
  
  geom_boxplot(
    outlier.shape = NA,
    width = 0.58,
    alpha = 0.80
  ) +
  
  geom_jitter(
    width = 0.12,
    alpha = 0.65,
    size = 1.8
  ) +
  
  facet_wrap(
    ~ Life_stage_display,
    nrow = 1
  ) +
  
  scale_fill_manual(
    values = plastic_cols
  ) +
  
  geom_text(
    
    data =
      ann_df,
    
    aes(
      x = 1.5,
      y = y,
      label = significance
    ),
    
    inherit.aes =
      FALSE,
    
    size =
      5.5
  ) +
  
  labs(
    
    title =
      "B. Microplastic effect on beta dispersion",
    
    subtitle =
      "Control vs MP within developmental stage",
    
    x =
      NULL,
    
    y =
      "Distance to group centroid"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    
    legend.position =
      "none",
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    strip.text =
      element_text(
        face = "bold",
        size = 12
      )
  )


#combine image A and B

fig_development <- p_stage + p_mp +
  
  plot_layout(
    widths = c(
      0.95,
      1.35
    )
  ) +
  
  plot_annotation(
    
    title =
      "Developmental-stage differences in microbiome beta dispersion",
    
    theme =
      theme(
        
        plot.title =
          element_text(
            face = "bold",
            size = 17
          )
      )
  )


print(
  fig_development
)


# SAVE COMBINED FIGURE


ggsave(
  
  "Figure_DevelopmentalStage_BetaDispersion.png",
  
  fig_development,
  
  width = 14,
  
  height = 6,
  
  dpi = 300
)



# SAVE PANELS SEPARATELY


ggsave(
  
  "Figure_DevelopmentalStage_BetaDispersion_PanelA.png",
  
  p_stage,
  
  width = 6.5,
  
  height = 5.5,
  
  dpi = 300
)


ggsave(
  
  "Figure_DevelopmentalStage_BetaDispersion_PanelB.png",
  
  p_mp,
  
  width = 8,
  
  height = 5.5,
  
  dpi = 300
)


#  PCoA — ADULT vs JUVENILE COMPOSITION


ord <- ordinate(
  ps,
  method = "PCoA",
  distance = d
)


ve <- round(
  100 *
    ord$values$Relative_eig[
      1:2
    ],
  1
)


odf <- as.data.frame(
  ord$vectors[
    ,
    1:2
  ]
)


colnames(
  odf
) <- c(
  "PCoA1",
  "PCoA2"
)


odf$Sample <- rownames(
  odf
)


odf <- left_join(
  
  odf,
  
  md %>%
    
    tibble::rownames_to_column(
      "Sample"
    ) %>%
    
    select(
      Sample,
      Stage
    ),
  
  by =
    "Sample"
)


stage_cols_full <- c(
  
  "23-day exposure – adult" =
    "#55A868",
  
  "23-day exposure – juvenile" =
    "#4C72B0"
)


p_pcoa <- ggplot(
  
  odf,
  
  aes(
    x = PCoA1,
    y = PCoA2,
    color = Stage
  )
  
) +
  
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  
  stat_ellipse(
    aes(
      group = Stage
    ),
    type = "norm",
    linewidth = 0.8
  ) +
  
  scale_color_manual(
    values = stage_cols_full
  ) +
  
  labs(
    
    title =
      "Adult and juvenile microbiome composition after 23 days of exposure",
    
    subtitle =
      sprintf(
        "Bray–Curtis PCoA · PERMANOVA R² = %.3f, p = %.4f",
        stage_permanova_R2,
        stage_permanova_p
      ),
    
    x =
      paste0(
        "PCoA 1 (",
        ve[1],
        "%)"
      ),
    
    y =
      paste0(
        "PCoA 2 (",
        ve[2],
        "%)"
      ),
    
    color =
      "Developmental stage"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold"
      )
  )


print(
  p_pcoa
)


ggsave(
  
  "AdultVsJuvenile_BrayCurtis_PCoA.png",
  
  p_pcoa,
  
  width = 7.5,
  
  height = 5.8,
  
  dpi = 300
)


# OPTIONAL BETADISPER ORDINATION

png(
  
  "AdultVsJuvenile_Betadisper_Ordination.png",
  
  width = 900,
  
  height = 700,
  
  res = 130
)


plot(
  
  bd_stage,
  
  hull = FALSE,
  
  ellipse = TRUE,
  
  seg.col = "grey60",
  
  main =
    "Beta dispersion: adult vs juvenile after 23 days of exposure"
)


dev.off()



cat(
  paste0(
    "\nFiles written:\n",
    "  Figure_DevelopmentalStage_BetaDispersion.png\n",
    "  Figure_DevelopmentalStage_BetaDispersion_PanelA.png\n",
    "  Figure_DevelopmentalStage_BetaDispersion_PanelB.png\n",
    "  AdultVsJuvenile_BrayCurtis_PCoA.png\n",
    "  AdultVsJuvenile_Betadisper_Ordination.png\n",
    "  DevelopmentalStage_BetaDispersion_Summary.csv\n",
    "  DevelopmentalStage_BetaDispersion_PerSample.csv\n",
    "  DevelopmentalStage_Statistical_Results.csv\n",
    "  DevelopmentalStage_MP_Control_Statistics.csv\n",
    "  DevelopmentalStage_MP_Control_PerSample.csv\n"
  )
)