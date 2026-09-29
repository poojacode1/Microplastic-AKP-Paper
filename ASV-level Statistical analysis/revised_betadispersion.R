library(phyloseq)
library(vegan)
library(tidyverse)

setwd("C://Users//u0175886/Documents//Microbiome_Microplastic//Exposure_Microplastic")
pseq <- readRDS("daphnia_phyloseq.rds")


# ASV-level community dispersion across ADULT exposure stages


# checking the data
cat("\n=== Exposure_Condition (before filter) ===\n")
print(table(as.character(sample_data(pseq)$Exposure_Condition), useNA = "ifany"))
cat("\n=== Plastic_content (before filter) ===\n")
print(table(as.character(sample_data(pseq)$Plastic_content), useNA = "ifany"))

# filtering
keep_adult <- c("Pre_exposure", "Post_exposure", "Post_long_exposure")

# display labels only: raw metadata values above are kept unchanged for filtering/statistics
group_labels <- c(
  "Pre_exposure"       = "Pre-exposure",
  "Post_exposure"      = "23-day exposure – adult",
  "Post_long_exposure" = "Long-term survivors"
)
display_levels <- unname(group_labels[keep_adult])
md0  <- data.frame(sample_data(pseq), stringsAsFactors = FALSE)
keep <- as.character(md0$Exposure_Condition) %in% keep_adult &
  as.character(md0$Plastic_content)    %in% c("YES", "NO")
pseq <- prune_samples(sample_names(pseq)[keep], pseq)

sample_data(pseq)$Exposure_Condition <- factor(
  sample_data(pseq)$Exposure_Condition,
  levels = keep_adult
)

cat("\n=== After filter (adults, known plastic status) ===\n")
print(table(sample_data(pseq)$Exposure_Condition))

# taxon filtering 
pseq_filt <- filter_taxa(
  pseq,
  function(x) sum(x > 0) >= 3 && sum(x) >= 10,
  prune = TRUE
)

#relative abundance for Bray-Curtis
pseq_rel <- transform_sample_counts(pseq_filt, function(x) x / sum(x))

#Bray-Curtis distance 
bray_dist <- phyloseq::distance(pseq_rel, method = "bray")

#aligned metadata 
meta_df <- data.frame(sample_data(pseq_rel), stringsAsFactors = FALSE)
meta_df <- meta_df[labels(bray_dist), , drop = FALSE]
meta_df$Exposure_Condition <- factor(
  meta_df$Exposure_Condition,
  levels = keep_adult
)

meta_df$Exposure_Label <- factor(
  unname(group_labels[as.character(meta_df$Exposure_Condition)]),
  levels = display_levels
)

# checks 
cat("\nGroup sizes for testing:\n"); print(table(meta_df$Exposure_Condition))
stopifnot(!anyNA(meta_df$Exposure_Condition))   # guard: no NA reaches adonis2
stopifnot(identical(rownames(meta_df), labels(bray_dist)))

# PERMANOVA

set.seed(1)
perm_res <- vegan::adonis2(
  bray_dist ~ Exposure_Condition,
  data = meta_df,
  permutations = 9999
)
cat("\n-- PERMANOVA (Bray-Curtis, adult exposure stages) --\n")
print(perm_res)


# PERMDISP (beta-dispersion)

bd <- betadisper(bray_dist, group = meta_df$Exposure_Condition, type = "centroid")
cat("\n-- betadisper ANOVA --\n");   print(anova(bd))
cat("\n-- betadisper permutest (pairwise) --\n")
print(permutest(bd, pairwise = TRUE, permutations = 9999))
cat("\nGroup mean distance to centroid:\n"); print(round(bd$group.distances, 4))

#per-sample dispersion table + summary 
disp_df <- data.frame(
  Sample = names(bd$distances),
  Exposure_Condition = meta_df[names(bd$distances), "Exposure_Condition"],
  Distance_to_Centroid = as.numeric(bd$distances)
)
disp_df$Exposure_Condition <- factor(disp_df$Exposure_Condition, levels = keep_adult)
disp_df$Exposure_Label <- factor(
  unname(group_labels[as.character(disp_df$Exposure_Condition)]),
  levels = display_levels
)

disp_summary <- disp_df %>%
  group_by(Exposure_Label) %>%
  summarise(
    n = n(),
    mean_distance   = mean(Distance_to_Centroid),
    median_distance = median(Distance_to_Centroid),
    sd_distance     = sd(Distance_to_Centroid),
    se_distance     = sd_distance / sqrt(n),
    .groups = "drop"
  )
print(disp_summary)
write.csv(disp_summary, "BetaDispersion_Summary_Adults.csv", row.names = FALSE)


# FIGURES

pal_cond <- c(
  "Pre-exposure"             = "#999999",
  "23-day exposure – adult"  = "#E69F00",
  "Long-term survivors"      = "#56B4E9"
)

#dispersion boxplot 
p_disp_box <- ggplot(disp_df, aes(Exposure_Label, Distance_to_Centroid, fill = Exposure_Label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.65, size = 2) +
  scale_fill_manual(values = pal_cond) +
  labs(title = "Community dispersion across adult exposure stages",
       subtitle = "Distance to group centroid (Bray-Curtis beta-dispersion)",
       x = NULL, y = "Distance to centroid") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))
print(p_disp_box)
ggsave("Adult_Dispersion_Boxplot.png", p_disp_box, width = 7.5, height = 5.5, dpi = 300)

#mean +/- SE bar 
p_disp_mean <- ggplot(disp_summary, aes(Exposure_Label, mean_distance, fill = Exposure_Label)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_errorbar(aes(ymin = mean_distance - se_distance, ymax = mean_distance + se_distance),
                width = 0.15, linewidth = 0.7) +
  geom_text(aes(label = round(mean_distance, 3)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = pal_cond) +
  labs(title = "Mean community dispersion by adult exposure stage",
       subtitle = "Higher values indicate greater between-sample heterogeneity",
       x = NULL, y = "Mean distance to centroid") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))
print(p_disp_mean)
ggsave("Adult_Dispersion_MeanSE.png", p_disp_mean, width = 7, height = 5.2, dpi = 300)

# PCoA 
ord <- ordinate(pseq_rel, method = "PCoA", distance = bray_dist)
ve  <- round(100 * ord$values$Relative_eig[1:2], 1)
ord_df <- as.data.frame(ord$vectors[, 1:2]) %>% rownames_to_column("Sample")
colnames(ord_df)[2:3] <- c("PCoA1", "PCoA2")
ord_df <- left_join(ord_df, rownames_to_column(meta_df, "Sample"), by = "Sample")

p_pcoa <- ggplot(ord_df, aes(PCoA1, PCoA2, color = Exposure_Label)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(aes(group = Exposure_Label), type = "norm", linewidth = 0.8) +
  scale_color_manual(values = pal_cond) +
  labs(title = "PCoA of adult Daphnia microbiome",
       subtitle = "Bray-Curtis distance",
       x = paste0("PCoA 1 (", ve[1], "%)"),
       y = paste0("PCoA 2 (", ve[2], "%)"),
       color = "Condition") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
print(p_pcoa)
ggsave("Adult_PCoA_clean.png", p_pcoa, width = 7.5, height = 5.8, dpi = 300)

# betadisper spider plot
bd_plot <- betadisper(
  bray_dist,
  group = meta_df$Exposure_Label,
  type = "centroid"
)

png("Adult_Betadisper_Spiderplot.png", width = 900, height = 700, res = 130)
plot(bd_plot, hull = FALSE, ellipse = TRUE, seg.col = "grey60",
     main = "Beta-dispersion by adult exposure stage")
dev.off()


# Pairwise PERMANOVA (composition)

pairwise_adonis <- function(dist_mat, grouping, nperm = 9999) {
  grouping <- factor(grouping)
  levs <- levels(grouping)
  out <- list()
  for (i in 1:(length(levs) - 1)) {
    for (j in (i + 1):length(levs)) {
      keep <- grouping %in% c(levs[i], levs[j])
      sub_dist <- as.dist(as.matrix(dist_mat)[keep, keep])
      sub_meta <- data.frame(group = droplevels(grouping[keep]))
      set.seed(1)
      res <- adonis2(sub_dist ~ group, data = sub_meta, permutations = nperm)
      out[[paste(levs[i], "vs", levs[j])]] <- res
    }
  }
  out
}

cat("\n-- Pairwise PERMANOVA (composition) --\n")
pw <- pairwise_adonis(bray_dist, meta_df$Exposure_Label, nperm = 9999)
print(pw)

cat("\nDone. Wrote: Adult_Dispersion_Boxplot.png, Adult_Dispersion_MeanSE.png, ",
    "Adult_PCoA_clean.png, Adult_Betadisper_Spiderplot.png, BetaDispersion_Summary_Adults.csv\n")



# Pairwise Bray-Curtis dissimilarity distributions


dist_mat <- as.matrix(bray_dist)
samples  <- rownames(meta_df)

# build the pairwise table 
pair_df <- expand.grid(Sample1 = samples, Sample2 = samples,
                       stringsAsFactors = FALSE) %>%
  filter(Sample1 < Sample2) %>%
  mutate(
    BrayCurtis = map2_dbl(Sample1, Sample2, ~ dist_mat[.x, .y]),
    Group1 = meta_df$Exposure_Condition[match(Sample1, rownames(meta_df))],
    Group2 = meta_df$Exposure_Condition[match(Sample2, rownames(meta_df))]
  ) %>%
  mutate(
    Comparison = case_when(
      Group1 == "Pre_exposure"       & Group2 == "Pre_exposure"       ~ "Pre-exposure vs Pre-exposure",
      Group1 == "Post_exposure"      & Group2 == "Post_exposure"      ~ "23-day exposure – adult vs 23-day exposure – adult",
      Group1 == "Post_long_exposure" & Group2 == "Post_long_exposure" ~ "Long-term survivors vs Long-term survivors",
      (Group1 == "Pre_exposure"       & Group2 == "Post_exposure")      |
        (Group1 == "Post_exposure"      & Group2 == "Pre_exposure")       ~ "Pre-exposure vs 23-day exposure – adult",
      (Group1 == "Pre_exposure"       & Group2 == "Post_long_exposure") |
        (Group1 == "Post_long_exposure" & Group2 == "Pre_exposure")       ~ "Pre-exposure vs Long-term survivors",
      (Group1 == "Post_exposure"      & Group2 == "Post_long_exposure") |
        (Group1 == "Post_long_exposure" & Group2 == "Post_exposure")      ~ "23-day exposure – adult vs Long-term survivors",
      TRUE ~ "Other"
    )
  )

pair_df$Comparison <- factor(
  pair_df$Comparison,
  levels = c(
    "Pre-exposure vs Pre-exposure",
    "23-day exposure – adult vs 23-day exposure – adult",
    "Long-term survivors vs Long-term survivors",
    "Pre-exposure vs 23-day exposure – adult",
    "Pre-exposure vs Long-term survivors",
    "23-day exposure – adult vs Long-term survivors"
  )
)

write.csv(pair_df, "BrayCurtis_pairwise_distributions_adults.csv", row.names = FALSE)

#summary
pair_summary <- pair_df %>%
  group_by(Comparison) %>%
  summarise(
    n = n(),
    mean_Bray   = mean(BrayCurtis, na.rm = TRUE),
    median_Bray = median(BrayCurtis, na.rm = TRUE),
    sd_Bray     = sd(BrayCurtis, na.rm = TRUE),
    .groups = "drop"
  )
print(pair_summary)
write.csv(pair_summary, "BrayCurtis_pairwise_summary_adults.csv", row.names = FALSE)

#between-stage figure 
pair_between <- pair_df %>%
  filter(Comparison %in% c(
    "Pre-exposure vs 23-day exposure – adult",
    "Pre-exposure vs Long-term survivors",
    "23-day exposure – adult vs Long-term survivors"
  )) %>%
  mutate(
    Comparison = factor(
      Comparison,
      levels = c(
        "Pre-exposure vs 23-day exposure – adult",
        "Pre-exposure vs Long-term survivors",
        "23-day exposure – adult vs Long-term survivors"
      )
    )
  )

pal_cmp <- c(
  "Pre-exposure vs 23-day exposure – adult" = "#E69F00",
  "Pre-exposure vs Long-term survivors" = "#56B4E9",
  "23-day exposure – adult vs Long-term survivors" = "#009E73"
)

p_bray_between <- ggplot(pair_between, aes(Comparison, BrayCurtis, fill = Comparison)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = pal_cmp) +
  labs(title = "Bray-Curtis dissimilarity between adult exposure stages",
       x = NULL, y = "Bray-Curtis dissimilarity") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))
print(p_bray_between)
ggsave("Adult_BrayCurtis_Between_Stages.png", p_bray_between, width = 8, height = 5.5, dpi = 300)





# Pairwise Bray-Curtis dissimilarity


pair_between_plot <- pair_df %>%
  mutate(
    Comparison_display = case_when(
      
      # Pre vs 23-day exposure
      (Group1 == "Pre_exposure" & Group2 == "Post_exposure") |
        (Group1 == "Post_exposure" & Group2 == "Pre_exposure") ~
        "Pre-exposure vs\n23-day exposure",
      
      # Pre vs long-term survivors
      (Group1 == "Pre_exposure" & Group2 == "Post_long_exposure") |
        (Group1 == "Post_long_exposure" & Group2 == "Pre_exposure") ~
        "Pre-exposure vs\nlong-term survivors",
      
      # 23-day exposure vs long-term survivors
      (Group1 == "Post_exposure" & Group2 == "Post_long_exposure") |
        (Group1 == "Post_long_exposure" & Group2 == "Post_exposure") ~
        "23-day exposure vs\nlong-term survivors",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Comparison_display))


# Fix order
pair_between_plot$Comparison_display <- factor(
  pair_between_plot$Comparison_display,
  levels = c(
    "Pre-exposure vs\n23-day exposure",
    "Pre-exposure vs\nlong-term survivors",
    "23-day exposure vs\nlong-term survivors"
  )
)


# Colours
pal_cmp_new <- c(
  "Pre-exposure vs\n23-day exposure" = "#E69F00",
  "Pre-exposure vs\nlong-term survivors" = "#56B4E9",
  "23-day exposure vs\nlong-term survivors" = "#009E73"
)


p_bray_panel <- ggplot(
  pair_between_plot,
  aes(
    x = Comparison_display,
    y = BrayCurtis,
    fill = Comparison_display
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.70
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.90
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3,
    fill = "white"
  ) +
  scale_fill_manual(values = pal_cmp_new) +
  labs(
    title = "B. Pairwise Bray–Curtis dissimilarity",
    x = NULL,
    y = "Bray–Curtis dissimilarity"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 10
    )
  )



# PCoA


p_pcoa_panel <- p_pcoa +
  labs(
    title = "A. Bray–Curtis PCoA",
    subtitle = NULL,
    color = "Exposure stage"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    legend.position = "bottom"
  )


# COMBINE

library(patchwork)

fig6 <- p_pcoa_panel + p_bray_panel +
  plot_layout(
    widths = c(1, 1.15)
  ) +
  plot_annotation(
    title = "Adult microbiome composition across exposure stages",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 17
      )
    )
  )

print(fig6)


ggsave(
  "Figure6_Adult_BrayCurtis_Combined.png",
  fig6,
  width = 14,
  height = 6.5,
  dpi = 300
)

# SENSITIVITY ANALYSIS — unequal sample sizes
# Repeated balanced subsampling to n = 9 per exposure stage


set.seed(123)

N_ITER   <- 1000
TARGET_N <- min(table(meta_df$Exposure_Condition))


cat("BALANCED SUBSAMPLING SENSITIVITY ANALYSIS\n")
cat("Iterations:", N_ITER, "\n")
cat("Samples per group:", TARGET_N, "\n")


# Sample IDs belonging to each exposure stage
ids_pre <- rownames(meta_df)[
  meta_df$Exposure_Condition == "Pre_exposure"
]

ids_post <- rownames(meta_df)[
  meta_df$Exposure_Condition == "Post_exposure"
]

ids_long <- rownames(meta_df)[
  meta_df$Exposure_Condition == "Post_long_exposure"
]

# Convert once, rather than repeatedly inside loop
bray_mat <- as.matrix(bray_dist)

# Store results
sens_results <- data.frame(
  iteration = seq_len(N_ITER),
  
  mean_pre  = NA_real_,
  mean_post = NA_real_,
  mean_long = NA_real_,
  
  post_minus_pre  = NA_real_,
  post_minus_long = NA_real_,
  long_minus_pre  = NA_real_,
  
  overall_F = NA_real_,
  overall_p = NA_real_
)


for (i in seq_len(N_ITER)) {
  
  
  # Random balanced sample: n = 9 from every group

  
  sel_pre <- sample(
    ids_pre,
    TARGET_N,
    replace = FALSE
  )
  
  sel_post <- sample(
    ids_post,
    TARGET_N,
    replace = FALSE
  )
  
  sel_long <- sample(
    ids_long,
    TARGET_N,
    replace = FALSE
  )
  
  sel_ids <- c(
    sel_pre,
    sel_post,
    sel_long
  )
  
 
  # Subset Bray-Curtis distance matrix
 
  
  sub_bray <- as.dist(
    bray_mat[
      sel_ids,
      sel_ids,
      drop = FALSE
    ]
  )
  
  # Matching metadata/group vector
  sub_group <- factor(
    meta_df[
      sel_ids,
      "Exposure_Condition"
    ],
    levels = keep_adult
  )
  

  # PERMDISP
  # bias.adjust = TRUE is useful here because n is small
  
  
  sub_bd <- vegan::betadisper(
    sub_bray,
    group = sub_group,
    type = "centroid",
    bias.adjust = TRUE
  )
  
  sub_test <- vegan::permutest(
    sub_bd,
    permutations = 999
  )
  
  # Mean distance to centroid
  means <- tapply(
    sub_bd$distances,
    sub_group,
    mean
  )
  
  sens_results$mean_pre[i] <-
    means["Pre_exposure"]
  
  sens_results$mean_post[i] <-
    means["Post_exposure"]
  
  sens_results$mean_long[i] <-
    means["Post_long_exposure"]
  
  # Difference in dispersion
  sens_results$post_minus_pre[i] <-
    means["Post_exposure"] -
    means["Pre_exposure"]
  
  sens_results$post_minus_long[i] <-
    means["Post_exposure"] -
    means["Post_long_exposure"]
  
  sens_results$long_minus_pre[i] <-
    means["Post_long_exposure"] -
    means["Pre_exposure"]
  
  # Overall PERMDISP
  sens_results$overall_F[i] <-
    sub_test$tab[1, "F"]
  
  sens_results$overall_p[i] <-
    sub_test$tab[1, "Pr(>F)"]
  
  if (i %% 100 == 0) {
    cat("Completed", i, "of", N_ITER, "iterations\n")
  }
}


=
# SUMMARY


cat("ROBUSTNESS RESULTS\n")


cat("\nMedian balanced mean distance to centroid:\n")

cat(
  "Pre-exposure:",
  round(median(sens_results$mean_pre), 3),
  "\n"
)

cat(
  "23-day exposure adults:",
  round(median(sens_results$mean_post), 3),
  "\n"
)

cat(
  "Long-term survivors:",
  round(median(sens_results$mean_long), 3),
  "\n"
)


# How consistently does post-exposure remain more dispersed?

post_gt_pre <- mean(
  sens_results$post_minus_pre > 0
) * 100

post_gt_long <- mean(
  sens_results$post_minus_long > 0
) * 100

cat("\nDirection consistency:\n")

cat(
  "23-day adults > Pre-exposure:",
  round(post_gt_pre, 1),
  "% of iterations\n"
)

cat(
  "23-day adults > Long-term survivors:",
  round(post_gt_long, 1),
  "% of iterations\n"
)


# 95% empirical intervals
cat("\nPost-exposure minus Pre-exposure:\n")

print(
  round(
    quantile(
      sens_results$post_minus_pre,
      c(0.025, 0.50, 0.975)
    ),
    4
  )
)


cat("\nPost-exposure minus Long-term survivors:\n")

print(
  round(
    quantile(
      sens_results$post_minus_long,
      c(0.025, 0.50, 0.975)
    ),
    4
  )
)


cat("\nOverall PERMDISP significant in:\n")

cat(
  round(
    mean(sens_results$overall_p < 0.05) * 100,
    1
  ),
  "% of balanced iterations\n"
)



# SAVE


write.csv(
  sens_results,
  "Adult_PERMDISP_Balanced_Sensitivity.csv",
  row.names = FALSE
)



# FIGURE — distribution of balanced dispersion estimates


sens_plot_df <- sens_results %>%
  
  select(
    iteration,
    mean_pre,
    mean_post,
    mean_long
  ) %>%
  
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "Group",
    values_to = "Mean_Distance"
  ) %>%
  
  mutate(
    Group = recode(
      Group,
      mean_pre  = "Pre-exposure",
      mean_post = "23-day exposure – adult",
      mean_long = "Long-term survivors"
    ),
    
    Group = factor(
      Group,
      levels = display_levels
    )
  )


p_sensitivity <- ggplot(
  sens_plot_df,
  aes(
    x = Group,
    y = Mean_Distance,
    fill = Group
  )
) +
  
  geom_violin(
    trim = FALSE,
    alpha = 0.65
  ) +
  
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.9
  ) +
  
  scale_fill_manual(
    values = pal_cond
  ) +
  
  labs(
    title = "Balanced-subsampling sensitivity analysis",
    subtitle = paste0(
      N_ITER,
      " iterations; n = ",
      TARGET_N,
      " per exposure stage"
    ),
    x = NULL,
    y = "Mean distance to centroid"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )


print(p_sensitivity)

ggsave(
  "Adult_PERMDISP_Balanced_Sensitivity.png",
  p_sensitivity,
  width = 7.5,
  height = 5.5,
  dpi = 300
)


cat("\nSaved:\n")
cat("Adult_PERMDISP_Balanced_Sensitivity.csv\n")
cat("Adult_PERMDISP_Balanced_Sensitivity.png\n")





# sensitivity-analysis figure
# Median + central 95% subsampling interval


sens_summary <- sens_results %>%
  select(mean_pre, mean_post, mean_long) %>%
  pivot_longer(
    everything(),
    names_to = "Group",
    values_to = "Mean_Distance"
  ) %>%
  mutate(
    Group = recode(
      Group,
      mean_pre  = "Pre-exposure",
      mean_post = "23-day exposure – adult",
      mean_long = "Long-term survivors"
    ),
    Group = factor(Group, levels = display_levels)
  ) %>%
  group_by(Group) %>%
  summarise(
    median = median(Mean_Distance),
    lower  = quantile(Mean_Distance, 0.025),
    upper  = quantile(Mean_Distance, 0.975),
    .groups = "drop"
  )

print(sens_summary)

p_sensitivity_clean <- ggplot(
  sens_summary,
  aes(x = Group, y = median, color = Group)
) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.10,
    linewidth = 1
  ) +
  geom_point(size = 4) +
  scale_color_manual(values = pal_cond) +
  labs(
    title = "Balanced-subsampling sensitivity analysis",
    subtitle = "1,000 iterations; n = 9 per exposure stage",
    x = NULL,
    y = "Mean distance to centroid"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

print(p_sensitivity_clean)

ggsave(
  "Adult_PERMDISP_Balanced_Sensitivity_Clean.png",
  p_sensitivity_clean,
  width = 7,
  height = 5,
  dpi = 300
)

bd_bias <- betadisper(
  bray_dist,
  group = meta_df$Exposure_Condition,
  type = "centroid",
  bias.adjust = TRUE
)

permutest(
  bd_bias,
  pairwise = TRUE,
  permutations = 9999
)

bd_bias$group.distances