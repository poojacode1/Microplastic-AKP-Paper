
# AKP SENSITIVITY ANALYSIS
#
# Question:
# Is microbiome dispersion still higher in MP-exposed adults when:
#   1. age is controlled (23-day adults only)
#   2. PET / Nylon / PLA are considered separately
#   3. sample numbers are equalized
#


library(phyloseq)
library(vegan)
library(tidyverse)




setwd("C://Users//u0175886/Documents//Microbiome_Microplastic//Exposure_Microplastic")

pseq <- readRDS("daphnia_phyloseq.rds")




meta0 <- data.frame(
  sample_data(pseq),
  stringsAsFactors = FALSE
)

print(colnames(meta0))


print(
  table(
    meta0$Exposure_Condition,
    useNA = "ifany"
  )
)


print(
  table(
    meta0$Plastic_content,
    useNA = "ifany"
  )
)


print(
  table(
    meta0$Plastic,
    useNA = "ifany"
  )
)



#filtering

adult_stages <- c(
  "Pre_exposure",
  "Post_exposure",
  "Post_long_exposure"
)

keep <- (
  meta0$Exposure_Condition %in% adult_stages &
    meta0$Plastic_content %in% c("YES", "NO")
)

pseq_adult <- prune_samples(
  sample_names(pseq)[keep],
  pseq
)





pseq_filt <- filter_taxa(
  pseq_adult,
  function(x) {
    sum(x > 0) >= 3 &&
      sum(x) >= 10
  },
  prune = TRUE
)

# RELATIVE ABUNDANCE


pseq_rel <- transform_sample_counts(
  pseq_filt,
  function(x) x / sum(x)
)





meta_rel <- data.frame(
  sample_data(pseq_rel),
  stringsAsFactors = FALSE
)

ids_23 <- rownames(meta_rel)[
  meta_rel$Exposure_Condition == "Post_exposure"
]

ps23 <- prune_samples(
  ids_23,
  pseq_rel
)

ps23 <- prune_taxa(
  taxa_sums(ps23) > 0,
  ps23
)

meta23 <- data.frame(
  sample_data(ps23),
  stringsAsFactors = FALSE
)


# CREATE CONTROL / PLASTIC-TYPE VARIABLE



meta23$Treatment <- NA_character_

# No-plastic controls
meta23$Treatment[
  meta23$Plastic_content == "NO"
] <- "Control"

# Plastic-exposed animals
meta23$Treatment[
  meta23$Plastic_content == "YES"
] <- as.character(
  meta23$Plastic[
    meta23$Plastic_content == "YES"
  ]
)


# Remove unresolved treatment values
keep <- !is.na(meta23$Treatment) &
  meta23$Treatment != ""

meta23 <- meta23[
  keep,
  ,
  drop = FALSE
]

ps23 <- prune_samples(
  rownames(meta23),
  ps23
)

ps23 <- prune_taxa(
  taxa_sums(ps23) > 0,
  ps23
)




c
print(
  table(meta23$Treatment)
)

cat("\nPlastic_content x Plastic:\n")

print(
  table(
    meta23$Plastic_content,
    meta23$Plastic,
    useNA = "ifany"
  )
)


#identify plastic types
plastic_types <- unique(
  meta23$Treatment[
    meta23$Treatment != "Control"
  ]
)

plastic_types <- sort(
  plastic_types
)

cat("\nPlastic types detected:\n")
print(plastic_types)


# BRAY-CURTIS DISTANCE


bray23 <- phyloseq::distance(
  ps23,
  method = "bray"
)

bray_mat <- as.matrix(
  bray23
)

# Align metadata
meta23 <- meta23[
  rownames(bray_mat),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(meta23),
    rownames(bray_mat)
  )
)


# CONTROL + EACH PLASTIC TYPE AS SEPARATE GROUPS


meta23$Treatment <- factor(
  meta23$Treatment,
  levels = c(
    "Control",
    plastic_types
  )
)

bd_all <- vegan::betadisper(
  bray23,
  group = meta23$Treatment,
  type = "centroid",
  bias.adjust = TRUE
)

test_all <- vegan::permutest(
  bd_all,
  pairwise = TRUE,
  permutations = 9999
)


print(test_all)

cat("\nMean distance to centroid:\n")

print(
  round(
    bd_all$group.distances,
    4
  )
)



# DISPERSION TABLE


disp_df <- data.frame(
  Sample = names(
    bd_all$distances
  ),
  
  Treatment = meta23[
    names(bd_all$distances),
    "Treatment"
  ],
  
  Distance_to_Centroid =
    as.numeric(
      bd_all$distances
    )
)


disp_summary <- disp_df %>%
  group_by(Treatment) %>%
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
    
    .groups = "drop"
  )


cat("\nDispersion summary:\n")
print(disp_summary)


write.csv(
  disp_summary,
  "AKP_23day_PlasticType_Dispersion.csv",
  row.names = FALSE
)



# 12. CONTROL VS EACH PLASTIC TYPE SEPARATELY
#
# For every comparison:
# Control has its own centroid
# PET has its own centroid
#
# then:
# Control vs Nylon
# Control vs PLA
# etc.


compare_one_plastic <- function(plastic_name) {
  
  ids <- rownames(meta23)[
    meta23$Treatment %in%
      c(
        "Control",
        plastic_name
      )
  ]
  
  d <- as.dist(
    bray_mat[
      ids,
      ids,
      drop = FALSE
    ]
  )
  
  grp <- factor(
    meta23[
      ids,
      "Treatment"
    ],
    levels = c(
      "Control",
      plastic_name
    )
  )
  
  
  bd <- vegan::betadisper(
    d,
    group = grp,
    type = "centroid",
    bias.adjust = TRUE
  )
  
  
  perm <- vegan::permutest(
    bd,
    permutations = 9999
  )
  
  
  means <- tapply(
    bd$distances,
    grp,
    mean
  )
  
  
  data.frame(
    
    Plastic = plastic_name,
    
    n_control =
      sum(
        grp == "Control"
      ),
    
    n_plastic =
      sum(
        grp == plastic_name
      ),
    
    mean_control =
      as.numeric(
        means["Control"]
      ),
    
    mean_plastic =
      as.numeric(
        means[plastic_name]
      ),
    
    difference =
      as.numeric(
        means[plastic_name] -
          means["Control"]
      ),
    
    F =
      as.numeric(
        perm$tab[1, "F"]
      ),
    
    p =
      as.numeric(
        perm$tab[
          1,
          "Pr(>F)"
        ]
      )
  )
}


pairwise_results <- map_dfr(
  plastic_types,
  compare_one_plastic
)



print(pairwise_results)


write.csv(
  pairwise_results,
  "AKP_Control_vs_Each_Plastic_PERMDISP.csv",
  row.names = FALSE
)


# =============================================================================
# 13. BALANCED SUBSAMPLING
# For each comparison:
# Control vs PET
# Control vs Nylon
# Control vs PLA
# choose:
# n = size of smaller group
# repeat 1000 times


set.seed(123)

N_ITER <- 1000


balanced_one_plastic <- function(plastic_name) {
  
  
  control_ids <- rownames(meta23)[
    meta23$Treatment == "Control"
  ]
  
  
  plastic_ids <- rownames(meta23)[
    meta23$Treatment == plastic_name
  ]
  
  
  TARGET_N <- min(
    length(control_ids),
    length(plastic_ids)
  )
  
  
  
  cat("Control vs", plastic_name, "\n")
  cat("Control n:", length(control_ids), "\n")
  cat(plastic_name, "n:", length(plastic_ids), "\n")
  cat("Balanced n:", TARGET_N, "per group\n")
 
  
  
  if (TARGET_N < 3) {
    
    warning(
      paste0(
        plastic_name,
        ": fewer than 3 samples."
      )
    )
    
    return(NULL)
  }
  
  
  result <- data.frame(
    
    iteration =
      1:N_ITER,
    
    Plastic =
      plastic_name,
    
    n_per_group =
      TARGET_N,
    
    control_dispersion =
      NA_real_,
    
    plastic_dispersion =
      NA_real_,
    
    difference =
      NA_real_
  )
  
  
  for (i in 1:N_ITER) {
    
    
    
    # Random equal-size groups
  
    
    c_ids <- sample(
      control_ids,
      TARGET_N,
      replace = FALSE
    )
    
    
    p_ids <- sample(
      plastic_ids,
      TARGET_N,
      replace = FALSE
    )
    
    
    ids <- c(
      c_ids,
      p_ids
    )
    
    
    
    # Subset Bray-Curtis
    
    
    d <- as.dist(
      bray_mat[
        ids,
        ids,
        drop = FALSE
      ]
    )
    
    
    
    # Treatment group
    
    
    grp <- factor(
      meta23[
        ids,
        "Treatment"
      ],
      levels = c(
        "Control",
        plastic_name
      )
    )
    
    
    
    # PERMDISP
    
    bd <- vegan::betadisper(
      d,
      group = grp,
      type = "centroid",
      bias.adjust = TRUE
    )
    
    
    means <- tapply(
      bd$distances,
      grp,
      mean
    )
    
    
    result$control_dispersion[i] <-
      means["Control"]
    
    
    result$plastic_dispersion[i] <-
      means[plastic_name]
    
    
    result$difference[i] <-
      means[plastic_name] -
      means["Control"]
    
    
    if (i %% 100 == 0) {
      
      cat(
        plastic_name,
        ":",
        i,
        "/",
        N_ITER,
        "\n"
      )
    }
  }
  
  
  return(result)
}


balanced_results <- map_dfr(
  plastic_types,
  balanced_one_plastic
)


# SUMMARY OF BALANCED ANALYSIS


balanced_summary <- balanced_results %>%
  
  group_by(
    Plastic,
    n_per_group
  ) %>%
  
  summarise(
    
    median_control =
      median(
        control_dispersion
      ),
    
    median_plastic =
      median(
        plastic_dispersion
      ),
    
    median_difference =
      median(
        difference
      ),
    
    lower_95 =
      quantile(
        difference,
        0.025
      ),
    
    upper_95 =
      quantile(
        difference,
        0.975
      ),
    
    percent_plastic_greater =
      mean(
        difference > 0
      ) * 100,
    
    .groups = "drop"
  )



print(
  balanced_summary,
  n = Inf
)


write.csv(
  balanced_results,
  "AKP_Balanced_All_Iterations.csv",
  row.names = FALSE
)


write.csv(
  balanced_summary,
  "AKP_Balanced_Summary.csv",
  row.names = FALSE
)



# RAW 23-DAY DISPERSION BY PLASTIC TYPE


p1 <- ggplot(
  disp_df,
  aes(
    x = Treatment,
    y = Distance_to_Centroid,
    fill = Treatment
  )
) +
  
  geom_boxplot(
    width = 0.65,
    alpha = 0.8,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.12,
    alpha = 0.65,
    size = 2
  ) +
  
  labs(
    title =
      "Microbiome beta dispersion after 23 days",
    subtitle =
      "Age-matched adults; plastic types analysed separately",
    x = NULL,
    y = "Distance to group centroid"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    legend.position = "none",
    plot.title =
      element_text(
        face = "bold"
      )
  )


print(p1)


ggsave(
  "AKP_23day_Dispersion_by_Plastic.png",
  p1,
  width = 8,
  height = 5.5,
  dpi = 300
)



# BALANCED SUBSAMPLING
#
# > 0 = plastic more dispersed than control
# < 0 = control more dispersed


p2 <- ggplot(
  balanced_results,
  aes(
    x = Plastic,
    y = difference,
    fill = Plastic
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
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  labs(
    title =
      "Balanced sensitivity analysis",
    subtitle =
      "23-day adults: each plastic type vs age-matched control",
    x = "Plastic type",
    y =
      "Difference in mean distance to centroid\n(plastic − control)"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    legend.position = "none",
    plot.title =
      element_text(
        face = "bold"
      )
  )


print(p2)


ggsave(
  "AKP_Balanced_Control_vs_Plastic.png",
  p2,
  width = 8,
  height = 5.5,
  dpi = 300
)



#  PERMANOVA:
# Is each plastic also shifting community composition?


permanova_one_plastic <- function(plastic_name) {
  
  
  ids <- rownames(meta23)[
    meta23$Treatment %in%
      c(
        "Control",
        plastic_name
      )
  ]
  
  
  d <- as.dist(
    bray_mat[
      ids,
      ids,
      drop = FALSE
    ]
  )
  
  
  m <- data.frame(
    Treatment = factor(
      meta23[
        ids,
        "Treatment"
      ],
      levels = c(
        "Control",
        plastic_name
      )
    )
  )
  
  
  set.seed(1)
  
  
  x <- vegan::adonis2(
    d ~ Treatment,
    data = m,
    permutations = 9999
  )
  
  
  data.frame(
    
    Plastic =
      plastic_name,
    
    R2 =
      x$R2[1],
    
    F =
      x$F[1],
    
    p =
      x$`Pr(>F)`[1]
  )
}


permanova_results <- map_dfr(
  plastic_types,
  permanova_one_plastic
)




print(
  permanova_results
)


write.csv(
  permanova_results,
  "AKP_Control_vs_Each_Plastic_PERMANOVA.csv",
  row.names = FALSE
)



cat("\n\nANALYSIS COMPLETE\n")


cat("\n1. Sample sizes:\n")
print(table(meta23$Treatment))

cat("\n2. Full PERMDISP comparisons:\n")
print(pairwise_results)

cat("\n3. Balanced subsampling:\n")
print(balanced_summary, n = Inf)

cat("\n4. PERMANOVA:\n")
print(permanova_results)