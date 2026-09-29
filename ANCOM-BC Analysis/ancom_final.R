library(phyloseq)
library(tidyverse)

setwd("C:/Users/u0175886/Documents/Microbiome_Microplastic/Exposure_Microplastic")


# Family-level differential abundance (ANCOM-BC), ADULT trajectory



#helper to normalize sample IDs uniformly 
norm_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- toupper(x)
  gsub("[^A-Z0-9._-]", "", x)
}

# Load raw files 
otu_raw  <- read.csv("otu_table.csv", row.names = 1, check.names = FALSE)
tax_raw  <- read.csv("taxonomy_assignments.csv", row.names = 1, check.names = FALSE)
meta_raw <- read.csv("meta-edited.csv", sep = ";", row.names = 1, check.names = FALSE)

rownames(otu_raw)  <- norm_id(rownames(otu_raw))
rownames(meta_raw) <- norm_id(rownames(meta_raw))
if (!is.null(meta_raw$Sample.nr)) {
  rownames(meta_raw) <- norm_id(meta_raw$Sample.nr); meta_raw$Sample.nr <- NULL
}

# Daphnia only
if (!"Type" %in% colnames(meta_raw)) stop("Metadata has no 'Type' column.")
meta_raw$Type <- trimws(meta_raw$Type)
daph_meta <- meta_raw %>% filter(grepl("^DAPHNIA$", toupper(Type)))
if (nrow(daph_meta) == 0) stop("No Daphnia samples found.")

#Exposure_Condition
if (!"Category" %in% colnames(daph_meta)) stop("Metadata has no 'Category' column.")
daph_meta$Category <- trimws(daph_meta$Category)
daph_meta <- daph_meta %>%
  mutate(Exposure_Condition = dplyr::case_when(
    grepl("PRE[- ]?EXPOSURE", toupper(Category)) ~ "Pre_exposure",
    grepl("^DAPHNIA:\\s*POST[- ]?EXPOSURE$", toupper(Category)) ~ "Post_exposure",
    grepl("POST[- ]?LONG[- ]?EXPOSURE", toupper(Category)) ~ "Post_long_exposure",
    grepl("JUVENILE", toupper(Category)) ~ "Juvenile",
    TRUE ~ "Other"
  ))

# Fitering adults only and known plastic status  
# drops juveniles/env AND the post-long "not sure" samples, so this

if (!"Plastic_content" %in% colnames(daph_meta))
  stop("Metadata has no 'Plastic_content' column.")

daph_meta_adults <- daph_meta %>%
  filter(Exposure_Condition %in% c("Pre_exposure", "Post_exposure", "Post_long_exposure"),
         Plastic_content %in% c("YES", "NO"))

cat("Group sizes after filter (metadata):\n")
print(table(daph_meta_adults$Exposure_Condition))
if (nrow(daph_meta_adults) == 0) stop("No adult samples remain after filtering.")

#Overlap with OTU 
common_ids <- intersect(rownames(daph_meta_adults), rownames(otu_raw))
cat("Common samples (OTU vs metadata):", length(common_ids), "\n")
if (length(common_ids) == 0) stop("No overlapping sample IDs.")

daph_meta_adults <- daph_meta_adults[common_ids, , drop = FALSE]
otu_adults       <- otu_raw[common_ids, , drop = FALSE]
stopifnot(identical(rownames(daph_meta_adults), rownames(otu_adults)))

#Build phyloseq
asvs_keep <- intersect(colnames(otu_adults), rownames(tax_raw))
if (length(asvs_keep) == 0) stop("No ASV overlap between OTU and taxonomy.")
otu_adults <- otu_adults[, asvs_keep, drop = FALSE]
tax_use    <- as.matrix(tax_raw[asvs_keep, , drop = FALSE])

OTU  <- otu_table(t(as.matrix(otu_adults)), taxa_are_rows = TRUE)
TAX  <- tax_table(tax_use)
META <- sample_data(daph_meta_adults)
pseq <- phyloseq(OTU, TAX, META)

otu_mat <- as.matrix(otu_table(pseq))
meta_df <- as.data.frame(sample_data(pseq))
meta_df <- meta_df[match(colnames(otu_mat), rownames(meta_df)), , drop = FALSE]
stopifnot(identical(rownames(meta_df), colnames(otu_mat)))
sample_data(pseq) <- sample_data(meta_df)

sample_data(pseq)$Exposure_Condition <- factor(
  sample_data(pseq)$Exposure_Condition,
  levels = c("Pre_exposure", "Post_exposure", "Post_long_exposure")
)

cat("nsamples =", nsamples(pseq), " | ntaxa =", ntaxa(pseq), "\n")


saveRDS(pseq, "daphnia_phyloseq_ancom.rds")   


# ANCOM-BC (Family level)  —  BH correction (was "holm")

set.seed(123)
library(ANCOMBC)

cat("Samples:", nsamples(pseq), " ASVs:", ntaxa(pseq), "\n")
print(table(sample_data(pseq)$Exposure_Condition))

out_ancombc <- ancombc(
  data = pseq,
  tax_level = "Family",
  formula = "Exposure_Condition",
  p_adj_method = "BH",          
  prv_cut = 0.10,
  lib_cut = 1000,
  group = "Exposure_Condition",
  struc_zero = TRUE,
  neg_lb = TRUE,
  tol = 1e-5,
  max_iter = 100,
  conserve = TRUE,
  alpha = 0.05,
  global = TRUE,
  n_cl = 1,
  verbose = TRUE
)
cat("\n\u2713 ANCOM-BC completed\n")

#how many samples survived ANCOM-BC's internal lib_cut? 
cat("Samples entering ANCOM-BC by group (check lib_cut didn't shrink groups):\n")
print(table(sample_data(pseq)$Exposure_Condition,
            sample_sums(pseq) >= 1000))

# save results 
saveRDS(out_ancombc, "ancombc_results.rds")   


# Build the summary table (Post-vs-Pre, Post-long-vs-Pre, Global)

library(stringr); library(tibble)

`%||%` <- function(x, y) if (is.null(x)) y else x
clean_rank <- function(x) {
  x <- as.character(x); x <- str_trim(x)
  x <- str_replace_all(x, "^[kpcogfs]__+", "")
  x <- ifelse(is.na(x) | x == "" | toupper(x) %in% c("NA","UNCLASSIFIED"), NA, x)
  toupper(x)
}
is_all_integer_strings <- function(x) all(grepl("^[0-9]+$", x %||% character()))

res <- out_ancombc$res
res_global <- out_ancombc$res_global
taxon_names <- rownames(res$lfc)

if (is_all_integer_strings(taxon_names)) {
  candidate_labels <- if (!is.null(out_ancombc$taxa_abn)) rownames(out_ancombc$taxa_abn)
  else if (!is.null(out_ancombc$feature_table)) rownames(out_ancombc$feature_table)
  else NULL
  if (!is.null(candidate_labels) && length(candidate_labels) == length(taxon_names)) {
    taxon_names <- candidate_labels
  } else {
    fam_levels <- sort(unique(na.omit(clean_rank(as.data.frame(tax_table(pseq))$Family))))
    if (length(fam_levels) == length(taxon_names)) taxon_names <- fam_levels
    else stop("Could not recover family labels; inspect rownames(out_ancombc$res$lfc).")
  }
}

summary_table <- tibble(
  Taxon = taxon_names,
  Post_vs_Pre_LFC  = res$lfc$`Exposure_ConditionPost_exposure`,
  Post_vs_Pre_SE   = res$se$`Exposure_ConditionPost_exposure`,
  Post_vs_Pre_W    = res$W$`Exposure_ConditionPost_exposure`,
  Post_vs_Pre_pval = res$p_val$`Exposure_ConditionPost_exposure`,
  Post_vs_Pre_qval = res$q_val$`Exposure_ConditionPost_exposure`,
  Post_vs_Pre_sig  = as.logical(res$diff_abn$`Exposure_ConditionPost_exposure`),
  PostLong_vs_Pre_LFC  = res$lfc$`Exposure_ConditionPost_long_exposure`,
  PostLong_vs_Pre_SE   = res$se$`Exposure_ConditionPost_long_exposure`,
  PostLong_vs_Pre_W    = res$W$`Exposure_ConditionPost_long_exposure`,
  PostLong_vs_Pre_pval = res$p_val$`Exposure_ConditionPost_long_exposure`,
  PostLong_vs_Pre_qval = res$q_val$`Exposure_ConditionPost_long_exposure`,
  PostLong_vs_Pre_sig  = as.logical(res$diff_abn$`Exposure_ConditionPost_long_exposure`)
) %>%
  mutate(.row_key = rownames(res$lfc)) %>%
  left_join(tibble(
    .row_key = rownames(res_global),
    Global_W = res_global$W, Global_pval = res_global$p_val,
    Global_qval = res_global$q_val, Global_sig = as.logical(res_global$diff_abn)
  ), by = ".row_key") %>%
  select(-.row_key)

# attach higher ranks
tax_family <- as.data.frame(tax_table(pseq)) %>%
  rownames_to_column("ASV") %>%
  mutate(Family_c = clean_rank(Family)) %>%
  group_by(Family_c) %>%
  summarise(Phylum = dplyr::first(na.omit(Phylum)),
            Class  = dplyr::first(na.omit(Class)),
            Order  = dplyr::first(na.omit(Order)),
            Genus  = dplyr::first(na.omit(Genus)), .groups = "drop")

summary_table <- summary_table %>%
  mutate(Taxon_clean = clean_rank(Taxon)) %>%
  left_join(tax_family, by = c("Taxon_clean" = "Family_c")) %>%
  mutate(Family = ifelse(is.na(Taxon) | Taxon == "" |
                           toupper(Taxon) %in% c("NA","UNCLASSIFIED"),
                         "Unclassified_Family", Taxon)) %>%
  select(Taxon, Family, Phylum, Class, Order, Genus,
         starts_with("Post_vs_Pre_"), starts_with("PostLong_vs_Pre_"),
         starts_with("Global_"))

write.csv(summary_table, "ancombc_summary_all_taxa.csv", row.names = FALSE)
write.csv(filter(summary_table, Global_sig), "ancombc_significant_taxa.csv", row.names = FALSE)

cat("\nSignificant families (Global):", sum(summary_table$Global_sig, na.rm = TRUE), "\n")
cat("Significant Post-vs-Pre:",  sum(summary_table$Post_vs_Pre_sig,  na.rm = TRUE), "\n")
cat("Significant Post-long-vs-Pre:", sum(summary_table$PostLong_vs_Pre_sig, na.rm = TRUE), "\n")


# Significant-family counts + volcano plots

suppressPackageStartupMessages({ library(ggplot2); library(gridExtra) })

pal_cond <- c("Pre_exposure"="#999999","Post_exposure"="#E69F00","Post_long_exposure"="#56B4E9")

sig_counts <- tibble(
  Comparison = factor(c("Post vs Pre","Post-long vs Pre","Global (Any)"),
                      levels = c("Post vs Pre","Post-long vs Pre","Global (Any)")),
  Count = c(sum(summary_table$Post_vs_Pre_sig, na.rm = TRUE),
            sum(summary_table$PostLong_vs_Pre_sig, na.rm = TRUE),
            sum(summary_table$Global_sig, na.rm = TRUE)))

p1 <- ggplot(sig_counts, aes(Comparison, Count, fill = Comparison)) +
  geom_col(width = 0.6) + geom_text(aes(label = Count), vjust = -0.4, size = 5) +
  scale_fill_manual(values = c("Post vs Pre"="#E69F00","Post-long vs Pre"="#56B4E9",
                               "Global (Any)"="#009E73")) +
  labs(title = "Differentially abundant families (ANCOM-BC, Family level)",
       x = NULL, y = "Number of families") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))
ggsave("Fig_ANCOM_Significant_Counts.png", p1, width = 8, height = 5.5, dpi = 300)

volcano_df <- function(lfc, q, sig) tibble(lfc = lfc, qval = q, sig = as.logical(sig)) %>%
  filter(is.finite(lfc), is.finite(qval))
plot_volcano <- function(df, title, col_sig)
  ggplot(df, aes(lfc, -log10(qval), color = sig)) +
  geom_point(alpha = 0.65, size = 2.4) +
  scale_color_manual(values = c("grey70", col_sig),
                     labels = c("Not significant","Significant"), name = NULL) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  labs(title = title, x = "Log fold change", y = "-log10(q)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5))

p2a <- plot_volcano(volcano_df(res$lfc$Exposure_ConditionPost_exposure,
                               res$q_val$Exposure_ConditionPost_exposure,
                               res$diff_abn$Exposure_ConditionPost_exposure),
                    "Post vs Pre", "#E69F00")
p2b <- plot_volcano(volcano_df(res$lfc$Exposure_ConditionPost_long_exposure,
                               res$q_val$Exposure_ConditionPost_long_exposure,
                               res$diff_abn$Exposure_ConditionPost_long_exposure),
                    "Post-long vs Pre", "#56B4E9")
ggsave("Fig_ANCOM_Volcano_Combined.png", grid.arrange(p2a, p2b, ncol = 2),
       width = 14, height = 5.2, dpi = 300)


# Boxplots: top significant families (Post vs Pre, by q) — bias-corrected

meta <- as.data.frame(sample_data(pseq))
meta$Exposure_Condition <- factor(meta$Exposure_Condition,
                                  levels = c("Pre_exposure","Post_exposure","Post_long_exposure"))

top_fams <- summary_table %>% filter(Post_vs_Pre_sig) %>%
  arrange(Post_vs_Pre_qval) %>% distinct(Taxon, .keep_all = TRUE) %>%
  slice_head(n = 9) %>% pull(Taxon)

if (length(top_fams) > 0) {
  asv_mat <- as.matrix(otu_table(pseq))
  asv2fam <- setNames(as.character(as.data.frame(tax_table(pseq))$Family), taxa_names(pseq))
  keep    <- !is.na(asv2fam[rownames(asv_mat)]) & asv2fam[rownames(asv_mat)] != ""
  fam_mat <- rowsum(asv_mat[keep, , drop = FALSE], group = asv2fam[rownames(asv_mat)[keep]])
  
  samp_frac <- out_ancombc$samp_frac; samp_frac[is.na(samp_frac)] <- 0
  if (!is.null(names(samp_frac))) samp_frac <- samp_frac[colnames(fam_mat)]
  fam_bc <- sweep(log(fam_mat + 1), 2, samp_frac, "-")
  
  common <- intersect(colnames(fam_bc), rownames(meta))
  fam_bc <- fam_bc[, common, drop = FALSE]; meta <- meta[common, , drop = FALSE]
  
  df_long <- bind_rows(lapply(intersect(top_fams, rownames(fam_bc)), function(f)
    data.frame(Family = f, Sample = colnames(fam_bc),
               Abundance = as.numeric(fam_bc[f, ]),
               Exposure_Condition = as.character(meta$Exposure_Condition),
               stringsAsFactors = FALSE))) %>%
    filter(!is.na(Exposure_Condition)) %>%
    mutate(Exposure_Condition = factor(Exposure_Condition,
                                       levels = c("Pre_exposure","Post_exposure","Post_long_exposure")))
  
  p4 <- ggplot(df_long, aes(Exposure_Condition, Abundance, fill = Exposure_Condition)) +
    geom_boxplot(width = 0.7, alpha = 0.85, outlier.shape = NA) +
    geom_jitter(width = 0.18, alpha = 0.5, size = 1.5) +
    scale_fill_manual(values = pal_cond) +
    facet_wrap(~ Family, scales = "free_y", ncol = 3) +
    labs(title = "Top differentially abundant families (Post vs Pre)",
         subtitle = "Bias-corrected log abundance (ANCOM-BC samp_frac)",
         x = NULL, y = "Bias-corrected log abundance") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave("Fig_ANCOM_Top_Family_Boxplots.png", p4, width = 12,
         height = max(6, ceiling(length(top_fams)/3) * 3), dpi = 300)
} else message("No significant Post-vs-Pre families — skipping boxplots.")

cat("\nDone. Wrote: ancombc_summary_all_taxa.csv, ancombc_significant_taxa.csv,\n",
    "  ancombc_results.rds, daphnia_phyloseq_ancom.rds, and ANCOM figures.\n")




#  ANCOM-BC bias-corrected abundance boxplots
#   for the FAMILIES that define the LDA microbial components 

suppressPackageStartupMessages({
  library(phyloseq); library(dplyr); library(tidyr)
  library(tibble);   library(ggplot2); library(readr)
})

if (!exists("pseq"))        pseq        <- readRDS("daphnia_phyloseq_ancom.rds")
if (!exists("out_ancombc")) out_ancombc <- readRDS("ancombc_results.rds")

# --- STEP 0: SEE your real family names first (run once, then fix list) ------
cat("\n=== Families present in tax_table(pseq) ===\n")
print(sort(unique(as.character(tax_table(pseq)[, "Family"]))))


component_families <- c(
  "Comamonadaceae",      # MC2/MC3/MC5: f_Comamonadaceae, g_Hydrogenophaga, g_Polynucleobacter
  "Rhodobacteraceae",    # MC0/MC3/MC5: f_Rhodobacteraceae, g_Gemmobacter, g_Pseudorhodobacter
  "Flavobacteriaceae",   # MC2: g_Flavobacterium
  "Spirosomaceae",       # MC2: g_Lacihabitans        (VERIFY: sometimes Cytophagaceae)
  "Pseudomonadaceae",    # MC2: g_Pseudomonas
  "Rubinisphaeraceae",   # MC3: f_Rubinisphaeraceae
  "Pirellulaceae",       # MC3: g_Rhodopirellula      (o_Planctomycetales member)
  "Verrucomicrobiaceae"  # MC0: g_Prosthecobacter     (VERIFY against taxonomy)

)


# Build Family x Sample bias-corrected matrix 

asv_mat <- as.matrix(otu_table(pseq))
asv2fam <- setNames(as.character(as.data.frame(tax_table(pseq))$Family), taxa_names(pseq))

keep    <- !is.na(asv2fam[rownames(asv_mat)]) & asv2fam[rownames(asv_mat)] != ""
fam_mat <- rowsum(asv_mat[keep, , drop = FALSE], group = asv2fam[rownames(asv_mat)[keep]])

samp_frac <- out_ancombc$samp_frac; samp_frac[is.na(samp_frac)] <- 0
if (!is.null(names(samp_frac))) samp_frac <- samp_frac[colnames(fam_mat)]
fam_bc <- sweep(log(fam_mat + 1), 2, samp_frac, "-")   # Family x Sample

meta <- as.data.frame(sample_data(pseq))
meta$Exposure_Condition <- factor(meta$Exposure_Condition,
                                  levels = c("Pre_exposure","Post_exposure","Post_long_exposure"))
common <- intersect(colnames(fam_bc), rownames(meta))
fam_bc <- fam_bc[, common, drop = FALSE]; meta <- meta[common, , drop = FALSE]
stopifnot(identical(colnames(fam_bc), rownames(meta)))


present <- intersect(component_families, rownames(fam_bc))
missing <- setdiff(component_families, rownames(fam_bc))
cat("\nPlotting component families:", paste(present, collapse = ", "), "\n")
if (length(missing))
  message("NOT found in taxonomy (skipped): ", paste(missing, collapse = ", "),
          "\n  -> fix names against the printed family list above.")
stopifnot(length(present) > 0)

#annotate each family with its ANCOM-BC significance (
if (exists("summary_table")) {
  sig_lookup <- summary_table %>%
    transmute(Family,
              PvP  = ifelse(Post_vs_Pre_sig %in% TRUE, "*", ""),
              PLvP = ifelse(PostLong_vs_Pre_sig %in% TRUE, "*", ""))
} else sig_lookup <- NULL


#figure

df_long <- bind_rows(lapply(present, function(f)
  data.frame(Family = f, Sample = colnames(fam_bc),
             Abundance = as.numeric(fam_bc[f, ]),
             Exposure_Condition = as.character(meta$Exposure_Condition),
             stringsAsFactors = FALSE))) %>%
  filter(!is.na(Exposure_Condition)) %>%
  mutate(Family = factor(Family, levels = present),
         Exposure_Condition = factor(Exposure_Condition,
                                     levels = c("Pre_exposure","Post_exposure","Post_long_exposure")))

pal_cond <- c("Pre_exposure"="#999999","Post_exposure"="#E69F00","Post_long_exposure"="#56B4E9")

p_comp <- ggplot(df_long, aes(Exposure_Condition, Abundance, fill = Exposure_Condition)) +
  geom_boxplot(width = 0.7, alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.18, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = pal_cond) +
  facet_wrap(~ Family, scales = "free_y", ncol = 3) +
  labs(title = "Microbial-component families across exposure stages",
       subtitle = "ANCOM-BC bias-corrected log abundance (adults; post-long = known plastic status)",
       x = NULL, y = "Bias-corrected log abundance") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("Fig_ANCOM_ComponentFamilies_Boxplots.png", p_comp,
       width = 12, height = max(6, ceiling(length(present)/3) * 3), dpi = 300)
cat("\u2713 Wrote: Fig_ANCOM_ComponentFamilies_Boxplots.png\n")

#summary + significance table for these families 
comp_summary <- df_long %>%
  group_by(Family, Exposure_Condition) %>%
  summarise(n = dplyr::n(), mean = mean(Abundance),
            median = median(Abundance), sd = sd(Abundance), .groups = "drop")
if (!is.null(sig_lookup))
  comp_summary <- left_join(comp_summary, sig_lookup, by = "Family")
write_csv(comp_summary, "ComponentFamilies_Summary.csv")
cat("\u2713 Wrote: ComponentFamilies_Summary.csv\n")

#print which component families are ANCOM-BC significant 
if (exists("summary_table")) {
  cat("\nANCOM-BC significance for component families:\n")
  summary_table %>%
    filter(Family %in% present) %>%
    select(Family, Post_vs_Pre_LFC, Post_vs_Pre_qval, Post_vs_Pre_sig,
           PostLong_vs_Pre_LFC, PostLong_vs_Pre_qval, PostLong_vs_Pre_sig) %>%
    arrange(Post_vs_Pre_qval) %>%
    print(n = Inf)
}