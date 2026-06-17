#!/usr/bin/env Rscript
# ============================================================
# Standalone taxonomy metagenomics pipeline
# Generated for GitHub/Zenodo reproducibility.
#
# Run from repository root, for example:
#   Rscript taxonomy/scripts/run_taxonomy_pipeline.R
#
# Edit PROJECT_DIR/NATCOM_DIR below only if you do not run from
# the repository root. You can also set environment variables:
#   PROJECT_DIR=/path/to/repo NATCOM_DIR=/path/to/output Rscript ...
# ============================================================

options(stringsAsFactors = FALSE)
set.seed(1)

get_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)
    return(normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

PROJECT_DIR <- Sys.getenv(
  "PROJECT_DIR",
  unset = get_repo_root()
)

# manuscript_figures root
NATCOM_DIR <- Sys.getenv(
  "NATCOM_DIR",
  unset = file.path(PROJECT_DIR, "manuscript_figures")
)

# Main manuscript figures (Figures 1–5)
MAIN_FIG_DIR <- file.path(
  NATCOM_DIR,
  "Main_Figures"
)

# Supplementary taxonomy figures
SUPP_TAX_DIR <- file.path(
  NATCOM_DIR,
  "Supplementary_Taxonomy_Figures"
)

# Create directories
dir.create(NATCOM_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MAIN_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_TAX_DIR, recursive = TRUE, showWarnings = FALSE)

message("PROJECT_DIR: ", PROJECT_DIR)
message("NATCOM_DIR:  ", NATCOM_DIR)
message("SUPP_TAX_DIR: ", SUPP_TAX_DIR)

run_section <- function(label, required = TRUE, code) {
  message("\n============================================================")
  message(label)
  message("============================================================")
  tryCatch(
    code(),
    error = function(e) {
      if (isTRUE(required)) stop(conditionMessage(e), call. = FALSE)
      message("Optional section skipped: ", conditionMessage(e))
    }
  )
}

run_section('Step 1 — Generate taxonomy PCoA panel for Figure 1B', required = FALSE, code = function() {
# ============================================================
# Nature Communications-style global taxonomy PCoA
# Figure 1B: Environmental filtering across sites
# Input: Bracken genus tables for root and rhizosphere
# Output: compact PDF/PNG suitable for merging with map + pathway PCoA
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggplot2)
})

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
rename <- dplyr::rename
arrange <- dplyr::arrange
summarise <- dplyr::summarise
left_join <- dplyr::left_join

# ----------------------------
# USER PATHS
# ----------------------------
PROJECT_DIR <- file.path(PROJECT_DIR, "taxonomy")
ROOT_DIR  <- file.path(PROJECT_DIR, "data", "root")
RHIZO_DIR <- file.path(PROJECT_DIR, "data", "rhizosphere")

root_meta_file  <- file.path(ROOT_DIR,  "Metadata.txt")
rhizo_meta_file <- file.path(RHIZO_DIR, "Metadata.txt")
root_genus_file  <- file.path(ROOT_DIR,  "Bracken_Genus.tsv")
rhizo_genus_file <- file.path(RHIZO_DIR, "Bracken_Genus.tsv")

out_dir <- file.path(PROJECT_DIR, "results", "Figure1_NatCom")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# HELPERS
# ----------------------------
norm_sampleid <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^X", "", x)
  x <- gsub("\\.sorted.*$", "", x)
  x <- gsub("\\.bam$", "", x)
  x <- gsub("\\.fastq.*$", "", x)
  x <- gsub("\\.fq.*$", "", x)
  x <- gsub("\\.gz$", "", x)
  x
}

read_tax_table <- function(file) {
  df <- read.delim(file, sep = "\t", check.names = FALSE)
  if (colnames(df)[1] != "name") df <- df %>% rename(name = all_of(colnames(df)[1]))
  df$name <- trimws(as.character(df$name))
  colnames(df) <- c("name", norm_sampleid(colnames(df)[-1]))
  df
}

read_meta <- function(meta_path, compartment_label) {
  meta <- read.delim(meta_path, sep = "\t", check.names = FALSE)
  if (!("Location" %in% names(meta))) {
    loc_col <- intersect(names(meta), c("Site", "SITE", "site", "SamplingSite", "FieldSite"))[1]
    if (!is.na(loc_col)) meta <- meta %>% rename(Location = all_of(loc_col))
  }
  stopifnot(all(c("SampleID", "Location", "Inoculation") %in% names(meta)))
  meta %>%
    mutate(
      SampleID = norm_sampleid(SampleID),
      Location = trimws(as.character(Location)),
      Inoculation = trimws(as.character(Inoculation)),
      Compartment = compartment_label
    ) %>%
    mutate(
      Location = case_when(
        Location %in% c("Eyrewell Forest", "Eyrewell_Forest", "Eyerell_Forest", "Eyerell Forest") ~ "Eyrewell_Forest",
        Location %in% c("West Coast", "West_Coast") ~ "West_Coast",
        TRUE ~ Location
      ),
      Location = factor(Location, levels = c("Eyrewell_Forest", "Kowhai", "LU_H8", "Rolleston", "West_Coast")),
      Inoculation = factor(Inoculation, levels = c("Control", "Panch")),
      Compartment = factor(Compartment, levels = c("Rhizosphere", "Root"))
    )
}

wide_to_props <- function(file, meta) {
  df <- read_tax_table(file)
  mat <- df %>% column_to_rownames("name") %>% as.matrix()
  colnames(mat) <- norm_sampleid(colnames(mat))
  mat[is.na(mat)] <- 0
  keep <- intersect(meta$SampleID, colnames(mat))
  mat <- mat[, keep, drop = FALSE]
  cs <- colSums(mat)
  props <- sweep(mat, 2, cs, "/")
  props[is.na(props)] <- 0
  props
}

expand_mat <- function(mat, taxa) {
  out <- matrix(0, nrow = length(taxa), ncol = ncol(mat), dimnames = list(taxa, colnames(mat)))
  out[rownames(mat), colnames(mat)] <- mat
  out
}

theme_natcom_pcoa <- function(base_size = 7, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(face = "bold", size = base_size + 1),
      axis.text = element_text(color = "black", size = base_size),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 0.2),
      legend.key.size = unit(0.32, "cm"),
      legend.spacing.y = unit(0.02, "cm"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.3),
      plot.margin = margin(3, 3, 3, 3)
    )
}

loc_cols <- c(
  "Eyrewell_Forest" = "#D55E00",
  "Kowhai"          = "#E69F00",
  "LU_H8"           = "#009E73",
  "Rolleston"       = "#0072B2",
  "West_Coast"      = "#CC79A7"
)

# ----------------------------
# PREP COMBINED DATA
# ----------------------------
meta_root  <- read_meta(root_meta_file, "Root")
meta_rhizo <- read_meta(rhizo_meta_file, "Rhizosphere")

props_root  <- wide_to_props(root_genus_file, meta_root)
props_rhizo <- wide_to_props(rhizo_genus_file, meta_rhizo)

all_taxa <- union(rownames(props_root), rownames(props_rhizo))
props_root2  <- expand_mat(props_root, all_taxa)
props_rhizo2 <- expand_mat(props_rhizo, all_taxa)
props_combined <- cbind(props_rhizo2, props_root2)

meta_combined <- bind_rows(meta_rhizo, meta_root) %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  filter(SampleID %in% colnames(props_combined)) %>%
  arrange(Location, Compartment, Inoculation, SampleID)

props_combined <- props_combined[, meta_combined$SampleID, drop = FALSE]
stopifnot(identical(colnames(props_combined), meta_combined$SampleID))

# Hellinger transformation
props_hell <- sqrt(props_combined)

# Bray-Curtis + PCoA
bc <- vegdist(t(props_hell), method = "bray")
pcoa <- cmdscale(bc, k = 2, eig = TRUE)

pcoa_df <- as.data.frame(pcoa$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$SampleID <- rownames(pcoa_df)
pcoa_df <- left_join(pcoa_df, meta_combined, by = "SampleID")

pos_eig <- pcoa$eig[pcoa$eig > 0]
var_expl <- pos_eig / sum(pos_eig)
xlab <- paste0("PCoA1 (", round(var_expl[1] * 100, 1), "%)")
ylab <- paste0("PCoA2 (", round(var_expl[2] * 100, 1), "%)")

# PERMANOVA for annotation
set.seed(1)
perm <- adonis2(bc ~ Location + Compartment + Inoculation,
                data = meta_combined,
                permutations = 999,
                by = "margin")
perm_tbl <- as.data.frame(perm) %>% tibble::rownames_to_column("Term")
write_tsv(perm_tbl, file.path(out_dir, "PERMANOVA_taxonomy_global.tsv"))

loc_r2 <- perm_tbl$R2[perm_tbl$Term == "Location"]
loc_p  <- perm_tbl$`Pr(>F)`[perm_tbl$Term == "Location"]
p_label <- ifelse(is.na(loc_p), "NA", ifelse(loc_p < 0.001, "<0.001", paste0("= ", signif(loc_p, 2))))

# ----------------------------
# NATURE COMMUNICATIONS STYLE GLOBAL PLOT
# ----------------------------
p_natcom <- ggplot(pcoa_df, aes(PCoA1, PCoA2)) +
  geom_point(aes(color = Location, shape = Compartment), size = 1.9, alpha = 0.92, stroke = 0.25) +
  scale_color_manual(values = loc_cols, drop = FALSE) +
  scale_shape_manual(values = c(Rhizosphere = 16, Root = 17), drop = FALSE) +
  labs(x = xlab, y = ylab, color = "Location", shape = "Compartment") +
  guides(color = guide_legend(override.aes = list(size = 2.1), order = 1),
         shape = guide_legend(override.aes = list(size = 2.1), order = 2)) +
  theme_natcom_pcoa(base_size = 7)

# Compact panel for merging into Fig. 1
# 85 x 70 mm is suitable for a 2-column layout with Panel A map above or beside it.
ggsave(file.path(out_dir, "Fig1B_Taxonomy_PCoA_NatCom_panel.pdf"), p_natcom,
       width = 85, height = 70, units = "mm", useDingbats = FALSE)
ggsave(file.path(out_dir, "Fig1B_Taxonomy_PCoA_NatCom_panel.png"), p_natcom,
       width = 85, height = 70, units = "mm", dpi = 600)

# Wider inspection version
ggsave(file.path(out_dir, "Fig1B_Taxonomy_PCoA_NatCom_wide.pdf"), p_natcom,
       width = 110, height = 80, units = "mm", useDingbats = FALSE)
ggsave(file.path(out_dir, "Fig1B_Taxonomy_PCoA_NatCom_wide.png"), p_natcom,
       width = 110, height = 80, units = "mm", dpi = 600)

message("Done. Saved Nature Communications-style taxonomy PCoA to: ", out_dir)

})


run_section('Step 2 — Run main taxonomy workflow for root and rhizosphere', required = TRUE, code = function() {
# ============================================================
# Chapter 2 — ROOT + RHIZOSPHERE metagenomics taxonomy pipeline
# Kowhai shown twice (Rainfed & Irrigated); other locations collapsed
# Publication-ready + reproducible + consistent palettes
#
# UPDATED:
# 1) Standardised site name to Eyrewell_Forest
# 2) Rhizosphere main PCoA removes visual outlier R_W7
# 3) Full PCoA retained as supplementary
# 4) Statistics remain based on the full dataset
# 5) Bray–Curtis is calculated on square-root-transformed relative abundance
#
# Run:
#   Rscript scripts/run_taxonomy_both_compartments.R
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(digest)
})

# Avoid function masking problems from other loaded packages.
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
rename    <- dplyr::rename
arrange   <- dplyr::arrange
summarise <- dplyr::summarise
left_join <- dplyr::left_join

if (!requireNamespace("ggsci", quietly = TRUE)) {
  stop("Please install ggsci: install.packages('ggsci')")
}

# ----------------------------
# REPRODUCIBILITY SETTINGS
# ----------------------------
set.seed(1)

# ----------------------------
# USER SETTINGS
# ----------------------------
TOPN <- list(Kingdom5 = 6, Phylum = 10, Family = 30, Genus = 50)

USE_GENUS_RARITY_BINS <- TRUE
GENUS_BINS <- list(
  "Low (0.4–0.5%)"         = c(0.4, 0.5),
  "Very Low (0.3–0.4%)"    = c(0.3, 0.4),
  "Rare (0.2–0.3%)"        = c(0.2, 0.3),
  "Ultra Rare (0.05–0.2%)" = c(0.05, 0.2),
  "Very Rare (0.01–0.05%)" = c(0.01, 0.05)
)
GENUS_OTHER_LABEL <- "Other (<0.01%)"
KEEP_VIRUSES_AS_BIN <- TRUE
REMOVE_VERTEBRATES <- TRUE

# UPDATED: Eyrewell_Forest standardised here
site_levels <- c("Eyrewell_Forest", "Kowhai", "LU_H8", "Rolleston", "West_Coast")
inoc_levels <- c("Control", "Panch")

# Facet order (Kowhai appears twice)
site_levels_k <- c("Kowhai (Rainfed)", "Kowhai (Irrigated)",
                   setdiff(site_levels, "Kowhai"))

# Legend columns
LEGEND_NCOL_PHYLUM  <- 1
LEGEND_NCOL_FAMILY  <- 2
LEGEND_NCOL_GENUS   <- 2
LEGEND_NCOL_KINGDOM <- 1

# ----------------------------
# PCoA main/supplementary settings
# ----------------------------
REMOVE_RHIZO_PCOA_OUTLIER_FOR_MAIN <- TRUE
RHIZO_PCOA_OUTLIER_IDS <- c("R_W7")

# ----------------------------
# PATHS (EDIT THESE)
# ----------------------------
library(here)
PROJECT_DIR <- PROJECT_DIR

DATASETS <- list(
  Root = list(
    data_dir = file.path(PROJECT_DIR, "taxonomy", "data", "root"),
    meta_file   = "Metadata.txt",
    phylum_file = "Bracken_Phylum.tsv",
    family_file = "Bracken_Family.tsv",
    genus_file  = "Bracken_Genus.tsv"
  ),
  Rhizosphere = list(
    data_dir = file.path(PROJECT_DIR, "taxonomy", "data", "rhizosphere"),
    meta_file   = "Metadata.txt",
    phylum_file = "Bracken_Phylum.tsv",
    family_file = "Bracken_Family.tsv",
    genus_file  = "Bracken_Genus.tsv"
  )
)

OUT_DIR <- file.path(PROJECT_DIR, "taxonomy", "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Extra output folder for manuscript-ready figures
NATCOM_DIR <- NATCOM_DIR
dir.create(NATCOM_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Utilities
# ----------------------------
norm_sampleid <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^X", "", x)
  x <- gsub("\\.sorted.*$", "", x)
  x <- gsub("\\.bam$", "", x)
  x <- gsub("\\.fastq.*$", "", x)
  x <- gsub("\\.fq.*$", "", x)
  x <- gsub("\\.gz$", "", x)
  x
}

sha1_file <- function(path) digest::digest(file = path, algo = "sha1")

theme_natureish <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_text(size = base_size, margin = margin(t = 7)),
      axis.title.y = element_text(size = base_size, margin = margin(r = 7)),
      axis.text = element_text(size = base_size),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      axis.ticks.length = unit(2, "pt"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      legend.title = element_text(size = base_size + 2, face = "bold"),
      legend.text  = element_text(size = base_size + 1),
      legend.key.size = unit(4, "mm"),
      legend.spacing.y = unit(2, "mm"),
      legend.box.spacing = unit(4, "mm"),
      panel.spacing = unit(4, "mm"),
      plot.margin = margin(10, 12, 10, 10)
    )
}

okabe_ito <- c(
  orange  = "#E69F00",
  sky     = "#56B4E9",
  green   = "#009E73",
  yellow  = "#F0E442",
  blue    = "#0072B2",
  verm    = "#D55E00",
  purple  = "#CC79A7",
  grey    = "#999999"
)

loc_cols <- c(
  "Eyrewell_Forest" = unname(okabe_ito["verm"]),
  "Kowhai"          = unname(okabe_ito["orange"]),
  "LU_H8"           = unname(okabe_ito["green"]),
  "Rolleston"       = unname(okabe_ito["blue"]),
  "West_Coast"      = unname(okabe_ito["purple"])
)

theme_nature_ordination <- function(base_size = 7, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title   = element_text(face = "bold", size = base_size + 1, hjust = 0),
      axis.title   = element_text(face = "bold", size = base_size + 1, colour = "black"),
      axis.text    = element_text(colour = "black", size = base_size),
      axis.line    = element_line(linewidth = 0.30, colour = "black"),
      axis.ticks   = element_line(linewidth = 0.30, colour = "black"),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text  = element_text(size = base_size),
      legend.key   = element_blank(),
      legend.key.size = unit(3.5, "mm"),
      legend.position = "right",
      panel.border = element_rect(fill = NA, linewidth = 0.35, colour = "black"),
      plot.margin  = margin(4, 4, 4, 4)
    )
}

save_pdf <- function(path, p, w = 235, h = 125) {
  ggsave(path, p, width = w, height = h, units = "mm",
         dpi = 600, useDingbats = FALSE, limitsize = FALSE)
}

save_png <- function(path, p, w = 235, h = 125, dpi = 600) {
  ggsave(path, p, width = w, height = h, units = "mm",
         dpi = dpi, device = "png", limitsize = FALSE)
}

# ----------------------------
# Filters
# ----------------------------
remove_human_contamination <- function(df) {
  df <- df %>% mutate(name = trimws(as.character(name)))
  
  human_like <- c(
    "^Homo$", "Homo\\s*sapiens", "^Hominidae$", "^Homininae$",
    "^Primates$", "^Mammalia$", "^Metazoa$", "^Animalia$"
  )
  vertebrate_like <- c("^Chordata$")
  
  pat <- paste(c(
    human_like,
    if (isTRUE(REMOVE_VERTEBRATES)) vertebrate_like else character(0)
  ), collapse = "|")
  
  df %>% filter(!stringr::str_detect(name, stringr::regex(pat, ignore_case = TRUE)))
}

remove_unclassified_like <- function(df, name_col = "name") {
  x <- trimws(as.character(df[[name_col]]))
  df %>%
    filter(!str_detect(x, regex("unclassified|unassigned|unknown|root", ignore_case = TRUE)))
}

# ----------------------------
# I/O
# ----------------------------
read_tax_table <- function(file) {
  df <- read.delim(file, sep = "\t", check.names = FALSE)
  if (colnames(df)[1] != "name") df <- df %>% rename(name = all_of(colnames(df)[1]))
  df$name <- trimws(as.character(df$name))
  colnames(df) <- c("name", norm_sampleid(colnames(df)[-1]))
  df
}

read_meta <- function(meta_path) {
  meta <- read.delim(meta_path, sep = "\t", check.names = FALSE)
  
  if (!("Location" %in% names(meta))) {
    loc_col <- intersect(names(meta), c("Site", "SITE", "site", "SamplingSite", "FieldSite"))[1]
    if (!is.na(loc_col)) meta <- meta %>% rename(Location = all_of(loc_col))
  }
  if (!("Water" %in% names(meta))) {
    water_col <- intersect(names(meta), c("Water", "WATER", "water", "Irrigation", "Regime", "TreatmentWater"))[1]
    if (!is.na(water_col)) meta <- meta %>% rename(Water = all_of(water_col))
  }
  
  stopifnot(all(c("SampleID", "Location", "Inoculation", "Water") %in% names(meta)))
  
  meta <- meta %>%
    mutate(
      SampleID    = norm_sampleid(SampleID),
      Location    = trimws(as.character(Location)),
      Inoculation = trimws(as.character(Inoculation)),
      Water       = trimws(as.character(Water))
    ) %>%
    mutate(
      Location = na_if(Location, ""),
      Location = na_if(Location, "NA"),
      Location = na_if(Location, "N/A"),
      Water    = na_if(Water, ""),
      Water    = na_if(Water, "NA"),
      Water    = na_if(Water, "N/A")
    ) %>%
    mutate(
      Water_plot = if_else(
        Location == "Kowhai",
        if_else(str_detect(Water, regex("irr", ignore_case = TRUE)), "Irrigated", "Rainfed"),
        "All"
      )
    )
  
  meta$Location    <- factor(meta$Location, levels = site_levels)
  meta$Inoculation <- factor(meta$Inoculation, levels = inoc_levels)
  meta$Water_plot  <- factor(meta$Water_plot, levels = c("All", "Rainfed", "Irrigated"))
  
  stopifnot(!any(is.na(meta$Location)))
  stopifnot(!any(is.na(meta$Inoculation)))
  stopifnot(!any(is.na(meta$Water_plot)))
  
  meta
}

tax_long_with_meta <- function(file, meta) {
  df <- read_tax_table(file) %>%
    remove_human_contamination()
  
  long <- df %>%
    pivot_longer(-name, names_to = "SampleID", values_to = "Value") %>%
    mutate(SampleID = norm_sampleid(SampleID)) %>%
    left_join(meta %>% mutate(SampleID = as.character(SampleID)), by = "SampleID")
  
  miss <- mean(is.na(long$Location))
  if (!is.na(miss) && miss > 0) {
    bad <- head(long$SampleID[is.na(long$Location)], 25)
    stop("❌ Join produced missing Location (", round(miss * 100, 2),
         "%). Example SampleIDs: ", paste(bad, collapse = ", "))
  }
  
  long$Value <- suppressWarnings(as.numeric(long$Value))
  long$Value[is.na(long$Value)] <- 0
  
  if (max(long$Value, na.rm = TRUE) <= 1.5) long <- long %>% mutate(Percent = Value * 100)
  else long <- long %>% mutate(Percent = Value)
  
  long %>%
    mutate(
      Location    = factor(as.character(Location), levels = site_levels),
      Inoculation = factor(as.character(Inoculation), levels = inoc_levels),
      Water_plot  = factor(as.character(Water_plot), levels = c("All", "Rainfed", "Irrigated")),
      LocationK = if_else(
        as.character(Location) == "Kowhai",
        paste0("Kowhai (", as.character(Water_plot), ")"),
        as.character(Location)
      ),
      LocationK = factor(LocationK, levels = site_levels_k)
    )
}

# ----------------------------
# Aggregation (uses LocationK)
# ----------------------------
prep_group_mean <- function(file, meta, top_n = 15, other_label = "Other", always_keep = character(0)) {
  long <- tax_long_with_meta(file, meta) %>% mutate(name = as.character(name))
  
  grp_all <- long %>%
    group_by(LocationK, Inoculation, name) %>%
    summarise(Percent = mean(Percent, na.rm = TRUE), .groups = "drop") %>%
    complete(LocationK, Inoculation, name, fill = list(Percent = 0)) %>%
    group_by(LocationK, Inoculation) %>%
    mutate(tot = sum(Percent, na.rm = TRUE),
           Percent = ifelse(tot > 0, Percent / tot * 100, 0)) %>%
    ungroup() %>% dplyr::select(-tot)
  
  top_taxa <- grp_all %>%
    group_by(name) %>%
    summarise(overall = mean(Percent, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(overall)) %>%
    slice_head(n = top_n) %>%
    pull(name)
  
  keep_extra <- character(0)
  if (length(always_keep) > 0) {
    keep_extra <- grp_all %>%
      filter(str_detect(name, regex(paste(always_keep, collapse = "|"), ignore_case = TRUE))) %>%
      pull(name) %>% unique()
  }
  keep_set <- unique(c(top_taxa, keep_extra))
  
  grp2 <- grp_all %>%
    mutate(Taxon = ifelse(name %in% keep_set, name, other_label)) %>%
    group_by(LocationK, Inoculation, Taxon) %>%
    summarise(Percent = sum(Percent, na.rm = TRUE), .groups = "drop") %>%
    complete(LocationK, Inoculation, Taxon, fill = list(Percent = 0)) %>%
    group_by(LocationK, Inoculation) %>%
    mutate(Percent = Percent / sum(Percent, na.rm = TRUE) * 100) %>%
    ungroup()
  
  grp2
}

prep_genus_with_rarity_bins <- function(genus_file, meta, top_n = 50, bins = GENUS_BINS,
                                        other_label = GENUS_OTHER_LABEL) {
  long <- tax_long_with_meta(genus_file, meta)
  
  grp_all <- long %>%
    group_by(LocationK, Inoculation, name) %>%
    summarise(Percent = mean(Percent, na.rm = TRUE), .groups = "drop") %>%
    complete(LocationK, Inoculation, name, fill = list(Percent = 0)) %>%
    group_by(LocationK, Inoculation) %>%
    mutate(tot = sum(Percent, na.rm = TRUE),
           Percent = ifelse(tot > 0, Percent / tot * 100, 0)) %>%
    ungroup() %>% dplyr::select(-tot)
  
  top_genera <- grp_all %>%
    group_by(name) %>%
    summarise(overall = mean(Percent, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(overall)) %>%
    slice_head(n = top_n) %>%
    pull(name)
  
  global_mean <- grp_all %>%
    group_by(name) %>%
    summarise(mu = mean(Percent, na.rm = TRUE), .groups = "drop")
  
  bin_one <- function(mu) {
    for (nm in names(bins)) {
      lo <- bins[[nm]][1]
      hi <- bins[[nm]][2]
      if (mu >= lo && mu < hi) return(nm)
    }
    if (mu >= 0.5) return(">=0.5% (not in Top50)")
    other_label
  }
  
  bin_map <- global_mean %>%
    mutate(Taxon = if_else(name %in% top_genera, name, map_chr(mu, bin_one))) %>%
    dplyr::select(name, Taxon)
  
  grp2 <- grp_all %>%
    left_join(bin_map, by = "name") %>%
    group_by(LocationK, Inoculation, Taxon) %>%
    summarise(Percent = sum(Percent, na.rm = TRUE), .groups = "drop") %>%
    complete(LocationK, Inoculation, Taxon, fill = list(Percent = 0)) %>%
    group_by(LocationK, Inoculation) %>%
    mutate(Percent = Percent / sum(Percent, na.rm = TRUE) * 100) %>%
    ungroup()
  
  grp2
}

check_bar_sums <- function(df, label = "") {
  s <- df %>% group_by(LocationK, Inoculation) %>% summarise(sumP = sum(Percent), .groups = "drop")
  bad <- s %>% filter(abs(sumP - 100) > 0.01)
  if (nrow(bad) > 0) {
    message("⚠ Bars not summing to 100 for: ", label)
    print(bad)
  } else {
    message("✅ All bars sum to 100 for: ", label)
  }
}

# ----------------------------
# Kingdom5 mapping (from Phylum)
# ----------------------------
map_phylum_to_kingdom5 <- function(phylum_name) {
  x <- trimws(as.character(phylum_name))
  
  if (str_detect(x, regex("Uroviricota|virus|viruses", ignore_case = TRUE))) {
    return(if (isTRUE(KEEP_VIRUSES_AS_BIN)) "Viruses" else "Other")
  }
  if (str_detect(x, regex("Ascomycota|Basidiomycota|Mucoromycota|Chytridiomycota|Zoopagomycota|Glomeromycota|Microsporidia",
                          ignore_case = TRUE))) return("Fungi")
  if (str_detect(x, regex("Streptophyta|Chlorophyta", ignore_case = TRUE))) return("Plantae")
  if (str_detect(x, regex("Chordata|Arthropoda|Nematoda|Annelida|Mollusca|Cnidaria|Echinodermata|Platyhelminthes|Rotifera",
                          ignore_case = TRUE))) return("Animalia")
  if (str_detect(x, regex("Actinomycetota|Pseudomonadota|Bacteroidota|Bacillota|Acidobacteriota|Myxococcota|Planctomycetota|Chloroflexota|Verrucomicrobiota|Gemmatimonadota|Cyanobacteriota|Deinococcota|Euryarchaeota|Crenarchaeota|Thaumarchaeota|Nanoarchaeota|DPANN",
                          ignore_case = TRUE))) return("Monera")
  "Protista"
}

prep_kingdom5_from_phylum <- function(phylum_file, meta) {
  long <- tax_long_with_meta(phylum_file, meta) %>%
    filter(!is.na(name)) %>%
    mutate(name = as.character(name)) %>%
    remove_unclassified_like("name")
  
  grp <- long %>%
    group_by(LocationK, Inoculation, name) %>%
    summarise(Percent = mean(Percent, na.rm = TRUE), .groups = "drop") %>%
    complete(LocationK, Inoculation, name, fill = list(Percent = 0)) %>%
    group_by(LocationK, Inoculation) %>%
    mutate(tot = sum(Percent, na.rm = TRUE),
           Percent = ifelse(tot > 0, Percent / tot * 100, 0)) %>%
    ungroup() %>% dplyr::select(-tot) %>%
    mutate(Kingdom5 = map_chr(name, map_phylum_to_kingdom5)) %>%
    group_by(LocationK, Inoculation, Kingdom5) %>%
    summarise(Percent = sum(Percent, na.rm = TRUE), .groups = "drop") %>%
    complete(LocationK, Inoculation, Kingdom5, fill = list(Percent = 0)) %>%
    group_by(LocationK, Inoculation) %>%
    mutate(tot = sum(Percent, na.rm = TRUE),
           Percent = ifelse(tot > 0, Percent / tot * 100, 0)) %>%
    ungroup() %>% dplyr::select(-tot) %>%
    rename(Taxon = Kingdom5)
  
  grp
}

# ----------------------------
# Ordering: stack bottom=Other/bins, top=dominant; legend big->small
# ----------------------------
apply_nature_ordering <- function(df, other_label = "Other", is_genus = FALSE) {
  df2 <- df %>% mutate(Taxon = as.character(Taxon))
  
  mu <- df2 %>%
    group_by(Taxon) %>%
    summarise(mu = mean(Percent, na.rm = TRUE), .groups = "drop")
  
  taxa_all <- mu$Taxon
  other_like <- taxa_all[taxa_all == other_label | str_detect(taxa_all, "^Other")]
  
  bin_like <- character(0)
  if (is_genus) {
    bin_order_small_to_large <- c(
      "Other (<0.01%)",
      "Very Rare (0.01–0.05%)",
      "Ultra Rare (0.05–0.2%)",
      "Rare (0.2–0.3%)",
      "Very Low (0.3–0.4%)",
      "Low (0.4–0.5%)",
      ">=0.5% (not in Top50)"
    )
    bin_like <- intersect(bin_order_small_to_large, taxa_all)
  }
  
  named <- setdiff(taxa_all, c(other_like, bin_like))
  
  named_small_to_large <- mu %>%
    filter(Taxon %in% named) %>%
    arrange(mu) %>%
    pull(Taxon)
  
  named_large_to_small <- rev(named_small_to_large)
  
  stack_levels  <- unique(c(other_like, bin_like, named_small_to_large))
  legend_levels <- unique(c(named_large_to_small, rev(bin_like), other_like))
  
  df2 %>%
    mutate(
      Taxon_stack = factor(Taxon, levels = stack_levels),
      Taxon_leg   = factor(Taxon, levels = legend_levels)
    )
}

plot_faceted_composition <- function(df, title, legend_title, palette,
                                     legend_ncol = 1, base_size = 8) {
  
  stopifnot(all(c("LocationK", "Inoculation", "Taxon_stack", "Taxon_leg", "Percent") %in% names(df)))
  
  df2 <- df %>%
    mutate(
      LocationK   = factor(as.character(LocationK), levels = site_levels_k),
      Inoculation = factor(as.character(Inoculation), levels = c("Control", "Panch")),
      Taxon_leg   = factor(as.character(Taxon_leg), levels = levels(df$Taxon_leg)),
      Taxon_stack = factor(as.character(Taxon_stack), levels = levels(df$Taxon_stack))
    )
  
  lev_leg <- levels(df2$Taxon_leg)
  pal <- palette[lev_leg]
  
  if (any(is.na(pal))) pal[is.na(pal)] <- "grey75"
  names(pal) <- lev_leg
  
  ggplot(df2, aes(x = Inoculation, y = Percent)) +
    geom_col(
      aes(fill = Taxon_leg, group = Taxon_stack),
      width = 0.85, color = NA,
      position = position_stack(reverse = TRUE)
    ) +
    facet_wrap(~ LocationK, nrow = 1, drop = FALSE) +
    scale_y_continuous(limits = c(0, 100.0001),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_x_discrete(labels = c(Control = "Control", Panch = "Panch")) +
    scale_fill_manual(values = pal, breaks = lev_leg, drop = FALSE) +
    guides(fill = guide_legend(ncol = legend_ncol)) +
    labs(title = title, x = NULL, y = "Mean relative abundance (%)", fill = legend_title) +
    theme_natureish(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      legend.position = "right"
    )
}

# ----------------------------
# Global palettes
# ----------------------------
npg10 <- ggsci::pal_npg("nrc")(10)
npg10_rot <- npg10[c(3:10, 1:2)]

okabe_ito_no_black <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

bin_cols <- c(
  ">=0.5% (not in Top50)"  = "#222222",
  "Low (0.4–0.5%)"         = "#4D4D4D",
  "Very Low (0.3–0.4%)"    = "#7A7A7A",
  "Rare (0.2–0.3%)"        = "#A6A6A6",
  "Ultra Rare (0.05–0.2%)" = "#C2C2C2",
  "Very Rare (0.01–0.05%)" = "#D7D7D7",
  "Other (<0.01%)"         = "#E6E6E6"
)

kingdom5_cols <- c(
  "Monera"   = "#1b9e77",
  "Fungi"    = "#e7298a",
  "Protista" = "#7570b3",
  "Viruses"  = "#1f78b4",
  "Plantae"  = "#66a61e",
  "Animalia" = "#d95f02",
  "Other"    = "grey75"
)

phylum_npg_fixed <- c(
  "Pseudomonadota"    = npg10_rot[1],
  "Actinomycetota"    = npg10_rot[2],
  "Bacteroidota"      = npg10_rot[3],
  "Mucoromycota"      = npg10_rot[4],
  "Myxococcota"       = npg10_rot[5],
  "Ascomycota"        = npg10_rot[6],
  "Bacillota"         = npg10_rot[7],
  "Planctomycetota"   = npg10_rot[8],
  "Verrucomicrobiota" = npg10_rot[9],
  "Acidobacteriota"   = npg10_rot[10],
  "Other"             = "grey75"
)

make_name_palette_top10_npg_rest_ramp <- function(taxa_big_to_small, other_like = c("Other"), reserved_cols = NULL) {
  taxa <- unique(taxa_big_to_small)
  
  other_set <- taxa[taxa %in% other_like | str_detect(taxa, "^Other")]
  reserved_names <- if (is.null(reserved_cols)) character(0) else names(reserved_cols)
  
  named <- setdiff(taxa, c(other_set, reserved_names))
  top10 <- head(named, 10)
  rest  <- setdiff(named, top10)
  
  pal <- c()
  if (!is.null(reserved_cols)) pal <- c(pal, reserved_cols[names(reserved_cols) %in% taxa])
  
  if (length(top10) > 0) pal <- c(pal, setNames(npg10_rot[seq_along(top10)], top10))
  
  if (length(rest) > 0) {
    ramp <- grDevices::colorRampPalette(okabe_ito_no_black)(max(7, length(rest)))
    pal <- c(pal, setNames(ramp[seq_along(rest)], rest))
  }
  
  if (length(other_set) > 0) pal <- c(pal, setNames(rep("grey75", length(other_set)), other_set))
  
  pal
}

make_phylum_palette <- function(taxa_leg_levels) {
  lev <- as.character(taxa_leg_levels)
  pal <- phylum_npg_fixed
  if ("Other" %in% lev && !("Other" %in% names(pal))) pal["Other"] <- "grey75"
  
  missing <- setdiff(setdiff(lev, "Other"), names(pal))
  if (length(missing) > 0) {
    used <- unname(pal[names(pal) != "Other"])
    remaining_npg <- setdiff(npg10_rot, used)
    
    add <- character(0)
    if (length(remaining_npg) > 0) {
      take <- min(length(remaining_npg), length(missing))
      add <- setNames(remaining_npg[seq_len(take)], missing[seq_len(take)])
      missing <- missing[-seq_len(take)]
    }
    if (length(missing) > 0) {
      extra <- grDevices::colorRampPalette(okabe_ito_no_black)(max(7, length(missing)))
      add <- c(add, setNames(extra[seq_along(missing)], missing))
    }
    pal <- c(pal, add)
  }
  pal
}

global_order_big_to_small <- function(dfA, dfB, taxon_col = "Taxon_leg") {
  a <- dfA %>%
    mutate(Taxon = as.character(.data[[taxon_col]])) %>%
    dplyr::select(Taxon, Percent)
  b <- dfB %>%
    mutate(Taxon = as.character(.data[[taxon_col]])) %>%
    dplyr::select(Taxon, Percent)
  
  bind_rows(a, b) %>%
    group_by(Taxon) %>%
    summarise(mu = mean(Percent, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mu)) %>%
    pull(Taxon)
}

# ----------------------------
# Diversity stats (Genus)
# ----------------------------
wide_to_props <- function(file, meta) {
  df <- read_tax_table(file) %>% remove_human_contamination()
  mat <- df %>% column_to_rownames("name") %>% as.matrix()
  colnames(mat) <- norm_sampleid(colnames(mat))
  mat[is.na(mat)] <- 0
  
  keep <- intersect(meta$SampleID, colnames(mat))
  mat <- mat[, keep, drop = FALSE]
  
  cs <- colSums(mat)
  props <- sweep(mat, 2, cs, "/")
  props[is.na(props)] <- 0
  props
}

alpha_table <- function(file, meta, rank_label) {
  props <- wide_to_props(file, meta)
  shannon <- vegan::diversity(t(props), "shannon")
  simpson <- vegan::diversity(t(props), "simpson")
  observed <- colSums(props > 0)
  
  tibble(
    SampleID = names(shannon),
    Rank = rank_label,
    Observed = as.numeric(observed),
    Shannon  = as.numeric(shannon),
    Simpson  = as.numeric(simpson)
  ) %>%
    left_join(meta, by = "SampleID")
}

plot_alpha <- function(df, metric, title) {
  ggplot(df, aes(x = Inoculation, y = .data[[metric]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Location), width = 0.15, size = 2.2, alpha = 0.9) +
    facet_grid(. ~ Location, scales = "free_y") +
    labs(title = title, x = NULL, y = metric, color = "Location") +
    theme_natureish(base_size = 9)
}

ensure_factor <- function(x) if (is.factor(x)) x else factor(x)

plot_permdisp <- function(dist_obj, meta_df, factor_var, title, out_pdf, w = 180, h = 120) {
  if (!(factor_var %in% names(meta_df))) return(invisible(NULL))
  f <- ensure_factor(meta_df[[factor_var]])
  f <- droplevels(f)
  if (nlevels(f) < 2) return(invisible(NULL))
  
  bd <- vegan::betadisper(dist_obj, f)
  d  <- tibble(Distance = bd$distances, Level = f)
  
  p <- ggplot(d, aes(x = Level, y = Distance)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.85) +
    labs(title = title, x = factor_var, y = "Distance to centroid") +
    theme_natureish(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(out_pdf, p, width = w, height = h, units = "mm",
         dpi = 600, useDingbats = FALSE, limitsize = FALSE)
  invisible(list(betadisper = bd, plot = p))
}

pairwise_two_level_adonis <- function(dist_obj, meta_df, group_var = "Inoculation", permutations = 999) {
  if (!(group_var %in% names(meta_df))) return(NULL)
  g <- ensure_factor(meta_df[[group_var]])
  g <- droplevels(g)
  if (nlevels(g) != 2) return(NULL)
  
  dd <- data.frame(group = g)
  perm <- vegan::adonis2(dist_obj ~ group, data = dd, permutations = permutations)
  
  out <- as.data.frame(perm) %>%
    tibble::rownames_to_column("Term") %>%
    filter(.data$Term == "group") %>%
    transmute(
      Contrast = paste(levels(g), collapse = " vs "),
      Term = group_var,
      R2 = .data$R2,
      F  = .data$F,
      p  = .data$`Pr(>F)`
    )
  
  if (nrow(out) == 0) return(NULL)
  out
}

# ----------------------------
# Safe ellipse helper
# ----------------------------
make_ellipse_df <- function(pcoa_df) {
  grp_counts <- pcoa_df %>%
    count(Location, Inoculation, name = "n")
  
  pcoa_df %>%
    left_join(grp_counts, by = c("Location", "Inoculation")) %>%
    filter(n >= 3)
}

# ----------------------------
# Beta / PCoA
# ----------------------------
beta_pack_genus <- function(file, meta, fig_dir, tab_dir, prefix,
                            dataset_label = NULL,
                            permutations = 999, min_n_site = 6,
                            remove_pcoa_outlier_for_main = FALSE,
                            pcoa_outlier_ids = character(0),
                            main_plot_suffix = "") {
  
  props <- wide_to_props(file, meta)
  props_sqrt <- sqrt(props)
  bc <- vegan::vegdist(t(props_sqrt), method = "bray")
  
  labs <- attr(bc, "Labels")
  m <- meta %>%
    filter(.data$SampleID %in% labs) %>%
    mutate(
      Location    = factor(.data$Location, levels = site_levels),
      Inoculation = factor(.data$Inoculation, levels = inoc_levels)
    ) %>%
    as.data.frame()
  rownames(m) <- m$SampleID
  m <- droplevels(m)
  if (is.null(dataset_label)) dataset_label <- prefix
  
  # ----------------------------
  # FULL PCoA (all samples)
  # ----------------------------
  pcoa <- stats::cmdscale(bc, k = 2, eig = TRUE)
  pcoa_df <- as.data.frame(pcoa$points)
  colnames(pcoa_df) <- c("Axis1", "Axis2")
  pcoa_df$SampleID <- rownames(pcoa_df)
  pcoa_df <- left_join(pcoa_df, meta, by = "SampleID")
  
  pos_eig <- pcoa$eig[pcoa$eig > 0]
  var_expl <- pos_eig / sum(pos_eig)
  xlab <- paste0("PCoA1 (", round(var_expl[1] * 100, 1), "%)")
  ylab <- paste0("PCoA2 (", round(var_expl[2] * 100, 1), "%)")
  
  p_pcoa_full <- ggplot(pcoa_df, aes(Axis1, Axis2)) +
    geom_point(aes(color = Location, shape = Inoculation),
               size = 2.4, alpha = 0.9, stroke = 0.30) +
    scale_color_manual(values = loc_cols, drop = FALSE) +
    scale_shape_manual(values = c(Control = 16, Panch = 17), drop = FALSE) +
    labs(
      title = paste0(dataset_label, " PCoA (Bray–Curtis) — full dataset"),
      x = xlab, y = ylab,
      color = "Location", shape = "Inoculation"
    ) +
    theme_nature_ordination(base_size = 7)
  
  
  ggsave(
    file.path(fig_dir, paste0(prefix, "_PCoA_Genus_FULL_SUPPLEMENTARY.pdf")),
    p_pcoa_full, width = 190, height = 140, units = "mm",
    dpi = 600, useDingbats = FALSE, limitsize = FALSE
  )
  
  ggsave(
    file.path(fig_dir, paste0(prefix, "_PCoA_Genus_FULL_SUPPLEMENTARY.png")),
    p_pcoa_full, width = 190, height = 140, units = "mm",
    dpi = 600, limitsize = FALSE
  )
  
  write.table(
    pcoa_df,
    file.path(tab_dir, paste0(prefix, "_PCoA_coords_Genus.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  # ----------------------------
  # MAIN PCoA (optional outlier removal for visualisation only)
  # ----------------------------
  pcoa_df_main <- pcoa_df
  removed_ids <- character(0)
  
  if (remove_pcoa_outlier_for_main && length(pcoa_outlier_ids) > 0) {
    removed_ids <- intersect(pcoa_outlier_ids, pcoa_df_main$SampleID)
    pcoa_df_main <- pcoa_df_main %>%
      filter(!SampleID %in% removed_ids)
  }
  
  main_title <- paste0(dataset_label, " PCoA (Bray–Curtis)")
  if (length(removed_ids) > 0) {
    main_title <- paste0(
      dataset_label, " PCoA (Bray–Curtis) — main plot (outlier removed: ",
      paste(removed_ids, collapse = ", "), ")"
    )
  }
  
  p_pcoa_main <- ggplot(pcoa_df_main, aes(Axis1, Axis2)) +
    geom_point(aes(color = Location, shape = Inoculation),
               size = 2.4, alpha = 0.9, stroke = 0.30) +
    scale_color_manual(values = loc_cols, drop = FALSE) +
    scale_shape_manual(values = c(Control = 16, Panch = 17), drop = FALSE) +
    labs(
      title = main_title,
      x = xlab, y = ylab,
      color = "Location", shape = "Inoculation"
    ) +
    theme_nature_ordination(base_size = 7)
  
  
  ggsave(
    file.path(fig_dir, paste0(prefix, "_PCoA_Genus_MAIN", main_plot_suffix, ".pdf")),
    p_pcoa_main, width = 190, height = 140, units = "mm",
    dpi = 600, useDingbats = FALSE, limitsize = FALSE
  )
  
  ggsave(
    file.path(fig_dir, paste0(prefix, "_PCoA_Genus_MAIN", main_plot_suffix, ".png")),
    p_pcoa_main, width = 190, height = 140, units = "mm",
    dpi = 600, limitsize = FALSE
  )
  
  if (length(removed_ids) > 0) {
    write.table(
      tibble(SampleID = removed_ids),
      file.path(fig_dir, paste0(prefix, "_PCoA_MAIN_removed_samples", main_plot_suffix, ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  
  # ----------------------------
  # Global PERMANOVA (full dataset)
  # ----------------------------
  # ----------------------------
  # Global PERMANOVA (full dataset)
  # ----------------------------
  perm_global <- vegan::adonis2(
    bc ~ Location + Inoculation,
    data = m,
    permutations = permutations,
    by = "margin"
  )
  
  perm_global_df <- as.data.frame(perm_global) %>%
    tibble::rownames_to_column("Term")
  
  write.table(
    perm_global_df,
    file.path(tab_dir, paste0(prefix, "_PERMANOVA_GLOBAL_Location_Inoculation.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  # ----------------------------
  # PERMDISP plots (full dataset)
  # ----------------------------
  plot_permdisp(
    bc, m, "Inoculation",
    "PERMDISP (Genus, Global): Inoculation",
    file.path(fig_dir, paste0(prefix, "_PERMDISP_Global_Inoculation.pdf"))
  )
  plot_permdisp(
    bc, m, "Location",
    "PERMDISP (Genus, Global): Location",
    file.path(fig_dir, paste0(prefix, "_PERMDISP_Global_Location.pdf")),
    w = 220, h = 120
  )
  
  # ----------------------------
  # Within-site PERMANOVA (full dataset)
  # ----------------------------
  loc_levels <- levels(droplevels(m$Location))
  
  per_site <- purrr::map_dfr(loc_levels, function(loc) {
    ids <- rownames(m)[m$Location == loc]
    if (length(ids) < min_n_site) return(NULL)
    
    bc_loc <- as.dist(as.matrix(bc)[ids, ids])
    m_loc  <- droplevels(m[ids, , drop = FALSE])
    if (nlevels(droplevels(m_loc$Inoculation)) < 2) return(NULL)
    
    perm <- vegan::adonis2(bc_loc ~ Inoculation, data = m_loc, permutations = permutations)
    
    as.data.frame(perm) %>%
      tibble::rownames_to_column("Term") %>%
      mutate(Location = loc, Model = "Inoculation_only", N = length(ids))
  })
  
  if (nrow(per_site) == 0) {
    per_site <- tibble(
      Term = character(), Df = numeric(), SumOfSqs = numeric(), R2 = numeric(),
      F = numeric(), `Pr(>F)` = numeric(), Location = character(),
      Model = character(), N = integer()
    )
  }
  
  write.table(
    per_site,
    file.path(tab_dir, paste0(prefix, "_PERMANOVA_WITHIN_Location_Inoculation.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  # ----------------------------
  # Pairwise within site (full dataset)
  # ----------------------------
  pairwise_site <- purrr::map_dfr(loc_levels, function(loc) {
    ids <- rownames(m)[m$Location == loc]
    if (length(ids) < min_n_site) return(NULL)
    
    bc_loc <- as.dist(as.matrix(bc)[ids, ids])
    m_loc  <- droplevels(m[ids, , drop = FALSE])
    
    pw <- pairwise_two_level_adonis(bc_loc, m_loc, "Inoculation", permutations = permutations)
    if (is.null(pw)) return(NULL)
    pw %>% mutate(Location = loc, N = length(ids))
  })
  
  if (nrow(pairwise_site) == 0) {
    pairwise_site <- tibble(
      Contrast = character(), Term = character(), R2 = numeric(),
      F = numeric(), p = numeric(), Location = character(), N = integer()
    )
  }
  
  write.table(
    pairwise_site,
    file.path(tab_dir, paste0(prefix, "_PAIRWISE_Inoculation_WITHIN_Location.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  invisible(list(
    pcoa_plot_main = p_pcoa_main,
    pcoa_plot_full = p_pcoa_full,
    pcoa_df = pcoa_df,
    perm_global = perm_global,
    per_site = per_site,
    pairwise = pairwise_site,
    removed_ids = removed_ids
  ))
}

# ----------------------------
# Run one dataset
# ----------------------------
DEBUG_PRINT_BAR_TOTALS <- FALSE

run_one_dataset <- function(dataset_name, cfg) {
  message("\n============================\nRunning: ", dataset_name, "\n============================")
  
  data_dir <- cfg$data_dir
  meta_path   <- file.path(data_dir, cfg$meta_file)
  phylum_path <- file.path(data_dir, cfg$phylum_file)
  family_path <- file.path(data_dir, cfg$family_file)
  genus_path  <- file.path(data_dir, cfg$genus_file)
  
  stopifnot(file.exists(meta_path), file.exists(phylum_path), file.exists(family_path), file.exists(genus_path))
  
  out_base <- file.path(OUT_DIR, tolower(dataset_name))
  fig_dir  <- file.path(out_base, "figures")
  tab_dir  <- file.path(out_base, "tables")
  log_dir  <- file.path(out_base, "logs")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  
  meta <- read_meta(meta_path)
  
  king_df <- prep_kingdom5_from_phylum(phylum_path, meta) %>%
    apply_nature_ordering(other_label = "Other", is_genus = FALSE)
  
  phy_df <- prep_group_mean(phylum_path, meta, top_n = TOPN$Phylum, other_label = "Other") %>%
    apply_nature_ordering(other_label = "Other", is_genus = FALSE)
  
  if (isTRUE(DEBUG_PRINT_BAR_TOTALS)) {
    print(
      phy_df %>%
        group_by(LocationK, Inoculation) %>%
        summarise(total = sum(Percent), .groups = "drop") %>%
        arrange(desc(total))
    )
  }
  
  fam_df <- prep_group_mean(family_path, meta, top_n = TOPN$Family, other_label = "Other") %>%
    apply_nature_ordering(other_label = "Other", is_genus = FALSE)
  
  gen_df <- if (USE_GENUS_RARITY_BINS) {
    prep_genus_with_rarity_bins(
      genus_path, meta,
      top_n = TOPN$Genus,
      bins = GENUS_BINS,
      other_label = GENUS_OTHER_LABEL
    ) %>%
      apply_nature_ordering(other_label = GENUS_OTHER_LABEL, is_genus = TRUE)
  } else {
    prep_group_mean(genus_path, meta, top_n = TOPN$Genus, other_label = "Other") %>%
      apply_nature_ordering(other_label = "Other", is_genus = FALSE)
  }
  
  check_bar_sums(king_df, paste(dataset_name, "Kingdom5"))
  check_bar_sums(phy_df,  paste(dataset_name, "Phylum"))
  check_bar_sums(fam_df,  paste(dataset_name, "Family"))
  check_bar_sums(gen_df,  paste(dataset_name, "Genus"))
  
  list(
    dataset_name = dataset_name,
    meta = meta,
    paths = list(meta = meta_path, phylum = phylum_path, family = family_path, genus = genus_path),
    out = list(base = out_base, fig = fig_dir, tab = tab_dir, log = log_dir),
    king_df = king_df,
    phy_df = phy_df,
    fam_df = fam_df,
    gen_df = gen_df
  )
}

# ----------------------------
# 1) Load both datasets first
# ----------------------------
obj_root <- run_one_dataset("Root", DATASETS$Root)
obj_rhiz <- run_one_dataset("Rhizosphere", DATASETS$Rhizosphere)

# ----------------------------
# 2) Build GLOBAL palettes by taxon name (union across both)
# ----------------------------
pal_king <- kingdom5_cols

phy_global_big <- global_order_big_to_small(obj_root$phy_df, obj_rhiz$phy_df, "Taxon_leg")
fam_global_big <- global_order_big_to_small(obj_root$fam_df, obj_rhiz$fam_df, "Taxon_leg")
gen_global_big <- global_order_big_to_small(obj_root$gen_df, obj_rhiz$gen_df, "Taxon_leg")

pal_phy <- make_phylum_palette(phy_global_big)
pal_fam <- make_name_palette_top10_npg_rest_ramp(fam_global_big, other_like = c("Other"))

pal_gen <- make_name_palette_top10_npg_rest_ramp(
  gen_global_big,
  other_like = c(GENUS_OTHER_LABEL, "Other"),
  reserved_cols = bin_cols
)
for (nm in names(pal_gen)[str_detect(names(pal_gen), "^Other")]) pal_gen[nm] <- "grey75"

# ----------------------------
# 3) Plot + save (each compartment) using GLOBAL palettes
# ----------------------------
plot_and_save_all <- function(obj, dataset_title_prefix) {
  fig_dir <- obj$out$fig
  tab_dir <- obj$out$tab
  log_dir <- obj$out$log
  
  # Kingdom
  pK <- plot_faceted_composition(
    obj$king_df,
    title = paste0(dataset_title_prefix, " community composition (Kingdom5 from Phylum)"),
    legend_title = "Kingdom",
    palette = pal_king,
    legend_ncol = LEGEND_NCOL_KINGDOM
  )
  save_pdf(file.path(fig_dir, "Fig_Kingdom5_FACET_Inoculation_by_LocationK.pdf"), pK, 260, 125)
  save_png(file.path(fig_dir, "Fig_Kingdom5_FACET_Inoculation_by_LocationK.png"), pK, 260, 125)
  
  # Phylum
  pP <- plot_faceted_composition(
    obj$phy_df,
    title = paste0(dataset_title_prefix, " community composition (Phylum) — Top ", TOPN$Phylum),
    legend_title = "Phylum",
    palette = pal_phy,
    legend_ncol = LEGEND_NCOL_PHYLUM
  )
  save_pdf(file.path(fig_dir, "Fig_Phylum_Top10_FACET_Inoculation_by_LocationK.pdf"), pP, 300, 140)
  save_png(file.path(fig_dir, "Fig_Phylum_Top10_FACET_Inoculation_by_LocationK.png"), pP, 300, 140)
  
  # Family
  pF <- plot_faceted_composition(
    obj$fam_df,
    title = paste0(dataset_title_prefix, " community composition (Family) — Top ", TOPN$Family),
    legend_title = "Family",
    palette = pal_fam,
    legend_ncol = LEGEND_NCOL_FAMILY
  )
  save_pdf(file.path(fig_dir, "Fig_Family_Top30_FACET_Inoculation_by_LocationK.pdf"), pF, 320, 170)
  save_png(file.path(fig_dir, "Fig_Family_Top30_FACET_Inoculation_by_LocationK.png"), pF, 320, 170)
  
  # Genus
  pG <- plot_faceted_composition(
    obj$gen_df,
    title = paste0(
      dataset_title_prefix, " community composition (Genus) — Top ", TOPN$Genus,
      if (USE_GENUS_RARITY_BINS) " + rarity bins" else ""
    ),
    legend_title = "Genus",
    palette = pal_gen,
    legend_ncol = LEGEND_NCOL_GENUS
  )
  save_pdf(file.path(fig_dir, "Fig_Genus_Top50_FACET_Inoculation_by_LocationK.pdf"), pG, 340, 190)
  save_png(file.path(fig_dir, "Fig_Genus_Top50_FACET_Inoculation_by_LocationK.png"), pG, 340, 190)
  
  # Alpha
  alpha_gen <- alpha_table(obj$paths$genus, obj$meta, "Genus")
  write.table(alpha_gen, file.path(tab_dir, "Table_AlphaDiversity_Genus.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  ggsave(
    file.path(fig_dir, "Fig_Alpha_Genus_Observed.pdf"),
    plot_alpha(alpha_gen, "Observed", "Alpha diversity (Genus) — Observed richness"),
    width = 220, height = 120, units = "mm", dpi = 600,
    useDingbats = FALSE, limitsize = FALSE
  )
  
  ggsave(
    file.path(fig_dir, "Fig_Alpha_Genus_Shannon.pdf"),
    plot_alpha(alpha_gen, "Shannon", "Alpha diversity (Genus) — Shannon"),
    width = 220, height = 120, units = "mm", dpi = 600,
    useDingbats = FALSE, limitsize = FALSE
  )
  
  ggsave(
    file.path(fig_dir, "Fig_Alpha_Genus_Simpson.pdf"),
    plot_alpha(alpha_gen, "Simpson", "Alpha diversity (Genus) — Simpson"),
    width = 220, height = 120, units = "mm", dpi = 600,
    useDingbats = FALSE, limitsize = FALSE
  )
  
  # Beta / PERMANOVA / Pairwise
  is_rhizo <- identical(dataset_title_prefix, "Rhizosphere")
  
  beta <- beta_pack_genus(
    obj$paths$genus,
    obj$meta,
    fig_dir,
    tab_dir,
    prefix = "Genus",
    dataset_label = dataset_title_prefix,
    remove_pcoa_outlier_for_main = is_rhizo && REMOVE_RHIZO_PCOA_OUTLIER_FOR_MAIN,
    pcoa_outlier_ids = if (is_rhizo) RHIZO_PCOA_OUTLIER_IDS else character(0),
    main_plot_suffix = if (is_rhizo) "_Rhizosphere" else "_Root"
  )
  
  # Logs
  settings <- list(
    TOPN = TOPN,
    USE_GENUS_RARITY_BINS = USE_GENUS_RARITY_BINS,
    GENUS_BINS = GENUS_BINS,
    GENUS_OTHER_LABEL = GENUS_OTHER_LABEL,
    KEEP_VIRUSES_AS_BIN = KEEP_VIRUSES_AS_BIN,
    REMOVE_VERTEBRATES = REMOVE_VERTEBRATES,
    site_levels = site_levels,
    site_levels_k = site_levels_k,
    inoc_levels = inoc_levels,
    REMOVE_RHIZO_PCOA_OUTLIER_FOR_MAIN = REMOVE_RHIZO_PCOA_OUTLIER_FOR_MAIN,
    RHIZO_PCOA_OUTLIER_IDS = RHIZO_PCOA_OUTLIER_IDS
  )
  writeLines(capture.output(str(settings)), file.path(log_dir, "Settings_snapshot.txt"))
  writeLines(capture.output(sessionInfo()), file.path(log_dir, "SessionInfo.txt"))
  
  hashes <- tibble(
    file = c(obj$paths$meta, obj$paths$phylum, obj$paths$family, obj$paths$genus),
    sha1 = map_chr(file, sha1_file)
  )
  write.table(hashes, file.path(log_dir, "Input_file_hashes.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  list(
    pK = pK, pP = pP, pF = pF, pG = pG,
    pPCOA_main = beta$pcoa_plot_main,
    pPCOA_full = beta$pcoa_plot_full,
    beta = beta
  )
}

plots_root <- plot_and_save_all(obj_root, "Root")
plots_rhiz <- plot_and_save_all(obj_rhiz, "Rhizosphere")

# ============================================================
# Supplementary taxonomy figures — NatCom style
# Merged panels:
#   Figure S1: Rhizosphere taxonomy composition
#      A Phylum-level composition
#      B Family-level composition
#   Figure S2: Root taxonomy composition
#      A Phylum-level composition
#      B Family-level composition
#
# Purpose:
#   Merge the previous separate phylum/family supplementary plots
#   into two multi-panel figures with larger fonts and cleaner spacing.
# ============================================================

if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
library(cowplot)

dir.create(SUPP_TAX_DIR, recursive = TRUE, showWarnings = FALSE)

theme_supp_taxonomy <- function(base_size = 11) {
  theme_natureish(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      strip.text = element_text(size = base_size, face = "bold"),
      axis.title.y = element_text(size = base_size + 1, face = "bold"),
      axis.text.y = element_text(size = base_size),
      axis.text.x = element_text(size = base_size, colour = "black"),
      legend.title = element_text(size = base_size + 1, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(4.2, "mm"),
      legend.spacing.y = unit(1.5, "mm"),
      legend.box.spacing = unit(3, "mm"),
      panel.spacing.x = unit(2.5, "mm"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

polish_supp_tax_panel <- function(p, panel_title, legend_ncol = 1, base_size = 11) {
  p +
    labs(
      title = panel_title,
      x = NULL,
      y = "Mean relative abundance (%)"
    ) +
    guides(fill = guide_legend(ncol = legend_ncol, byrow = TRUE, title.position = "top")) +
    theme_supp_taxonomy(base_size = base_size) +
    theme(
      legend.position = "right"
    )
}

# Rhizosphere supplementary Figure S1: A Phylum, B Family
p_supp_rhiz_phylum <- polish_supp_tax_panel(
  plots_rhiz$pP,
  panel_title = "Phylum-level composition",
  legend_ncol = 1,
  base_size = 11
)

p_supp_rhiz_family <- polish_supp_tax_panel(
  plots_rhiz$pF,
  panel_title = "Family-level composition",
  legend_ncol = 2,
  base_size = 11
)

figS1_rhiz_taxonomy <- cowplot::plot_grid(
  p_supp_rhiz_phylum,
  p_supp_rhiz_family,
  ncol = 1,
  labels = c("A", "B"),
  label_size = 18,
  label_fontface = "bold",
  rel_heights = c(1.0, 1.25),
  align = "v",
  axis = "lr"
)

ggsave(
  file.path(SUPP_TAX_DIR, "FigureS1_Rhizosphere_Taxonomic_Composition_Merged_NatCom.pdf"),
  figS1_rhiz_taxonomy,
  width = 320,
  height = 250,
  units = "mm",
  dpi = 600,
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  file.path(SUPP_TAX_DIR, "FigureS1_Rhizosphere_Taxonomic_Composition_Merged_NatCom.png"),
  figS1_rhiz_taxonomy,
  width = 320,
  height = 250,
  units = "mm",
  dpi = 600,
  limitsize = FALSE
)

ggsave(
  file.path(SUPP_TAX_DIR, "FigureS1_Rhizosphere_Taxonomic_Composition_Merged_NatCom.svg"),
  figS1_rhiz_taxonomy,
  width = 320,
  height = 250,
  units = "mm",
  limitsize = FALSE
)

# Root supplementary Figure S2: A Phylum, B Family
p_supp_root_phylum <- polish_supp_tax_panel(
  plots_root$pP,
  panel_title = "Phylum-level composition",
  legend_ncol = 1,
  base_size = 11
)

p_supp_root_family <- polish_supp_tax_panel(
  plots_root$pF,
  panel_title = "Family-level composition",
  legend_ncol = 2,
  base_size = 11
)

figS2_root_taxonomy <- cowplot::plot_grid(
  p_supp_root_phylum,
  p_supp_root_family,
  ncol = 1,
  labels = c("A", "B"),
  label_size = 18,
  label_fontface = "bold",
  rel_heights = c(1.0, 1.25),
  align = "v",
  axis = "lr"
)

ggsave(
  file.path(SUPP_TAX_DIR, "FigureS2_Root_Taxonomic_Composition_Merged_NatCom.pdf"),
  figS2_root_taxonomy,
  width = 320,
  height = 250,
  units = "mm",
  dpi = 600,
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  file.path(SUPP_TAX_DIR, "FigureS2_Root_Taxonomic_Composition_Merged_NatCom.png"),
  figS2_root_taxonomy,
  width = 320,
  height = 250,
  units = "mm",
  dpi = 600,
  limitsize = FALSE
)

ggsave(
  file.path(SUPP_TAX_DIR, "FigureS2_Root_Taxonomic_Composition_Merged_NatCom.svg"),
  figS2_root_taxonomy,
  width = 320,
  height = 250,
  units = "mm",
  limitsize = FALSE
)

saveRDS(figS1_rhiz_taxonomy, file.path(SUPP_TAX_DIR, "FigureS1_Rhizosphere_Taxonomic_Composition_Merged_NatCom.rds"))
saveRDS(figS2_root_taxonomy, file.path(SUPP_TAX_DIR, "FigureS2_Root_Taxonomic_Composition_Merged_NatCom.rds"))

message("NatCom supplementary taxonomy figures written to: ", SUPP_TAX_DIR)


# ============================================================
# Manuscript Figure 2 — taxonomic response + genus composition
# NatCom style
# Panels:
# A Rhizosphere PCoA
# B Root PCoA
# C PERMANOVA inoculation effect
# D Rhizosphere genus composition
# E Root genus composition
# ============================================================

if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
library(cowplot)

NATCOM_DIR <- NATCOM_DIR
dir.create(NATCOM_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# A–B: PCoA panels
# ----------------------------
p_fig2_rhiz <- plots_rhiz$pPCOA_main +
  labs(title = "Rhizosphere", colour = "Location", shape = "Inoculation") +
  theme_nature_ordination(base_size = 7) +
  theme(legend.position = "none")

p_fig2_root <- plots_root$pPCOA_main +
  labs(title = "Root", colour = "Location", shape = "Inoculation") +
  theme_nature_ordination(base_size = 7)

# ----------------------------
# C: PERMANOVA treatment R2
# ----------------------------
perm_fig2 <- tibble::tibble(
  Compartment = factor(c("Rhizosphere", "Root"), levels = c("Rhizosphere", "Root")),
  Treatment_R2 = c(0.013, 0.066),
  Label = c("R² = 0.013\nP = 0.246\nns",
            "R² = 0.066\nP = 0.002\n**")
)

p_fig2_perm <- ggplot(perm_fig2, aes(x = Compartment, y = Treatment_R2)) +
  geom_col(width = 0.55, fill = c("grey75", "#0072B2"),
           colour = "black", linewidth = 0.25) +
  geom_text(aes(label = Label), vjust = -0.20, size = 2.2, lineheight = 0.9) +
  scale_y_continuous(limits = c(0, 0.085), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "Inoculation effect", x = NULL, y = "Treatment R²") +
  theme_nature_ordination(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),
    legend.position = "none"
  )

# ----------------------------
# D–E: genus stacked composition panels
# ----------------------------
p_fig2D_rhiz_genus <- plots_rhiz$pG +
  labs(
    title = "Rhizosphere genus composition",
    x = NULL,
    y = "Mean relative abundance (%)",
    fill = "Genus"
  ) +
  theme_natureish(base_size = 7) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 9,
      hjust = 0
    ),
    strip.text = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(size = 6, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 7, face = "bold"),
    legend.position = "none",
    panel.spacing.x = unit(1, "mm"),
    plot.margin = margin(3, 3, 3, 3)
  )

p_fig2E_root_genus <- plots_root$pG +
  labs(
    title = "Root genus composition",
    x = NULL,
    y = "Mean relative abundance (%)",
    fill = "Genus"
  ) +
  theme_natureish(base_size = 7) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 9,
      hjust = 0
    ),
    strip.text = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(size = 6, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 7, face = "bold"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 6),
    legend.key.height = unit(1.8, "mm"),
    legend.key.width = unit(1.8, "mm"),
    legend.position = "right",
    panel.spacing.x = unit(1, "mm"),
    plot.margin = margin(3, 3, 3, 3)
  )

# ----------------------------
# Top row: A–C
# ----------------------------
fig2_top <- cowplot::plot_grid(
  p_fig2_rhiz,
  p_fig2_root,
  p_fig2_perm,
  nrow = 1,
  labels = c("A", "B", "C"),
  label_size = 12,
  label_fontface = "bold",
  rel_widths = c(1, 1.25, 0.65),
  align = "h",
  axis = "tb"
)

# ----------------------------
# Bottom row: D–E
# ----------------------------
fig2_bottom <- cowplot::plot_grid(
  p_fig2D_rhiz_genus,
  p_fig2E_root_genus,
  nrow = 1,
  labels = c("D", "E"),
  label_size = 12,
  label_fontface = "bold",
  rel_widths = c(1, 1.4),
  align = "h",
  axis = "tb"
)

# ----------------------------
# Final merged Figure 2
# ----------------------------
fig2_taxonomy_response <- cowplot::plot_grid(
  fig2_top,
  fig2_bottom,
  ncol = 1,
  rel_heights = c(1.0, 0.75),
  align = "v"
)

ggsave(
  file.path(MAIN_FIG_DIR, "Figure2_Taxonomic_Response_NatCom.png"),
  fig2_taxonomy_response,
  width = 350,
  height = 200,
  units = "mm",
  dpi = 600,
  limitsize = FALSE
)

ggsave(
  file.path(MAIN_FIG_DIR, "Figure2_Taxonomic_Response_NatCom.svg"),
  fig2_taxonomy_response,
  width = 350,
  height = 200,
  units = "mm",
  limitsize = FALSE
)

ggsave(
  file.path(MAIN_FIG_DIR, "Figure2_Taxonomic_Response_NatCom.pdf"),
  fig2_taxonomy_response,
  width = 350,
  height = 200,
  units = "mm",
  useDingbats = FALSE,
  limitsize = FALSE
)

# Save RDS panels for later editing
saveRDS(p_fig2_rhiz, file.path(NATCOM_DIR, "p_Fig2A_Rhizosphere_PCoA.rds"))
saveRDS(p_fig2_root, file.path(NATCOM_DIR, "p_Fig2B_Root_PCoA.rds"))
saveRDS(p_fig2_perm, file.path(NATCOM_DIR, "p_Fig2C_PERMANOVA.rds"))
saveRDS(p_fig2D_rhiz_genus, file.path(NATCOM_DIR, "p_Fig2D_Rhizosphere_Genus_Composition.rds"))
saveRDS(p_fig2E_root_genus, file.path(NATCOM_DIR, "p_Fig2E_Root_Genus_Composition.rds"))

message("Manuscript-ready Figure 2 written to: ", NATCOM_DIR)

message("\nDONE ✅")
message("Outputs written to: ", OUT_DIR)
message("Manuscript-ready Figure 2 written to: ", NATCOM_DIR)

})


run_section('Step 3 — Merge Figure 1 panels after map, taxonomy PCoA, and pathway PCoA panels exist', required = FALSE, code = function() {
# ============================================================
# Merge Figure 1 panels: A map, B taxonomy PCoA, C pathway PCoA
# ============================================================

# install.packages(c("magick", "cowplot", "ggplot2"))

library(magick)
library(cowplot)
library(ggplot2)
library(grid)

out_dir <- NATCOM_DIR

# Input files
map_file <- file.path(out_dir, "Fig1A_Study_Sites_Map_with_Canterbury_Inset.png")
tax_file <- file.path(out_dir, "Fig1B_Taxonomy_PCoA_NatCom_panel.png")
path_file <- file.path(out_dir, "Fig1C_Pathway_PCoA_NatCom_panel.png")

# Read images
map_img  <- image_read(map_file)
tax_img  <- image_read(tax_file)
path_img <- image_read(path_file)

# Convert to grobs
map_grob  <- rasterGrob(as.raster(map_img),  interpolate = TRUE)
tax_grob  <- rasterGrob(as.raster(tax_img),  interpolate = TRUE)
path_grob <- rasterGrob(as.raster(path_img), interpolate = TRUE)

# Convert grobs to ggdraw objects
pA <- ggdraw() + draw_grob(map_grob)
pB <- ggdraw() + draw_grob(tax_grob)
pC <- ggdraw() + draw_grob(path_grob)

# Combine panels
fig1 <- plot_grid(
  pA,
  plot_grid(pB, pC, nrow = 1, labels = c("B", "C"), label_size = 14),
  ncol = 1,
  labels = c("A", ""),
  label_size = 14,
  rel_heights = c(0.8, 1.2)
)

# Save final figure
ggsave(
  file.path(out_dir, "Figure1_NatCom_merged.png"),
  fig1,
  width = 180,
  height = 180,
  units = "mm",
  dpi = 600
)

ggsave(
  file.path(out_dir, "Figure1_NatCom_merged.pdf"),
  fig1,
  width = 180,
  height = 180,
  units = "mm",
  useDingbats = FALSE
)

})

message("\nTaxonomy standalone pipeline complete.")
