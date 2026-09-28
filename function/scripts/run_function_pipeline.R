# ============================================================
# Chapter 2 metagenomics — Functional potential pipeline
# STANDALONE NatCom/GitHub-ready script
#
# Run from repository root:
#   setwd("C:/Users/soths/OneDrive - Lincoln University/Writing/NatCom/Github")
#   source("function/scripts/run_function_pipeline.R")
# ============================================================

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(repo_root, "function", "data"))) {
  stop("Run this script from the repository root containing function/data", call. = FALSE)
}

STAGING_DIR <- file.path(repo_root, "function", "results", "merged")
MAIN_FIG_DIR <- file.path(repo_root, "manuscript_figures", "Main_Figures")
SUPP_FUNC_DIR <- file.path(repo_root, "manuscript_figures", "Supplementary_Function_Figures")
dir.create(STAGING_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MAIN_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_FUNC_DIR, recursive = TRUE, showWarnings = FALSE)

# Text sizes for panels that are assembled into multi-panel manuscript figures.
# These values improve readability after reduction without changing the data,
# scales, panel arrangement, or final export dimensions.
MERGED_BASE_SIZE         <- 9
MERGED_TITLE_SIZE        <- 11
MERGED_LEGEND_TITLE_SIZE <- 10
MERGED_LEGEND_TEXT_SIZE  <- 9
MERGED_TAG_SIZE          <- 14

# Figure 1C is exported without a legend because Figure 1 uses one shared
# legend to the right of Panels B and C. Keep these dimensions and text sizes
# identical to the Figure 1B settings in run_taxonomy_pipeline.R.
FIG1_PCOA_BASE_SIZE         <- 11
FIG1_PCOA_AXIS_TITLE_SIZE   <- 13
FIG1_PCOA_AXIS_TEXT_SIZE    <- 11
FIG1_PCOA_EXPORT_WIDTH_MM   <- 85
FIG1_PCOA_EXPORT_HEIGHT_MM  <- 70
FIG1_PANEL_VERSION          <- "Functional pipeline v20: location and environment interaction models"
message(FIG1_PANEL_VERSION)

# Pathway screens are exploratory: nominal ALDEx2 we.ep < 0.10 and
# |effect| > 0.25. BH-adjusted we.eBH < 0.05 is reported separately.
# Do not describe pathways selected by the nominal screen as FDR significant.
PATHWAY_NOMINAL_P <- 0.10
PATHWAY_EFFECT_MIN <- 0.25
PATHWAY_ADJUSTED_Q <- 0.05

# Figure 1C and the site-specific volcanoes show six field environments. The
# global PERMANOVA and cross-location consistency analyses retain five sites.
FUNCTION_ENV_LEVELS <- c(
  "Eyrewell_Forest", "Kowhai_Irrigated", "Kowhai_Rainfed",
  "LU_H8", "Rolleston", "West_Coast"
)

# Figures 1C, 3 and 4 display six field environments. The five-location
# loc_cols palette below remains available for site-level analyses.
functional_environment_cols <- c(
  "Eyrewell_Forest"   = "#D55E00",
  "Kowhai_Irrigated"  = "#8C510A",
  "Kowhai_Rainfed"    = "#E69F00",
  "LU_H8"             = "#009E73",
  "Rolleston"         = "#0072B2",
  "West_Coast"        = "#CC79A7"
)
pretty_environment_label <- function(x) {
  dplyr::recode(
    as.character(x),
    "Eyrewell_Forest" = "Eyrewell Forest",
    "Kowhai_Irrigated" = "Kowhai (Irrigated)",
    "Kowhai_Rainfed" = "Kowhai (Rainfed)",
    "West_Coast" = "West Coast",
    .default = as.character(x)
  )
}
fig4_environment_axis_label <- function(x) {
  y <- pretty_environment_label(x)
  y[y == "Eyrewell Forest"] <- "Eyrewell"
  y[y == "Kowhai (Irrigated)"] <- "Kowhai Irr."
  y[y == "Kowhai (Rainfed)"] <- "Kowhai Rain."
  y
}

# Import code accepts the original workbook treatment name and standardises it
# to Trichoderma before modelling. Plotmath italicises the genus in figures.
trichoderma_inoculation_labels <- function(x) {
  parse(text = ifelse(x == "Trichoderma", "italic(Trichoderma)", "'Control'"))
}
trichoderma_direction_labels <- function(x) {
  parse(text = ifelse(
    x == "Higher in Trichoderma",
    "'Higher in ' * italic(Trichoderma)",
    paste0("'", x, "'")
  ))
}
trichoderma_title <- function(x) {
  parts <- strsplit(x, "Trichoderma", fixed = TRUE)[[1]]
  if (length(parts) == 1L) return(x)
  out <- parts[[1]]
  for (part in parts[-1]) {
    out <- call("*", call("*", out, quote(italic(Trichoderma))), part)
  }
  out
}

# Figure 3 matches the compact PCoA geometry and typography used in the top
# row of taxonomy Figure 2. The A and B plotting regions are taller than wide,
# leaving enough room for readable labels and the shared legend in Panel B.
FIG3_PCOA_BASE_SIZE         <- 15
FIG3_PCOA_TITLE_SIZE        <- 19
FIG3_PCOA_AXIS_TITLE_SIZE   <- 17
FIG3_PCOA_AXIS_TEXT_SIZE    <- 15
FIG3_PCOA_LEGEND_TITLE_SIZE <- 16
FIG3_PCOA_LEGEND_TEXT_SIZE  <- 14
FIG3_PCOA_ASPECT_RATIO      <- 1.18
FIG3_TAG_SIZE               <- 20
FIG3_WIDTH_MM               <- 350
FIG3_HEIGHT_MM              <- 115

# Figure 4 is exported on a wide 530-mm canvas and is subsequently reduced for
# the manuscript. These larger source sizes remain readable after reduction.
FIG4_SOURCE_BASE_SIZE    <- 11
FIG4_PATHWAY_TEXT_SIZE   <- 17
FIG4_SITE_TEXT_SIZE      <- 19
FIG4_TITLE_SIZE          <- 20
FIG4_LEGEND_TITLE_SIZE   <- 19
FIG4_LEGEND_TEXT_SIZE    <- 18
FIG4_TAG_SIZE            <- 30
FIG4_FIELD_TEXT_SIZE     <- 16

# Figure 5 typography is set at the source-figure scale so the labels remain
# readable when the two contributor panels are reduced for the manuscript.
FIG5_SOURCE_BASE_SIZE    <- 12
FIG5_PATHWAY_TEXT_SIZE   <- 13
FIG5_PATHWAY_TEXT_ANGLE  <- 55
FIG5_TAXON_TEXT_SIZE     <- 13
FIG5_TITLE_SIZE          <- 18
FIG5_LEGEND_TITLE_SIZE   <- 14
FIG5_LEGEND_TEXT_SIZE    <- 12
FIG5_TAG_SIZE            <- 24
FIG5_EXPORT_WIDTH_MM     <- 380
FIG5_EXPORT_HEIGHT_MM    <- 170

# Supplementary Figure S5 keeps the quantitative colour and size encodings,
# but uses a visible minimum point size, a darker blue endpoint and outlines.
FIGS5_POINT_SIZE_RANGE   <- c(2.6, 8)
FIGS5_LOW_COLOUR         <- "#08306B"
FIGS5_OUTLINE_COLOUR     <- "grey25"
FIGS5_OUTLINE_STROKE     <- 0.35

message("Repository root: ", repo_root)
message("Staging/merged output: ", STAGING_DIR)
message("Main manuscript figures: ", MAIN_FIG_DIR)
message("Supplementary function figures: ", SUPP_FUNC_DIR)

# Avoid namespace masking conflicts (e.g., MASS::select vs dplyr::select)
if (requireNamespace("dplyr", quietly = TRUE)) {
  select    <- dplyr::select
  filter    <- dplyr::filter
  mutate    <- dplyr::mutate
  rename    <- dplyr::rename
  arrange   <- dplyr::arrange
  summarise <- dplyr::summarise
  left_join <- dplyr::left_join
  bind_rows <- dplyr::bind_rows
}

run_section <- function(label, code) {
  message("\n============================================================")
  message(label)
  message("============================================================")
  force(code)
}

copy_matching <- function(from_dir, patterns, to_dir) {
  if (!dir.exists(from_dir)) return(invisible(character()))
  files <- list.files(from_dir, recursive = TRUE, full.names = TRUE)
  keep <- Reduce(`|`, lapply(patterns, function(p) grepl(p, basename(files), ignore.case = FALSE)))
  files <- files[keep]
  if (length(files)) {
    dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(files, file.path(to_dir, basename(files)), overwrite = TRUE)
    message("Copied ", sum(ok), " file(s) to ", to_dir)
  }
  invisible(files)
}

copy_outputs <- function() {
  main_patterns <- c(
    "^Fig1C_Pathway_PCoA_NatCom_panel\\.(pdf|png|rds)$",
    "^Figure3_Functional_Response_NatCom\\.(pdf|png|svg)$",
    "^Figure4_Functional_Pathway_Consistency_NatCom\\.(pdf|png|svg)$",
    "^Figure5_Stratified_Pathway_Contributors_NatCom\\.(pdf|png|svg)$"
  )
  supp_patterns <- c(
    "^FigureS[0-9].*\\.(pdf|png|svg)$",
    "^FigS[0-9].*\\.(pdf|png|svg)$",
    "^FigureS5_WITH_Unclassified.*\\.(pdf|png|svg)$"
  )
  copy_matching(STAGING_DIR, main_patterns, MAIN_FIG_DIR)
  copy_matching(STAGING_DIR, supp_patterns, SUPP_FUNC_DIR)
  copy_matching(file.path(repo_root, "function", "results"), supp_patterns, SUPP_FUNC_DIR)
}

# ============================================================
# Figure 1C pathway PCoA
# ============================================================
run_section("Figure 1C pathway PCoA", {
local({
# ============================================================
# Nature Communications-style global pathway PCoA
# Figure 1C: Environmental filtering across sites
# Input: HUMAnN4 unstratified CPM data for root and rhizosphere
# Output: compact PDF/PNG suitable for merging with map + taxonomy PCoA
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ape)
  library(ggplot2)
})

# ----------------------------
# Paths
# ----------------------------
proj_dir <- repo_root
rhizo_dir <- file.path(repo_root, "function", "data", "rhizosphere")
root_dir  <- file.path(repo_root, "function", "data", "root")

rhizo_cpm_fp  <- file.path(rhizo_dir, "Rhizosphere_pathabundance_CPM_unstratified_Rprefix.tsv.gz")
rhizo_meta_fp <- file.path(rhizo_dir, "Rhizo_metadata.tsv")

root_cpm_fp   <- file.path(root_dir, "Root_pathabundance_UNSTRATIFIED_CPM.tsv.gz")
root_meta_fp  <- file.path(root_dir, "Root_metadata.tsv.gz")

out_dir <- STAGING_DIR
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Colours and theme
# ----------------------------
loc_cols <- c(
  "Eyrewell_Forest" = "#D55E00",
  "Kowhai"          = "#E69F00",
  "LU_H8"           = "#009E73",
  "Rolleston"       = "#0072B2",
  "West_Coast"      = "#CC79A7"
)
fig1_environment_cols <- c(
  "Eyrewell_Forest"    = "#D55E00",
  "Kowhai (Irrigated)" = "#8C510A",
  "Kowhai (Rainfed)"   = "#E69F00",
  "LU_H8"              = "#009E73",
  "Rolleston"          = "#0072B2",
  "West_Coast"         = "#CC79A7"
)

theme_natcom_pcoa <- function(base_size = 7, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(face = "bold", size = FIG1_PCOA_AXIS_TITLE_SIZE),
      axis.text = element_text(color = "black", size = FIG1_PCOA_AXIS_TEXT_SIZE),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(0.32, "cm"),
      legend.spacing.y = unit(0.02, "cm"),
      panel.border = element_rect(fill = NA, linewidth = 0.3, color = "black"),
      legend.position = "none",
      plot.margin = margin(4, 4, 4, 4)
    )
}

# ----------------------------
# Readers
# ----------------------------
read_humann_table <- function(fp) {
  tbl <- readr::read_tsv(fp, show_col_types = FALSE, comment = "")
  names(tbl)[1] <- "Pathway"
  mat <- as.data.frame(tbl)
  rownames(mat) <- mat$Pathway
  mat$Pathway <- NULL
  mat[] <- lapply(mat, as.numeric)
  mat[is.na(mat)] <- 0
  mat
}

read_meta_rhizo <- function(fp) {
  readr::read_tsv(fp, show_col_types = FALSE) %>%
    mutate(
      SampleID    = trimws(as.character(SampleID)),
      Location    = trimws(as.character(Location)),
      Inoculation = trimws(as.character(Inoculation)),
      Water       = if ("Water" %in% names(.)) trimws(as.character(Water)) else NA_character_,
      Compartment = "Rhizosphere"
    ) %>%
    mutate(
      Location = case_when(
        Location %in% c("Eyerell_Forest", "Eyerell Forest", "Eyrewell Forest") ~ "Eyrewell_Forest",
        Location %in% c("West Coast", "West_Coast") ~ "West_Coast",
        TRUE ~ Location
      ),
      Inoculation = case_when(
        str_to_lower(Inoculation) == "control" ~ "Control",
        str_to_lower(Inoculation) %in% c("panch", "trichoderma")   ~ "Trichoderma",
        TRUE ~ Inoculation
      )
    ) %>%
    select(SampleID, Location, Inoculation, Water, Compartment)
}

read_meta_root <- function(fp) {
  readr::read_tsv(fp, show_col_types = FALSE) %>%
    mutate(
      SampleID    = trimws(as.character(SampleID)),
      Location    = trimws(as.character(Location)),
      Inoculation = trimws(as.character(Inoculation)),
      Water       = if ("Water" %in% names(.)) trimws(as.character(Water)) else NA_character_,
      Compartment = "Root"
    ) %>%
    mutate(
      Location = case_when(
        Location %in% c("Eyerell_Forest", "Eyerell Forest", "Eyrewell Forest") ~ "Eyrewell_Forest",
        Location %in% c("West Coast", "West_Coast") ~ "West_Coast",
        TRUE ~ Location
      ),
      Inoculation = case_when(
        str_to_lower(Inoculation) == "control" ~ "Control",
        str_to_lower(Inoculation) %in% c("panch", "trichoderma")   ~ "Trichoderma",
        TRUE ~ Inoculation
      )
    ) %>%
    select(SampleID, Location, Inoculation, Water, Compartment)
}

expand_mat <- function(mat, features) {
  out <- matrix(0, nrow = length(features), ncol = ncol(mat), dimnames = list(features, colnames(mat)))
  out[rownames(mat), colnames(mat)] <- as.matrix(mat)
  out
}

# ----------------------------
# Load data
# ----------------------------
rhizo_meta <- read_meta_rhizo(rhizo_meta_fp)
root_meta  <- read_meta_root(root_meta_fp)

rhizo_cpm <- read_humann_table(rhizo_cpm_fp)
root_cpm  <- read_humann_table(root_cpm_fp)

# ----------------------------
# Align pathways and samples
# ----------------------------
all_paths <- union(rownames(rhizo_cpm), rownames(root_cpm))
rhizo_cpm2 <- expand_mat(rhizo_cpm, all_paths)
root_cpm2  <- expand_mat(root_cpm,  all_paths)

meta_combined <- bind_rows(rhizo_meta, root_meta) %>% distinct(SampleID, .keep_all = TRUE)
X <- cbind(rhizo_cpm2, root_cpm2)
common_samples <- intersect(meta_combined$SampleID, colnames(X))

meta_combined <- meta_combined %>%
  filter(SampleID %in% common_samples) %>%
  mutate(
    Location = factor(Location, levels = names(loc_cols)),
    Inoculation = factor(Inoculation, levels = c("Control", "Trichoderma")),
    Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")),
    Environment = case_when(
      as.character(Location) == "Kowhai" &
        stringr::str_detect(Water, stringr::regex("^irr", ignore_case = TRUE)) ~
        "Kowhai (Irrigated)",
      as.character(Location) == "Kowhai" &
        stringr::str_detect(Water, stringr::regex("^rain", ignore_case = TRUE)) ~
        "Kowhai (Rainfed)",
      as.character(Location) == "Kowhai" ~ NA_character_,
      TRUE ~ as.character(Location)
    )
  ) %>%
  arrange(Location, Compartment, Inoculation, SampleID)

X <- X[, meta_combined$SampleID, drop = FALSE]
stopifnot(identical(colnames(X), meta_combined$SampleID))

# The model pools Kowhai irrigated and rainfed samples into one Location.
# Environment records six field conditions in the data, not a sixth site.
model_scope <- tibble::tibble(
  Model = "Figure 1C pathway global PERMANOVA",
  Location_levels_in_model = n_distinct(meta_combined$Location, na.rm = TRUE),
  Field_environments_in_dataset = n_distinct(meta_combined$Environment),
  Definition = paste(
    "Location has five sites; Kowhai irrigated and rainfed are pooled.",
    "The ordination colours show six environments; the PERMANOVA model uses five locations."
  )
)
readr::write_tsv(model_scope, file.path(out_dir, "Figure1C_model_scope.tsv"))
if (model_scope$Location_levels_in_model != 5 ||
    model_scope$Field_environments_in_dataset != 6) {
  warning("Expected five locations and six field environments; check Water labels.")
}
if (anyNA(meta_combined$Environment) ||
    !setequal(unique(meta_combined$Environment), names(fig1_environment_cols))) {
  stop("Figure 1C needs all six field environments with valid Kowhai Water labels.")
}
meta_combined$Environment <- factor(
  meta_combined$Environment, levels = names(fig1_environment_cols)
)

X_path <- t(X)
storage.mode(X_path) <- "numeric"
X_path[is.na(X_path)] <- 0
X_path <- X_path[, colSums(X_path) > 0, drop = FALSE]

# ----------------------------
# Bray-Curtis and PERMANOVA
# No Hellinger transform for HUMAnN CPM table
# ----------------------------
set.seed(1)
d_bray <- vegan::vegdist(X_path, method = "bray")

# Five-location model; Kowhai water regimes are not separate Location levels.
perm_global <- vegan::adonis2(
  d_bray ~ Location + Compartment + Inoculation,
  data = meta_combined,
  permutations = 999,
  by = "margin"
)

perm_global_tbl <- as.data.frame(perm_global) %>% tibble::rownames_to_column("Term")
readr::write_tsv(perm_global_tbl,
  file.path(out_dir, "PERMANOVA_pathway_global_additive.tsv"))

# Five-geographic-location x inoculation test adjusted for compartment.
# Kowhai water regimes are pooled here. Root and rhizosphere share plots, but
# plot IDs are unavailable in this metadata; interpret the combined test with
# this pairing limitation in mind. Treatment labels are permuted within site.
if (!identical(attr(d_bray, "Labels"), as.character(meta_combined$SampleID))) {
  stop("Combined pathway metadata do not match Bray-Curtis sample order.")
}
location_permutations <- permute::how(nperm = 999, blocks = meta_combined$Location)
set.seed(1)
perm_location_interaction <- vegan::adonis2(
  d_bray ~ Compartment + Location * Inoculation,
  data = meta_combined, permutations = location_permutations, by = "margin"
)
perm_location_tbl <- as.data.frame(perm_location_interaction) %>%
  tibble::rownames_to_column("Term")
perm_location_row <- perm_location_tbl %>%
  dplyr::filter(Term %in% c("Location:Inoculation", "Inoculation:Location"))
if (nrow(perm_location_row) != 1L) {
  stop("Location x Inoculation term missing from pathway model: ",
       paste(perm_location_tbl$Term, collapse = ", "))
}
readr::write_tsv(perm_location_row,
  file.path(out_dir, "PERMANOVA_pathway_Location_by_Inoculation.tsv"))
# Keep the original additive output separately; the global file identifies
# which rows are drawn from distinct model specifications.
readr::write_tsv(dplyr::bind_rows(
  dplyr::mutate(perm_global_tbl, Model = "Five-location additive", .before = 1),
  dplyr::mutate(perm_location_row,
    Model = "Five-location interaction; restricted permutations within location",
    .before = 1)
), file.path(out_dir, "PERMANOVA_pathway_global.tsv"))
message("Pathway Location x Inoculation interaction saved (five geographic locations).")

# Direct root-versus-rhizosphere comparison on the Figure 1C pathway matrix.
# Keep the six field environments as permutation blocks (Kowhai irrigated and
# rainfed separate) even though the original global Location factor has five.
comp_meta <- meta_combined %>%
  mutate(Environment = factor(Environment))
if (anyNA(comp_meta[, c("Environment", "Compartment", "Inoculation")]) ||
    dplyr::n_distinct(comp_meta$Compartment) != 2L ||
    dplyr::n_distinct(comp_meta$Inoculation) != 2L) {
  stop("Compartment comparison requires both compartments and treatments with complete metadata.")
}
comp_counts <- comp_meta %>%
  dplyr::count(Environment, Compartment, Inoculation, name = "N")
readr::write_tsv(comp_counts, file.path(out_dir, "PERMANOVA_pathway_compartment_group_counts.tsv"))
comp_permutations <- permute::how(
  nperm = 999, blocks = comp_meta$Environment
)
set.seed(1)
comp_main <- vegan::adonis2(
  d_bray ~ Environment + Inoculation + Compartment,
  data = comp_meta, permutations = comp_permutations, by = "margin"
)
set.seed(1)
comp_interaction <- vegan::adonis2(
  d_bray ~ Environment + Inoculation * Compartment,
  data = comp_meta, permutations = comp_permutations, by = "margin"
)
comp_main_tbl <- as.data.frame(comp_main) %>%
  tibble::rownames_to_column("Term")
comp_interaction_tbl <- as.data.frame(comp_interaction) %>%
  tibble::rownames_to_column("Term")
readr::write_tsv(comp_main_tbl, file.path(out_dir, "PERMANOVA_pathway_compartment_additive.tsv"))
readr::write_tsv(comp_interaction_tbl, file.path(out_dir, "PERMANOVA_pathway_compartment_interaction.tsv"))
comp_summary <- dplyr::bind_rows(
  comp_main_tbl %>% dplyr::filter(Term == "Compartment") %>%
    dplyr::mutate(Comparison = "Root versus rhizosphere, adjusted for environment and inoculation"),
  comp_interaction_tbl %>%
    dplyr::filter(grepl("Inoculation", Term) & grepl("Compartment", Term) & grepl(":", Term)) %>%
    dplyr::mutate(Comparison = "Difference in inoculation response between compartments")
)
if (nrow(comp_summary) != 2L) stop("Could not extract both compartment PERMANOVA terms.")
comp_summary <- comp_summary %>%
  dplyr::mutate(
    N_total = nrow(comp_meta),
    N_rhizosphere = sum(comp_meta$Compartment == "Rhizosphere"),
    N_root = sum(comp_meta$Compartment == "Root"),
    Permutations = 999L,
    Permutation_blocks = "Six field environments (Kowhai Irrigated and Rainfed separate)",
    Distance = "Bray-Curtis on untransformed HUMAnN pathway CPM"
  ) %>%
  dplyr::select(Comparison, dplyr::everything())
readr::write_tsv(comp_summary, file.path(out_dir, "PERMANOVA_pathway_compartment_comparison.tsv"))

# Check compartment dispersion before interpreting PERMANOVA separation.
comp_disp <- vegan::betadisper(d_bray, comp_meta$Compartment)
set.seed(1)
comp_disp_test <- vegan::permutest(comp_disp, permutations = comp_permutations)
comp_disp_tbl <- as.data.frame(comp_disp_test$tab) %>%
  tibble::rownames_to_column("Term")
readr::write_tsv(comp_disp_tbl, file.path(out_dir, "PERMDISP_pathway_compartment.tsv"))

loc_r2 <- perm_global_tbl$R2[perm_global_tbl$Term == "Location"]
loc_p  <- perm_global_tbl$`Pr(>F)`[perm_global_tbl$Term == "Location"]
p_label <- ifelse(is.na(loc_p), "NA", ifelse(loc_p < 0.001, "<0.001", paste0("= ", signif(loc_p, 2))))

# ----------------------------
# PCoA
# ----------------------------
pcoa <- ape::pcoa(as.matrix(d_bray))

scores <- as.data.frame(pcoa$vectors[, 1:2]) %>%
  rownames_to_column("SampleID") %>%
  left_join(meta_combined, by = "SampleID") %>%
  rename(PCoA1 = Axis.1, PCoA2 = Axis.2)

readr::write_tsv(scores, file.path(out_dir, "PCoA_scores_pathway_global.tsv"))

xlab <- paste0("PCoA1 (", round(pcoa$values$Relative_eig[1] * 100, 1), "%)")
ylab <- paste0("PCoA2 (", round(pcoa$values$Relative_eig[2] * 100, 1), "%)")

# ----------------------------
# NATURE COMMUNICATIONS STYLE GLOBAL PLOT
# ----------------------------
p_natcom <- ggplot(scores, aes(PCoA1, PCoA2)) +
  geom_point(aes(color = Environment, shape = Compartment), size = 1.9, alpha = 0.92, stroke = 0.25) +
  scale_color_manual(values = fig1_environment_cols, drop = FALSE) +
  scale_shape_manual(values = c(Rhizosphere = 16, Root = 17), drop = FALSE) +
  labs(x = xlab, y = ylab, color = "Field environment", shape = "Compartment") +
  guides(color = guide_legend(override.aes = list(size = 2.1), order = 1),
         shape = guide_legend(override.aes = list(size = 2.1), order = 2)) +
  theme_natcom_pcoa(base_size = FIG1_PCOA_BASE_SIZE)

# Compact, legend-free panel for merging into Fig. 1. Its dimensions exactly
# match Figure 1B so the two ordination plotting regions remain the same size.
saveRDS(
  p_natcom,
  file.path(out_dir, "Fig1C_Pathway_PCoA_NatCom_panel.rds")
)
ggsave(file.path(out_dir, "Fig1C_Pathway_PCoA_NatCom_panel.pdf"), p_natcom,
       width = FIG1_PCOA_EXPORT_WIDTH_MM, height = FIG1_PCOA_EXPORT_HEIGHT_MM,
       units = "mm", useDingbats = FALSE)
ggsave(file.path(out_dir, "Fig1C_Pathway_PCoA_NatCom_panel.png"), p_natcom,
       width = FIG1_PCOA_EXPORT_WIDTH_MM, height = FIG1_PCOA_EXPORT_HEIGHT_MM,
       units = "mm", dpi = 600)

# Wider inspection version
ggsave(file.path(out_dir, "Fig1C_Pathway_PCoA_NatCom_wide.pdf"), p_natcom,
       width = 110, height = 80, units = "mm", useDingbats = FALSE)
ggsave(file.path(out_dir, "Fig1C_Pathway_PCoA_NatCom_wide.png"), p_natcom,
       width = 110, height = 80, units = "mm", dpi = 600)

message("Done. Saved Nature Communications-style pathway PCoA to: ", out_dir)
message("Editable Figure 1C RDS: ", file.path(out_dir, "Fig1C_Pathway_PCoA_NatCom_panel.rds"))

})
})


# ============================================================
# Rhizosphere unstratified functional analysis (Figures 3A, 4A, S3)
# ============================================================
run_section("Rhizosphere unstratified functional analysis (Figures 3A, 4A, S3)", {
local({
# ============================================================
# RHIZOSPHERE shotgun metagenomics (48) — HUMAnN4 UNSTRATIFIED
# Focus: Inoculation (Trichoderma vs Control), controlling for Location
# Outputs:
#   - PERMANOVA, PCoA, NMDS
#   - ALDEx2 (global + within-Location)
#   - Volcanoes (Nature-style)
#   - Consistency table + headline pathways
#   - Heatmap A: top variable pathways
#   - Heatmap B: inoculation-associated pathways
#   - Section 11: Top pathways responding to Trichoderma inoculation
# UPDATED:
#   - standardised Eyrewell_Forest spelling
#   - Section 11 safely integrated without overwriting out_dir
# ============================================================

# ----------------------------
# 0) Packages (reproducible)
# ----------------------------
pkgs_cran <- c("tidyverse","vegan","ape","pheatmap","matrixStats","ggrepel")
to_install <- pkgs_cran[!sapply(pkgs_cran, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs_cran, library, character.only = TRUE))

if (!requireNamespace("ALDEx2", quietly = TRUE)) {
  message("Installing ALDEx2 from Bioconductor ...")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ALDEx2")
}
library(ALDEx2)

# ----------------------------
# 1) Paths
# ----------------------------
library(here)
project_dir <- repo_root
data_dir <- file.path(project_dir, "function", "data", "rhizosphere")
out_dir  <- file.path(project_dir, "function", "results", "rhizosphere", "unstratified")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cpm_fp  <- file.path(data_dir, "Rhizosphere_pathabundance_CPM_unstratified_Rprefix.tsv.gz")
raw_fp  <- file.path(data_dir, "Rhizosphere_pathabundance_RAW_unstratified_Rprefix.tsv.gz")
meta_fp <- file.path(data_dir, "Rhizo_metadata.tsv")

fig_dir <- file.path(out_dir, "Figures")
tab_dir <- file.path(out_dir, "Tables")
log_dir <- file.path(out_dir, "Logs")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 2) Nature Microbiology-ish helpers
# ----------------------------
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

dir_cols <- c(
  "Higher in Trichoderma"   = unname(okabe_ito["blue"]),
  "Higher in Control" = unname(okabe_ito["verm"])
)

theme_nature <- function(base_size = 8, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title   = element_text(face = "bold", size = base_size + 1),
      axis.title   = element_text(face = "bold", size = base_size + 1),
      axis.text    = element_text(color = "black", size = base_size),
      axis.line    = element_line(linewidth = 0.35),
      axis.ticks   = element_line(linewidth = 0.35),
      legend.title = element_text(face = "bold", size = base_size + 0.5),
      legend.text  = element_text(size = base_size),
      legend.key   = element_blank(),
      legend.position = "right",
      panel.border = element_rect(fill = NA, linewidth = 0.35, color = "black"),
      plot.margin  = margin(6, 8, 6, 6)
    )
}

theme_nature_matrix <- function(base_size = 9, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title   = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      axis.title   = element_text(face = "bold"),
      axis.text.x  = element_text(color = "black", size = base_size, face = "bold"),
      axis.text.y  = element_text(color = "black", size = base_size - 1),
      axis.line    = element_line(linewidth = 0.35),
      axis.ticks   = element_line(linewidth = 0.35),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text  = element_text(size = base_size - 1),
      legend.key   = element_blank(),
      legend.position = "right",
      plot.margin  = margin(10, 40, 10, 10)
    )
}

save_fig <- function(p, name, w = 6.2, h = 4.6) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")),
         p, width = w, height = h, useDingbats = FALSE)
  ggsave(file.path(fig_dir, paste0(name, ".png")),
         p, width = w, height = h, dpi = 600)
}

write_tsv2 <- function(df, filename) {
  readr::write_tsv(df, file.path(tab_dir, filename))
}

read_humann_table <- function(fp) {
  tbl <- readr::read_tsv(fp, show_col_types = FALSE, comment = "")
  names(tbl)[1] <- "Pathway"
  mat <- as.data.frame(tbl)
  rownames(mat) <- mat$Pathway
  mat$Pathway <- NULL
  mat[] <- lapply(mat, function(x) as.numeric(x))
  mat[is.na(mat)] <- 0
  mat
}

save_pheatmap <- function(mat, anno, name, main,
                          w_pdf = 14, h_pdf = 7,
                          w_png = 6000, h_png = 3200, res = 600,
                          cluster_cols = TRUE, cluster_rows = TRUE,
                          show_colnames = FALSE, show_rownames = TRUE,
                          fontsize = 9, fontsize_row = 7, fontsize_col = 7,
                          gaps_col = NULL, gaps_row = NULL,
                          annotation_names_col = TRUE,
                          annotation_fontsize = 8) {
  
  pdf(file.path(fig_dir, paste0(name, ".pdf")),
      width = w_pdf, height = h_pdf, useDingbats = FALSE)
  pheatmap::pheatmap(
    mat,
    annotation_col = anno,
    cluster_cols = cluster_cols,
    cluster_rows = cluster_rows,
    show_colnames = show_colnames,
    show_rownames = show_rownames,
    border_color = NA,
    gaps_col = gaps_col,
    gaps_row = gaps_row,
    fontsize = fontsize,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    annotation_names_col = annotation_names_col,
    annotation_names_row = TRUE,
    annotation_fontsize = annotation_fontsize,
    main = main
  )
  dev.off()
  
  png(file.path(fig_dir, paste0(name, ".png")),
      width = w_png, height = h_png, res = res)
  pheatmap::pheatmap(
    mat,
    annotation_col = anno,
    cluster_cols = cluster_cols,
    cluster_rows = cluster_rows,
    show_colnames = show_colnames,
    show_rownames = show_rownames,
    border_color = NA,
    gaps_col = gaps_col,
    gaps_row = gaps_row,
    fontsize = fontsize,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    annotation_names_col = annotation_names_col,
    annotation_names_row = TRUE,
    annotation_fontsize = annotation_fontsize,
    main = main
  )
  dev.off()
}

plot_volcano_splitlabels <- function(res_tbl, title, fname_base,
                                     p_cut = 0.10, effect_cut = 0.25,
                                     label_n_right = 10, label_n_left = 6) {
  
  df <- res_tbl %>%
    mutate(
      neglog10p = -log10(we.ep + 1e-300),
      Sig = (we.ep < p_cut) & (abs(effect) > effect_cut),
      Direction = case_when(
        Sig & effect > 0 ~ "Higher in Control",
        Sig & effect < 0 ~ "Higher in Trichoderma",
        TRUE ~ "Not highlighted"
      )
    )
  
  p <- ggplot(df, aes(effect, neglog10p)) +
    geom_point(aes(colour = Direction), size = 1.7, alpha = 0.75) +
    scale_colour_manual(
      values = c(
        "Higher in Trichoderma" = "#0072B2",
        "Higher in Control" = "#D55E00",
        "Not highlighted" = "grey70"
      ),
      labels = trichoderma_direction_labels
    ) +
    geom_vline(xintercept = c(-effect_cut, effect_cut),
               linetype = 2, linewidth = 0.35) +
    geom_hline(yintercept = -log10(p_cut),
               linetype = 2, linewidth = 0.35) +
    theme_nature(base_size = 8) +
    labs(
      title = trichoderma_title(title),
      x = "Effect size (ALDEx2)",
      y = expression(-log[10](p)),
      colour = NULL
    )
  
  save_fig(p, fname_base, w = 7.2, h = 4.8)
}

# ----------------------------
# 3) Read metadata
# ----------------------------
meta <- readr::read_tsv(meta_fp, show_col_types = FALSE) %>%
  mutate(
    SampleID    = trimws(as.character(SampleID)),
    Location    = trimws(as.character(Location)),
    Inoculation = trimws(as.character(Inoculation)),
    Water       = trimws(as.character(Water))
  ) %>%
  mutate(
    Location = case_when(
      Location %in% c("Eyerell_Forest", "Eyerell Forest", "Eyrewell Forest") ~ "Eyrewell_Forest",
      TRUE ~ Location
    ),
    Inoculation = case_when(
      str_to_lower(Inoculation) %in% c("control") ~ "Control",
      str_to_lower(Inoculation) %in% c("panch", "trichoderma")   ~ "Trichoderma",
      TRUE ~ Inoculation
    )
  )

stopifnot(all(c("SampleID","Location","Inoculation","Water") %in% names(meta)))
if (anyDuplicated(meta$SampleID)) stop("Duplicate SampleID in metadata.")

# ----------------------------
# 4) Read HUMAnN tables
# ----------------------------
cpm_mat <- read_humann_table(cpm_fp)
raw_mat <- read_humann_table(raw_fp)

# ----------------------------
# 5) Align samples
# ----------------------------
common_samples <- Reduce(intersect, list(meta$SampleID, colnames(cpm_mat), colnames(raw_mat)))
if (length(common_samples) < 6) stop("Too few shared samples. Check SampleID naming/columns.")

meta2 <- meta %>%
  filter(SampleID %in% common_samples) %>%
  arrange(SampleID) %>%
  mutate(
    Inoculation = factor(Inoculation, levels = c("Control","Trichoderma")),
    Water = stringr::str_to_title(trimws(as.character(Water))),
    Environment = factor(
      ifelse(Location == "Kowhai", paste0("Kowhai_", Water), Location),
      levels = FUNCTION_ENV_LEVELS
    )
  )

if (anyNA(meta2$Environment)) {
  stop("Cannot assign all samples to six environments; check Kowhai Water metadata.")
}

bad_loc <- setdiff(sort(unique(meta2$Location)), names(loc_cols))
if (length(bad_loc) > 0) {
  stop("These Location values are not in loc_cols: ", paste(bad_loc, collapse = ", "))
}
meta2 <- meta2 %>%
  mutate(Location = factor(Location, levels = names(loc_cols)))

cpm_mat2 <- cpm_mat[, meta2$SampleID, drop = FALSE]
raw_mat2 <- raw_mat[, meta2$SampleID, drop = FALSE]
stopifnot(identical(colnames(cpm_mat2), meta2$SampleID))
stopifnot(identical(colnames(raw_mat2), meta2$SampleID))

X_cpm <- t(as.matrix(cpm_mat2))
storage.mode(X_cpm) <- "numeric"
X_cpm[is.na(X_cpm)] <- 0
X_cpm <- X_cpm[, colSums(X_cpm) > 0, drop = FALSE]

X_raw <- t(as.matrix(raw_mat2))
storage.mode(X_raw) <- "numeric"
X_raw[is.na(X_raw)] <- 0
X_raw <- X_raw[, colSums(X_raw) > 0, drop = FALSE]

message("Samples used: ", nrow(X_cpm))
message("CPM pathways (nonzero): ", ncol(X_cpm))
message("RAW pathways (nonzero): ", ncol(X_raw))

write_tsv2(meta2, "Metadata_aligned.tsv")

# ----------------------------
# 6) Bray–Curtis + PERMANOVA + dispersion
# ----------------------------
meta2_perm <- meta2 %>% filter(!is.na(Location), !is.na(Inoculation), !is.na(Environment))
env_counts <- meta2_perm %>% dplyr::count(Environment, Inoculation, name = "n")
write_tsv2(env_counts, "DIAGNOSTIC_Environment_by_Inoculation_counts.tsv")
if (nrow(env_counts) != 2L * length(FUNCTION_ENV_LEVELS) || any(env_counts$n < 2L)) {
  stop("Six environments with at least two samples per treatment are required for site-specific ALDEx2 volcanoes.")
}
X_cpm_perm <- X_cpm[meta2_perm$SampleID, , drop = FALSE]

set.seed(1)
d_bray <- vegan::vegdist(X_cpm_perm, method = "bray")

# Five-location global model: Kowhai irrigated and rainfed share one Location.
# Field_environment is descriptive and is not a sixth level in PERMANOVA.
field_environment <- ifelse(
  as.character(meta2_perm$Location) == "Kowhai",
  paste("Kowhai", as.character(meta2_perm$Water)),
  as.character(meta2_perm$Location)
)
model_scope <- tibble::tibble(
  Model = "Pathway PERMANOVA: Location + Inoculation",
  Location_levels_in_model = dplyr::n_distinct(meta2_perm$Location),
  Field_environments_in_dataset = dplyr::n_distinct(field_environment),
  Definition = paste(
    "Five locations in PERMANOVA; Kowhai water regimes pooled.",
    "Within-environment PERMANOVA and ALDEx2 analyse Kowhai Irrigated and Rainfed separately."
  )
)
write_tsv2(model_scope, "PERMANOVA_model_scope.tsv")
if (model_scope$Location_levels_in_model != 5 ||
    model_scope$Field_environments_in_dataset != 6) {
  warning("Expected five locations and six field environments; check Water metadata.")
}

# ----------------------------
# 6A) Global five-location PERMANOVA
# ----------------------------
perm_main <- vegan::adonis2(
  d_bray ~ Location + Inoculation,
  data = meta2_perm,
  permutations = 999,
  by = "margin"
)

perm_main_tbl <- as.data.frame(perm_main) %>%
  tibble::rownames_to_column("Term")

write_tsv2(
  perm_main_tbl,
  "PERMANOVA_main_Location_Inoculation.tsv"
)

# ----------------------------
# 6A2) Does the inoculation response vary among six field environments?
# ----------------------------
# This interaction compares treatment-associated shifts between environments
# within this compartment. It is distinct from six separate treatment tests.
# Treatment labels are permuted within environment. The metadata currently
# lack a plot/block ID, so the test does not represent paired plot blocking.
if (!identical(attr(d_bray, "Labels"), as.character(meta2_perm$SampleID))) {
  stop("PERMANOVA metadata order does not match pathway distance matrix.")
}
env_permutations <- permute::how(nperm = 999, blocks = meta2_perm$Environment)
set.seed(1)
perm_env_interaction <- vegan::adonis2(
  d_bray ~ Environment * Inoculation,
  data = meta2_perm, permutations = env_permutations, by = "margin"
)
perm_env_interaction_df <- as.data.frame(perm_env_interaction) %>%
  tibble::rownames_to_column("Term")
perm_env_interaction_row <- perm_env_interaction_df %>%
  dplyr::filter(Term %in% c("Environment:Inoculation", "Inoculation:Environment"))
if (nrow(perm_env_interaction_row) != 1L) {
  stop("Could not extract pathway Environment x Inoculation interaction: ",
       paste(perm_env_interaction_df$Term, collapse = ", "))
}
perm_env_interaction_row <- perm_env_interaction_row %>%
  dplyr::mutate(
    N = nrow(meta2_perm), Field_environments = length(FUNCTION_ENV_LEVELS),
    Permutations = 999L,
    Permutation_blocks = "Environment (six fields; Kowhai split)",
    .before = 1
  )
write_tsv2(perm_env_interaction_row,
           "PERMANOVA_Environment_by_Inoculation.tsv")
message("Functional Environment x Inoculation PERMANOVA saved in ", tab_dir)

# ----------------------------
# 6B) Inoculation effect with permutations constrained within Environment
# ----------------------------
perm_inoc_strata <- vegan::adonis2(
  d_bray ~ Inoculation,
  data = meta2_perm,
  permutations = 999,
  strata = meta2_perm$Environment
)

perm_inoc_strata_tbl <- as.data.frame(perm_inoc_strata) %>%
  tibble::rownames_to_column("Term")

write_tsv2(
  perm_inoc_strata_tbl,
  "PERMANOVA_Inoculation_stratifiedByEnvironment.tsv"
)
unlink(file.path(tab_dir, "PERMANOVA_Inoculation_stratifiedByLocation.tsv"))

# 6C) Within-environment PERMANOVA: six separate field environments
# ----------------------------
# The five-location global model above is retained. Field-specific treatment
# tests below split Kowhai into Irrigated and Rainfed (four samples per group).
perm_by_environment <- purrr::map_dfr(FUNCTION_ENV_LEVELS, function(env) {
  meta_env <- meta2_perm %>%
    dplyr::filter(as.character(Environment) == env) %>%
    droplevels()
  group_counts <- table(meta_env$Inoculation)
  if (nrow(meta_env) < 4L || length(group_counts) != 2L ||
      any(group_counts < 2L)) {
    stop("Expected both inoculation groups in field environment ", env)
  }
  X_env <- X_cpm_perm[meta_env$SampleID, , drop = FALSE]
  d_env <- vegan::vegdist(X_env, method = "bray")
  fit <- vegan::adonis2(
    d_env ~ Inoculation, data = meta_env, permutations = 999
  )
  fit_df <- as.data.frame(fit)
  if (nrow(fit_df) < 1L ||
      !rownames(fit_df)[1] %in% c("Model", "Inoculation")) {
    stop("Unexpected adonis2 output in ", env, ": ",
         paste(rownames(fit_df), collapse = ", "))
  }
  fit_df[1, , drop = FALSE] %>%
    dplyr::mutate(Environment = env, N = nrow(meta_env),
                  Term = "Inoculation", .before = 1) %>%
    dplyr::select(Environment, N, Term, Df, SumOfSqs, R2, F, `Pr(>F)`)
})
if (nrow(perm_by_environment) != 6L ||
    !setequal(perm_by_environment$Environment, FUNCTION_ENV_LEVELS)) {
  stop("Expected one inoculation PERMANOVA row per six field environments.")
}
write_tsv2(perm_by_environment, "PERMANOVA_withinEnvironment_Inoculation.tsv")
# Remove pooled five-location outputs produced by older script versions.
unlink(file.path(tab_dir, "PERMANOVA_withinLocation_Inoculation.tsv"))
message("Functional within-environment PERMANOVA: six separate field results saved.")

# ----------------------------
# 6D) Combined PERMANOVA summary table (global + six environments)
# ----------------------------
perm_main_combined <- perm_main_tbl %>%
  filter(Term %in% c("Location", "Inoculation")) %>%
  mutate(
    Model = "Global",
    Location_label = "All sites",
    N = nrow(meta2_perm)
  ) %>%
  dplyr::select(Model, Location_label, N, Term, Df, SumOfSqs, R2, F, `Pr(>F)`)

perm_by_environment_combined <- perm_by_environment %>%
  dplyr::mutate(
    Model = "Within-environment",
    Location_label = pretty_environment_label(Environment)
  ) %>%
  dplyr::select(Model, Location_label, N, Term, Df, SumOfSqs, R2, F, `Pr(>F)`)

perm_combined <- dplyr::bind_rows(
  perm_main_combined,
  perm_env_interaction_row %>%
    dplyr::mutate(Model = "Environment interaction", Location_label = "Six environments") %>%
    dplyr::select(Model, Location_label, N, Term, Df, SumOfSqs, R2, F, `Pr(>F)`),
  perm_by_environment_combined
)

write_tsv2(
  perm_combined,
  "TABLE_PERMANOVA_Function_Global_and_WithinEnvironment.tsv"
)
unlink(file.path(tab_dir, "TABLE_PERMANOVA_Function_Global_and_WithinSite.tsv"))

# ----------------------------
# 6E) Dispersion test
# ----------------------------
bd <- vegan::betadisper(d_bray, meta2_perm$Inoculation)
capture.output(anova(bd), file = file.path(log_dir, "Betadisper_Inoculation_ANOVA.txt"))
capture.output(
  vegan::permutest(bd, permutations = 999),
  file = file.path(log_dir, "Betadisper_Inoculation_permutest.txt")
)

# ----------------------------
# 7) Ordinations
# ----------------------------
pcoa <- ape::pcoa(as.matrix(d_bray))
scores_pcoa <- as.data.frame(pcoa$vectors[, 1:2]) %>%
  rownames_to_column("SampleID") %>%
  left_join(meta2_perm, by = "SampleID") %>%
  rename(PCoA1 = Axis.1, PCoA2 = Axis.2)
write_tsv2(scores_pcoa, "PCoA_scores.tsv")

p_pcoa <- ggplot(scores_pcoa, aes(PCoA1, PCoA2)) +
  geom_point(aes(color = Environment, shape = Inoculation),
             size = 2.2, alpha = 0.9, stroke = 0.35) +
  scale_color_manual(
    values = functional_environment_cols,
    labels = pretty_environment_label,
    drop = FALSE
  ) +
  scale_shape_manual(values = c(Control = 16, Trichoderma = 17),
                     labels = trichoderma_inoculation_labels, drop = FALSE) +
  theme_nature(base_size = 8) +
  labs(
    title = "Rhizosphere PCoA (Bray–Curtis)",
    x = paste0("PCoA1 (", round(pcoa$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa$values$Relative_eig[2] * 100, 1), "%)")
  )
save_fig(p_pcoa, "Fig_PCoA_Location_Inoculation_Nature", w = 6.5, h = 4.8)

# ============================================================
# NatCom Figure 3A — Rhizosphere functional pathway PCoA
# Consistent style with taxonomy Figure 2
# ============================================================
NATCOM_DIR <- STAGING_DIR
dir.create(NATCOM_DIR, recursive = TRUE, showWarnings = FALSE)

p_fig3_rhiz_func <- p_pcoa +
  labs(
    title = "Rhizosphere",
    colour = "Field environment",
    shape = "Inoculation"
  ) +
  theme_nature(base_size = FIG3_PCOA_BASE_SIZE) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0, size = FIG3_PCOA_TITLE_SIZE),
    axis.title = element_text(face = "bold", size = FIG3_PCOA_AXIS_TITLE_SIZE),
    axis.text = element_text(colour = "black", size = FIG3_PCOA_AXIS_TEXT_SIZE),
    aspect.ratio = FIG3_PCOA_ASPECT_RATIO,
    plot.margin = margin(4, 4, 4, 4)
  ) +
  guides(
    colour = guide_legend(override.aes = list(size = 3.2), order = 1),
    shape  = guide_legend(override.aes = list(size = 3.2), order = 2)
  )

saveRDS(
  p_fig3_rhiz_func,
  file.path(NATCOM_DIR, "p_fig3_rhiz_func.rds")
)

ggsave(
  file.path(NATCOM_DIR, "Fig3A_Rhizosphere_Function_PCoA_NatCom_panel.pdf"),
  p_fig3_rhiz_func,
  width = 95,
  height = FIG3_HEIGHT_MM,
  units = "mm",
  useDingbats = FALSE
)

ggsave(
  file.path(NATCOM_DIR, "Fig3A_Rhizosphere_Function_PCoA_NatCom_panel.png"),
  p_fig3_rhiz_func,
  width = 95,
  height = FIG3_HEIGHT_MM,
  units = "mm",
  dpi = 600
)

set.seed(1)
nmds <- vegan::metaMDS(X_cpm_perm, distance = "bray", k = 2, trymax = 200,
                       autotransform = FALSE, trace = FALSE)
scores_nmds <- as.data.frame(vegan::scores(nmds, display = "sites")) %>%
  rownames_to_column("SampleID") %>%
  left_join(meta2_perm, by = "SampleID")
write_tsv2(scores_nmds, "NMDS_scores.tsv")

p_nmds <- ggplot(scores_nmds, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = Location, shape = Inoculation),
             size = 2.2, alpha = 0.9, stroke = 0.35) +
  scale_color_manual(values = loc_cols, drop = FALSE) +
  scale_shape_manual(values = c(Control = 16, Trichoderma = 17),
                     labels = trichoderma_inoculation_labels, drop = FALSE) +
  theme_nature(base_size = 8) +
  labs(title = paste0("Rhizosphere NMDS (stress=", round(nmds$stress, 3), ")"))
save_fig(p_nmds, "Fig_NMDS_Location_Inoculation_Nature", w = 6.5, h = 4.8)

# ----------------------------
# 8) Heatmap A
# ----------------------------
N_var <- 30
vars <- matrixStats::colVars(X_cpm_perm)
top_var <- names(sort(vars, decreasing = TRUE))[1:min(N_var, length(vars))]

Z <- scale(X_cpm_perm[, top_var, drop = FALSE])
hmA <- t(Z)

annoA <- meta2_perm %>%
  dplyr::select(SampleID, Location, Inoculation) %>%
  as.data.frame()
rownames(annoA) <- annoA$SampleID
annoA$SampleID <- NULL
annoA <- annoA[colnames(hmA), , drop = FALSE]

save_pheatmap(
  hmA, annoA,
  name = "Fig_Heatmap_TopVariablePathways_CPM_GLOBAL",
  main = paste0("Top ", N_var, " variable pathways (z-scored CPM)"),
  w_pdf = 16, h_pdf = 7.5,
  cluster_cols = TRUE, cluster_rows = TRUE,
  show_colnames = FALSE, show_rownames = TRUE,
  fontsize = 9, fontsize_row = 7, fontsize_col = 7,
  annotation_fontsize = 8,
  annotation_names_col = TRUE
)

col_orderA <- meta2_perm %>%
  arrange(Inoculation, Location, SampleID) %>%
  pull(SampleID)

hmA_blk <- hmA[, col_orderA, drop = FALSE]

annoA_blk <- meta2_perm %>%
  dplyr::select(SampleID, Location, Inoculation) %>%
  as.data.frame()
rownames(annoA_blk) <- annoA_blk$SampleID
annoA_blk$SampleID <- NULL
annoA_blk <- annoA_blk[col_orderA, , drop = FALSE]

n_control <- sum(meta2_perm$Inoculation == "Control", na.rm = TRUE)

save_pheatmap(
  hmA_blk, annoA_blk,
  name = "Fig_Heatmap_TopVariablePathways_CPM_BLOCKED_ControlVsTrichoderma",
  main = paste0("Top ", N_var, " variable pathways (z-scored CPM) — blocked by inoculation"),
  w_pdf = 16, h_pdf = 7.5,
  cluster_cols = FALSE, cluster_rows = TRUE,
  show_colnames = FALSE, show_rownames = TRUE,
  gaps_col = n_control,
  fontsize = 9, fontsize_row = 7, fontsize_col = 7,
  annotation_fontsize = 8,
  annotation_names_col = TRUE
)

# ----------------------------
# 9) ALDEx2 + volcanoes + consistency
# ----------------------------
X_raw_int <- round(X_raw[meta2_perm$SampleID, , drop = FALSE])
X_raw_int[X_raw_int < 0] <- 0
storage.mode(X_raw_int) <- "integer"

run_aldex2_one <- function(X_sub, group, label, mc_samples = 128, denom = "all") {
  group_f <- droplevels(factor(group))
  if (nlevels(group_f) != 2) stop("ALDEx2 requires exactly 2 groups: ", label)
  
  X_sub <- as.matrix(X_sub)
  storage.mode(X_sub) <- "numeric"
  X_sub[is.na(X_sub)] <- 0
  X_sub <- round(X_sub); X_sub[X_sub < 0] <- 0
  
  counts <- t(X_sub)
  storage.mode(counts) <- "integer"
  
  ald <- ALDEx2::aldex.clr(
    reads = counts,
    conds = as.character(group_f),
    mc.samples = as.integer(mc_samples),
    denom = denom,
    verbose = FALSE
  )
  tt  <- ALDEx2::aldex.ttest(ald, paired.test = FALSE, verbose = FALSE)
  eff <- ALDEx2::aldex.effect(ald, verbose = FALSE)
  
  tibble::tibble(
    Pathway  = rownames(tt),
    we.eBH   = tt$we.eBH,
    we.ep    = tt$we.ep,
    effect   = eff$effect,
    diff.btw = eff$diff.btw,
    overlap  = eff$overlap,
    Contrast = label,
    Group1   = levels(group_f)[1],
    Group2   = levels(group_f)[2]
  ) %>% arrange(we.ep, desc(abs(effect)))
}

aldex_global <- run_aldex2_one(X_raw_int, meta2_perm$Inoculation, "Inoculation_Global")
write_tsv2(aldex_global, "ALDEx2_Inoculation_Global.tsv")

plot_volcano_splitlabels(
  aldex_global,
  "Differential pathways: Trichoderma vs Control (Global; RAW)",
  "Fig_ALDEx2_Volcano_Inoculation_Global_Nature",
  p_cut = 0.10, effect_cut = 0.25,
  label_n_right = 10, label_n_left = 6
)

aldex_by_loc <- purrr::map_dfr(levels(meta2_perm$Location), function(loc) {
  idx <- meta2_perm$Location == loc
  g <- droplevels(meta2_perm$Inoculation[idx])
  if (nlevels(g) != 2) return(tibble())
  
  run_aldex2_one(
    X_sub = X_raw_int[idx, , drop = FALSE],
    group = g,
    label = paste0("Inoculation_", loc),
    mc_samples = 128,
    denom = "all"
  ) %>% mutate(Location = loc)
})
write_tsv2(aldex_by_loc, "ALDEx2_Inoculation_byLocation.tsv")

# Separate field-environment contrasts for Supplementary Figures S3/S4.
# Existing five-location results and cross-location pathway screens stay pooled.
aldex_by_env <- purrr::map_dfr(FUNCTION_ENV_LEVELS, function(env) {
  idx <- as.character(meta2_perm$Environment) == env
  g <- droplevels(meta2_perm$Inoculation[idx])
  if (nlevels(g) != 2L || any(table(g) < 2L)) {
    stop("Insufficient treatment replicates for environment: ", env)
  }
  run_aldex2_one(
    X_sub = X_raw_int[idx, , drop = FALSE],
    group = g,
    label = paste0("Inoculation_", env),
    mc_samples = 128,
    denom = "all"
  ) %>% mutate(Environment = env, Location = as.character(meta2_perm$Location[which(idx)[1]]))
})
write_tsv2(aldex_by_env, "ALDEx2_Inoculation_byEnvironment.tsv")

for (loc in levels(meta2_perm$Location)) {
  sub <- aldex_by_loc %>% filter(Location == loc)
  if (nrow(sub) == 0) next
  
  plot_volcano_splitlabels(
    sub,
    paste0("Differential pathways: Trichoderma vs Control — ", loc, " (RAW)"),
    paste0("Fig_ALDEx2_Volcano_Inoculation_", loc, "_Nature"),
    p_cut = 0.10, effect_cut = 0.25,
    label_n_right = 10, label_n_left = 6
  )
}

# Screening for cross-location consistency is nominal, not BH significant.
MATRIX_EFFECT_MIN <- 0.20  # Rhizosphere site matrix only
support_ep  <- PATHWAY_NOMINAL_P
support_eff <- PATHWAY_EFFECT_MIN
write_tsv2(tibble::tibble(
  Selection = c("Exploratory pathway screen", "BH-adjusted pathway evidence"),
  P_column = c("we.ep", "we.eBH"),
  P_cutoff = c(PATHWAY_NOMINAL_P, PATHWAY_ADJUSTED_Q),
  Minimum_absolute_effect = PATHWAY_EFFECT_MIN,
  Matrix_absolute_effect = MATRIX_EFFECT_MIN,
  Used_for_headline_selection = c(TRUE, FALSE)
), "TABLE_PathwayThresholds.tsv")

consistency <- aldex_by_loc %>%
  mutate(
    supported = (we.ep < support_ep) & (abs(effect) > support_eff),
    adjusted_supported = (we.eBH < PATHWAY_ADJUSTED_Q) &
      (abs(effect) > support_eff),
    dir = case_when(
      effect > 0 ~ "Higher_in_Control",
      effect < 0 ~ "Higher_in_Trichoderma",
      TRUE ~ "Zero"
    )
  ) %>%
  group_by(Pathway) %>%
  summarise(
    n_locations_tested    = n_distinct(Location),
    n_locations_supported = sum(supported, na.rm = TRUE),
    n_locations_adjusted = sum(adjusted_supported, na.rm = TRUE),
    dir_supported         = paste(sort(unique(dir[supported])), collapse = ";"),
    median_effect         = median(effect, na.rm = TRUE),
    min_weep              = min(we.ep, na.rm = TRUE),
    min_weebh             = min(we.eBH, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_locations_supported), min_weep, desc(abs(median_effect)))

write_tsv2(consistency, "TABLE_PathwayConsistency_acrossLocations.tsv")
write_tsv2(
  consistency %>% filter(n_locations_adjusted >= 1),
  "TABLE_PathwayConsistency_BHAdjusted.tsv"
)

headline <- consistency %>%
  filter(n_locations_supported >= 2) %>%
  slice_head(n = 40) %>%
  mutate(Selection_basis = "Exploratory: we.ep < 0.10 and |effect| > 0.25 in at least two locations; not BH adjusted")
write_tsv2(headline, "TABLE_HeadlinePathways_sigIn2plusLocations.tsv")
message("Exploratory nominal-screen headline pathways (>=2 locations; not BH adjusted): ", nrow(headline))

headline1 <- tibble::tibble(Pathway = character())
if (nrow(headline) == 0) {
  headline1 <- consistency %>%
    filter(n_locations_supported >= 1) %>%
    slice_head(n = 40)
  write_tsv2(headline1, "TABLE_HeadlinePathways_sigIn1plusLocations_FALLBACK.tsv")
  message("Fallback headline (>=1 location supported): ", nrow(headline1))
}

# ----------------------------
# 10) Heatmap B
# ----------------------------
hm_paths <- if (nrow(headline) > 0) headline$Pathway else headline1$Pathway
hm_paths <- intersect(unique(hm_paths), colnames(X_cpm_perm))
hm_paths <- hm_paths[1:min(30, length(hm_paths))]
message("Inoculation-associated pathways available for heatmap: ", length(hm_paths))

if (length(hm_paths) >= 10) {
  
  ZB <- scale(X_cpm_perm[, hm_paths, drop = FALSE])
  hmB <- t(ZB)
  
  col_orderB <- meta2_perm %>% arrange(Inoculation, Location, SampleID) %>% pull(SampleID)
  hmB_blk <- hmB[, col_orderB, drop = FALSE]
  
  annoB_blk <- meta2_perm %>%
    dplyr::select(SampleID, Location, Inoculation) %>%
    as.data.frame()
  rownames(annoB_blk) <- annoB_blk$SampleID
  annoB_blk$SampleID <- NULL
  annoB_blk <- annoB_blk[col_orderB, , drop = FALSE]
  
  save_pheatmap(
    hmB_blk, annoB_blk,
    name = "Fig_Heatmap_HeadlineInoculationPathways_CPM_BLOCKED_ControlVsTrichoderma",
    main = "Inoculation-associated pathways (z-scored CPM)",
    w_pdf = 16, h_pdf = 7.0,
    cluster_cols = FALSE, cluster_rows = TRUE,
    show_colnames = FALSE, show_rownames = TRUE,
    gaps_col = n_control,
    fontsize = 9, fontsize_row = 7, fontsize_col = 7,
    annotation_fontsize = 8,
    annotation_names_col = TRUE
  )
  
  write_tsv2(
    tibble(Pathway = hm_paths),
    "TABLE_HeatmapPathways_InoculationAssociated.tsv"
  )
  
} else {
  message("Skipping inoculation-associated heatmap (need >=10 pathways). Found: ", length(hm_paths))
}

# =========================
# 11) RHIZOSPHERE top pathways responding to Trichoderma inoculation
#     Colour-fixed version, matched to root plotting style
# =========================
# Exploratory rhizosphere matrix: nominal P and a 0.20 effect cutoff.
# This cutoff is more permissive than the 0.25 headline screen.
p_cut_matrix   <- PATHWAY_NOMINAL_P
eff_cut_matrix <- MATRIX_EFFECT_MIN
min_environments_supported <- 1
top_n_matrix <- 15
focus_matrix <- "Both"   # "Trichoderma" or "Both"

aldex_plot_fp <- file.path(tab_dir, "ALDEx2_Inoculation_byEnvironment.tsv")
if (file.exists(aldex_plot_fp)) {
  
  df_matrix <- readr::read_tsv(aldex_plot_fp, show_col_types = FALSE) %>%
    mutate(
      supported = (we.ep < p_cut_matrix) & (abs(effect) > eff_cut_matrix),
      direction = case_when(
        effect < 0 ~ "Higher in Trichoderma",
        effect > 0 ~ "Higher in Control",
        TRUE ~ "No change"
      ),
      direction = factor(direction, levels = c("Higher in Trichoderma", "Higher in Control", "No change"))
    )
  
  rank_tbl <- df_matrix %>%
    filter(direction != "No change") %>%
    group_by(Pathway) %>%
    summarise(
      n_environments_supported = sum(supported, na.rm = TRUE),
      median_abs_effect = median(abs(effect), na.rm = TRUE),
      min_p = min(we.ep, na.rm = TRUE),
      n_trichoderma_supported = sum(supported & direction == "Higher in Trichoderma", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_environments_supported >= min_environments_supported)
  
  if (focus_matrix == "Trichoderma") {
    rank_tbl <- rank_tbl %>% filter(n_trichoderma_supported >= 1)
  }
  
  rank_tbl <- rank_tbl %>%
    arrange(desc(n_environments_supported), min_p, desc(median_abs_effect)) %>%
    slice_head(n = top_n_matrix)
  
  top_paths_matrix <- rank_tbl$Pathway
  
  if (length(top_paths_matrix) > 0) {
    env_levels <- FUNCTION_ENV_LEVELS
    
    plot_df_matrix <- df_matrix %>%
      filter(Pathway %in% top_paths_matrix) %>%
      mutate(
        Pathway  = factor(Pathway, levels = rev(top_paths_matrix)),
        Environment = factor(Environment, levels = env_levels),
        alpha_val = ifelse(supported, 0.95, 0.20),
        color_key = case_when(
          direction == "Higher in Trichoderma" ~ "Higher in Trichoderma",
          direction == "Higher in Control" ~ "Higher in Control",
          TRUE ~ "No change"
        ),
        color_key = factor(color_key, levels = c("Higher in Trichoderma", "Higher in Control", "No change"))
      )
    
    title_txt <- paste0(
      "Top ", top_n_matrix,
      ifelse(focus_matrix == "Trichoderma", " rhizosphere pathways with Trichoderma signal", " rhizosphere pathways responding"),
      " supported across \u2265", min_environments_supported, " environment"
    )
    
    p_matrix <- ggplot(plot_df_matrix, aes(x = Environment, y = Pathway)) +
      geom_point(
        aes(
          size  = abs(effect),
          shape = direction,
          color = color_key,
          alpha = alpha_val
        ),
        stroke = 0.25
      ) +
      scale_x_discrete(labels = fig4_environment_axis_label, drop = FALSE) +
      scale_alpha_identity(guide = "none") +
      scale_size_continuous(name = "Effect size |effect|", range = c(1.5, 8)) +
      scale_shape_manual(
        name = "Direction",
        values = c("Higher in Trichoderma" = 16, "Higher in Control" = 17, "No change" = 1),
        drop = FALSE,
        guide = "none"
      ) +
      scale_color_manual(
        name = "Direction",
        values = c(dir_cols, "No change" = "grey80"),
        breaks = c("Higher in Trichoderma", "Higher in Control"),
        labels = trichoderma_direction_labels,
        guide = guide_legend(override.aes = list(alpha = 1, size = 3,
                                                 shape = c(16, 17)))
      ) +
      labs(title = trichoderma_title(title_txt), x = NULL, y = NULL) +
      theme_nature_matrix(base_size = FIG4_SOURCE_BASE_SIZE) +
      theme(
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.35),
        panel.grid.major.y = element_blank()
      )
    
    ggsave(file.path(fig_dir, "Fig_Rhizo_Top15_Pathways_SiteMatrix_Nature.pdf"),
           p_matrix, width = 10.8, height = 6.4, useDingbats = FALSE)
    ggsave(file.path(fig_dir, "Fig_Rhizo_Top15_Pathways_SiteMatrix_Nature.png"),
           p_matrix, width = 10.8, height = 6.4, dpi = 600)
    
    saveRDS(
      p_matrix,
      file.path(fig_dir, "p_Fig4A_Rhizo_Top15_Pathways_SiteMatrix.rds")
    )
    
    write_tsv2(rank_tbl, "TABLE_Rhizo_Top15_Pathways_SiteMatrix_Ranking.tsv")
    message("Section 11 figure saved to: ", file.path(fig_dir, "Fig_Rhizo_Top15_Pathways_SiteMatrix_Nature.pdf"))
  } else {
    message("Section 11 skipped: no pathways met the ranking criteria.")
  }
} else {
  message("Section 11 skipped: missing file ", aldex_plot_fp)
}


# ============================================================
# Supplementary Figure S3/S4 — merged site-specific volcano plots
# NatCom style
# ============================================================

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)
NATCOM_DIR <- STAGING_DIR
dir.create(NATCOM_DIR, recursive = TRUE, showWarnings = FALSE)

pretty_location <- function(x) {
  dplyr::recode(
    as.character(x),
    "Eyrewell_Forest" = "Eyrewell Forest",
    "Kowhai_Irrigated" = "Kowhai (Irrigated)",
    "Kowhai_Rainfed" = "Kowhai (Rainfed)",
    "West_Coast" = "West Coast",
    .default = as.character(x)
  )
}

plot_volcano_panel <- function(res_tbl, panel_title,
                               p_cut = 0.10,
                               effect_cut = 0.25,
                               label_n_right = 4,
                               label_n_left = 4,
                               xlim = c(-3, 3),
                               ylim = c(0, 3)) {
  
  df <- res_tbl %>%
    mutate(
      neglog10p = -log10(we.ep + 1e-300),
      Sig = (we.ep < p_cut) & (abs(effect) > effect_cut),
      Direction = case_when(
        Sig & effect > 0 ~ "Higher in Control",
        Sig & effect < 0 ~ "Higher in Trichoderma",
        TRUE ~ "Not highlighted"
      )
    )
  
  lab_right <- df %>%
    filter(Sig, effect > 0) %>%
    arrange(we.ep) %>%
    slice_head(n = label_n_right)
  
  lab_left <- df %>%
    filter(Sig, effect < 0) %>%
    arrange(we.ep) %>%
    slice_head(n = label_n_left)
  
  ggplot(df, aes(effect, neglog10p)) +
    geom_point(
      aes(colour = Direction),
      size = 1.2,
      alpha = 0.75
    ) +
    scale_colour_manual(
      values = c(
        "Higher in Trichoderma" = "#0072B2",
        "Higher in Control" = "#D55E00",
        "Not highlighted" = "grey70"
      ),
      breaks = c("Higher in Trichoderma", "Higher in Control", "Not highlighted"),
      labels = trichoderma_direction_labels
    ) +
    geom_vline(
      xintercept = c(-effect_cut, effect_cut),
      linetype = "dashed",
      linewidth = 0.3
    ) +
    geom_hline(
      yintercept = -log10(p_cut),
      linetype = "dashed",
      linewidth = 0.3
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim, clip = "off") +
    labs(
      title = panel_title,
      x = "Effect size (ALDEx2)",
      y = expression(-log[10](p)),
      colour = NULL
    ) +
    theme_nature(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0),
      axis.title = element_text(face = "bold", size = 9),
      axis.text = element_text(size = 8, colour = "black"),
      legend.position = "none",
      plot.margin = margin(6, 10, 6, 6)
    )
}

# Automatically set common y-axis limit from all site-level results
common_ylim <- c(
  0,
  ceiling(max(-log10(aldex_by_env$we.ep + 1e-300), na.rm = TRUE))
)

# Optional: restrict very tall p-value axis if labels become compressed
common_ylim[2] <- min(common_ylim[2], 4)

common_xlim <- c(
  floor(min(aldex_by_env$effect, na.rm = TRUE)),
  ceiling(max(aldex_by_env$effect, na.rm = TRUE))
)

# Optional: make all panels visually comparable
common_xlim[1] <- max(common_xlim[1], -4)
common_xlim[2] <- min(common_xlim[2], 4)

site_order <- FUNCTION_ENV_LEVELS

volcano_panels <- lapply(site_order, function(loc) {
  sub <- aldex_by_env %>% filter(Environment == loc)
  if (nrow(sub) == 0) return(NULL)
  
  plot_volcano_panel(
    sub,
    panel_title = pretty_location(loc),
    p_cut = 0.10,
    effect_cut = 0.25,
    xlim = common_xlim,
    ylim = common_ylim
  )
})

volcano_panels <- volcano_panels[!vapply(volcano_panels, is.null, logical(1))]

supp_volcano <- wrap_plots(
  volcano_panels,
  ncol = 3
) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 16
    )
  )

# Change filename depending on script:
# Rhizosphere script = FigureS3
# Root script = FigureS4

ggsave(
  file.path(NATCOM_DIR, "FigureS3_Rhizosphere_Volcano_Merged_NatCom.pdf"),
  supp_volcano,
  width = 240,
  height = 160,
  units = "mm",
  dpi = 600,
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  file.path(NATCOM_DIR, "FigureS3_Rhizosphere_Volcano_Merged_NatCom.png"),
  supp_volcano,
  width = 240,
  height = 160,
  units = "mm",
  dpi = 600,
  limitsize = FALSE
)

ggsave(
  file.path(NATCOM_DIR, "FigureS3_Rhizosphere_Volcano_Merged_NatCom.svg"),
  supp_volcano,
  width = 240,
  height = 160,
  units = "mm",
  limitsize = FALSE
)

# ----------------------------
# 12) Reproducibility snapshot
# ----------------------------
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))
message("DONE. Outputs saved to: ", out_dir)


})
})


# ============================================================
# Root unstratified functional analysis (Figures 3B, 4B, S4)
# ============================================================
run_section("Root unstratified functional analysis (Figures 3B, 4B, S4)", {
local({
# ============================================================
# ROOT shotgun metagenomics (48) — HUMAnN4 UNSTRATIFIED
# Focus: Inoculation (Trichoderma vs Control), controlling for Location
# Outputs:
#   - PERMANOVA, PCoA, NMDS
#   - ALDEx2 (global + within-Location)
#   - Volcanoes (Nature-style)
#   - Consistency table + headline pathways
#   - Heatmap A: top variable pathways
#   - Heatmap B: inoculation-associated pathways
#   - Section 11: Top pathways responding to Trichoderma inoculation
# UPDATED:
#   - standardised Eyrewell_Forest spelling
#   - Section 11 safely integrated without overwriting out_dir
# ============================================================

# ----------------------------
# 0) Packages (reproducible)
# ----------------------------
pkgs_cran <- c("tidyverse","vegan","ape","pheatmap","matrixStats","ggrepel")
to_install <- pkgs_cran[!sapply(pkgs_cran, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs_cran, library, character.only = TRUE))

if (!requireNamespace("ALDEx2", quietly = TRUE)) {
  message("Installing ALDEx2 from Bioconductor ...")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ALDEx2")
}
library(ALDEx2)

# ----------------------------
# 1) Paths
# ----------------------------
library(here)
project_dir <- repo_root
data_dir <- file.path(project_dir, "function", "data", "root")
out_dir  <- file.path(project_dir, "function", "results", "root", "unstratified")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cpm_fp  <- file.path(data_dir, "Root_pathabundance_UNSTRATIFIED_CPM.tsv.gz")
raw_fp  <- file.path(data_dir, "Root_pathabundance_UNSTRATIFIED.tsv.gz")
meta_fp <- file.path(data_dir, "Root_metadata.tsv.gz")

fig_dir <- file.path(out_dir, "Figures")
tab_dir <- file.path(out_dir, "Tables")
log_dir <- file.path(out_dir, "Logs")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 2) Nature Microbiology-ish helpers
# ----------------------------
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

# MUST match metadata values exactly
loc_cols <- c(
  "Eyrewell_Forest" = unname(okabe_ito["verm"]),
  "Kowhai"          = unname(okabe_ito["orange"]),
  "LU_H8"           = unname(okabe_ito["green"]),
  "Rolleston"       = unname(okabe_ito["blue"]),
  "West_Coast"      = unname(okabe_ito["purple"])
)

dir_cols <- c(
  "Higher in Trichoderma"   = unname(okabe_ito["blue"]),
  "Higher in Control" = unname(okabe_ito["verm"])
)

theme_nature <- function(base_size = 8, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title   = element_text(face = "bold", size = base_size + 1),
      axis.title   = element_text(face = "bold", size = base_size + 1),
      axis.text    = element_text(color = "black", size = base_size),
      axis.line    = element_line(linewidth = 0.35),
      axis.ticks   = element_line(linewidth = 0.35),
      legend.title = element_text(face = "bold", size = base_size + 0.5),
      legend.text  = element_text(size = base_size),
      legend.key   = element_blank(),
      legend.position = "right",
      panel.border = element_rect(fill = NA, linewidth = 0.35, color = "black"),
      plot.margin  = margin(6, 8, 6, 6)
    )
}

theme_nature_matrix <- function(base_size = 9, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title   = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      axis.title   = element_text(face = "bold"),
      axis.text.x  = element_text(color = "black", size = base_size, face = "bold"),
      axis.text.y  = element_text(color = "black", size = base_size - 1),
      axis.line    = element_line(linewidth = 0.35),
      axis.ticks   = element_line(linewidth = 0.35),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text  = element_text(size = base_size - 1),
      legend.key   = element_blank(),
      legend.position = "right",
      plot.margin  = margin(10, 40, 10, 10)
    )
}

save_fig <- function(p, name, w = 6.5, h = 4.8) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")),
         p, width = w, height = h, useDingbats = FALSE)
  ggsave(file.path(fig_dir, paste0(name, ".png")),
         p, width = w, height = h, dpi = 600)
}

write_tsv2 <- function(df, filename) {
  readr::write_tsv(df, file.path(tab_dir, filename))
}

read_humann_table <- function(fp) {
  tbl <- readr::read_tsv(fp, show_col_types = FALSE, comment = "")
  names(tbl)[1] <- "Pathway"
  mat <- as.data.frame(tbl)
  rownames(mat) <- mat$Pathway
  mat$Pathway <- NULL
  mat[] <- lapply(mat, function(x) as.numeric(x))
  mat[is.na(mat)] <- 0
  mat
}

save_pheatmap <- function(mat, anno, name, main,
                          w_pdf = 16, h_pdf = 7.5,
                          w_png = 6000, h_png = 3200, res = 600,
                          cluster_cols = TRUE, cluster_rows = TRUE,
                          show_colnames = FALSE, show_rownames = TRUE,
                          fontsize = 9, fontsize_row = 7, fontsize_col = 7,
                          gaps_col = NULL, gaps_row = NULL,
                          annotation_names_col = TRUE,
                          annotation_fontsize = 8) {
  
  pdf(file.path(fig_dir, paste0(name, ".pdf")),
      width = w_pdf, height = h_pdf, useDingbats = FALSE)
  pheatmap::pheatmap(
    mat,
    annotation_col = anno,
    cluster_cols = cluster_cols,
    cluster_rows = cluster_rows,
    show_colnames = show_colnames,
    show_rownames = show_rownames,
    border_color = NA,
    gaps_col = gaps_col,
    gaps_row = gaps_row,
    fontsize = fontsize,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    annotation_names_col = annotation_names_col,
    annotation_names_row = TRUE,
    annotation_fontsize = annotation_fontsize,
    main = main
  )
  dev.off()
  
  png(file.path(fig_dir, paste0(name, ".png")),
      width = w_png, height = h_png, res = res)
  pheatmap::pheatmap(
    mat,
    annotation_col = anno,
    cluster_cols = cluster_cols,
    cluster_rows = cluster_rows,
    show_colnames = show_colnames,
    show_rownames = show_rownames,
    border_color = NA,
    gaps_col = gaps_col,
    gaps_row = gaps_row,
    fontsize = fontsize,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    annotation_names_col = annotation_names_col,
    annotation_names_row = TRUE,
    annotation_fontsize = annotation_fontsize,
    main = main
  )
  dev.off()
}

plot_volcano_splitlabels <- function(res_tbl, title, fname_base,
                                     p_cut = 0.10, effect_cut = 0.25,
                                     label_n_right = 10, label_n_left = 6) {
  
  df <- res_tbl %>%
    mutate(
      neglog10p = -log10(we.ep + 1e-300),
      Sig = (we.ep < p_cut) & (abs(effect) > effect_cut),
      Direction = case_when(
        Sig & effect > 0 ~ "Higher in Control",
        Sig & effect < 0 ~ "Higher in Trichoderma",
        TRUE ~ "Not highlighted"
      )
    )
  
  p <- ggplot(df, aes(effect, neglog10p)) +
    geom_point(aes(colour = Direction), size = 1.7, alpha = 0.75) +
    scale_colour_manual(
      values = c(
        "Higher in Trichoderma" = "#0072B2",
        "Higher in Control" = "#D55E00",
        "Not highlighted" = "grey70"
      ),
      labels = trichoderma_direction_labels
    ) +
    geom_vline(xintercept = c(-effect_cut, effect_cut),
               linetype = 2, linewidth = 0.35) +
    geom_hline(yintercept = -log10(p_cut),
               linetype = 2, linewidth = 0.35) +
    theme_nature(base_size = 8) +
    labs(
      title = trichoderma_title(title),
      x = "Effect size (ALDEx2)",
      y = expression(-log[10](p)),
      colour = NULL
    )
  
  save_fig(p, fname_base, w = 7.2, h = 4.8)
}

# ----------------------------
# 3) Read metadata (SAFE)
# ----------------------------
meta <- readr::read_tsv(meta_fp, show_col_types = FALSE) %>%
  mutate(
    SampleID    = trimws(as.character(SampleID)),
    Location    = trimws(as.character(Location)),
    Inoculation = trimws(as.character(Inoculation)),
    Water       = trimws(as.character(Water))
  ) %>%
  mutate(
    Location = case_when(
      Location %in% c("Eyerell_Forest", "Eyerell Forest", "Eyrewell Forest") ~ "Eyrewell_Forest",
      TRUE ~ Location
    ),
    Inoculation = case_when(
      str_to_lower(Inoculation) %in% c("control") ~ "Control",
      str_to_lower(Inoculation) %in% c("panch", "trichoderma")   ~ "Trichoderma",
      TRUE ~ Inoculation
    )
  )

stopifnot(all(c("SampleID","Location","Inoculation","Water") %in% names(meta)))
if (anyDuplicated(meta$SampleID)) stop("Duplicate SampleID in metadata.")

# ----------------------------
# 4) Read HUMAnN tables (CPM + RAW)
# ----------------------------
cpm_mat <- read_humann_table(cpm_fp)
raw_mat <- read_humann_table(raw_fp)

# ----------------------------
# 5) Align samples
# ----------------------------
common_samples <- Reduce(intersect, list(meta$SampleID, colnames(cpm_mat), colnames(raw_mat)))
if (length(common_samples) < 6) stop("Too few shared samples. Check SampleID naming/columns.")

meta2 <- meta %>%
  filter(SampleID %in% common_samples) %>%
  arrange(SampleID) %>%
  mutate(
    Inoculation = factor(Inoculation, levels = c("Control","Trichoderma")),
    Water = stringr::str_to_title(trimws(as.character(Water))),
    Environment = factor(
      ifelse(Location == "Kowhai", paste0("Kowhai_", Water), Location),
      levels = FUNCTION_ENV_LEVELS
    )
  )

if (anyNA(meta2$Environment)) {
  stop("Cannot assign all samples to six environments; check Kowhai Water metadata.")
}

bad_loc <- setdiff(sort(unique(meta2$Location)), names(loc_cols))
if (length(bad_loc) > 0) {
  stop("These Location values are not in loc_cols: ", paste(bad_loc, collapse = ", "))
}
meta2 <- meta2 %>% mutate(Location = factor(Location, levels = names(loc_cols)))

cpm_mat2 <- cpm_mat[, meta2$SampleID, drop = FALSE]
raw_mat2 <- raw_mat[, meta2$SampleID, drop = FALSE]
stopifnot(identical(colnames(cpm_mat2), meta2$SampleID))
stopifnot(identical(colnames(raw_mat2), meta2$SampleID))

X_cpm <- t(as.matrix(cpm_mat2))
storage.mode(X_cpm) <- "numeric"
X_cpm[is.na(X_cpm)] <- 0
X_cpm <- X_cpm[, colSums(X_cpm) > 0, drop = FALSE]

X_raw <- t(as.matrix(raw_mat2))
storage.mode(X_raw) <- "numeric"
X_raw[is.na(X_raw)] <- 0
X_raw <- X_raw[, colSums(X_raw) > 0, drop = FALSE]

message("Samples used: ", nrow(X_cpm))
message("CPM pathways (nonzero): ", ncol(X_cpm))
message("RAW pathways (nonzero): ", ncol(X_raw))

write_tsv2(meta2, "Metadata_aligned.tsv")

# ----------------------------
# 6) Bray–Curtis + PERMANOVA + dispersion
# ----------------------------
meta2_perm <- meta2 %>% dplyr::filter(!is.na(Location), !is.na(Inoculation), !is.na(Environment))
env_counts <- meta2_perm %>% dplyr::count(Environment, Inoculation, name = "n")
write_tsv2(env_counts, "DIAGNOSTIC_Environment_by_Inoculation_counts.tsv")
if (nrow(env_counts) != 2L * length(FUNCTION_ENV_LEVELS) || any(env_counts$n < 2L)) {
  stop("Six environments with at least two samples per treatment are required for site-specific ALDEx2 volcanoes.")
}
X_cpm_perm <- X_cpm[meta2_perm$SampleID, , drop = FALSE]

set.seed(1)
d_bray <- vegan::vegdist(X_cpm_perm, method = "bray")

# Five-location global model: Kowhai irrigated and rainfed share one Location.
# Field_environment is descriptive and is not a sixth level in PERMANOVA.
field_environment <- ifelse(
  as.character(meta2_perm$Location) == "Kowhai",
  paste("Kowhai", as.character(meta2_perm$Water)),
  as.character(meta2_perm$Location)
)
model_scope <- tibble::tibble(
  Model = "Pathway PERMANOVA: Location + Inoculation",
  Location_levels_in_model = dplyr::n_distinct(meta2_perm$Location),
  Field_environments_in_dataset = dplyr::n_distinct(field_environment),
  Definition = paste(
    "Five locations in PERMANOVA; Kowhai water regimes pooled.",
    "Within-environment PERMANOVA and ALDEx2 analyse Kowhai Irrigated and Rainfed separately."
  )
)
write_tsv2(model_scope, "PERMANOVA_model_scope.tsv")
if (model_scope$Location_levels_in_model != 5 ||
    model_scope$Field_environments_in_dataset != 6) {
  warning("Expected five locations and six field environments; check Water metadata.")
}

# ----------------------------
# 6A) Global five-location PERMANOVA
# ----------------------------
perm_main <- vegan::adonis2(
  d_bray ~ Location + Inoculation,
  data = meta2_perm,
  permutations = 999,
  by = "margin"
)

perm_main_tbl <- as.data.frame(perm_main) %>%
  tibble::rownames_to_column("Term")

write_tsv2(
  perm_main_tbl,
  "PERMANOVA_main_Location_Inoculation.tsv"
)

# ----------------------------
# 6A2) Does the inoculation response vary among six field environments?
# ----------------------------
# This interaction compares treatment-associated shifts between environments
# within this compartment. It is distinct from six separate treatment tests.
# Treatment labels are permuted within environment. The metadata currently
# lack a plot/block ID, so the test does not represent paired plot blocking.
if (!identical(attr(d_bray, "Labels"), as.character(meta2_perm$SampleID))) {
  stop("PERMANOVA metadata order does not match pathway distance matrix.")
}
env_permutations <- permute::how(nperm = 999, blocks = meta2_perm$Environment)
set.seed(1)
perm_env_interaction <- vegan::adonis2(
  d_bray ~ Environment * Inoculation,
  data = meta2_perm, permutations = env_permutations, by = "margin"
)
perm_env_interaction_df <- as.data.frame(perm_env_interaction) %>%
  tibble::rownames_to_column("Term")
perm_env_interaction_row <- perm_env_interaction_df %>%
  dplyr::filter(Term %in% c("Environment:Inoculation", "Inoculation:Environment"))
if (nrow(perm_env_interaction_row) != 1L) {
  stop("Could not extract pathway Environment x Inoculation interaction: ",
       paste(perm_env_interaction_df$Term, collapse = ", "))
}
perm_env_interaction_row <- perm_env_interaction_row %>%
  dplyr::mutate(
    N = nrow(meta2_perm), Field_environments = length(FUNCTION_ENV_LEVELS),
    Permutations = 999L,
    Permutation_blocks = "Environment (six fields; Kowhai split)",
    .before = 1
  )
write_tsv2(perm_env_interaction_row,
           "PERMANOVA_Environment_by_Inoculation.tsv")
message("Functional Environment x Inoculation PERMANOVA saved in ", tab_dir)

# ----------------------------
# 6B) Inoculation effect with permutations constrained within Environment
# ----------------------------
perm_inoc_strata <- vegan::adonis2(
  d_bray ~ Inoculation,
  data = meta2_perm,
  permutations = 999,
  strata = meta2_perm$Environment
)

perm_inoc_strata_tbl <- as.data.frame(perm_inoc_strata) %>%
  tibble::rownames_to_column("Term")

write_tsv2(
  perm_inoc_strata_tbl,
  "PERMANOVA_Inoculation_stratifiedByEnvironment.tsv"
)
unlink(file.path(tab_dir, "PERMANOVA_Inoculation_stratifiedByLocation.tsv"))

# 6C) Within-environment PERMANOVA: six separate field environments
# ----------------------------
# The five-location global model above is retained. Field-specific treatment
# tests below split Kowhai into Irrigated and Rainfed (four samples per group).
perm_by_environment <- purrr::map_dfr(FUNCTION_ENV_LEVELS, function(env) {
  meta_env <- meta2_perm %>%
    dplyr::filter(as.character(Environment) == env) %>%
    droplevels()
  group_counts <- table(meta_env$Inoculation)
  if (nrow(meta_env) < 4L || length(group_counts) != 2L ||
      any(group_counts < 2L)) {
    stop("Expected both inoculation groups in field environment ", env)
  }
  X_env <- X_cpm_perm[meta_env$SampleID, , drop = FALSE]
  d_env <- vegan::vegdist(X_env, method = "bray")
  fit <- vegan::adonis2(
    d_env ~ Inoculation, data = meta_env, permutations = 999
  )
  fit_df <- as.data.frame(fit)
  if (nrow(fit_df) < 1L ||
      !rownames(fit_df)[1] %in% c("Model", "Inoculation")) {
    stop("Unexpected adonis2 output in ", env, ": ",
         paste(rownames(fit_df), collapse = ", "))
  }
  fit_df[1, , drop = FALSE] %>%
    dplyr::mutate(Environment = env, N = nrow(meta_env),
                  Term = "Inoculation", .before = 1) %>%
    dplyr::select(Environment, N, Term, Df, SumOfSqs, R2, F, `Pr(>F)`)
})
if (nrow(perm_by_environment) != 6L ||
    !setequal(perm_by_environment$Environment, FUNCTION_ENV_LEVELS)) {
  stop("Expected one inoculation PERMANOVA row per six field environments.")
}
write_tsv2(perm_by_environment, "PERMANOVA_withinEnvironment_Inoculation.tsv")
# Remove pooled five-location outputs produced by older script versions.
unlink(file.path(tab_dir, "PERMANOVA_withinLocation_Inoculation.tsv"))
message("Functional within-environment PERMANOVA: six separate field results saved.")

# ----------------------------
# 6D) Combined summary table: global + six environments
# ----------------------------
perm_main_combined <- perm_main_tbl %>%
  dplyr::filter(Term %in% c("Location", "Inoculation")) %>%
  dplyr::mutate(
    Model = "Global model",
    Location_label = "All sites",
    N = nrow(meta2_perm)
  ) %>%
  dplyr::select(Model, Location_label, N, Term, Df, SumOfSqs, R2, F, `Pr(>F)`)

perm_by_environment_combined <- perm_by_environment %>%
  dplyr::mutate(
    Model = "Within-environment",
    Location_label = pretty_environment_label(Environment)
  ) %>%
  dplyr::select(Model, Location_label, N, Term, Df, SumOfSqs, R2, F, `Pr(>F)`)

perm_combined <- dplyr::bind_rows(
  perm_main_combined,
  perm_env_interaction_row %>%
    dplyr::mutate(Model = "Environment interaction", Location_label = "Six environments") %>%
    dplyr::select(Model, Location_label, N, Term, Df, SumOfSqs, R2, F, `Pr(>F)`),
  perm_by_environment_combined
)

write_tsv2(
  perm_combined,
  "TABLE_PERMANOVA_ROOT_Function_Global_and_WithinEnvironment.tsv"
)
unlink(file.path(tab_dir, "TABLE_PERMANOVA_ROOT_Function_Global_and_WithinSite.tsv"))

# ----------------------------
# 6E) Dispersion tests
# ----------------------------
bd_inoc <- vegan::betadisper(d_bray, meta2_perm$Inoculation)
bd_loc  <- vegan::betadisper(d_bray, meta2_perm$Location)

capture.output(anova(bd_inoc), file = file.path(log_dir, "Betadisper_Inoculation_ANOVA.txt"))
capture.output(
  vegan::permutest(bd_inoc, permutations = 999),
  file = file.path(log_dir, "Betadisper_Inoculation_permutest.txt")
)

capture.output(anova(bd_loc), file = file.path(log_dir, "Betadisper_Location_ANOVA.txt"))
capture.output(
  vegan::permutest(bd_loc, permutations = 999),
  file = file.path(log_dir, "Betadisper_Location_permutest.txt")
)

# ----------------------------
# 7) Ordinations
# ----------------------------
pcoa <- ape::pcoa(as.matrix(d_bray))
scores_pcoa <- as.data.frame(pcoa$vectors[, 1:2]) %>%
  rownames_to_column("SampleID") %>%
  left_join(meta2_perm, by = "SampleID") %>%
  rename(PCoA1 = Axis.1, PCoA2 = Axis.2)
write_tsv2(scores_pcoa, "PCoA_scores.tsv")

p_pcoa <- ggplot(scores_pcoa, aes(PCoA1, PCoA2)) +
  geom_point(aes(color = Environment, shape = Inoculation),
             size = 2.2, alpha = 0.9, stroke = 0.35) +
  scale_color_manual(
    values = functional_environment_cols,
    labels = pretty_environment_label,
    drop = FALSE
  ) +
  scale_shape_manual(values = c(Control = 16, Trichoderma = 17),
                     labels = trichoderma_inoculation_labels, drop = FALSE) +
  theme_nature(base_size = 8) +
  labs(
    title = "Root PCoA (Bray–Curtis)",
    x = paste0("PCoA1 (", round(pcoa$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa$values$Relative_eig[2] * 100, 1), "%)")
  )
save_fig(p_pcoa, "Fig_PCoA_Location_Inoculation_Nature", w = 6.5, h = 4.8)

# ============================================================
# NatCom Figure 3B — Root functional pathway PCoA
# NatCom Figure 3C — PERMANOVA treatment-effect summary
# Combined Figure 3 — functional response to inoculation
# Consistent style with taxonomy Figure 2
# ============================================================
NATCOM_DIR <- STAGING_DIR
dir.create(NATCOM_DIR, recursive = TRUE, showWarnings = FALSE)

# Root pathway PCoA panel; keep the only legend here
p_fig3_root_func <- p_pcoa +
  labs(
    title = "Root",
    colour = "Field environment",
    shape = "Inoculation"
  ) +
  theme_nature(base_size = FIG3_PCOA_BASE_SIZE) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0, size = FIG3_PCOA_TITLE_SIZE),
    axis.title = element_text(face = "bold", size = FIG3_PCOA_AXIS_TITLE_SIZE),
    axis.text = element_text(colour = "black", size = FIG3_PCOA_AXIS_TEXT_SIZE),
    legend.title = element_text(face = "bold", size = FIG3_PCOA_LEGEND_TITLE_SIZE),
    legend.text = element_text(size = FIG3_PCOA_LEGEND_TEXT_SIZE),
    legend.key.size = grid::unit(4.2, "mm"),
    legend.spacing.y = grid::unit(1.2, "mm"),
    aspect.ratio = FIG3_PCOA_ASPECT_RATIO,
    plot.margin = margin(4, 4, 4, 4)
  ) +
  guides(
    colour = guide_legend(override.aes = list(size = 3.2), order = 1),
    shape  = guide_legend(override.aes = list(size = 3.2), order = 2)
  )

saveRDS(
  p_fig3_root_func,
  file.path(NATCOM_DIR, "p_fig3_root_func.rds")
)

ggsave(
  file.path(NATCOM_DIR, "Fig3B_Root_Function_PCoA_NatCom_panel.pdf"),
  p_fig3_root_func,
  width = 155,
  height = FIG3_HEIGHT_MM,
  units = "mm",
  useDingbats = FALSE
)

ggsave(
  file.path(NATCOM_DIR, "Fig3B_Root_Function_PCoA_NatCom_panel.png"),
  p_fig3_root_func,
  width = 155,
  height = FIG3_HEIGHT_MM,
  units = "mm",
  dpi = 600
)

# Panel C: global pathway PERMANOVA treatment effect by compartment
perm_fig3 <- tibble::tibble(
  Compartment = factor(
    c("Rhizosphere", "Root"),
    levels = c("Rhizosphere", "Root")
  ),
  Treatment_R2 = c(0.0068, 0.036),
  P_value = c(0.567, 0.009),
  Significance = c("Not significant", "Significant"),
  Label = c(
    "R² = 0.0068\nP = 0.567\nns",
    "R² = 0.036\nP = 0.009\n**"
  )
)

p_fig3_perm <- ggplot(
  perm_fig3,
  aes(x = Compartment, y = Treatment_R2, fill = Significance)
) +
  geom_col(
    width = 0.55,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = Label),
    vjust = -0.25,
    size = 4.6,
    lineheight = 0.9
  ) +
  scale_fill_manual(
    values = c(
      "Significant" = "#0072B2",
      "Not significant" = "grey75"
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 0.05),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    title = "Inoculation effect",
    x = NULL,
    y = expression(Treatment~R^2)
  ) +
  theme_nature(base_size = FIG3_PCOA_BASE_SIZE) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0, size = FIG3_PCOA_TITLE_SIZE),
    axis.title = element_text(face = "bold", size = FIG3_PCOA_AXIS_TITLE_SIZE),
    axis.text = element_text(colour = "black", size = FIG3_PCOA_AXIS_TEXT_SIZE),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.margin = margin(4, 4, 4, 4)
  )

saveRDS(
  p_fig3_perm,
  file.path(NATCOM_DIR, "p_fig3_perm.rds")
)

ggsave(
  file.path(NATCOM_DIR, "Fig3C_PERMANOVA_Function_Treatment_Effect_NatCom_panel.pdf"),
  p_fig3_perm,
  width = 100,
  height = FIG3_HEIGHT_MM,
  units = "mm",
  useDingbats = FALSE
)

ggsave(
  file.path(NATCOM_DIR, "Fig3C_PERMANOVA_Function_Treatment_Effect_NatCom_panel.png"),
  p_fig3_perm,
  width = 100,
  height = FIG3_HEIGHT_MM,
  units = "mm",
  dpi = 600
)

# Combined Figure 3. Run the rhizosphere script first so Panel A exists.
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("png", quietly = TRUE)) install.packages("png")
library(cowplot)
library(png)
library(grid)

rhiz_panel_png <- file.path(NATCOM_DIR, "Fig3A_Rhizosphere_Function_PCoA_NatCom_panel.png")

if (file.exists(rhiz_panel_png)) {
  p_fig3_rhiz_func_img <- cowplot::ggdraw() +
    cowplot::draw_grob(
      grid::rasterGrob(
        png::readPNG(rhiz_panel_png),
        interpolate = TRUE
      )
    )
  
  fig3_function_response <- cowplot::plot_grid(
    p_fig3_rhiz_func_img,
    p_fig3_root_func,
    p_fig3_perm,
    nrow = 1,
    labels = c("A", "B", "C"),
    label_size = FIG3_TAG_SIZE,
    label_fontface = "bold",
    rel_widths = c(0.82, 1.33, 0.88),
    align = "h",
    axis = "tb"
  )
  
  ggsave(
    file.path(NATCOM_DIR, "Figure3_Functional_Response_NatCom.pdf"),
    fig3_function_response,
    width = FIG3_WIDTH_MM,
    height = FIG3_HEIGHT_MM,
    units = "mm",
    dpi = 600,
    useDingbats = FALSE
  )
  
  ggsave(
    file.path(NATCOM_DIR, "Figure3_Functional_Response_NatCom.png"),
    fig3_function_response,
    width = FIG3_WIDTH_MM,
    height = FIG3_HEIGHT_MM,
    units = "mm",
    dpi = 600
  )
  
  ggsave(
    file.path(NATCOM_DIR, "Figure3_Functional_Response_NatCom.svg"),
    fig3_function_response,
    width = FIG3_WIDTH_MM,
    height = FIG3_HEIGHT_MM,
    units = "mm"
  )
} else {
  message(
    "Rhizosphere panel not found: ", rhiz_panel_png,
    "\nRun 03_rhizosphere_unstratified_NatCom_Fig3.R first, then rerun this root script to create merged Figure 3."
  )
}

set.seed(1)
nmds <- vegan::metaMDS(X_cpm_perm, distance = "bray", k = 2, trymax = 200,
                       autotransform = FALSE, trace = FALSE)
scores_nmds <- as.data.frame(vegan::scores(nmds, display = "sites")) %>%
  rownames_to_column("SampleID") %>%
  left_join(meta2_perm, by = "SampleID")
write_tsv2(scores_nmds, "NMDS_scores.tsv")

p_nmds <- ggplot(scores_nmds, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = Location, shape = Inoculation),
             size = 2.2, alpha = 0.9, stroke = 0.35) +
  scale_color_manual(values = loc_cols, drop = FALSE) +
  scale_shape_manual(values = c(Control = 16, Trichoderma = 17),
                     labels = trichoderma_inoculation_labels, drop = FALSE) +
  theme_nature(base_size = 8) +
  labs(title = paste0("Root NMDS (stress=", round(nmds$stress, 3), ")"))
save_fig(p_nmds, "Fig_NMDS_Location_Inoculation_Nature", w = 6.5, h = 4.8)

# ----------------------------
# 8) Heatmap A: Top variable pathways
# ----------------------------
N_var <- 30
vars <- matrixStats::colVars(X_cpm_perm)
top_var <- names(sort(vars, decreasing = TRUE))[1:min(N_var, length(vars))]

ZA <- scale(X_cpm_perm[, top_var, drop = FALSE])
hmA <- t(ZA)

annoA <- meta2_perm %>% dplyr::select(SampleID, Location, Inoculation) %>% as.data.frame()
rownames(annoA) <- annoA$SampleID
annoA$SampleID <- NULL
annoA <- annoA[colnames(hmA), , drop = FALSE]

save_pheatmap(
  hmA, annoA,
  name = "Fig_Heatmap_TopVariablePathways_CPM_GLOBAL",
  main = paste0("Top ", N_var, " variable pathways (z-scored CPM)"),
  w_pdf = 16, h_pdf = 7.0,
  cluster_cols = TRUE, cluster_rows = TRUE,
  show_colnames = FALSE, show_rownames = TRUE,
  fontsize = 9, fontsize_row = 7, fontsize_col = 7,
  annotation_fontsize = 8,
  annotation_names_col = TRUE
)

col_orderA <- meta2_perm %>% arrange(Inoculation, Location, SampleID) %>% pull(SampleID)
hmA_blk <- hmA[, col_orderA, drop = FALSE]
annoA_blk <- meta2_perm %>% dplyr::select(SampleID, Location, Inoculation) %>% as.data.frame()
rownames(annoA_blk) <- annoA_blk$SampleID
annoA_blk$SampleID <- NULL
annoA_blk <- annoA_blk[col_orderA, , drop = FALSE]
n_control <- sum(meta2_perm$Inoculation == "Control", na.rm = TRUE)

save_pheatmap(
  hmA_blk, annoA_blk,
  name = "Fig_Heatmap_TopVariablePathways_CPM_BLOCKED_ControlVsTrichoderma",
  main = paste0("Top ", N_var, " variable pathways (z-scored CPM) — blocked by inoculation"),
  w_pdf = 16, h_pdf = 7.0,
  cluster_cols = FALSE, cluster_rows = TRUE,
  show_colnames = FALSE, show_rownames = TRUE,
  gaps_col = n_control,
  fontsize = 9, fontsize_row = 7, fontsize_col = 7,
  annotation_fontsize = 8,
  annotation_names_col = TRUE
)

# ----------------------------
# 9) ALDEx2 + volcanoes + consistency
# ----------------------------
X_raw_int <- round(X_raw[meta2_perm$SampleID, , drop = FALSE])
X_raw_int[X_raw_int < 0] <- 0
storage.mode(X_raw_int) <- "integer"

run_aldex2_one <- function(X_sub, group, label, mc_samples = 128, denom = "all") {
  group_f <- droplevels(factor(group))
  if (nlevels(group_f) != 2) stop("ALDEx2 requires exactly 2 groups: ", label)
  
  X_sub <- as.matrix(X_sub)
  storage.mode(X_sub) <- "numeric"
  X_sub[is.na(X_sub)] <- 0
  X_sub <- round(X_sub)
  X_sub[X_sub < 0] <- 0
  
  counts <- t(X_sub)
  storage.mode(counts) <- "integer"
  
  ald <- ALDEx2::aldex.clr(
    reads = counts,
    conds = as.character(group_f),
    mc.samples = as.integer(mc_samples),
    denom = denom,
    verbose = FALSE
  )
  tt  <- ALDEx2::aldex.ttest(ald, paired.test = FALSE, verbose = FALSE)
  eff <- ALDEx2::aldex.effect(ald, verbose = FALSE)
  
  tibble::tibble(
    Pathway  = rownames(tt),
    we.eBH   = tt$we.eBH,
    we.ep    = tt$we.ep,
    effect   = eff$effect,
    diff.btw = eff$diff.btw,
    overlap  = eff$overlap,
    Contrast = label,
    Group1   = levels(group_f)[1],
    Group2   = levels(group_f)[2]
  ) %>% arrange(we.ep, desc(abs(effect)))
}

aldex_global <- run_aldex2_one(X_raw_int, meta2_perm$Inoculation, "Inoculation_Global")
write_tsv2(aldex_global, "ALDEx2_Inoculation_Global.tsv")

plot_volcano_splitlabels(
  aldex_global,
  "Differential pathways: Trichoderma vs Control (Global; RAW)",
  "Fig_ALDEx2_Volcano_Inoculation_Global_Nature",
  p_cut = 0.10, effect_cut = 0.25,
  label_n_right = 10, label_n_left = 10
)

aldex_by_loc <- purrr::map_dfr(levels(meta2_perm$Location), function(loc) {
  idx <- meta2_perm$Location == loc
  g <- droplevels(meta2_perm$Inoculation[idx])
  if (nlevels(g) != 2) return(tibble())
  
  run_aldex2_one(
    X_sub = X_raw_int[idx, , drop = FALSE],
    group = g,
    label = paste0("Inoculation_", loc),
    mc_samples = 128,
    denom = "all"
  ) %>% mutate(Location = loc)
})
write_tsv2(aldex_by_loc, "ALDEx2_Inoculation_byLocation.tsv")

# Separate field-environment contrasts for Supplementary Figures S3/S4.
# Existing five-location results and cross-location pathway screens stay pooled.
aldex_by_env <- purrr::map_dfr(FUNCTION_ENV_LEVELS, function(env) {
  idx <- as.character(meta2_perm$Environment) == env
  g <- droplevels(meta2_perm$Inoculation[idx])
  if (nlevels(g) != 2L || any(table(g) < 2L)) {
    stop("Insufficient treatment replicates for environment: ", env)
  }
  run_aldex2_one(
    X_sub = X_raw_int[idx, , drop = FALSE],
    group = g,
    label = paste0("Inoculation_", env),
    mc_samples = 128,
    denom = "all"
  ) %>% mutate(Environment = env, Location = as.character(meta2_perm$Location[which(idx)[1]]))
})
write_tsv2(aldex_by_env, "ALDEx2_Inoculation_byEnvironment.tsv")

for (loc in levels(meta2_perm$Location)) {
  sub <- aldex_by_loc %>% filter(Location == loc)
  if (nrow(sub) == 0) next
  
  plot_volcano_splitlabels(
    sub,
    paste0("Differential pathways: Trichoderma vs Control — ", loc, " (RAW)"),
    paste0("Fig_ALDEx2_Volcano_Inoculation_", loc, "_Nature"),
    p_cut = 0.10, effect_cut = 0.25,
    label_n_right = 10, label_n_left = 10
  )
}

# Screening for cross-location consistency is nominal, not BH significant.
MATRIX_EFFECT_MIN <- PATHWAY_EFFECT_MIN  # Root site matrix
support_ep  <- PATHWAY_NOMINAL_P
support_eff <- PATHWAY_EFFECT_MIN
write_tsv2(tibble::tibble(
  Selection = c("Exploratory pathway screen", "BH-adjusted pathway evidence"),
  P_column = c("we.ep", "we.eBH"),
  P_cutoff = c(PATHWAY_NOMINAL_P, PATHWAY_ADJUSTED_Q),
  Minimum_absolute_effect = PATHWAY_EFFECT_MIN,
  Matrix_absolute_effect = MATRIX_EFFECT_MIN,
  Used_for_headline_selection = c(TRUE, FALSE)
), "TABLE_PathwayThresholds.tsv")

consistency <- aldex_by_loc %>%
  mutate(
    supported = (we.ep < support_ep) & (abs(effect) > support_eff),
    adjusted_supported = (we.eBH < PATHWAY_ADJUSTED_Q) &
      (abs(effect) > support_eff),
    dir = case_when(
      effect > 0 ~ "Higher_in_Control",
      effect < 0 ~ "Higher_in_Trichoderma",
      TRUE ~ "Zero"
    )
  ) %>%
  group_by(Pathway) %>%
  summarise(
    n_locations_tested    = n_distinct(Location),
    n_locations_supported = sum(supported, na.rm = TRUE),
    n_locations_adjusted = sum(adjusted_supported, na.rm = TRUE),
    dir_supported         = paste(sort(unique(dir[supported])), collapse = ";"),
    median_effect         = median(effect, na.rm = TRUE),
    min_weep              = min(we.ep, na.rm = TRUE),
    min_weebh             = min(we.eBH, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_locations_supported), min_weep, desc(abs(median_effect)))

write_tsv2(consistency, "TABLE_PathwayConsistency_acrossLocations.tsv")
write_tsv2(
  consistency %>% filter(n_locations_adjusted >= 1),
  "TABLE_PathwayConsistency_BHAdjusted.tsv"
)

headline <- consistency %>%
  filter(n_locations_supported >= 2) %>%
  slice_head(n = 40) %>%
  mutate(Selection_basis = "Exploratory: we.ep < 0.10 and |effect| > 0.25 in at least two locations; not BH adjusted")
write_tsv2(headline, "TABLE_HeadlinePathways_sigIn2plusLocations.tsv")
message("Exploratory nominal-screen headline pathways (>=2 locations; not BH adjusted): ", nrow(headline))

headline1 <- tibble::tibble(Pathway = character())
if (nrow(headline) == 0) {
  headline1 <- consistency %>%
    filter(n_locations_supported >= 1) %>%
    slice_head(n = 40)
  write_tsv2(headline1, "TABLE_HeadlinePathways_sigIn1plusLocations_FALLBACK.tsv")
  message("Fallback headline (>=1 location supported): ", nrow(headline1))
}

# ----------------------------
# 10) Heatmap B: Inoculation-associated pathways
# ----------------------------
hm_paths <- if (nrow(headline) > 0) headline$Pathway else headline1$Pathway
hm_paths <- intersect(unique(hm_paths), colnames(X_cpm_perm))
hm_paths <- hm_paths[1:min(30, length(hm_paths))]
message("Inoculation-associated pathways available for heatmap: ", length(hm_paths))

if (length(hm_paths) >= 10) {
  
  ZB <- scale(X_cpm_perm[, hm_paths, drop = FALSE])
  hmB <- t(ZB)
  
  col_orderB <- meta2_perm %>% arrange(Inoculation, Location, SampleID) %>% pull(SampleID)
  hmB_blk <- hmB[, col_orderB, drop = FALSE]
  
  annoB_blk <- meta2_perm %>%
    dplyr::select(SampleID, Location, Inoculation) %>%
    as.data.frame()
  rownames(annoB_blk) <- annoB_blk$SampleID
  annoB_blk$SampleID <- NULL
  annoB_blk <- annoB_blk[col_orderB, , drop = FALSE]
  
  save_pheatmap(
    hmB_blk, annoB_blk,
    name = "Fig_Heatmap_HeadlineInoculationPathways_CPM_BLOCKED_ControlVsTrichoderma",
    main = "Inoculation-associated pathways (z-scored CPM)",
    w_pdf = 16, h_pdf = 7.0,
    cluster_cols = FALSE, cluster_rows = TRUE,
    show_colnames = FALSE, show_rownames = TRUE,
    gaps_col = n_control,
    fontsize = 9, fontsize_row = 7, fontsize_col = 7,
    annotation_fontsize = 8,
    annotation_names_col = TRUE
  )
  
  write_tsv2(tibble::tibble(Pathway = hm_paths),
             "TABLE_HeatmapPathways_InoculationAssociated.tsv")
  
} else {
  message("Skipping inoculation-associated heatmap (need >=10 pathways). Found: ", length(hm_paths))
}

# =========================
# 11) Top pathways responding to Trichoderma inoculation
#     safely integrated
# =========================
# Exploratory root matrix: same nominal and effect cutoffs as headline.
p_cut_matrix   <- PATHWAY_NOMINAL_P
eff_cut_matrix <- MATRIX_EFFECT_MIN
min_environments_supported <- 2
top_n_matrix <- 15
focus_matrix <- "Trichoderma"   # "Trichoderma" or "Both"

aldex_plot_fp <- file.path(tab_dir, "ALDEx2_Inoculation_byEnvironment.tsv")
if (file.exists(aldex_plot_fp)) {
  
  df_matrix <- readr::read_tsv(aldex_plot_fp, show_col_types = FALSE) %>%
    mutate(
      supported = (we.ep < p_cut_matrix) & (abs(effect) > eff_cut_matrix),
      direction = case_when(
        effect < 0 ~ "Higher in Trichoderma",
        effect > 0 ~ "Higher in Control",
        TRUE ~ "No change"
      ),
      direction = factor(direction, levels = c("Higher in Trichoderma", "Higher in Control", "No change"))
    )
  
  rank_tbl <- df_matrix %>%
    filter(direction != "No change") %>%
    group_by(Pathway) %>%
    summarise(
      n_environments_supported = sum(supported, na.rm = TRUE),
      median_abs_effect = median(abs(effect), na.rm = TRUE),
      min_p = min(we.ep, na.rm = TRUE),
      n_trichoderma_supported = sum(supported & direction == "Higher in Trichoderma", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_environments_supported >= min_environments_supported)
  
  if (focus_matrix == "Trichoderma") {
    rank_tbl <- rank_tbl %>% filter(n_trichoderma_supported >= 1)
  }
  
  rank_tbl <- rank_tbl %>%
    arrange(desc(n_environments_supported), min_p, desc(median_abs_effect)) %>%
    slice_head(n = top_n_matrix)
  
  top_paths_matrix <- rank_tbl$Pathway
  
  if (length(top_paths_matrix) > 0) {
    env_levels <- FUNCTION_ENV_LEVELS
    
    plot_df_matrix <- df_matrix %>%
      filter(Pathway %in% top_paths_matrix) %>%
      mutate(
        Pathway  = factor(Pathway, levels = rev(top_paths_matrix)),
        Environment = factor(Environment, levels = env_levels),
        alpha_val = ifelse(supported, 0.95, 0.20),
        color_key = case_when(
          direction == "Higher in Trichoderma" ~ "Higher in Trichoderma",
          direction == "Higher in Control" ~ "Higher in Control",
          TRUE ~ "No change"
        )
      )
    
    title_txt <- paste0(
      "Top ", top_n_matrix,
      ifelse(focus_matrix == "Trichoderma", " pathways with Trichoderma signal", " pathways responding"),
      " supported across \u2265", min_environments_supported, " environments"
    )
    
    p_matrix <- ggplot(plot_df_matrix, aes(x = Environment, y = Pathway)) +
      geom_point(
        aes(
          size  = abs(effect),
          shape = direction,
          color = color_key,
          alpha = alpha_val
        ),
        stroke = 0.25
      ) +
      scale_x_discrete(labels = fig4_environment_axis_label, drop = FALSE) +
      scale_alpha_identity(guide = "none") +
      scale_size_continuous(name = "Effect size |effect|", range = c(1.5, 8)) +
      scale_shape_manual(
        name = "Direction",
        values = c("Higher in Trichoderma" = 16, "Higher in Control" = 17, "No change" = 1),
        drop = FALSE,
        guide = "none"
      ) +
      scale_color_manual(
        name = "Direction",
        values = c(dir_cols, "No change" = "grey80"),
        breaks = c("Higher in Trichoderma", "Higher in Control"),
        labels = trichoderma_direction_labels,
        guide = guide_legend(override.aes = list(alpha = 1, size = 3,
                                                 shape = c(16, 17)))
      ) +
      labs(title = trichoderma_title(title_txt), x = NULL, y = NULL) +
      theme_nature_matrix(base_size = FIG4_SOURCE_BASE_SIZE) +
      theme(
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.35),
        panel.grid.major.y = element_blank()
      )
    
    ggsave(file.path(fig_dir, "Fig_Top15_Pathways_SiteMatrix_ROOT_Nature.pdf"),
           p_matrix, width = 10.8, height = 6.4, useDingbats = FALSE)
    ggsave(file.path(fig_dir, "Fig_Top15_Pathways_SiteMatrix_ROOT_Nature.png"),
           p_matrix, width = 10.8, height = 6.4, dpi = 600)
    
    saveRDS(
      p_matrix,
      file.path(fig_dir, "p_Fig4B_Root_Top15_Pathways_SiteMatrix.rds")
    )
    
    write_tsv2(rank_tbl, "TABLE_Top15_Pathways_SiteMatrix_ROOT_Ranking.tsv")
    message("Section 11 figure saved to: ", file.path(fig_dir, "Fig_Top15_Pathways_SiteMatrix_ROOT_Nature.pdf"))
  } else {
    message("Section 11 skipped: no pathways met the ranking criteria.")
  }
} else {
  message("Section 11 skipped: missing file ", aldex_plot_fp)
}

# ============================================================
# Supplementary Figure S3/S4 — merged site-specific volcano plots
# NatCom style
# ============================================================

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)
NATCOM_DIR <- STAGING_DIR
dir.create(NATCOM_DIR, recursive = TRUE, showWarnings = FALSE)

pretty_location <- function(x) {
  dplyr::recode(
    as.character(x),
    "Eyrewell_Forest" = "Eyrewell Forest",
    "Kowhai_Irrigated" = "Kowhai (Irrigated)",
    "Kowhai_Rainfed" = "Kowhai (Rainfed)",
    "West_Coast" = "West Coast",
    .default = as.character(x)
  )
}

plot_volcano_panel <- function(res_tbl, panel_title,
                               p_cut = 0.10,
                               effect_cut = 0.25,
                               xlim = c(-3, 3),
                               ylim = c(0, 3)) {
  
  df <- res_tbl %>%
    mutate(
      neglog10p = -log10(we.ep + 1e-300),
      Sig = (we.ep < p_cut) & (abs(effect) > effect_cut),
      Direction = case_when(
        Sig & effect > 0 ~ "Higher in Control",
        Sig & effect < 0 ~ "Higher in Trichoderma",
        TRUE ~ "Not highlighted"
      )
    )
  ggplot(df, aes(effect, neglog10p)) +
    geom_point(
      aes(colour = Direction),
      size = 1.2,
      alpha = 0.75
    ) +
    scale_colour_manual(
      values = c(
        "Higher in Trichoderma" = "#0072B2",
        "Higher in Control" = "#D55E00",
        "Not highlighted" = "grey70"
      ),
      breaks = c("Higher in Trichoderma", "Higher in Control", "Not highlighted"),
      labels = trichoderma_direction_labels
    ) +
    geom_vline(
      xintercept = c(-effect_cut, effect_cut),
      linetype = "dashed",
      linewidth = 0.3
    ) +
    geom_hline(
      yintercept = -log10(p_cut),
      linetype = "dashed",
      linewidth = 0.3
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim, clip = "off") +
    labs(
      title = panel_title,
      x = "Effect size (ALDEx2)",
      y = expression(-log[10](p)),
      colour = NULL
    ) +
    theme_nature(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0),
      axis.title = element_text(face = "bold", size = 9),
      axis.text = element_text(size = 8, colour = "black"),
      legend.position = "none",
      plot.margin = margin(6, 10, 6, 6)
    )
}

# Automatically set common y-axis limit from all site-level results
common_ylim <- c(
  0,
  ceiling(max(-log10(aldex_by_env$we.ep + 1e-300), na.rm = TRUE))
)

# Optional: restrict very tall p-value axis if labels become compressed
common_ylim[2] <- min(common_ylim[2], 4)

common_xlim <- c(
  floor(min(aldex_by_env$effect, na.rm = TRUE)),
  ceiling(max(aldex_by_env$effect, na.rm = TRUE))
)

# Optional: make all panels visually comparable
common_xlim[1] <- max(common_xlim[1], -4)
common_xlim[2] <- min(common_xlim[2], 4)

site_order <- FUNCTION_ENV_LEVELS

volcano_panels <- lapply(site_order, function(loc) {
  sub <- aldex_by_env %>% filter(Environment == loc)
  if (nrow(sub) == 0) return(NULL)
  
  plot_volcano_panel(
    sub,
    panel_title = pretty_location(loc),
    p_cut = 0.10,
    effect_cut = 0.25,
    xlim = common_xlim,
    ylim = common_ylim
  )
})

volcano_panels <- volcano_panels[!vapply(volcano_panels, is.null, logical(1))]

supp_volcano <- wrap_plots(
  volcano_panels,
  ncol = 3
) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 16
    )
  )

# Change filename depending on script:
# Rhizosphere script = FigureS3
# Root script = FigureS4

ggsave(
  file.path(NATCOM_DIR, "FigureS4_Root_Volcano_Merged_NatCom.pdf"),
  supp_volcano,
  width = 240,
  height = 160,
  units = "mm",
  dpi = 600,
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  file.path(NATCOM_DIR, "FigureS4_Root_Volcano_Merged_NatCom.png"),
  supp_volcano,
  width = 240,
  height = 160,
  units = "mm",
  dpi = 600,
  limitsize = FALSE
)

ggsave(
  file.path(NATCOM_DIR, "FigureS4_Root_Volcano_Merged_NatCom.svg"),
  supp_volcano,
  width = 240,
  height = 160,
  units = "mm",
  limitsize = FALSE
)

# ----------------------------
# 12) Reproducibility snapshot
# ----------------------------
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))
message("DONE. Outputs saved to: ", out_dir)


})
})


# ============================================================
# Merge Figure 3 functional response
# ============================================================
run_section("Merge Figure 3 functional response", {
local({
# ============================================================
# Figure 3 merge script — functional pathway response
# Run order:
#   1) 03_rhizosphere_unstratified_NatCom_Fig3_size_consistent.R
#   2) 04_root_unstratified_NatCom_Fig3_size_consistent.R
#   3) this script
# ============================================================

library(cowplot)
library(ggplot2)
NATCOM_DIR <- STAGING_DIR
# ----------------------------
# Load ggplot objects
# ----------------------------
p_fig3_rhiz_func <- readRDS(
  file.path(NATCOM_DIR, "p_fig3_rhiz_func.rds")
)

p_fig3_root_func <- readRDS(
  file.path(NATCOM_DIR, "p_fig3_root_func.rds")
)

p_fig3_perm <- readRDS(
  file.path(NATCOM_DIR, "p_fig3_perm.rds")
)

# ----------------------------
# Match taxonomy Figure 2 layout
# Panel A: no legend
# Panel B: one shared legend on right
# Panel C: compact PERMANOVA summary
# ----------------------------
p_fig3_rhiz_clean <- p_fig3_rhiz_func +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0, size = FIG3_PCOA_TITLE_SIZE),
    axis.title = element_text(face = "bold", size = FIG3_PCOA_AXIS_TITLE_SIZE),
    axis.text = element_text(colour = "black", size = FIG3_PCOA_AXIS_TEXT_SIZE),
    aspect.ratio = FIG3_PCOA_ASPECT_RATIO,
    plot.margin = margin(4, 4, 4, 4)
  )

p_fig3_root_clean <- p_fig3_root_func +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0, size = FIG3_PCOA_TITLE_SIZE),
    axis.title = element_text(face = "bold", size = FIG3_PCOA_AXIS_TITLE_SIZE),
    axis.text = element_text(colour = "black", size = FIG3_PCOA_AXIS_TEXT_SIZE),
    legend.title = element_text(face = "bold", size = FIG3_PCOA_LEGEND_TITLE_SIZE),
    legend.text = element_text(size = FIG3_PCOA_LEGEND_TEXT_SIZE),
    legend.key.size = grid::unit(4.2, "mm"),
    legend.spacing.y = grid::unit(1.2, "mm"),
    aspect.ratio = FIG3_PCOA_ASPECT_RATIO,
    plot.margin = margin(4, 4, 4, 4)
  ) +
  guides(
    colour = guide_legend(override.aes = list(size = 3.2), order = 1),
    shape = guide_legend(override.aes = list(size = 3.2), order = 2)
  )

p_fig3_perm_clean <- p_fig3_perm +
  labs(title = "Inoculation effect") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = FIG3_PCOA_TITLE_SIZE),
    axis.title = element_text(face = "bold", size = FIG3_PCOA_AXIS_TITLE_SIZE),
    axis.text = element_text(colour = "black", size = FIG3_PCOA_AXIS_TEXT_SIZE),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.margin = margin(4, 4, 4, 4)
  )

fig3_function_response <- cowplot::plot_grid(
  p_fig3_rhiz_clean,
  p_fig3_root_clean,
  p_fig3_perm_clean,
  nrow = 1,
  labels = c("A", "B", "C"),
  label_size = FIG3_TAG_SIZE,
  label_fontface = "bold",
  rel_widths = c(0.82, 1.33, 0.88),
  align = "h",
  axis = "tb"
)

# ----------------------------
# Save outputs
# ----------------------------
ggsave(
  file.path(NATCOM_DIR, "Figure3_Functional_Response_NatCom.pdf"),
  fig3_function_response,
  width = FIG3_WIDTH_MM,
  height = FIG3_HEIGHT_MM,
  units = "mm",
  dpi = 600,
  useDingbats = FALSE
)

ggsave(
  file.path(NATCOM_DIR, "Figure3_Functional_Response_NatCom.png"),
  fig3_function_response,
  width = FIG3_WIDTH_MM,
  height = FIG3_HEIGHT_MM,
  units = "mm",
  dpi = 600
)

ggsave(
  file.path(NATCOM_DIR, "Figure3_Functional_Response_NatCom.svg"),
  fig3_function_response,
  width = FIG3_WIDTH_MM,
  height = FIG3_HEIGHT_MM,
  units = "mm"
)

})
})


# ============================================================
# Merge Figure 4 site-matrix bubble plots
# ============================================================
run_section("Merge Figure 4 site-matrix bubble plots", {
local({
# ============================================================
# Figure 4. Functional pathway consistency bubble plots
# Merge rhizosphere and root ggplot RDS objects
# Nature Communications style
# ============================================================

# ----------------------------
# 0) Packages
# ----------------------------
pkgs <- c("ggplot2", "cowplot")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

# ----------------------------
# 1) Paths
# ----------------------------
project_dir <- repo_root
NATCOM_DIR <- STAGING_DIR
dir.create(NATCOM_DIR, recursive = TRUE, showWarnings = FALSE)

rhiz_rds <- file.path(
  project_dir,
  "function", "results", "rhizosphere", "unstratified", "Figures",
  "p_Fig4A_Rhizo_Top15_Pathways_SiteMatrix.rds"
)

root_rds <- file.path(
  project_dir,
  "function", "results", "root", "unstratified", "Figures",
  "p_Fig4B_Root_Top15_Pathways_SiteMatrix.rds"
)


stopifnot(file.exists(rhiz_rds))
stopifnot(file.exists(root_rds))

# ----------------------------
# 2) Read ggplot objects
# ----------------------------
p_rhiz <- readRDS(rhiz_rds)
p_root <- readRDS(root_rds)

# Short pathway labels for Figure 4 site-matrix panels
short_pathway_labels_fig4 <- function(x) {
  y <- as.character(x)
  y <- gsub("^.*?:\\s*", "", y)              # remove PWY-ID prefix
  y <- gsub("&alpha;|α|@", "alpha", y)
  y <- gsub("&beta;|β", "beta", y)
  y <- gsub("\\([^)]*\\)", "", y)          # remove parenthetical qualifiers
  y <- gsub("superpathway of the ", "", y, ignore.case = TRUE)
  y <- gsub("superpathway of ", "", y, ignore.case = TRUE)
  y <- gsub(" biosynthesis I$", " biosynthesis", y, ignore.case = TRUE)
  y <- gsub(" biosynthesis II$", " biosynthesis II", y, ignore.case = TRUE)
  y <- gsub(" and salvage$", "", y, ignore.case = TRUE)
  y <- gsub(" salvage I$", " salvage", y, ignore.case = TRUE)
  y <- gsub("^S-adenosyl-L-methionine salvage$", "SAM Salvage", y, ignore.case = TRUE)
  y <- gsub("^3,8-divinyl-chlorophyllide a biosynthesis I$", "Chlorophyllide Biosynthesis", y, ignore.case = TRUE)
  y <- gsub("^dTDP-alpha-D-ravidosamine and dTDP-4-acetyl-alpha-D-ravidosamine biosynthesis$", "dTDP-Ravidosamine Biosynthesis", y, ignore.case = TRUE)
  y <- gsub("^cytosolic NADPH production$", "Cytosolic NADPH Production", y, ignore.case = TRUE)
  y <- gsub("^D-galactose degradation I$", "D-Galactose Degradation I", y, ignore.case = TRUE)
  y <- gsub("^L-lysine fermentation to acetate and butanoate$", "L-Lysine Fermentation", y, ignore.case = TRUE)
  y <- gsub("^hexitol fermentation to lactate, formate, ethanol and acetate$", "Hexitol Fermentation", y, ignore.case = TRUE)
  y <- gsub("^glycogen biosynthesis II.*$", "Glycogen Biosynthesis II", y, ignore.case = TRUE)
  y <- gsub("^CMP-legionaminate biosynthesis I$", "CMP-Legionaminate Biosynthesis", y, ignore.case = TRUE)
  y <- gsub("^queuosine biosynthesis I.*$", "Queuosine Biosynthesis", y, ignore.case = TRUE)
  y <- gsub("^fatty acid biosynthesis I.*$", "Fatty Acid Biosynthesis I", y, ignore.case = TRUE)
  y <- gsub("^fatty acid biosynthesis II.*$", "Fatty Acid Biosynthesis II", y, ignore.case = TRUE)
  y <- gsub("^geranylgeranyldiphosphate biosynthesis( I.*)?$", "GGPP Biosynthesis", y, ignore.case = TRUE)
  y <- gsub("^3[-−]hydroxypropanoate/4[-−]hydroxybutanate cycle$", "3-HP/4-HB Cycle", y, ignore.case = TRUE)
  y <- gsub("^homocysteine and cysteine interconversion$", "Cys–Homocys Interconversion", y, ignore.case = TRUE)
  y <- gsub("^tetrahydrofolate biosynthesis$", "Tetrahydrofolate Biosynthesis", y, ignore.case = TRUE)
  y <- gsub("^pyridoxal 5'-phosphate biosynthesis$", "Pyridoxal 5'-Phosphate Biosynthesis", y, ignore.case = TRUE)
  y <- gsub("^taurine degradation$", "Taurine Degradation", y, ignore.case = TRUE)
  y <- gsub("^stachyose degradation$", "Stachyose Degradation", y, ignore.case = TRUE)
  y <- gsub("^mannitol cycle$", "Mannitol Cycle", y, ignore.case = TRUE)
  y <- gsub("^petroselinate biosynthesis$", "Petroselinate Biosynthesis", y, ignore.case = TRUE)
  y <- gsub("^sulfide oxidation$", "Sulfide Oxidation", y, ignore.case = TRUE)
  y <- gsub("\\s+", " ", y)
  y <- trimws(y)
  y <- tools::toTitleCase(y)
  # Display abbreviations only. Full MetaCyc descriptions are saved in the
  # adjacent Figure4_Pathway_Label_Key.tsv for unambiguous interpretation.
  y <- gsub("^Hentriaconta-.*-Nonaene Biosynthesis$",
            "Hentriaconta-Nonaene Biosynth.", y, ignore.case = TRUE)
  y <- gsub("^4-Amino-2-Methyl-5-Diphosphomethylpyrimidine Biosynthesis II$",
            "5-DPMP Biosynth. II", y, ignore.case = TRUE)
  y <- gsub("^D-Myo-Inositol.?Trisphosphate Biosynthesis$",
            "Inositol Trisphosphate Biosynth.", y, ignore.case = TRUE)
  y <- gsub("^Nucleoside and Nucleotide Degradation$",
            "Nucleoside/Nucleotide Degrad.", y, ignore.case = TRUE)
  y <- gsub("^Pyridoxal 5'-Phosphate Biosynthesis$",
            "Pyridoxal 5'-P Biosynth.", y, ignore.case = TRUE)
  y <- gsub("^Demethylmenaquinol-8 Biosynthesis$",
            "Demethylmenaquinol-8 Biosynth.", y, ignore.case = TRUE)
  y <- gsub("^Fatty Acid Biosynthesis (I|II)$",
            "Fatty Acid Biosynth. \\1", y, ignore.case = TRUE)
  y <- gsub("^CMP-Legionaminate Biosynthesis$",
            "CMP-Legionaminate Biosynth.", y, ignore.case = TRUE)
  y <- gsub("^Glycogen Biosynthesis II$",
            "Glycogen Biosynth. II", y, ignore.case = TRUE)
  # Match the displayed long names regardless of pathway prefixes or case.
  y[grepl("Geranylgeranyldiphosphate", y, ignore.case = TRUE)] <-
    "GGPP Biosynth. I"
  y[grepl("^Hexitol Fermentation", y, ignore.case = TRUE)] <-
    "Hexitol Fermentation"
  y[grepl("Hydroxypropanoate/4-Hydroxybutanate", y, ignore.case = TRUE)] <-
    "3-HP/4-HB Cycle"
  y[grepl("^D-Myo-Inositol.*Trisphosphate", y, ignore.case = TRUE)] <-
    "Inositol Trisphosphate"
  y[grepl("^Hentriaconta.*Nonaene", y, ignore.case = TRUE)] <-
    "Hentriaconta-Nonaene"
  y[grepl("^Demethylmenaquinol-8", y, ignore.case = TRUE)] <-
    "Demethylmenaquinol-8"
  y[grepl("^Nucleoside/Nucleotide", y, ignore.case = TRUE)] <-
    "Nucleoside/Nucleotide"
  y[grepl("^CMP-Legionaminate", y, ignore.case = TRUE)] <-
    "CMP-Legionaminate"
  y <- gsub(" Biosynthesis", " Biosynth.", y, fixed = TRUE)
  y <- gsub(" Degradation", " Degrad.", y, fixed = TRUE)
  y <- sub("^-", "", y)
  # Compact Figure 4 labels; preserve pathway distinctions (such as I/II).
  # Figure4_Pathway_Label_Key.tsv retains the full MetaCyc descriptions.
  replacements <- c(
    "^GGPP Biosynth\\. I$" = "GGPP biosyn. I",
    "^Hexitol Fermentation$" = "Hexitol ferment.",
    "^Glycogen Biosynth\\. II$" = "Glycogen biosyn. II",
    "^CMP-Legionaminate.*$" = "CMP-legionaminate",
    "^Nucleoside/Nucleotide.*$" = "Nucleoside/nucleotide",
    "^Ethanolamine Utilization$" = "Ethanolamine use",
    "^Acetylene Degrad\\.$" = "Acetylene degrad.",
    "^Fatty Acid Biosynth\\. (I|II)$" = "Fatty acid biosyn. \\1",
    "^Petroselinate Biosynth\\.$" = "Petroselinate biosyn.",
    "^Mevalonate Pathway I$" = "Mevalonate I",
    "^Polyamine Biosynth\\. II$" = "Polyamine biosyn. II",
    "^Queuosine Biosynth\\.( I)?$" = "Queuosine biosyn. I",
    "^Hentriaconta-Nonaene.*$" = "Hentriaconta-nonaene",
    "^Cis-Alkene Biosynth\\.$" = "Cis-alkene biosyn.",
    "^D-Cycloserine Biosynth\\.$" = "D-Cycloserine biosyn.",
    "^Pyridoxal 5'-P Biosynth\\.$" = "Pyridoxal 5'-P",
    "^Inositol Trisphosphate.*$" = "Inositol trisphos.",
    "^Glyoxylate Assimilation$" = "Glyoxylate assim.",
    "^1,3-Propanediol Biosynth\\.$" = "1,3-Propanediol",
    "^Chitin Degrad\\. I$" = "Chitin degrad. I",
    "^3-HP/4-HB Cycle$" = "3-HP/4-HB cycle",
    "^Dodecenoate Biosynth\\. II$" = "Dodecenoate biosyn. II",
    "^Nitrate Reduction I$" = "Nitrate reduction I",
    "^5-DPMP Biosynth\\. II$" = "5-DPMP biosyn. II",
    "^Demethylmenaquinol-8$" = "Demethylmenaquinol-8",
    "^Stachyose Degrad\\.$" = "Stachyose degrad.",
    "^Camphor Degrad\\.$" = "Camphor degrad."
  )
  for (pattern in names(replacements)) {
    y <- gsub(pattern, replacements[[pattern]], y, ignore.case = TRUE)
  }
  # Bound unfamiliar future names at a word boundary rather than allowing
  # a long y-axis label to consume the width of both dot grids.
  long <- nchar(y) > 23L
  y[long] <- paste0(sub("[[:space:][:punct:]]+$", "", substr(y[long], 1L, 20L)), "…")
  y
}

clean_fig4_plot <- function(p) {
  if (!is.null(p$data) && "Pathway" %in% names(p$data)) {
    old <- as.character(p$data$Pathway)
    new <- short_pathway_labels_fig4(old)
    old_levels <- if (is.factor(p$data$Pathway)) levels(p$data$Pathway) else unique(old)
    new_levels <- short_pathway_labels_fig4(old_levels)
    if (anyDuplicated(new_levels)) {
      stop("Figure 4 has duplicate short pathway names; revise the label map.")
    }
    p$data$Pathway <- factor(new, levels = new_levels)
  }
  p +
    scale_y_discrete(labels = identity) +
    guides(shape = "none")
}

fig4_label_key <- dplyr::bind_rows(
  tibble::tibble(Panel = "A: Rhizosphere",
                 Full_pathway = unique(as.character(p_rhiz$data$Pathway))),
  tibble::tibble(Panel = "B: Root",
                 Full_pathway = unique(as.character(p_root$data$Pathway)))
) %>%
  dplyr::mutate(Figure_label = short_pathway_labels_fig4(Full_pathway))
readr::write_tsv(
  fig4_label_key,
  file.path(NATCOM_DIR, "Figure4_Pathway_Label_Key.tsv")
)

# Replace titles
p_rhiz <- clean_fig4_plot(p_rhiz) +
  labs(title = trichoderma_title("Top 15 rhizosphere pathways with Trichoderma signal"))

p_root <- clean_fig4_plot(p_root) +
  labs(title = trichoderma_title("Top 15 root pathways with Trichoderma signal"))

# ----------------------------
# 3) Clean panels for merging
# ----------------------------

# Panel A: remove legend and increase readability
p_rhiz_clean <- p_rhiz +
  guides(shape = "none") +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = FIG4_PATHWAY_TEXT_SIZE, colour = "black"),
    axis.text.x = element_text(
      angle = 25,
      hjust = 1,
      vjust = 1,
      size = FIG4_FIELD_TEXT_SIZE,
      face = "bold",
      colour = "black"
    ),
    axis.title = element_text(size = FIG4_SITE_TEXT_SIZE, face = "bold"),
    plot.title = element_text(
      size = FIG4_TITLE_SIZE,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 4)
    ),
    plot.title.position = "plot",
    plot.margin = margin(6, 8, 14, 8)
  )

# Panel B: keep one shared legend and increase readability
p_root_clean <- p_root +
  guides(
    shape = "none",
    colour = guide_legend(
      title = "Direction",
      order = 1,
      override.aes = list(alpha = 1, size = 5, shape = c(16, 17))
    ),
    size = guide_legend(
      title = "Effect size |effect|",
      order = 2
    )
  ) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = grid::unit(0.35, "cm"),
    legend.margin = margin(4, 4, 4, 4),
    legend.title = element_text(size = FIG4_LEGEND_TITLE_SIZE, face = "bold"),
    legend.text = element_text(size = FIG4_LEGEND_TEXT_SIZE),
    legend.key.height = grid::unit(7, "mm"),
    axis.text.y = element_text(size = FIG4_PATHWAY_TEXT_SIZE, colour = "black"),
    axis.text.x = element_text(
      angle = 25,
      hjust = 1,
      vjust = 1,
      size = FIG4_FIELD_TEXT_SIZE,
      face = "bold",
      colour = "black"
    ),
    axis.title = element_text(size = FIG4_SITE_TEXT_SIZE, face = "bold"),
    plot.title = element_text(size = FIG4_TITLE_SIZE, face = "bold", hjust = 0.5),
    plot.margin = margin(6, 8, 14, 8)
  )

# ----------------------------
# 4) Combine panels
# ----------------------------

# Put the shared legend in its own column. Give the two legend-free plots
# identical gtable widths so their six-environment dot grids have the same width,
# even when their pathway labels have different lengths.
fig4_legend <- cowplot::get_legend(p_root_clean)
if (is.null(fig4_legend)) stop("Could not extract the Figure 4 legend.")
fig4_rhiz_grob <- ggplot2::ggplotGrob(p_rhiz_clean)
fig4_root_grob <- ggplot2::ggplotGrob(p_root_clean + theme(legend.position = "none"))
if (length(fig4_rhiz_grob$widths) != length(fig4_root_grob$widths)) {
  stop("Figure 4 panels have incompatible layouts; cannot equalise plot widths.")
}
fig4_shared_widths <- grid::unit.pmax(fig4_rhiz_grob$widths, fig4_root_grob$widths)
fig4_rhiz_grob$widths <- fig4_shared_widths
fig4_root_grob$widths <- fig4_shared_widths

fig4_legend_width <- 0.62
fig4_body <- cowplot::plot_grid(
  fig4_rhiz_grob,
  fig4_root_grob,
  fig4_legend,
  nrow = 1,
  rel_widths = c(1, 1, fig4_legend_width)
)

fig4_function_consistency <- cowplot::ggdraw() +
  cowplot::draw_plot(
    fig4_body,
    x = 0,
    y = 0,
    width = 1,
    height = 0.96
  ) +
  cowplot::draw_label(
    "A",
    x = 0.005,
    y = 0.995,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = FIG4_TAG_SIZE
  ) +
  cowplot::draw_label(
    "B",
    x = 1 / (2 + fig4_legend_width) + 0.005,
    y = 0.995,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = FIG4_TAG_SIZE
  )

# ----------------------------
# 5) Save outputs
# ----------------------------
ggsave(
  file.path(NATCOM_DIR, "Figure4_Functional_Pathway_Consistency_NatCom.pdf"),
  fig4_function_consistency,
  width = 530,
  height = 205,
  units = "mm",
  useDingbats = FALSE,
  bg = "white"
)

ggsave(
  file.path(NATCOM_DIR, "Figure4_Functional_Pathway_Consistency_NatCom.png"),
  fig4_function_consistency,
  width = 530,
  height = 205,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(NATCOM_DIR, "Figure4_Functional_Pathway_Consistency_NatCom.svg"),
  fig4_function_consistency,
  width = 530,
  height = 205,
  units = "mm",
  bg = "white"
)
message("Figure 4 saved to: ", NATCOM_DIR)


})
})


# ============================================================
# Rhizosphere stratified contributors (Figure 5 and S5 panels)
# ============================================================
run_section("Rhizosphere stratified contributors (Figure 5 and S5 panels)", {
local({
# ============================================================
# RHIZOSPHERE shotgun metagenomics (48) — HUMAnN4 STRATIFIED
# Goal: "Who does what?" taxon-resolved contributors to pathways
# Output:
#   - Per-pathway stacked bars (sample + group means)
#   - Heatmap (log10 proportion)
#   - Dotplots (two options: incl Unclassified vs classified-only)
# Reproducibility:
#   - Run log + session info + input MD5 + output manifest
# UPDATED: standardised Eyrewell_Forest spelling
# ============================================================

# ----------------------------
# 0) Packages
# ----------------------------
pkgs <- c("tidyverse","pheatmap","matrixStats","digest","sessioninfo","grid")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

options(error = NULL)
set.seed(1)

# ----------------------------
# 1) Paths
# ----------------------------
library(here)
project_dir <- repo_root
data_dir <- file.path(project_dir, "function", "data", "rhizosphere")
out_dir  <- file.path(project_dir, "function", "results", "rhizosphere", "stratified")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

strat_fp <- file.path(data_dir, "Rhizosphere_pathabundance_CPM_stratified_Rprefix.tsv.gz")
meta_fp  <- file.path(data_dir, "Rhizo_metadata.tsv")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

subdirs <- c("ByPathway_SampleStacks","ByPathway_GroupMeans","Heatmaps","DotPlots","Tables","RunInfo")
for (d in subdirs) dir.create(file.path(out_dir, d), showWarnings = FALSE, recursive = TRUE)

headline_fp <- file.path(
  project_dir,
  "function", "results", "rhizosphere", "unstratified", "Tables",
  "TABLE_HeadlinePathways_sigIn2plusLocations.tsv"
)

# ----------------------------
# 1b) Reproducibility: run log, session info, hashes
# ----------------------------
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_fp <- file.path(out_dir, "RunInfo", paste0("RUNLOG_", run_id, ".txt"))

zz <- file(log_fp, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")

.onExit_sinks <- function() {
  try(sink(type="message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(zz), silent = TRUE)
}
on.exit(.onExit_sinks(), add = TRUE)

assert_exists <- function(fp, label = fp) {
  if (!file.exists(fp)) stop("Missing file: ", label, "\nPath: ", fp, call. = FALSE)
}

md5_or_na <- function(fp) if (file.exists(fp)) as.character(tools::md5sum(fp)) else NA_character_

message("=== HUMAnN4 STRATIFIED WhoDoesWhat (RHIZO) ===")
message("Run ID: ", run_id)
message("Timestamp: ", Sys.time())
message("R version: ", R.version.string)
message("Platform: ", R.version$platform)
message("project_dir: ", project_dir)
message("data_dir:    ", data_dir)
message("strat_fp:    ", strat_fp, "  MD5=", md5_or_na(strat_fp))
message("meta_fp:     ", meta_fp,  "  MD5=", md5_or_na(meta_fp))
message("headline_fp: ", headline_fp, "  MD5=", md5_or_na(headline_fp))
message("out_dir:     ", out_dir)

si_fp <- file.path(out_dir, "RunInfo", paste0("SESSIONINFO_", run_id, ".txt"))
writeLines(capture.output(sessioninfo::session_info()), si_fp)
message("Session info saved: ", si_fp)

# ----------------------------
# 2) Helpers
# ----------------------------
write_tsv2 <- function(df, filename) {
  readr::write_tsv(df, file.path(out_dir, "Tables", filename))
}

save_gg <- function(p, outbase, w=12, h=7) {
  ggsave(paste0(outbase, ".pdf"), p, width=w, height=h, useDingbats = FALSE)
  ggsave(paste0(outbase, ".png"), p, width=w, height=h, dpi=600)
}

theme_nm <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "white"),
      legend.title = element_text(face = "bold")
    )
}

parse_strat <- function(x) {
  parts <- strsplit(x, "\\|")[[1]]
  if (length(parts) == 1) return(list(PathwayBase = parts[1], Taxon = "UNSTRATIFIED"))
  list(PathwayBase = parts[1], Taxon = paste(parts[-1], collapse="|"))
}

clean_taxon <- function(taxon) {
  taxon %>%
    stringr::str_replace("^s__", "") %>%
    stringr::str_replace("\\.t__SGB\\d+$", "") %>%
    stringr::str_replace("^unclassified$", "Unclassified")
}

drop_pwy_id <- function(x) sub("^.*?[:：∶]\\s*", "", as.character(x), perl = TRUE)

wrap_text <- function(x, width = 28) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

collapse_other <- function(df, top_n = 12, min_prop = 0.01) {
  tot <- sum(df$Value, na.rm = TRUE)
  df <- df %>% mutate(Prop = ifelse(tot > 0, Value / tot, 0))
  keep <- df %>% arrange(desc(Value)) %>% slice_head(n = top_n) %>% pull(Taxon)
  keep <- union(keep, df %>% filter(Prop >= min_prop) %>% pull(Taxon))
  
  df %>%
    mutate(Taxon2 = ifelse(Taxon %in% keep, Taxon, "Other")) %>%
    group_by(Taxon2) %>%
    summarise(Value = sum(Value, na.rm = TRUE), .groups="drop") %>%
    rename(Taxon = Taxon2)
}

# ----------------------------
# 3) Read metadata
# ----------------------------
assert_exists(meta_fp, "meta_fp")
meta <- readr::read_tsv(meta_fp, show_col_types = FALSE) %>%
  mutate(
    SampleID    = trimws(as.character(SampleID)),
    Location    = trimws(as.character(Location)),
    Inoculation = trimws(as.character(Inoculation)),
    Water       = trimws(as.character(Water))
  ) %>%
  mutate(
    # safety standardisation
    Location = case_when(
      Location %in% c("Eyerell_Forest", "Eyerell Forest", "Eyrewell Forest") ~ "Eyrewell_Forest",
      TRUE ~ Location
    ),
    Inoculation = case_when(
      str_to_lower(Inoculation) == "control" ~ "Control",
      str_to_lower(Inoculation) %in% c("panch", "trichoderma") ~ "Trichoderma",
      TRUE ~ Inoculation
    )
  ) %>%
  mutate(
    Location    = factor(Location),
    Inoculation = factor(Inoculation, levels = c("Control", "Trichoderma")),
    Water       = factor(Water)
  )

stopifnot(all(c("SampleID","Location","Inoculation") %in% names(meta)))
stopifnot(!anyDuplicated(meta$SampleID))
message("Metadata rows: ", nrow(meta))

# ----------------------------
# 4) Read stratified CPM table
# ----------------------------
assert_exists(strat_fp, "strat_fp")
strat_tbl <- readr::read_tsv(strat_fp, show_col_types = FALSE)

nm1 <- names(strat_tbl)[1]
if (grepl("^#\\s*Pathway", nm1)) names(strat_tbl)[1] <- "Pathway"
if (!"Pathway" %in% names(strat_tbl)) strat_tbl <- strat_tbl %>% dplyr::rename(Pathway = 1)

mat <- as.data.frame(strat_tbl)
rownames(mat) <- mat$Pathway
mat$Pathway <- NULL
mat[] <- lapply(mat, function(x) suppressWarnings(as.numeric(x)))
mat[is.na(mat)] <- 0

message("Stratified table: rows=", nrow(mat), " cols=", ncol(mat))

meta2 <- meta %>% filter(SampleID %in% colnames(mat))
stopifnot(nrow(meta2) > 0)
mat2 <- mat[, meta2$SampleID, drop = FALSE]
stopifnot(identical(colnames(mat2), meta2$SampleID))

keep_rows <- !grepl("^(UNMAPPED|UNINTEGRATED)", rownames(mat2))
mat2 <- mat2[keep_rows, , drop = FALSE]
message("Rows after dropping UNMAPPED/UNINTEGRATED: ", nrow(mat2))

# ----------------------------
# 5) Select target pathways
# ----------------------------
targets <- character()

if (file.exists(headline_fp)) {
  message("Using headline targets from: ", headline_fp)
  head_tbl <- readr::read_tsv(headline_fp, show_col_types = FALSE)
  
  if ("Pathway" %in% names(head_tbl)) {
    targets <- head_tbl %>% pull(Pathway)
  } else if ("PathwayBase" %in% names(head_tbl)) {
    targets <- head_tbl %>% pull(PathwayBase)
  } else {
    stop("Headline table missing Pathway/PathwayBase column.", call. = FALSE)
  }
  
  targets <- targets %>% as.character() %>% discard(is.na) %>% unique()
  
} else {
  message("No headline table found; selecting top variable pathways (base totals).")
  rowinfo <- tibble(Row = rownames(mat2)) %>%
    mutate(PathwayBase = purrr::map_chr(Row, ~parse_strat(.x)$PathwayBase))
  
  long0 <- mat2 %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Row") %>%
    tidyr::pivot_longer(-Row, names_to="SampleID", values_to="Value") %>%
    left_join(rowinfo, by="Row") %>%
    group_by(PathwayBase, SampleID) %>%
    summarise(Value = sum(Value, na.rm=TRUE), .groups="drop")
  
  wide0 <- long0 %>% pivot_wider(names_from=SampleID, values_from=Value, values_fill=0)
  X0 <- wide0 %>% column_to_rownames("PathwayBase") %>% as.matrix()
  vars <- matrixStats::rowVars(X0)
  
  targets <- names(sort(vars, decreasing=TRUE))[1:30] %>% as.character() %>% unique()
}

write_tsv2(tibble(PathwayBase = targets), "Selected_TargetPathways.tsv")
message("N target pathways: ", length(targets))

# ----------------------------
# 6) Build tidy long table for targets
# ----------------------------
row_map <- tibble(Row = rownames(mat2)) %>%
  mutate(
    PathwayBase = purrr::map_chr(Row, ~parse_strat(.x)$PathwayBase),
    TaxonRaw    = purrr::map_chr(Row, ~parse_strat(.x)$Taxon),
    Taxon       = clean_taxon(TaxonRaw)
  )

rows_keep <- row_map$PathwayBase %in% targets
mat_tgt <- mat2[rows_keep, , drop = FALSE]
row_map_tgt <- row_map[rows_keep, , drop = FALSE]

df_long <- mat_tgt %>%
  as.data.frame() %>%
  rownames_to_column("Row") %>%
  pivot_longer(-Row, names_to="SampleID", values_to="Value") %>%
  left_join(row_map_tgt, by="Row") %>%
  select(PathwayBase, Taxon, SampleID, Value) %>%
  left_join(meta2, by="SampleID")

write_tsv2(df_long, "STRATIFIED_TargetPathways_Long.tsv")
message("df_long rows: ", nrow(df_long))

# ----------------------------
# 7) Contributor tables
# ----------------------------
top_contrib_overall <- df_long %>%
  group_by(PathwayBase, Taxon) %>%
  summarise(Total = sum(Value, na.rm=TRUE), .groups="drop") %>%
  group_by(PathwayBase) %>%
  mutate(Prop = ifelse(sum(Total) > 0, Total / sum(Total), 0)) %>%
  arrange(PathwayBase, desc(Total)) %>%
  slice_head(n = 25) %>%
  ungroup()
write_tsv2(top_contrib_overall, "TABLE_TopContributors_overall_top25.tsv")

top_contrib_by_loc <- df_long %>%
  group_by(Location, PathwayBase, Taxon) %>%
  summarise(Total = sum(Value, na.rm=TRUE), .groups="drop") %>%
  group_by(Location, PathwayBase) %>%
  mutate(Prop = ifelse(sum(Total) > 0, Total / sum(Total), 0)) %>%
  arrange(Location, PathwayBase, desc(Total)) %>%
  slice_head(n = 25) %>%
  ungroup()
write_tsv2(top_contrib_by_loc, "TABLE_TopContributors_byLocation_top25.tsv")

# ----------------------------
# 8) Per-pathway plots
# ----------------------------
make_sample_stack_plot <- function(pwy, top_n_taxa = 12, min_prop = 0.01) {
  d <- df_long %>% filter(PathwayBase == pwy)
  if (nrow(d) == 0) return(invisible(NULL))
  
  d2 <- d %>%
    group_by(SampleID) %>%
    group_modify(~{
      x <- .x %>% select(Taxon, Value)
      collapse_other(x, top_n = top_n_taxa, min_prop = min_prop)
    }) %>%
    ungroup() %>%
    left_join(meta2, by="SampleID") %>%
    mutate(SampleID = factor(SampleID, levels = meta2 %>% arrange(Location, Inoculation, SampleID) %>% pull(SampleID)))
  
  tax_order <- d2 %>%
    group_by(Taxon) %>%
    summarise(Total = sum(Value, na.rm=TRUE), .groups="drop") %>%
    arrange(desc(Total)) %>% pull(Taxon)
  
  d2 <- d2 %>% mutate(Taxon = factor(Taxon, levels = tax_order))
  
  p <- ggplot(d2, aes(x = SampleID, y = Value, fill = Taxon)) +
    geom_col(width = 0.95) +
    facet_grid(~Location, scales="free_x", space="free_x") +
    theme_nm(10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(title = paste0("Taxon-resolved contributions — ", pwy, " (Rhizosphere)"),
         subtitle = "Stacked bars per sample (taxa collapsed into Other per sample)",
         x = NULL, y = "Pathway abundance (CPM)")
  
  outbase <- file.path(out_dir, "ByPathway_SampleStacks",
                       paste0("Fig_STRATIFIED_SampleStack_", make.names(substr(pwy,1,80))))
  save_gg(p, outbase, w=14, h=6)
}

make_group_mean_plot <- function(pwy, top_n_taxa = 12, min_prop = 0.01) {
  d <- df_long %>% filter(PathwayBase == pwy)
  if (nrow(d) == 0) return(invisible(NULL))
  
  d_mean <- d %>%
    group_by(Location, Inoculation, Taxon) %>%
    summarise(Value = mean(Value, na.rm=TRUE), .groups="drop")
  
  d2 <- d_mean %>%
    group_by(Location, Inoculation) %>%
    group_modify(~collapse_other(.x %>% select(Taxon, Value), top_n = top_n_taxa, min_prop = min_prop)) %>%
    ungroup()
  
  tax_order <- d2 %>%
    group_by(Taxon) %>%
    summarise(Total = sum(Value, na.rm=TRUE), .groups="drop") %>%
    arrange(desc(Total)) %>% pull(Taxon)
  
  d2 <- d2 %>% mutate(Taxon = factor(Taxon, levels = tax_order))
  
  p <- ggplot(d2, aes(x = Inoculation, y = Value, fill = Taxon)) +
    geom_col(width = 0.8) +
    facet_wrap(~Location, nrow = 1, scales="free_y") +
    theme_nm(11) +
    labs(title = paste0("Mean taxon contributions — ", pwy, " (Rhizosphere)"),
         subtitle = "Means per Location × Inoculation (taxa collapsed into Other per group)",
         x = NULL, y = "Mean pathway abundance (CPM)")
  
  outbase <- file.path(out_dir, "ByPathway_GroupMeans",
                       paste0("Fig_STRATIFIED_GroupMeans_", make.names(substr(pwy,1,80))))
  save_gg(p, outbase, w=14, h=5)
}

for (pwy in targets) {
  message("Plotting pathway: ", pwy)
  make_sample_stack_plot(pwy, top_n_taxa = 12, min_prop = 0.01)
  make_group_mean_plot(pwy, top_n_taxa = 12, min_prop = 0.01)
}

# ============================================================
# 9) HEATMAP: Taxa rows (RIGHT labels), Pathways cols (BOTTOM)
#     log10 proportion within pathway
# ============================================================
UNCL_LABEL <- "Unclassified"
TOP_PWY    <- 19
TOP_TAXA   <- 15
PSEUDO     <- 1e-6
LOG_RANGE  <- c(-4, 0)

pwy_global <- df_long %>%
  group_by(PathwayBase) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = TOP_PWY) %>%
  pull(PathwayBase)

tax_totals <- df_long %>%
  filter(PathwayBase %in% pwy_global) %>%
  group_by(Taxon) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total))

taxa_top <- tax_totals %>%
  filter(Taxon != UNCL_LABEL) %>%
  slice_head(n = TOP_TAXA) %>%
  pull(Taxon)

taxa_global <- c(taxa_top, UNCL_LABEL)

tax_pwy_hm <- df_long %>%
  filter(PathwayBase %in% pwy_global, Taxon %in% taxa_global) %>%
  group_by(PathwayBase, Taxon) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups="drop") %>%
  group_by(PathwayBase) %>%
  mutate(PathwayTotal = sum(Total, na.rm = TRUE),
         Prop = ifelse(PathwayTotal > 0, Total / PathwayTotal, 0)) %>%
  ungroup()

write_tsv2(tax_pwy_hm, "TABLE_Heatmap_TopTaxa_byTopPathways_Long.tsv")

hm <- tax_pwy_hm %>%
  select(PathwayBase, Taxon, Prop) %>%
  pivot_wider(names_from = Taxon, values_from = Prop, values_fill = 0) %>%
  column_to_rownames("PathwayBase") %>%
  as.matrix()

hm <- hm[, taxa_global, drop = FALSE]
hm_log <- log10(hm + PSEUDO)
hm_log[hm_log < LOG_RANGE[1]] <- LOG_RANGE[1]
hm_log[hm_log > LOG_RANGE[2]] <- LOG_RANGE[2]

hm_log2 <- t(hm_log)

rownames(hm_log2) <- gsub("_"," ", rownames(hm_log2))
colnames(hm_log2) <- wrap_text(drop_pwy_id(colnames(hm_log2)), width = 26)

cols <- colorRampPalette(c("#2C7BB6","#ABD9E9","#FFFFBF","#FDAE61","#D7191C"))(100)
bk   <- seq(LOG_RANGE[1], LOG_RANGE[2], length.out = length(cols) + 1)

ph <- pheatmap(
  hm_log2,
  color = cols,
  breaks = bk,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  fontsize_row = 9,
  fontsize_col = 7,
  angle_col = 45,
  cellwidth = 24,
  cellheight = 10,
  silent = TRUE
)

gt <- ph$gtable

nm_gt <- as.character(gt$layout$name)
cn_idx <- grep("col.*name", nm_gt, ignore.case = TRUE)
if (length(cn_idx) > 0) {
  cn_row <- unique(gt$layout$t[cn_idx])[1]
  gt$heights[cn_row] <- grid::unit(4.5, "cm")
}

rn_idx <- grep("row.*name", nm_gt, ignore.case = TRUE)
if (length(rn_idx) > 0) {
  rn_col <- unique(gt$layout$l[rn_idx])[1]
  gt$widths[rn_col] <- grid::unit(5.0, "cm")
}

save_gtable_safe <- function(gt, pdf_out, png_out, title,
                             extra_w_in = 1.5, extra_h_in = 0.2,
                             inset_mm = 6, title_in = 0.35, res = 600) {
  
  w0 <- grid::convertWidth(sum(gt$widths),  "in", valueOnly = TRUE)
  h0 <- grid::convertHeight(sum(gt$heights),"in", valueOnly = TRUE)
  w_in <- w0 + extra_w_in
  h_in <- h0 + extra_h_in + title_in
  
  draw_all <- function() {
    grid::grid.newpage()
    lay <- grid::grid.layout(nrow = 2, ncol = 1,
                             heights = grid::unit.c(grid::unit(title_in,"in"), grid::unit(1,"null")))
    grid::pushViewport(grid::viewport(layout = lay))
    
    grid::pushViewport(grid::viewport(layout.pos.row = 1))
    grid::grid.text(title, gp = grid::gpar(fontsize = 18, fontface = "bold"))
    grid::popViewport()
    
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    grid::pushViewport(grid::viewport(
      width  = grid::unit(1, "npc") - grid::unit(2*inset_mm, "mm"),
      height = grid::unit(1, "npc") - grid::unit(2*inset_mm, "mm"),
      clip="off"
    ))
    grid::grid.draw(gt)
    grid::popViewport(2)
    grid::popViewport()
  }
  
  grDevices::pdf(pdf_out, width = w_in, height = h_in, useDingbats = FALSE)
  draw_all(); grDevices::dev.off()
  
  grDevices::png(png_out, width = w_in, height = h_in, units="in", res=res)
  draw_all(); grDevices::dev.off()
}

hm_pdf <- file.path(out_dir, "Heatmaps", "Fig_STRATIFIED_Heatmap_TaxaRows_PathwaysCols_log10Prop.pdf")
hm_png <- file.path(out_dir, "Heatmaps", "Fig_STRATIFIED_Heatmap_TaxaRows_PathwaysCols_log10Prop.png")
save_gtable_safe(gt, hm_pdf, hm_png, title = "Major microbial contributors to pathways")
message("Heatmap written: ", hm_pdf)

# ============================================================
# 9b) DOTPLOTS (two options: keep both, reproducible)
# ============================================================
out_dotdir <- file.path(out_dir, "DotPlots")

make_dotplots <- function(df_long, TOP_PWY=19, TOP_TAXA=15, UNCL_LABEL="Unclassified",
                          PSEUDO=1e-6, option=c("within_pathway","classified_only")) {
  
  option <- match.arg(option)
  
  pwy_global <- df_long %>%
    group_by(PathwayBase) %>%
    summarise(Total = sum(Value, na.rm = TRUE), .groups="drop") %>%
    arrange(desc(Total)) %>%
    slice_head(n = TOP_PWY) %>%
    pull(PathwayBase)
  
  tax_totals <- df_long %>%
    filter(PathwayBase %in% pwy_global) %>%
    group_by(Taxon) %>%
    summarise(Total = sum(Value, na.rm = TRUE), .groups="drop") %>%
    arrange(desc(Total))
  
  taxa_top <- tax_totals %>%
    filter(Taxon != UNCL_LABEL) %>%
    slice_head(n = TOP_TAXA) %>%
    pull(Taxon)
  
  taxa_global <- c(taxa_top, UNCL_LABEL)
  
  tax_pwy <- df_long %>%
    filter(PathwayBase %in% pwy_global, Taxon %in% taxa_global) %>%
    group_by(PathwayBase, Taxon) %>%
    summarise(Total = sum(Value, na.rm = TRUE), .groups="drop") %>%
    group_by(PathwayBase) %>%
    mutate(
      PathwayTotal   = sum(Total, na.rm = TRUE),
      ClassifiedTotal = sum(Total[Taxon != UNCL_LABEL], na.rm = TRUE),
      Prop_within = ifelse(PathwayTotal > 0, Total / PathwayTotal, 0),
      Prop_classified = if_else(Taxon == UNCL_LABEL, NA_real_,
                                if_else(ClassifiedTotal > 0, Total / ClassifiedTotal, 0))
    ) %>%
    ungroup() %>%
    mutate(
      PathwayLabel = factor(drop_pwy_id(PathwayBase), levels = rev(unique(drop_pwy_id(pwy_global)))),
      TaxonLabel   = gsub("_"," ", Taxon),
      log10_within = log10(Prop_within + PSEUDO),
      log10_class  = log10(Prop_classified + PSEUDO)
    )
  
  tax_order <- tax_pwy %>%
    filter(Taxon != UNCL_LABEL) %>%
    group_by(TaxonLabel) %>%
    summarise(Total = sum(Total, na.rm=TRUE), .groups="drop") %>%
    arrange(desc(Total)) %>% pull(TaxonLabel)
  
  tax_pwy <- tax_pwy %>%
    mutate(TaxonLabel = factor(TaxonLabel, levels = c(tax_order, UNCL_LABEL)))
  
  if (option == "within_pathway") {
    plot_df <- tax_pwy %>% filter(Prop_within > 0)
    size_var <- "Prop_within"
    col_var  <- "log10_within"
    subtitle <- "Dot size = proportion within pathway (includes Unclassified); colour = log10(proportion)"
    out_tag  <- "WITH_Unclassified"
  } else {
    plot_df <- tax_pwy %>% filter(!is.na(Prop_classified), Prop_classified > 0)
    size_var <- "Prop_classified"
    col_var  <- "log10_class"
    subtitle <- "Dot size = proportion among classified taxa; colour = log10(proportion)"
    out_tag  <- "CLASSIFIED_only"
  }
  
  nm_cols <- c("#1b3a8a", "#3b7ddd", "#bfe3ff", "#f7f7f7", "#ffd08a", "#f07c2b")
  nm_vals <- scales::rescale(c(-4,-3,-2,-1,-0.5,0))
  
  p <- ggplot(plot_df, aes(x = TaxonLabel, y = PathwayLabel)) +
    geom_point(aes(size = .data[[size_var]], colour = .data[[col_var]]),
               alpha = 0.95, stroke = 0.2) +
    scale_size_continuous(range = c(0.6, 11),
                          labels = scales::percent_format(accuracy = 1),
                          name = "Proportion") +
    scale_colour_gradientn(colours = nm_cols, values = nm_vals,
                           limits = c(-4,0), name = "log10(prop)") +
    labs(title = "Major microbial contributors to pathways",
         subtitle = subtitle, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_line(color="grey92", linewidth=0.3),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(face="bold", size=14),
      legend.title = element_text(face="bold")
    )
  
  pdf_out <- file.path(out_dotdir, paste0("Fig_DOT_TaxaByPathway_", out_tag, ".pdf"))
  png_out <- file.path(out_dotdir, paste0("Fig_DOT_TaxaByPathway_", out_tag, ".png"))
  ggsave(pdf_out, p, width = 14, height = 8, useDingbats = FALSE)
  ggsave(png_out, p, width = 14, height = 8, dpi = 600)
  
  write_tsv2(tax_pwy, paste0("TABLE_DotPlot_TaxaByPathway_", out_tag, ".tsv"))
  
  message("Dotplot written: ", pdf_out)
  invisible(p)
}

make_dotplots(df_long, option = "within_pathway")
make_dotplots(df_long, option = "classified_only")

# ----------------------------
# 10) Output manifest (MD5 of all outputs)
# ----------------------------
outs <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
manifest <- tibble(
  file  = outs,
  bytes = file.info(outs)$size,
  md5   = vapply(outs, function(fp) as.character(tools::md5sum(fp)), FUN.VALUE = character(1))
)
readr::write_tsv(manifest, file.path(out_dir, "RunInfo", paste0("MANIFEST_outputs_", run_id, ".tsv")))

# ============================================================
# 11) NatCom Figure 5 panel: stratified pathway contributors
#     Taxon-resolved contributors to pathway responses
# ============================================================

COMPARTMENT <- "Rhizosphere"
PANEL_TAG <- "B"
NATCOM_DIR <- STAGING_DIR
dir.create(NATCOM_DIR, recursive = TRUE, showWarnings = FALSE)

TOP_PWY_NATCOM  <- 12
TOP_TAXA_NATCOM <- 12
UNCL_LABEL <- "Unclassified"
PSEUDO <- 1e-6

# ----------------------------
# Helper: italicise scientific names
# ----------------------------
italic_taxa <- function(x) {
  
  sapply(x, function(z) {
    
    # Genus sp. Strain
    if (grepl(" sp ", z)) {
      
      genus  <- sub(" sp .*", "", z)
      strain <- sub(".* sp ", "", z)
      
      return(
        paste0(
          "italic('", genus, "')~'sp.'~'",
          strain,
          "'"
        )
      )
    }
    
    # Genus species
    words <- strsplit(z, " ")[[1]]
    
    if (length(words) >= 2) {
      
      return(
        paste0(
          "italic('",
          paste(words[1:2], collapse = " "),
          "')"
        )
      )
    }
    
    # fallback
    paste0("'", z, "'")
  })
}

# ----------------------------
# Helper: clean pathway names
# ----------------------------
clean_pathway_natcom <- function(x) {
  x <- drop_pwy_id(x)
  x <- stringr::str_replace_all(x, "superpathway of ", "")
  x <- stringr::str_replace_all(x, "&alpha;", "α")
  x <- stringr::str_replace_all(x, "&beta;", "β")
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_squish(x)
  stringr::str_wrap(x, width = 24)
}

# ----------------------------
# Helper: sentence-case pathway names
# ----------------------------
sentence_case_pathway <- function(x) {
  
  x <- stringr::str_squish(x)
  
  # Capitalize first character only
  x <- paste0(
    toupper(substr(x, 1, 1)),
    substr(x, 2, nchar(x))
  )
  
  # Restore important abbreviations
  x <- x %>%
    stringr::str_replace_all("\\bNadph\\b", "NADPH") %>%
    stringr::str_replace_all("\\bTca\\b", "TCA") %>%
    stringr::str_replace_all("\\bIva\\b", "IVA") %>%
    stringr::str_replace_all("\\bAtp\\b", "ATP") %>%
    stringr::str_replace_all("\\bNad\\b", "NAD") %>%
    stringr::str_replace_all("\\bNadh\\b", "NADH") %>%
    stringr::str_replace_all("\\bUdp\\b", "UDP") %>%
    stringr::str_replace_all("\\bDtdp\\b", "dTDP")
  
  x
}

# ----------------------------
# Select major pathways
# ----------------------------
pwy_natcom <- df_long %>%
  group_by(PathwayBase) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = TOP_PWY_NATCOM) %>%
  pull(PathwayBase)

# ----------------------------
# Select major classified taxa
# ----------------------------
taxa_natcom <- df_long %>%
  filter(
    PathwayBase %in% pwy_natcom,
    Taxon != UNCL_LABEL
  ) %>%
  group_by(Taxon) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = TOP_TAXA_NATCOM) %>%
  pull(Taxon)

# ----------------------------
# Prepare plotting dataframe
# ----------------------------
plot_df_natcom <- df_long %>%
  filter(
    PathwayBase %in% pwy_natcom,
    Taxon %in% taxa_natcom
  ) %>%
  group_by(PathwayBase, Taxon) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  group_by(PathwayBase) %>%
  mutate(
    PathwayTotal = sum(Total, na.rm = TRUE),
    Prop = ifelse(PathwayTotal > 0, Total / PathwayTotal, 0),
    log10Prop = log10(Prop + PSEUDO)
  ) %>%
  ungroup() %>%
  mutate(
    PathwayLabel = clean_pathway_natcom(PathwayBase),
    PathwayLabel = sentence_case_pathway(PathwayLabel),
    TaxonLabel = stringr::str_replace_all(Taxon, "_", " ")
  )

# Shorten and wrap pathway labels safely
plot_df_natcom$PathwayLabel <- stringr::str_replace(
  plot_df_natcom$PathwayLabel,
  "dTDP-α-D-ravidosamine and dTDP-4-acetyl-α-D-ravidosamine biosynthesis",
  "dTDP-ravidosamine biosynthesis"
)

plot_df_natcom$PathwayLabel <- stringr::str_wrap(
  plot_df_natcom$PathwayLabel,
  width = 24
)

# Keep pathway order by total abundance
path_order <- plot_df_natcom %>%
  group_by(PathwayLabel) %>%
  summarise(Total = sum(Total, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(PathwayLabel)

# Keep taxon order by total contribution
tax_order <- plot_df_natcom %>%
  group_by(TaxonLabel) %>%
  summarise(Total = sum(Total, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(TaxonLabel)

plot_df_natcom <- plot_df_natcom %>%
  mutate(
    PathwayLabel = factor(PathwayLabel, levels = rev(path_order)),
    TaxonLabel = factor(TaxonLabel, levels = rev(tax_order))
  )

# ----------------------------
# Create Figure 5 panel
# ----------------------------
p_fig5_panel <- ggplot(
  plot_df_natcom,
  aes(x = PathwayLabel, y = TaxonLabel)
) +
  geom_point(
    aes(size = Prop, colour = log10Prop),
    alpha = 0.95
  ) +
  scale_y_discrete(
    labels = function(x)
      parse(text = italic_taxa(as.character(x)))
  ) +
  scale_size_continuous(
    name = "Contribution",
    range = c(1.2, 7),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_colour_gradientn(
    colours = c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C"),
    limits = c(-4, 0),
    name = "log10\nproportion"
  ) +
  labs(
    title = COMPARTMENT,
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black"),
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

# ----------------------------
# Save panel for later Figure 5 merge
# ----------------------------
ggsave(
  file.path(
    NATCOM_DIR,
    paste0("Fig5", PANEL_TAG, "_", COMPARTMENT, "_Stratified_Pathway_Contributors.pdf")
  ),
  p_fig5_panel,
  width = 140,
  height = 120,
  units = "mm",
  useDingbats = FALSE
)

ggsave(
  file.path(
    NATCOM_DIR,
    paste0("Fig5", PANEL_TAG, "_", COMPARTMENT, "_Stratified_Pathway_Contributors.png")
  ),
  p_fig5_panel,
  width = 140,
  height = 120,
  units = "mm",
  dpi = 600
)

saveRDS(
  p_fig5_panel,
  file.path(
    NATCOM_DIR,
    paste0("p_Fig5", PANEL_TAG, "_", COMPARTMENT, "_Stratified_Pathway_Contributors.rds")
  )
)

write_tsv2(
  plot_df_natcom,
  paste0("TABLE_Fig5", PANEL_TAG, "_", COMPARTMENT, "_Stratified_Pathway_Contributors.tsv")
)




# ============================================================
# 12) Supplementary Figure S5 panel: same style as Figure 5, WITH Unclassified
#     Top classified taxa plus Unclassified; proportions use all taxa as denominator
# ============================================================

SUPP_PANEL_TAG <- "A"
TOP_PWY_SUPP  <- TOP_PWY_NATCOM
TOP_TAXA_SUPP <- TOP_TAXA_NATCOM

# Select the same major pathways as the NatCom contributor panel
pwy_supp <- pwy_natcom

# Select top classified taxa, then force inclusion of Unclassified
taxa_supp_classified <- df_long %>%
  filter(
    PathwayBase %in% pwy_supp,
    Taxon != UNCL_LABEL
  ) %>%
  group_by(Taxon) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = TOP_TAXA_SUPP) %>%
  pull(Taxon)

taxa_supp <- c(taxa_supp_classified, UNCL_LABEL)

# Denominator: total pathway abundance across all taxa, including Unclassified
pwy_total_all <- df_long %>%
  filter(PathwayBase %in% pwy_supp) %>%
  group_by(PathwayBase) %>%
  summarise(PathwayTotalAll = sum(Value, na.rm = TRUE), .groups = "drop")

plot_df_supp <- df_long %>%
  filter(
    PathwayBase %in% pwy_supp,
    Taxon %in% taxa_supp
  ) %>%
  group_by(PathwayBase, Taxon) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  left_join(pwy_total_all, by = "PathwayBase") %>%
  mutate(
    Prop = ifelse(PathwayTotalAll > 0, Total / PathwayTotalAll, 0),
    log10Prop = log10(Prop + PSEUDO)
  ) %>%
  ungroup() %>%
  mutate(
    PathwayLabel = clean_pathway_natcom(PathwayBase),
    PathwayLabel = sentence_case_pathway(PathwayLabel),
    TaxonLabel = stringr::str_replace_all(Taxon, "_", " ")
  )

# Shorten and wrap pathway labels safely
plot_df_supp$PathwayLabel <- stringr::str_replace(
  plot_df_supp$PathwayLabel,
  "dTDP-α-D-ravidosamine and dTDP-4-acetyl-α-D-ravidosamine biosynthesis",
  "dTDP-ravidosamine biosynthesis"
)

plot_df_supp$PathwayLabel <- stringr::str_wrap(
  plot_df_supp$PathwayLabel,
  width = 24
)

# Keep pathway order consistent with the main contributor panel
path_order_supp <- plot_df_supp %>%
  group_by(PathwayLabel) %>%
  summarise(Total = sum(Total, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(PathwayLabel)

# Keep top classified taxa ordered by contribution, with Unclassified placed last
tax_order_supp_classified <- plot_df_supp %>%
  filter(TaxonLabel != UNCL_LABEL) %>%
  group_by(TaxonLabel) %>%
  summarise(Total = sum(Total, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(TaxonLabel)

tax_order_supp <- c(tax_order_supp_classified, UNCL_LABEL)

plot_df_supp <- plot_df_supp %>%
  mutate(
    PathwayLabel = factor(PathwayLabel, levels = rev(path_order_supp)),
    TaxonLabel = factor(TaxonLabel, levels = rev(tax_order_supp))
  )

p_figS5_panel <- ggplot(
  plot_df_supp,
  aes(x = PathwayLabel, y = TaxonLabel)
) +
  geom_point(
    aes(size = Prop, colour = log10Prop),
    alpha = 0.95
  ) +
  scale_y_discrete(
    labels = function(x)
      parse(text = italic_taxa(as.character(x)))
  ) +
  scale_size_continuous(
    name = "Contribution",
    range = c(1.2, 7),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_colour_gradientn(
    colours = c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C"),
    limits = c(-4, 0),
    name = "log10\nproportion"
  ) +
  labs(
    title = COMPARTMENT,
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black"),
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(
  file.path(
    NATCOM_DIR,
    paste0("FigS5", SUPP_PANEL_TAG, "_", COMPARTMENT, "_WITH_Unclassified_panel.pdf")
  ),
  p_figS5_panel,
  width = 140,
  height = 120,
  units = "mm",
  useDingbats = FALSE
)

ggsave(
  file.path(
    NATCOM_DIR,
    paste0("FigS5", SUPP_PANEL_TAG, "_", COMPARTMENT, "_WITH_Unclassified_panel.png")
  ),
  p_figS5_panel,
  width = 140,
  height = 120,
  units = "mm",
  dpi = 600
)

saveRDS(
  p_figS5_panel,
  file.path(
    NATCOM_DIR,
    paste0("p_FigS5", SUPP_PANEL_TAG, "_", COMPARTMENT, "_WITH_Unclassified.rds")
  )
)

write_tsv2(
  plot_df_supp,
  paste0("TABLE_FigS5", SUPP_PANEL_TAG, "_", COMPARTMENT, "_WITH_Unclassified.tsv")
)


message("DONE. Outputs saved: ", out_dir)
message("Run log: ", log_fp)


})
})


# ============================================================
# Root stratified contributors (Figure 5 and S5 panels)
# ============================================================
run_section("Root stratified contributors (Figure 5 and S5 panels)", {
local({
# ============================================================
# ROOT shotgun metagenomes (n=48) — HUMAnN4 STRATIFIED outputs
# Aim: "Who does what?" (taxon-resolved contributors to pathways)
#
# Inputs:
#   - Root_pathabundance_STRATIFIED.tsv.gz
#   - Root_metadata.tsv.gz
#
# Optional (recommended; from UNSTRATIFIED pipeline):
#   - TABLE_HeadlinePathways_sigIn2plusLocations.tsv
#   - TABLE_PathwayConsistency_acrossLocations.tsv
#   - ALDEx2_Inoculation_byLocation.tsv
#
# Outputs (PDF+PNG + TSV):
#   - Per-pathway taxon contribution stacked bars (per-sample)
#   - Per-pathway taxon contribution stacked bars (group means)
#   - Who-does-what Heatmap: Taxa(rows) × Pathway(cols), log10(prop within pathway)
#   - Dotplots (two versions): WITH Unclassified + CLASSIFIED-only
#   - Tables: top contributors overall and by Location
#   - RunInfo: run log, session info, MD5 inputs, outputs manifest
#
# UPDATED: standardised Eyrewell_Forest spelling
# ============================================================

# ----------------------------
# 0) Packages
# ----------------------------
pkgs <- c("tidyverse","pheatmap","matrixStats","digest","sessioninfo","grid")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

options(error = NULL)
set.seed(1)

select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
summarise <- dplyr::summarise

# ----------------------------
# 1) Paths + output folders
# ----------------------------
library(here)
project_dir <- repo_root
data_dir <- file.path(project_dir, "function", "data", "root")
out_dir  <- file.path(project_dir, "function", "results", "root", "stratified")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

strat_fp <- file.path(data_dir, "Root_pathabundance_STRATIFIED.tsv.gz")
meta_fp  <- file.path(data_dir, "Root_metadata.tsv.gz")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

subdirs <- c("ByPathway_SampleStacks","ByPathway_GroupMeans","Heatmaps","DotPlots","Tables","RunInfo")
for (d in subdirs) dir.create(file.path(out_dir, d), showWarnings = FALSE, recursive = TRUE)

headline_dir <- file.path(project_dir, "function", "results", "root", "unstratified", "Tables")
headline_fp1 <- file.path(headline_dir, "TABLE_HeadlinePathways_sigIn2plusLocations.tsv")
headline_fp2 <- file.path(headline_dir, "Root_TABLE_PathwayConsistency_acrossLocations.tsv")
aldex_fp     <- file.path(headline_dir, "ALDEx2_Inoculation_byLocation.tsv")

# ----------------------------
# 1b) Reproducibility: run log + session info + MD5 inputs
# ----------------------------
assert_exists <- function(fp, label = fp) {
  if (!file.exists(fp)) stop("Missing file: ", label, "\nPath: ", fp, call. = FALSE)
}

md5_or_na <- function(fp) if (file.exists(fp)) as.character(tools::md5sum(fp)) else NA_character_

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_fp <- file.path(out_dir, "RunInfo", paste0("RUNLOG_", run_id, ".txt"))

zz <- file(log_fp, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")

.onExit_sinks <- function() {
  try(sink(type="message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(zz), silent = TRUE)
}
on.exit(.onExit_sinks(), add = TRUE)

message("=== HUMAnN4 STRATIFIED WhoDoesWhat (ROOT) ===")
message("Run ID: ", run_id)
message("Timestamp: ", Sys.time())
message("R version: ", R.version.string)
message("Platform: ", R.version$platform)

message("project_dir: ", project_dir)
message("data_dir:    ", data_dir)
message("strat_fp:   ", strat_fp, "  MD5=", md5_or_na(strat_fp))
message("meta_fp:    ", meta_fp,  "  MD5=", md5_or_na(meta_fp))
message("headline1:  ", headline_fp1, "  MD5=", md5_or_na(headline_fp1))
message("headline2:  ", headline_fp2, "  MD5=", md5_or_na(headline_fp2))
message("aldex_fp:   ", aldex_fp, "  MD5=", md5_or_na(aldex_fp))
message("out_dir:    ", out_dir)

si_fp <- file.path(out_dir, "RunInfo", paste0("SESSIONINFO_", run_id, ".txt"))
writeLines(capture.output(sessioninfo::session_info()), si_fp)
message("Session info saved: ", si_fp)

# ----------------------------
# 2) Helpers
# ----------------------------
theme_nm <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "white"),
      legend.title = element_text(face = "bold")
    )
}

safe_fname <- function(x, max_len = 60) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (nchar(x) > max_len) x <- substr(x, 1, max_len)
  x
}

pathway_id <- function(pwy) {
  id <- sub(":.*$", "", pwy)
  id <- trimws(id)
  safe_fname(id, max_len = 25)
}

save_gg <- function(p, outpath_base, w = 12, h = 7, dpi = 600) {
  dir.create(dirname(outpath_base), showWarnings = FALSE, recursive = TRUE)
  
  pdf_fp <- paste0(outpath_base, ".pdf")
  png_fp <- paste0(outpath_base, ".png")
  
  path_too_long <- function(fp) {
    nchar(normalizePath(fp, winslash="\\", mustWork = FALSE)) > 235
  }
  
  if (path_too_long(pdf_fp) || path_too_long(png_fp)) {
    base_dir <- dirname(outpath_base)
    base_nm  <- basename(outpath_base)
    
    base_nm2 <- safe_fname(base_nm, max_len = 40)
    pdf_fp2 <- file.path(base_dir, paste0(base_nm2, ".pdf"))
    png_fp2 <- file.path(base_dir, paste0(base_nm2, ".png"))
    
    if (path_too_long(pdf_fp2)) {
      tok <- tail(strsplit(base_nm2, "_")[[1]], 1)
      tok <- safe_fname(tok, max_len = 30)
      pdf_fp2 <- file.path(base_dir, paste0(tok, ".pdf"))
      png_fp2 <- file.path(base_dir, paste0(tok, ".png"))
    }
    
    pdf_fp <- pdf_fp2
    png_fp <- png_fp2
  }
  
  ggplot2::ggsave(pdf_fp, p, width = w, height = h, useDingbats = FALSE)
  ggplot2::ggsave(png_fp, p, width = w, height = h, dpi = dpi)
}

write_tsv2 <- function(df, filename) {
  readr::write_tsv(df, file.path(out_dir, "Tables", filename))
}

parse_strat <- function(x) {
  parts <- strsplit(x, "\\|")[[1]]
  if (length(parts) == 1) return(list(PathwayBase = parts[1], Taxon = "UNSTRATIFIED"))
  list(PathwayBase = parts[1], Taxon = paste(parts[-1], collapse="|"))
}

clean_taxon <- function(taxon) {
  taxon %>%
    stringr::str_replace("^s__", "") %>%
    stringr::str_replace("\\.t__SGB\\d+$", "") %>%
    stringr::str_replace("^unclassified$", "Unclassified")
}

collapse_other <- function(df, top_n = 12, min_prop = 0.01) {
  tot <- sum(df$Value, na.rm = TRUE)
  df <- df %>% mutate(Prop = ifelse(tot > 0, Value / tot, 0))
  keep <- df %>% arrange(desc(Value)) %>% slice_head(n = top_n) %>% pull(Taxon)
  keep <- union(keep, df %>% filter(Prop >= min_prop) %>% pull(Taxon))
  
  df %>%
    mutate(Taxon2 = ifelse(Taxon %in% keep, Taxon, "Other")) %>%
    group_by(Taxon2) %>%
    summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop") %>%
    rename(Taxon = Taxon2)
}

drop_pwy_id <- function(x) sub("^.*?[:：∶]\\s*", "", as.character(x), perl = TRUE)

wrap_text <- function(x, width = 26) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

# ----------------------------
# 3) Import metadata and harmonise SampleID
# ----------------------------
assert_exists(meta_fp, "meta_fp")
meta <- readr::read_tsv(meta_fp, show_col_types = FALSE) %>%
  mutate(
    SampleID    = as.character(SampleID),
    SampleID    = trimws(SampleID),
    SampleID    = gsub("^R_", "", SampleID),
    Location    = trimws(as.character(Location)),
    Inoculation = trimws(as.character(Inoculation))
  ) %>%
  mutate(
    # safety standardisation
    Location = case_when(
      Location %in% c("Eyerell_Forest", "Eyerell Forest", "Eyrewell Forest") ~ "Eyrewell_Forest",
      TRUE ~ Location
    ),
    Inoculation = case_when(
      str_to_lower(Inoculation) == "control" ~ "Control",
      str_to_lower(Inoculation) %in% c("panch", "trichoderma")   ~ "Trichoderma",
      TRUE ~ Inoculation
    )
  ) %>%
  mutate(
    Location    = factor(Location),
    Inoculation = factor(Inoculation, levels = c("Control", "Trichoderma"))
  )

stopifnot(all(c("SampleID","Location","Inoculation") %in% names(meta)))
if (anyDuplicated(meta$SampleID)) stop("Duplicate SampleID in metadata.")
message("Metadata rows: ", nrow(meta))

# ----------------------------
# 4) Import stratified table and align samples
# ----------------------------
assert_exists(strat_fp, "strat_fp")
strat_tbl <- readr::read_tsv(strat_fp, show_col_types = FALSE, comment = "")
names(strat_tbl)[1] <- "Pathway"
stopifnot("Pathway" %in% names(strat_tbl))

mat <- as.data.frame(strat_tbl)
rownames(mat) <- mat$Pathway
mat$Pathway <- NULL

mat[] <- lapply(mat, function(x) suppressWarnings(as.numeric(x)))
mat[is.na(mat)] <- 0

colnames(mat) <- gsub("^R_", "", colnames(mat))
colnames(mat) <- trimws(colnames(mat))

meta2 <- meta %>% filter(SampleID %in% colnames(mat))
if (nrow(meta2) == 0) stop("No metadata SampleID matches stratified table columns. Check prefixes.")

meta2 <- meta2 %>% arrange(Location, Inoculation, SampleID)
mat2 <- mat[, meta2$SampleID, drop = FALSE]
stopifnot(identical(colnames(mat2), meta2$SampleID))

message("Samples aligned: ", ncol(mat2))
message("Rows (stratified): ", nrow(mat2))

drop_pat <- "^(UNMAPPED|UNINTEGRATED)"
keep_rows <- !grepl(drop_pat, rownames(mat2))
mat2 <- mat2[keep_rows, , drop = FALSE]
message("Rows after dropping unmapped/unintegrated: ", nrow(mat2))

# ----------------------------
# 5) Define target pathways
# ----------------------------
targets <- character()

if (file.exists(headline_fp1)) {
  message("Using headline pathways from: ", headline_fp1)
  targets <- readr::read_tsv(headline_fp1, show_col_types = FALSE) %>%
    pull(Pathway) %>% unique()
  
} else if (file.exists(headline_fp2)) {
  message("Using consistency table from: ", headline_fp2)
  tmp <- readr::read_tsv(headline_fp2, show_col_types = FALSE)
  
  if (all(c("Pathway","n_locations_sig","median_effect") %in% names(tmp))) {
    targets <- tmp %>%
      arrange(desc(n_locations_sig), desc(abs(median_effect))) %>%
      slice_head(n = 30) %>%
      pull(Pathway) %>% unique()
  } else if ("Pathway" %in% names(tmp)) {
    targets <- tmp %>% pull(Pathway) %>% unique()
    targets <- targets[1:min(30, length(targets))]
  }
  
} else if (file.exists(aldex_fp)) {
  message("Using ALDEx2 (BH<0.05) from: ", aldex_fp)
  targets <- readr::read_tsv(aldex_fp, show_col_types = FALSE) %>%
    filter(we.eBH < 0.05) %>%
    arrange(we.eBH, desc(abs(effect))) %>%
    slice_head(n = 30) %>%
    pull(Pathway) %>% unique()
  
} else {
  message("No unstratified selection found; selecting targets by variance (summed over taxa).")
  
  rowinfo <- tibble(Row = rownames(mat2)) %>%
    mutate(PathwayBase = purrr::map_chr(Row, ~parse_strat(.x)$PathwayBase))
  
  long0 <- mat2 %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Row") %>%
    pivot_longer(-Row, names_to="SampleID", values_to="Value") %>%
    left_join(rowinfo, by="Row") %>%
    group_by(PathwayBase, SampleID) %>%
    summarise(Value = sum(Value, na.rm=TRUE), .groups="drop")
  
  wide0 <- long0 %>% pivot_wider(names_from=SampleID, values_from=Value, values_fill=0)
  X0 <- wide0 %>% column_to_rownames("PathwayBase") %>% as.matrix()
  
  vars <- matrixStats::rowVars(X0)
  targets <- names(sort(vars, decreasing=TRUE))[1:30]
}

targets <- unique(targets[!is.na(targets)])
message("N target pathways: ", length(targets))
write_tsv2(tibble(PathwayBase = targets), "Selected_TargetPathways.tsv")

# ----------------------------
# 6) Tidy long format for targets
# ----------------------------
row_map <- tibble(Row = rownames(mat2)) %>%
  mutate(
    PathwayBase = purrr::map_chr(Row, ~parse_strat(.x)$PathwayBase),
    TaxonRaw    = purrr::map_chr(Row, ~parse_strat(.x)$Taxon),
    Taxon       = clean_taxon(TaxonRaw)
  )

rows_keep <- row_map$PathwayBase %in% targets
mat_tgt <- mat2[rows_keep, , drop = FALSE]
row_map_tgt <- row_map[rows_keep, , drop = FALSE]

df_long <- mat_tgt %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Row") %>%
  tidyr::pivot_longer(-Row, names_to="SampleID", values_to="Value") %>%
  left_join(row_map_tgt, by="Row") %>%
  select(PathwayBase, Taxon, SampleID, Value) %>%
  left_join(meta2, by="SampleID")

write_tsv2(df_long, "STRATIFIED_TargetPathways_Long.tsv")
message("df_long rows: ", nrow(df_long))

# ----------------------------
# 7) Contributor tables (overall + by Location)
# ----------------------------
top_contrib_overall <- df_long %>%
  group_by(PathwayBase, Taxon) %>%
  summarise(Total = sum(Value, na.rm=TRUE), .groups="drop") %>%
  group_by(PathwayBase) %>%
  mutate(Prop = ifelse(sum(Total) > 0, Total / sum(Total), 0)) %>%
  arrange(PathwayBase, desc(Total)) %>%
  slice_head(n = 25) %>%
  ungroup()
write_tsv2(top_contrib_overall, "TABLE_TopContributors_overall_top25.tsv")

top_contrib_by_loc <- df_long %>%
  group_by(Location, PathwayBase, Taxon) %>%
  summarise(Total = sum(Value, na.rm=TRUE), .groups="drop") %>%
  group_by(Location, PathwayBase) %>%
  mutate(Prop = ifelse(sum(Total) > 0, Total / sum(Total), 0)) %>%
  arrange(Location, PathwayBase, desc(Total)) %>%
  slice_head(n = 25) %>%
  ungroup()
write_tsv2(top_contrib_by_loc, "TABLE_TopContributors_byLocation_top25.tsv")

# ----------------------------
# 8) Per-pathway plots (sample stacks + group means)
# ----------------------------
make_sample_stack_plot <- function(pwy, top_n_taxa = 12, min_prop = 0.01) {
  d <- df_long %>% filter(PathwayBase == pwy)
  if (nrow(d) == 0) return(invisible(NULL))
  
  d2 <- d %>%
    group_by(SampleID) %>%
    group_modify(~{
      x <- .x %>% select(Taxon, Value)
      collapse_other(x, top_n = top_n_taxa, min_prop = min_prop)
    }) %>%
    ungroup() %>%
    left_join(meta2, by="SampleID") %>%
    mutate(SampleID = factor(SampleID, levels = meta2 %>% arrange(Location, Inoculation, SampleID) %>% pull(SampleID)))
  
  tax_order <- d2 %>%
    group_by(Taxon) %>%
    summarise(Total = sum(Value, na.rm=TRUE), .groups="drop") %>%
    arrange(desc(Total)) %>%
    pull(Taxon)
  d2 <- d2 %>% mutate(Taxon = factor(Taxon, levels = tax_order))
  
  p <- ggplot(d2, aes(x = SampleID, y = Value, fill = Taxon)) +
    geom_col(width = 0.95) +
    facet_grid(~Location, scales="free_x", space="free_x") +
    theme_nm(10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "right") +
    labs(title = paste0("Taxon-resolved contributions — ", pwy, " (Root)"),
         subtitle = "Stacked bars per sample; taxa collapsed to Other per sample",
         x = NULL, y = "Pathway abundance (CPM)")
  
  id   <- pathway_id(pwy)
  slug <- safe_fname(sub("^.*?:\\s*", "", pwy), max_len = 35)
  outbase <- file.path(out_dir, "ByPathway_SampleStacks", paste0("Fig_Strat_Sample_", id, "_", slug))
  save_gg(p, outbase, w = 14, h = 6)
  invisible(TRUE)
}

make_group_mean_plot <- function(pwy, top_n_taxa = 12, min_prop = 0.01) {
  d <- df_long %>% filter(PathwayBase == pwy)
  if (nrow(d) == 0) return(invisible(NULL))
  
  d_mean <- d %>%
    group_by(Location, Inoculation, Taxon) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  d2 <- d_mean %>%
    group_by(Location, Inoculation) %>%
    group_modify(~collapse_other(.x %>% select(Taxon, Value), top_n = top_n_taxa, min_prop = min_prop)) %>%
    ungroup()
  
  tax_order <- d2 %>%
    group_by(Taxon) %>%
    summarise(Total = sum(Value, na.rm=TRUE), .groups="drop") %>%
    arrange(desc(Total)) %>%
    pull(Taxon)
  d2 <- d2 %>% mutate(Taxon = factor(Taxon, levels = tax_order))
  
  p <- ggplot(d2, aes(x = Inoculation, y = Value, fill = Taxon)) +
    geom_col(width = 0.8) +
    facet_wrap(~Location, nrow = 1, scales = "free_y") +
    theme_nm(11) +
    labs(title = paste0("Mean taxon contributions — ", pwy, " (Root)"),
         subtitle = "Means per Location × Inoculation (taxa collapsed to Other per group)",
         x = NULL, y = "Mean pathway abundance (CPM)")
  
  id   <- pathway_id(pwy)
  slug <- safe_fname(sub("^.*?:\\s*", "", pwy), max_len = 35)
  outbase <- file.path(out_dir, "ByPathway_GroupMeans", paste0("Fig_Strat_GroupMean_", id, "_", slug))
  save_gg(p, outbase, w = 14, h = 5)
  invisible(TRUE)
}

for (pwy in targets) {
  message("Plotting pathway: ", pwy)
  make_sample_stack_plot(pwy, top_n_taxa = 12, min_prop = 0.01)
  make_group_mean_plot(pwy, top_n_taxa = 12, min_prop = 0.01)
}

# ============================================================
# 9) WHO-DOES-WHAT HEATMAP
#    Taxa (rows) × Pathways (cols)
#    Values = log10(proportion within pathway)
# ============================================================
UNCL_LABEL <- "Unclassified"
TOP_PWY    <- 13
TOP_TAXA   <- 15
PSEUDO     <- 1e-6
LOG_RANGE  <- c(-4, 0)

pwy_global <- df_long %>%
  group_by(PathwayBase) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups="drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = TOP_PWY) %>%
  pull(PathwayBase)

tax_totals <- df_long %>%
  filter(PathwayBase %in% pwy_global) %>%
  group_by(Taxon) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups="drop") %>%
  arrange(desc(Total))

taxa_top <- tax_totals %>%
  filter(Taxon != UNCL_LABEL) %>%
  slice_head(n = TOP_TAXA) %>%
  pull(Taxon)

taxa_global <- c(taxa_top, UNCL_LABEL)

tax_pwy_hm <- df_long %>%
  filter(PathwayBase %in% pwy_global, Taxon %in% taxa_global) %>%
  group_by(PathwayBase, Taxon) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups="drop") %>%
  group_by(PathwayBase) %>%
  mutate(PathwayTotal = sum(Total, na.rm = TRUE),
         Prop = ifelse(PathwayTotal > 0, Total / PathwayTotal, 0)) %>%
  ungroup()

write_tsv2(tax_pwy_hm, "TABLE_Heatmap_TaxaByPathway_Long.tsv")

hm <- tax_pwy_hm %>%
  select(PathwayBase, Taxon, Prop) %>%
  pivot_wider(names_from = Taxon, values_from = Prop, values_fill = 0) %>%
  column_to_rownames("PathwayBase") %>%
  as.matrix()

hm <- hm[, taxa_global, drop = FALSE]
hm_log <- log10(hm + PSEUDO)
hm_log[hm_log < LOG_RANGE[1]] <- LOG_RANGE[1]
hm_log[hm_log > LOG_RANGE[2]] <- LOG_RANGE[2]

hm_log2 <- t(hm_log)
rownames(hm_log2) <- gsub("_"," ", rownames(hm_log2))
colnames(hm_log2) <- wrap_text(drop_pwy_id(colnames(hm_log2)), width = 26)

cols <- colorRampPalette(c("#2C7BB6","#ABD9E9","#FFFFBF","#FDAE61","#D7191C"))(100)
bk   <- seq(LOG_RANGE[1], LOG_RANGE[2], length.out = length(cols) + 1)

ph <- pheatmap(
  hm_log2,
  color = cols,
  breaks = bk,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  fontsize_row = 9,
  fontsize_col = 7,
  angle_col = 45,
  cellwidth = 24,
  cellheight = 10,
  silent = TRUE
)

gt <- ph$gtable
nm_gt <- as.character(gt$layout$name)

cn_idx <- grep("col.*name", nm_gt, ignore.case = TRUE)
if (length(cn_idx) > 0) {
  cn_row <- unique(gt$layout$t[cn_idx])[1]
  gt$heights[cn_row] <- grid::unit(4.5, "cm")
}

rn_idx <- grep("row.*name", nm_gt, ignore.case = TRUE)
if (length(rn_idx) > 0) {
  rn_col <- unique(gt$layout$l[rn_idx])[1]
  gt$widths[rn_col] <- grid::unit(5.0, "cm")
}

save_gtable_safe <- function(gt, pdf_out, png_out, title,
                             extra_w_in = 1.5, extra_h_in = 0.2,
                             inset_mm = 6, title_in = 0.35, res = 600) {
  
  w0 <- grid::convertWidth(sum(gt$widths),  "in", valueOnly = TRUE)
  h0 <- grid::convertHeight(sum(gt$heights),"in", valueOnly = TRUE)
  w_in <- w0 + extra_w_in
  h_in <- h0 + extra_h_in + title_in
  
  draw_all <- function() {
    grid::grid.newpage()
    lay <- grid::grid.layout(nrow = 2, ncol = 1,
                             heights = grid::unit.c(grid::unit(title_in,"in"), grid::unit(1,"null")))
    grid::pushViewport(grid::viewport(layout = lay))
    
    grid::pushViewport(grid::viewport(layout.pos.row = 1))
    grid::grid.text(title, gp = grid::gpar(fontsize = 18, fontface = "bold"))
    grid::popViewport()
    
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    grid::pushViewport(grid::viewport(
      width  = grid::unit(1, "npc") - grid::unit(2*inset_mm, "mm"),
      height = grid::unit(1, "npc") - grid::unit(2*inset_mm, "mm"),
      clip="off"
    ))
    grid::grid.draw(gt)
    grid::popViewport(2)
    grid::popViewport()
  }
  
  grDevices::pdf(pdf_out, width = w_in, height = h_in, useDingbats = FALSE)
  draw_all(); grDevices::dev.off()
  
  grDevices::png(png_out, width = w_in, height = h_in, units="in", res=res)
  draw_all(); grDevices::dev.off()
}

hm_pdf <- file.path(out_dir, "Heatmaps", "Fig_STRATIFIED_Heatmap_TaxaRows_PathwaysCols_log10Prop.pdf")
hm_png <- file.path(out_dir, "Heatmaps", "Fig_STRATIFIED_Heatmap_TaxaRows_PathwaysCols_log10Prop.png")
save_gtable_safe(gt, hm_pdf, hm_png, title = "Major microbial contributors to pathways")
message("Who-does-what heatmap written: ", hm_pdf)

# ============================================================
# 9b) DOTPLOTS (BOTH versions)
# ============================================================
make_dotplots <- function(df_long, TOP_PWY=19, TOP_TAXA=15, UNCL_LABEL="Unclassified",
                          PSEUDO=1e-6, option=c("within_pathway","classified_only")) {
  
  option <- match.arg(option)
  
  pwy_global <- df_long %>%
    group_by(PathwayBase) %>%
    summarise(Total = sum(Value, na.rm = TRUE), .groups="drop") %>%
    arrange(desc(Total)) %>%
    slice_head(n = TOP_PWY) %>%
    pull(PathwayBase)
  
  tax_totals <- df_long %>%
    filter(PathwayBase %in% pwy_global) %>%
    group_by(Taxon) %>%
    summarise(Total = sum(Value, na.rm = TRUE), .groups="drop") %>%
    arrange(desc(Total))
  
  taxa_top <- tax_totals %>%
    filter(Taxon != UNCL_LABEL) %>%
    slice_head(n = TOP_TAXA) %>%
    pull(Taxon)
  
  taxa_global <- c(taxa_top, UNCL_LABEL)
  
  tax_pwy <- df_long %>%
    filter(PathwayBase %in% pwy_global, Taxon %in% taxa_global) %>%
    group_by(PathwayBase, Taxon) %>%
    summarise(Total = sum(Value, na.rm = TRUE), .groups="drop") %>%
    group_by(PathwayBase) %>%
    mutate(
      PathwayTotal    = sum(Total, na.rm = TRUE),
      ClassifiedTotal = sum(Total[Taxon != UNCL_LABEL], na.rm = TRUE),
      Prop_within     = ifelse(PathwayTotal > 0, Total / PathwayTotal, 0),
      Prop_classified = if_else(Taxon == UNCL_LABEL, NA_real_,
                                if_else(ClassifiedTotal > 0, Total / ClassifiedTotal, 0))
    ) %>%
    ungroup() %>%
    mutate(
      PathwayLabel = drop_pwy_id(PathwayBase),
      TaxonLabel   = gsub("_"," ", Taxon),
      log10_within = log10(Prop_within + PSEUDO),
      log10_class  = log10(Prop_classified + PSEUDO)
    )
  
  tax_order <- tax_pwy %>%
    filter(Taxon != UNCL_LABEL) %>%
    group_by(TaxonLabel) %>%
    summarise(Total = sum(Total, na.rm=TRUE), .groups="drop") %>%
    arrange(desc(Total)) %>%
    pull(TaxonLabel)
  
  tax_pwy <- tax_pwy %>%
    mutate(
      PathwayLabel = factor(PathwayLabel, levels = rev(drop_pwy_id(pwy_global))),
      TaxonLabel   = factor(TaxonLabel, levels = c(tax_order, UNCL_LABEL))
    )
  
  if (option == "within_pathway") {
    plot_df <- tax_pwy %>% filter(Prop_within > 0)
    size_var <- "Prop_within"
    col_var  <- "log10_within"
    subtitle <- "Dot size = proportion within pathway (includes Unclassified); colour = log10(proportion)"
    out_tag  <- "WITH_Unclassified"
    legend_size <- "Proportion\n(within pathway)"
  } else {
    plot_df <- tax_pwy %>% filter(!is.na(Prop_classified), Prop_classified > 0)
    size_var <- "Prop_classified"
    col_var  <- "log10_class"
    subtitle <- "Dot size = proportion among classified taxa; colour = log10(proportion)"
    out_tag  <- "CLASSIFIED_only"
    legend_size <- "Proportion\n(classified taxa)"
  }
  
  nm_cols <- c("#1b3a8a", "#3b7ddd", "#bfe3ff", "#f7f7f7", "#ffd08a", "#f07c2b")
  nm_vals <- scales::rescale(c(-4,-3,-2,-1,-0.5,0))
  
  p <- ggplot(plot_df, aes(x = TaxonLabel, y = PathwayLabel)) +
    geom_point(aes(size = .data[[size_var]], colour = .data[[col_var]]),
               alpha = 0.95, stroke = 0.2) +
    scale_size_continuous(range = c(0.6, 11),
                          labels = scales::percent_format(accuracy = 1),
                          name = legend_size) +
    scale_colour_gradientn(colours = nm_cols, values = nm_vals,
                           limits = c(-4,0), name = "log10(prop)") +
    labs(title = "Major microbial contributors to pathways",
         subtitle = subtitle, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_line(color="grey92", linewidth=0.3),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(face="bold", size=14),
      legend.title = element_text(face="bold")
    )
  
  pdf_out <- file.path(out_dir, "DotPlots", paste0("Fig_DOT_TaxaByPathway_", out_tag, ".pdf"))
  png_out <- file.path(out_dir, "DotPlots", paste0("Fig_DOT_TaxaByPathway_", out_tag, ".png"))
  ggsave(pdf_out, p, width = 14, height = 8, useDingbats = FALSE)
  ggsave(png_out, p, width = 14, height = 8, dpi = 600)
  
  write_tsv2(tax_pwy, paste0("TABLE_DotPlot_TaxaByPathway_", out_tag, ".tsv"))
  message("Dotplot written: ", pdf_out)
  
  invisible(p)
}

make_dotplots(df_long, option = "within_pathway")
make_dotplots(df_long, option = "classified_only")


# ----------------------------
# 10) Output manifest (MD5 of all outputs)
# ----------------------------
outs <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
manifest <- tibble::tibble(
  file  = outs,
  bytes = file.info(outs)$size,
  md5   = vapply(outs, function(fp) as.character(tools::md5sum(fp)), FUN.VALUE = character(1))
)
readr::write_tsv(manifest, file.path(out_dir, "RunInfo", paste0("MANIFEST_outputs_", run_id, ".tsv")))

# ============================================================
# 11) NatCom Figure 5 panel: stratified pathway contributors
#     Taxon-resolved contributors to pathway responses
# ============================================================

COMPARTMENT <- "Root"
PANEL_TAG <- "A"
NATCOM_DIR <- STAGING_DIR
dir.create(NATCOM_DIR, recursive = TRUE, showWarnings = FALSE)

TOP_PWY_NATCOM  <- 12
TOP_TAXA_NATCOM <- 12
UNCL_LABEL <- "Unclassified"
PSEUDO <- 1e-6

# ----------------------------
# Helper: italicise scientific names
# ----------------------------
italic_taxa <- function(x) {
  
  sapply(x, function(z) {
    
    # Genus sp. Strain
    if (grepl(" sp ", z)) {
      
      genus  <- sub(" sp .*", "", z)
      strain <- sub(".* sp ", "", z)
      
      return(
        paste0(
          "italic('", genus, "')~'sp.'~'",
          strain,
          "'"
        )
      )
    }
    
    # Genus species
    words <- strsplit(z, " ")[[1]]
    
    if (length(words) >= 2) {
      
      return(
        paste0(
          "italic('",
          paste(words[1:2], collapse = " "),
          "')"
        )
      )
    }
    
    # fallback
    paste0("'", z, "'")
  })
}

# ----------------------------
# Helper: clean pathway names
# ----------------------------
clean_pathway_natcom <- function(x) {
  x <- drop_pwy_id(x)
  x <- stringr::str_replace_all(x, "superpathway of ", "")
  x <- stringr::str_replace_all(x, "&alpha;", "α")
  x <- stringr::str_replace_all(x, "&beta;", "β")
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_squish(x)
  stringr::str_wrap(x, width = 24)
}

# ----------------------------
# Helper: sentence-case pathway names
# ----------------------------
sentence_case_pathway <- function(x) {
  
  x <- stringr::str_squish(x)
  
  # Capitalize first character only
  x <- paste0(
    toupper(substr(x, 1, 1)),
    substr(x, 2, nchar(x))
  )
  
  # Restore important abbreviations
  x <- x %>%
    stringr::str_replace_all("\\bNadph\\b", "NADPH") %>%
    stringr::str_replace_all("\\bTca\\b", "TCA") %>%
    stringr::str_replace_all("\\bIva\\b", "IVA") %>%
    stringr::str_replace_all("\\bAtp\\b", "ATP") %>%
    stringr::str_replace_all("\\bNad\\b", "NAD") %>%
    stringr::str_replace_all("\\bNadh\\b", "NADH") %>%
    stringr::str_replace_all("\\bUdp\\b", "UDP") %>%
    stringr::str_replace_all("\\bDtdp\\b", "dTDP")
  
  x
}

# ----------------------------
# Select major pathways
# ----------------------------
pwy_natcom <- df_long %>%
  group_by(PathwayBase) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = TOP_PWY_NATCOM) %>%
  pull(PathwayBase)

# ----------------------------
# Select major classified taxa
# ----------------------------
taxa_natcom <- df_long %>%
  filter(
    PathwayBase %in% pwy_natcom,
    Taxon != UNCL_LABEL
  ) %>%
  group_by(Taxon) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = TOP_TAXA_NATCOM) %>%
  pull(Taxon)

# ----------------------------
# Prepare plotting dataframe
# ----------------------------
plot_df_natcom <- df_long %>%
  filter(
    PathwayBase %in% pwy_natcom,
    Taxon %in% taxa_natcom
  ) %>%
  group_by(PathwayBase, Taxon) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  group_by(PathwayBase) %>%
  mutate(
    PathwayTotal = sum(Total, na.rm = TRUE),
    Prop = ifelse(PathwayTotal > 0, Total / PathwayTotal, 0),
    log10Prop = log10(Prop + PSEUDO)
  ) %>%
  ungroup() %>%
  mutate(
    PathwayLabel = clean_pathway_natcom(PathwayBase),
    PathwayLabel = sentence_case_pathway(PathwayLabel),
    TaxonLabel = stringr::str_replace_all(Taxon, "_", " ")
  )

# Shorten and wrap pathway labels safely
plot_df_natcom$PathwayLabel <- stringr::str_replace(
  plot_df_natcom$PathwayLabel,
  "dTDP-α-D-ravidosamine and dTDP-4-acetyl-α-D-ravidosamine biosynthesis",
  "dTDP-ravidosamine biosynthesis"
)

plot_df_natcom$PathwayLabel <- stringr::str_wrap(
  plot_df_natcom$PathwayLabel,
  width = 24
)

# Keep pathway order by total abundance
path_order <- plot_df_natcom %>%
  group_by(PathwayLabel) %>%
  summarise(Total = sum(Total, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(PathwayLabel)

# Keep taxon order by total contribution
tax_order <- plot_df_natcom %>%
  group_by(TaxonLabel) %>%
  summarise(Total = sum(Total, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(TaxonLabel)

plot_df_natcom <- plot_df_natcom %>%
  mutate(
    PathwayLabel = factor(PathwayLabel, levels = rev(path_order)),
    TaxonLabel = factor(TaxonLabel, levels = rev(tax_order))
  )

# ----------------------------
# Create Figure 5 panel
# ----------------------------
p_fig5_panel <- ggplot(
  plot_df_natcom,
  aes(x = PathwayLabel, y = TaxonLabel)
) +
  geom_point(
    aes(size = Prop, colour = log10Prop),
    alpha = 0.95
  ) +
  scale_y_discrete(
    labels = function(x)
      parse(text = italic_taxa(as.character(x)))
  ) +
  scale_size_continuous(
    name = "Contribution",
    range = c(1.2, 7),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_colour_gradientn(
    colours = c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C"),
    limits = c(-4, 0),
    name = "log10\nproportion"
  ) +
  labs(
    title = COMPARTMENT,
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black"),
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

# ----------------------------
# Save panel for later Figure 5 merge
# ----------------------------
ggsave(
  file.path(
    NATCOM_DIR,
    paste0("Fig5", PANEL_TAG, "_", COMPARTMENT, "_Stratified_Pathway_Contributors.pdf")
  ),
  p_fig5_panel,
  width = 140,
  height = 120,
  units = "mm",
  useDingbats = FALSE
)

ggsave(
  file.path(
    NATCOM_DIR,
    paste0("Fig5", PANEL_TAG, "_", COMPARTMENT, "_Stratified_Pathway_Contributors.png")
  ),
  p_fig5_panel,
  width = 140,
  height = 120,
  units = "mm",
  dpi = 600
)

saveRDS(
  p_fig5_panel,
  file.path(
    NATCOM_DIR,
    paste0("p_Fig5", PANEL_TAG, "_", COMPARTMENT, "_Stratified_Pathway_Contributors.rds")
  )
)

write_tsv2(
  plot_df_natcom,
  paste0("TABLE_Fig5", PANEL_TAG, "_", COMPARTMENT, "_Stratified_Pathway_Contributors.tsv")
)




# ============================================================
# 12) Supplementary Figure S5 panel: same style as Figure 5, WITH Unclassified
#     Top classified taxa plus Unclassified; proportions use all taxa as denominator
# ============================================================

SUPP_PANEL_TAG <- "B"
TOP_PWY_SUPP  <- TOP_PWY_NATCOM
TOP_TAXA_SUPP <- TOP_TAXA_NATCOM

# Select the same major pathways as the NatCom contributor panel
pwy_supp <- pwy_natcom

# Select top classified taxa, then force inclusion of Unclassified
taxa_supp_classified <- df_long %>%
  filter(
    PathwayBase %in% pwy_supp,
    Taxon != UNCL_LABEL
  ) %>%
  group_by(Taxon) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = TOP_TAXA_SUPP) %>%
  pull(Taxon)

taxa_supp <- c(taxa_supp_classified, UNCL_LABEL)

# Denominator: total pathway abundance across all taxa, including Unclassified
pwy_total_all <- df_long %>%
  filter(PathwayBase %in% pwy_supp) %>%
  group_by(PathwayBase) %>%
  summarise(PathwayTotalAll = sum(Value, na.rm = TRUE), .groups = "drop")

plot_df_supp <- df_long %>%
  filter(
    PathwayBase %in% pwy_supp,
    Taxon %in% taxa_supp
  ) %>%
  group_by(PathwayBase, Taxon) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  left_join(pwy_total_all, by = "PathwayBase") %>%
  mutate(
    Prop = ifelse(PathwayTotalAll > 0, Total / PathwayTotalAll, 0),
    log10Prop = log10(Prop + PSEUDO)
  ) %>%
  ungroup() %>%
  mutate(
    PathwayLabel = clean_pathway_natcom(PathwayBase),
    PathwayLabel = sentence_case_pathway(PathwayLabel),
    TaxonLabel = stringr::str_replace_all(Taxon, "_", " ")
  )

# Shorten and wrap pathway labels safely
plot_df_supp$PathwayLabel <- stringr::str_replace(
  plot_df_supp$PathwayLabel,
  "dTDP-α-D-ravidosamine and dTDP-4-acetyl-α-D-ravidosamine biosynthesis",
  "dTDP-ravidosamine biosynthesis"
)

plot_df_supp$PathwayLabel <- stringr::str_wrap(
  plot_df_supp$PathwayLabel,
  width = 24
)

# Keep pathway order consistent with the main contributor panel
path_order_supp <- plot_df_supp %>%
  group_by(PathwayLabel) %>%
  summarise(Total = sum(Total, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(PathwayLabel)

# Keep top classified taxa ordered by contribution, with Unclassified placed last
tax_order_supp_classified <- plot_df_supp %>%
  filter(TaxonLabel != UNCL_LABEL) %>%
  group_by(TaxonLabel) %>%
  summarise(Total = sum(Total, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(TaxonLabel)

tax_order_supp <- c(tax_order_supp_classified, UNCL_LABEL)

plot_df_supp <- plot_df_supp %>%
  mutate(
    PathwayLabel = factor(PathwayLabel, levels = rev(path_order_supp)),
    TaxonLabel = factor(TaxonLabel, levels = rev(tax_order_supp))
  )

p_figS5_panel <- ggplot(
  plot_df_supp,
  aes(x = PathwayLabel, y = TaxonLabel)
) +
  geom_point(
    aes(size = Prop, colour = log10Prop),
    alpha = 0.95
  ) +
  scale_y_discrete(
    labels = function(x)
      parse(text = italic_taxa(as.character(x)))
  ) +
  scale_size_continuous(
    name = "Contribution",
    range = c(1.2, 7),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_colour_gradientn(
    colours = c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C"),
    limits = c(-4, 0),
    name = "log10\nproportion"
  ) +
  labs(
    title = COMPARTMENT,
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black"),
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(
  file.path(
    NATCOM_DIR,
    paste0("FigS5", SUPP_PANEL_TAG, "_", COMPARTMENT, "_WITH_Unclassified_panel.pdf")
  ),
  p_figS5_panel,
  width = 140,
  height = 120,
  units = "mm",
  useDingbats = FALSE
)

ggsave(
  file.path(
    NATCOM_DIR,
    paste0("FigS5", SUPP_PANEL_TAG, "_", COMPARTMENT, "_WITH_Unclassified_panel.png")
  ),
  p_figS5_panel,
  width = 140,
  height = 120,
  units = "mm",
  dpi = 600
)

saveRDS(
  p_figS5_panel,
  file.path(
    NATCOM_DIR,
    paste0("p_FigS5", SUPP_PANEL_TAG, "_", COMPARTMENT, "_WITH_Unclassified.rds")
  )
)

write_tsv2(
  plot_df_supp,
  paste0("TABLE_FigS5", SUPP_PANEL_TAG, "_", COMPARTMENT, "_WITH_Unclassified.tsv")
)


message("DONE. Outputs saved to: ", out_dir)
message("Run log: ", log_fp)

})
})


# ============================================================
# Merge Figure 5 stratified pathway contributors
# ============================================================
run_section("Merge Figure 5 stratified pathway contributors", {
local({
# ============================================================
# Figure 5. Stratified pathway contributors
# Merge Rhizosphere and Root RDS panels
# Legend shown on Root panel only
# ============================================================

pkgs <- c("ggplot2", "cowplot")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})
NATCOM_DIR <- STAGING_DIR
dir.create(NATCOM_DIR, recursive = TRUE, showWarnings = FALSE)

root_rds <- file.path(
  NATCOM_DIR,
  "p_Fig5A_Root_Stratified_Pathway_Contributors.rds"
)

rhiz_rds <- file.path(
  NATCOM_DIR,
  "p_Fig5B_Rhizosphere_Stratified_Pathway_Contributors.rds"
)

stopifnot(file.exists(root_rds))
stopifnot(file.exists(rhiz_rds))

italic_taxa <- function(x) {
  
  sapply(x, function(z) {
    
    if (grepl(" sp\\. | sp ", z)) {
      genus <- sub(" sp\\. .*| sp .*", "", z)
      strain <- sub(".* sp\\. |.* sp ", "", z)
      
      return(
        paste0(
          "italic('", genus, "')~'sp.'~'",
          strain,
          "'"
        )
      )
    }
    
    words <- strsplit(z, " ")[[1]]
    
    if (length(words) >= 2) {
      return(
        paste0(
          "italic('",
          paste(words[1:2], collapse = " "),
          "')"
        )
      )
    }
    
    paste0("'", z, "'")
  })
}

p_root <- readRDS(root_rds)
p_rhiz <- readRDS(rhiz_rds)

# Compact display labels for the merged Figure 5 only. The complete pathway
# names remain unchanged in the data, RDS objects and exported tables.
short_pathway_labels_fig5 <- function(x) {
  y <- gsub("\n", " ", as.character(x))
  y <- gsub("\\s+", " ", trimws(y))

  vapply(y, function(z) {
    if (grepl("Hexitol fermentation", z, ignore.case = TRUE)) {
      return("Hexitol\nfermentation")
    }
    if (grepl("S-adenosyl|methionine salvage", z, ignore.case = TRUE)) {
      return("SAM\nsalvage")
    }
    if (grepl("Lipid IVA.*P\\. putida", z, ignore.case = TRUE)) {
      return("Lipid IVA\nbiosynthesis\n(P. putida)")
    }
    if (grepl("Lipid IVA.*E\\. coli", z, ignore.case = TRUE)) {
      return("Lipid IVA\nbiosynthesis\n(E. coli)")
    }
    if (grepl("Cytosolic NADPH", z, ignore.case = TRUE)) {
      return("Cytosolic\nNADPH\nproduction")
    }
    if (grepl("Queuosine", z, ignore.case = TRUE)) {
      return("Queuosine\nbiosynthesis")
    }
    if (grepl("Fatty acid biosynthesis", z, ignore.case = TRUE)) {
      suffix <- if (grepl("E\\. coli", z, ignore.case = TRUE)) "\n(E. coli)" else ""
      return(paste0("Fatty acid\nbiosynthesis", suffix))
    }
    if (grepl("TCA cycle VII", z, ignore.case = TRUE)) {
      return("TCA cycle VII")
    }
    if (grepl("D-galactose", z, ignore.case = TRUE)) {
      return("D-galactose\ndegradation")
    }
    if (grepl("Stachyose", z, ignore.case = TRUE)) {
      return("Stachyose\ndegradation")
    }
    if (grepl("Homocysteine|Cys.*Homocys", z, ignore.case = TRUE)) {
      return("Cys–Homocys\ninterconversion")
    }
    if (grepl("ravidosamine", z, ignore.case = TRUE)) {
      return("dTDP-\nravidosamine\nbiosynthesis")
    }

    stringr::str_wrap(z, width = 12)
  }, character(1))
}

# Panel A: Rhizosphere, remove legend
p_rhiz_clean <- p_rhiz +
  scale_x_discrete(labels = short_pathway_labels_fig5) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      size = FIG5_PATHWAY_TEXT_SIZE,
      angle = FIG5_PATHWAY_TEXT_ANGLE,
      hjust = 1,
      vjust = 1,
      lineheight = 0.9,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = FIG5_TAXON_TEXT_SIZE,
      colour = "black"
    ),
    axis.title = element_text(size = FIG5_SOURCE_BASE_SIZE + 1, face = "bold"),
    plot.title = element_text(
      size = FIG5_TITLE_SIZE,
      face = "bold",
      hjust = 0
    ),
    plot.margin = margin(8, 10, 8, 8)
  )

# Panel B: Root, keep legend on the right
p_root_clean <- p_root +
  scale_x_discrete(labels = short_pathway_labels_fig5) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    axis.text.x = element_text(
      size = FIG5_PATHWAY_TEXT_SIZE,
      angle = FIG5_PATHWAY_TEXT_ANGLE,
      hjust = 1,
      vjust = 1,
      lineheight = 0.9,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = FIG5_TAXON_TEXT_SIZE,
      colour = "black"
    ),
    axis.title = element_text(size = FIG5_SOURCE_BASE_SIZE + 1, face = "bold"),
    plot.title = element_text(
      size = FIG5_TITLE_SIZE,
      face = "bold",
      hjust = 0
    ),
    legend.title = element_text(size = FIG5_LEGEND_TITLE_SIZE, face = "bold"),
    legend.text = element_text(size = FIG5_LEGEND_TEXT_SIZE),
    legend.key.height = grid::unit(5, "mm"),
    legend.spacing.y = grid::unit(1.5, "mm"),
    plot.margin = margin(8, 8, 8, 8)
  )

fig5_stratified_contributors <- cowplot::plot_grid(
  p_rhiz_clean,
  p_root_clean,
  nrow = 1,
  labels = c("A", "B"),
  label_size = FIG5_TAG_SIZE,
  label_fontface = "bold",
  label_x = c(0.005, 0.005),
  label_y = c(0.995, 0.995),
  hjust = 0,
  vjust = 1,
  rel_widths = c(1.15, 1.22),
  align = "h",
  axis = "tb"
)

ggsave(
  file.path(NATCOM_DIR, "Figure5_Stratified_Pathway_Contributors_NatCom.pdf"),
  fig5_stratified_contributors,
  width = FIG5_EXPORT_WIDTH_MM,
  height = FIG5_EXPORT_HEIGHT_MM,
  units = "mm",
  useDingbats = FALSE
)

ggsave(
  file.path(NATCOM_DIR, "Figure5_Stratified_Pathway_Contributors_NatCom.png"),
  fig5_stratified_contributors,
  width = FIG5_EXPORT_WIDTH_MM,
  height = FIG5_EXPORT_HEIGHT_MM,
  units = "mm",
  dpi = 600
)

ggsave(
  file.path(NATCOM_DIR, "Figure5_Stratified_Pathway_Contributors_NatCom.svg"),
  fig5_stratified_contributors,
  width = FIG5_EXPORT_WIDTH_MM,
  height = FIG5_EXPORT_HEIGHT_MM,
  units = "mm"
)

message("Figure 5 saved to: ", NATCOM_DIR)
})
})


# ============================================================
# Merge Supplementary Figure S5 with unclassified
# ============================================================
run_section("Merge Supplementary Figure S5 with unclassified", {
local({
# ============================================================
# Supplementary Figure S5 — Taxonomic contributors WITH Unclassified
# Merge Rhizosphere + Root RDS panels, left-right layout
# ============================================================

pkgs <- c("ggplot2", "cowplot")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})
NATCOM_DIR <- STAGING_DIR
dir.create(NATCOM_DIR, recursive = TRUE, showWarnings = FALSE)

rhiz_rds <- file.path(
  NATCOM_DIR,
  "p_FigS5A_Rhizosphere_WITH_Unclassified.rds"
)

root_rds <- file.path(
  NATCOM_DIR,
  "p_FigS5B_Root_WITH_Unclassified.rds"
)

stopifnot(file.exists(rhiz_rds))
stopifnot(file.exists(root_rds))

short_pathway_labels <- function(x) {
  x <- as.character(x)
  
  x <- stringr::str_replace_all(x, "&alpha;|α|@", "alpha")
  x <- stringr::str_replace_all(x, "\\([^\\)]*\\)", "")
  x <- stringr::str_squish(x)
  
  x <- dplyr::case_when(
    stringr::str_detect(x, "ravidosamine") ~ "Ravidosamine\nbiosynthesis",
    stringr::str_detect(x, "Mannitol cycle") ~ "Mannitol\ncycle",
    stringr::str_detect(x, "Hexitol fermentation") ~ "Hexitol\nfermentation",
    stringr::str_detect(x, "Polyamine biosynthesis") ~ "Polyamine\nbiosynthesis",
    stringr::str_detect(x, "CMP-legionaminate") ~ "CMP-legionaminate\nbiosynthesis",
    stringr::str_detect(x, "tryptophan degradation") ~ "L-tryptophan\ndegradation",
    stringr::str_detect(x, "Fatty acid biosynthesis II") ~ "Fatty acid\nbiosynthesis II",
    stringr::str_detect(x, "Fatty acid biosynthesis I") ~ "Fatty acid\nbiosynthesis I",
    stringr::str_detect(x, "S-adenosyl") ~ "Methionine\nsalvage",
    stringr::str_detect(x, "NADPH") ~ "NADPH\nproduction",
    stringr::str_detect(x, "Lipid IVA") ~ "Lipid IVA\nbiosynthesis",
    stringr::str_detect(x, "Queuosine") ~ "Queuosine\nbiosynthesis",
    stringr::str_detect(x, "lysine fermentation") ~ "L-lysine\nfermentation",
    stringr::str_detect(x, "propanediol") ~ "1,3-propanediol\nbiosynthesis",
    stringr::str_detect(x, "lysine degradation") ~ "L-lysine\ndegradation",
    stringr::str_detect(x, "chlorophyllide") ~ "Chlorophyllide\nbiosynthesis",
    stringr::str_detect(x, "Tetrahydrofolate") ~ "Tetrahydrofolate\nbiosynthesis",
    stringr::str_detect(x, "Pyridoxal") ~ "Pyridoxal phosphate\nbiosynthesis",
    stringr::str_detect(x, "Homocysteine") ~ "Homocysteine-\ncysteine\ninterconversion",
    stringr::str_detect(x, "D-galactose") ~ "D-galactose\ndegradation",
    stringr::str_detect(x, "TCA cycle") ~ "TCA cycle VII",
    stringr::str_detect(x, "Stachyose") ~ "Stachyose\ndegradation",
    TRUE ~ stringr::str_wrap(x, width = 8)
  )
  
  x
}

# Needed because the saved ggplot objects use parse(text = italic_taxa(...))
italic_taxa <- function(x) {
  sapply(x, function(z) {
    if (z == "Unclassified") return("'Unclassified'")
    
    if (grepl(" sp\\. | sp ", z)) {
      genus <- sub(" sp\\. .*| sp .*", "", z)
      strain <- sub(".* sp\\. |.* sp ", "", z)
      return(paste0("italic('", genus, "')~'sp.'~'", strain, "'"))
    }
    
    words <- strsplit(z, " ")[[1]]
    if (length(words) >= 2) {
      return(paste0("italic('", paste(words[1:2], collapse = " "), "')"))
    }
    
    paste0("'", z, "'")
  })
}

p_rhiz <- readRDS(rhiz_rds) +
  scale_x_discrete(
    labels = function(x) gsub("\n", " ", short_pathway_labels(x))
  ) +
  scale_size_continuous(
    name = "Contribution",
    range = FIGS5_POINT_SIZE_RANGE,
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_colour_gradientn(
    colours = c(FIGS5_LOW_COLOUR, "#4292C6", "#FFFFBF", "#FDAE61", "#D7191C"),
    limits = c(-4, 0),
    oob = scales::squish,
    name = "log10\nproportion"
  ) +
  geom_point(
    aes(size = Prop),
    shape = 21,
    fill = NA,
    colour = FIGS5_OUTLINE_COLOUR,
    stroke = FIGS5_OUTLINE_STROKE,
    show.legend = FALSE
  ) +
  labs(title = "Rhizosphere") +
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 11,
      face = "plain",
      margin = margin(t = 0)
    ),
    
    axis.text.y = element_text(
      size = 12,
      face = "italic"
    ),
    
    plot.title = element_text(
      size = 20,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.ticks.length.x = unit(2, "pt"),
    plot.margin = margin(8, 10, 0, 8)
  )

p_root <- readRDS(root_rds) +
  scale_x_discrete(
    labels = function(x) gsub("\n", " ", short_pathway_labels(x))
  ) +
  scale_size_continuous(
    name = "Contribution",
    range = FIGS5_POINT_SIZE_RANGE,
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_colour_gradientn(
    colours = c(FIGS5_LOW_COLOUR, "#4292C6", "#FFFFBF", "#FDAE61", "#D7191C"),
    limits = c(-4, 0),
    oob = scales::squish,
    name = "log10\nproportion"
  ) +
  geom_point(
    aes(size = Prop),
    shape = 21,
    fill = NA,
    colour = FIGS5_OUTLINE_COLOUR,
    stroke = FIGS5_OUTLINE_STROKE,
    show.legend = FALSE
  ) +
  labs(title = "Root") +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 11,
      face = "plain",
      margin = margin(t = 0)
    ),
    
    axis.text.y = element_text(
      size = 12,
      face = "italic"
    ),
    
    plot.title = element_text(
      size = 20,
      face = "bold",
      hjust = 0.5
    ),
    
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 12
    ),
    
    axis.ticks.length.x = unit(2, "pt"),
    plot.margin = margin(8, 10, 0, 8)
  )

figS5_with_unclassified <- cowplot::plot_grid(
  p_rhiz,
  p_root,
  nrow = 1,
  labels = c("A", "B"),
  label_size = 24,
  label_fontface = "bold",
  label_x = c(0.005, 0.005),
  label_y = c(0.995, 0.995),
  hjust = 0,
  vjust = 1,
  rel_widths = c(1.05, 1.15),
  align = "h",
  axis = "tb"
)

ggsave(
  file.path(NATCOM_DIR, "FigureS5_Taxonomic_Contributors_WITH_Unclassified_NatCom.pdf"),
  figS5_with_unclassified,
  width = 460,
  height = 200,
  units = "mm",
  useDingbats = FALSE
)

ggsave(
  file.path(NATCOM_DIR, "FigureS5_Taxonomic_Contributors_WITH_Unclassified_NatCom.png"),
  figS5_with_unclassified,
  width = 460,
  height = 200,
  units = "mm",
  dpi = 600
)

ggsave(
  file.path(NATCOM_DIR, "FigureS5_Taxonomic_Contributors_WITH_Unclassified_NatCom.svg"),
  figS5_with_unclassified,
  width = 460,
  height = 200,
  units = "mm"
)

message("Supplementary Figure S5 with Unclassified saved to: ", NATCOM_DIR)

})
})

# ============================================================
# Copy manuscript-ready outputs
# ============================================================
copy_outputs()
message("\nDONE. Functional pipeline complete.")
