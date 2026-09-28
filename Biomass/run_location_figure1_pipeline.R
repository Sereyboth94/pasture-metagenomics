# ============================================================
# Figure 1 standalone pipeline
# Plant biomass (Figure 1A) + taxonomy and pathway PCoA panels
# Nature Communications style
#
# Repository location:
#   Biomass/run_location_figure1_pipeline.R
#
# Run from repository root in RStudio:
#   source("Biomass/run_location_figure1_pipeline.R")
#
# Or from command line:
#   Rscript Biomass/run_location_figure1_pipeline.R
#
# Expected repository structure:
#   Github/
#   ├── location/
#   │   └── results/
#   ├── Biomass/
#   │   ├── run_location_figure1_pipeline.R
#   │   └── Plant_Biomass_ANOVA_results/
#   ├── taxonomy/
#   │   └── results/
#   └── function/
#       └── results/
#
# Outputs:
#   location/results/Figure1_NatCom_merged_v24.pdf
#   location/results/Figure1_NatCom_merged_v24.png
#   location/results/Figure1_NatCom_merged_v24.rds
#   location/results/Figure1_NatCom_merged_v24_components.rds
#
# Manuscript figure copies:
#   manuscript_figures/Main_Figures/10_plant_biomass_panel_600dpi.png
#   manuscript_figures/Main_Figures/Fig1B_Taxonomy_PCoA_NatCom_panel.png
#   manuscript_figures/Main_Figures/Fig1C_Pathway_PCoA_NatCom_panel.png
#   manuscript_figures/Main_Figures/Figure1_NatCom_merged_v24.pdf
#   manuscript_figures/Main_Figures/Figure1_NatCom_merged_v24.png
#   manuscript_figures/Main_Figures/Figure1_NatCom_merged_v24.rds
#   manuscript_figures/Main_Figures/Figure1_NatCom_merged_v24_components.rds
# ============================================================

# ----------------------------
# 0) Packages
# ----------------------------
required_pkgs <- c(
  "ggplot2",
  "cowplot"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Please install missing packages before running this script:\n",
    paste0("install.packages(c(", paste(sprintf('\"%s\"', missing_pkgs), collapse = ", "), "))"),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

pipeline_version <- "Figure 1 pipeline v24: recolour stored B and C by six field environments"
message(pipeline_version)

# ----------------------------
# 1) Project paths
# ----------------------------
# Fixed project root for this analysis. All input and output paths below are
# derived from this directory, regardless of the RStudio working directory.
repo_root <- paste0(
  "C:/Users/soths/OneDrive - Lincoln University/Writing/NatCom/",
  "Github"
)

if (!dir.exists(repo_root)) {
  stop("Repository root not found: ", repo_root, call. = FALSE)
}

repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
location_dir <- file.path(repo_root, "location")

# Create a recoverable backup of this pipeline before producing the figure.
script_file <- file.path(
  repo_root,
  "Biomass",
  "run_location_figure1_pipeline.R"
)
backup_file <- paste0(script_file, ".backup")

if (!file.exists(script_file)) {
  stop("Pipeline script not found: ", script_file, call. = FALSE)
}

backup_ok <- file.copy(script_file, backup_file, overwrite = TRUE)
if (!isTRUE(backup_ok)) {
  warning("Could not create script backup: ", backup_file, call. = FALSE)
} else {
  message("Created script backup: ", backup_file)
}

out_dir <- file.path(location_dir, "results")
main_fig_dir <- file.path(repo_root, "manuscript_figures", "Main_Figures")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(main_fig_dir, recursive = TRUE, showWarnings = FALSE)

copy_to_main_figures <- function(files, overwrite = TRUE) {
  files <- files[file.exists(files)]
  if (length(files) == 0) return(invisible(character()))
  dest <- file.path(main_fig_dir, basename(files))
  ok <- file.copy(files, dest, overwrite = overwrite)
  if (any(!ok)) {
    warning(
      "Some files could not be copied to manuscript_figures/Main_Figures. ",
      "Close any open output files and rerun if needed: ",
      paste(basename(files)[!ok], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(dest[ok])
}

copy_matching_sidecars <- function(primary_file) {
  # Copy matching PDF/PNG/SVG/RDS versions when they exist.
  base <- tools::file_path_sans_ext(primary_file)
  sidecars <- paste0(base, c(".pdf", ".png", ".svg", ".rds"))
  copy_to_main_figures(sidecars)
}

message("Repository root: ", repo_root)
message("Location output directory: ", out_dir)
message("Main manuscript figure directory: ", main_fig_dir)

# ----------------------------
# 2) Regenerate Figure 1A, then merge Figure 1A, 1B, and 1C
# ----------------------------
find_first_existing <- function(candidates, label) {
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop(
      "Could not find ", label, ".\n",
      "Checked:\n",
      paste(candidates, collapse = "\n"),
      call. = FALSE
    )
  }
  normalizePath(existing[1], winslash = "/", mustWork = TRUE)
}

# Panel A is a raster image, so its embedded text cannot be relabelled during
# merging. Run the biomass source script to create a fresh PNG with the
# Trichoderma treatment and percentage labels before loading that image.
biomass_script <- find_first_existing(
  file.path(repo_root, "Biomass", "Plant_Biomass_ANOVA.R"),
  "updated plant biomass source script"
)
biomass_source_lines <- readLines(biomass_script, warn = FALSE)
if (!all(vapply(
  c("bolditalic(Trichoderma)", "control_sd = sd(.x$Control)",
    "trichoderma_sd = sd(.x$Trichoderma)"),
  function(marker) any(grepl(marker, biomass_source_lines, fixed = TRUE)),
  logical(1)
))) {
  stop(
    "Biomass/Plant_Biomass_ANOVA.R does not contain the updated ",
    "treatment mean ± SD (kg ha-1) annotations. Install the supplied ",
    "Plant_Biomass_ANOVA.R before running this merger.",
    call. = FALSE
  )
}
message("Regenerating Figure 1A with Trichoderma labels: ", biomass_script)
source(biomass_script, local = new.env(parent = globalenv()))

biomass_file <- find_first_existing(
  c(
    file.path(
      repo_root,
      "Biomass",
      "Plant_Biomass_ANOVA_results",
      "10_plant_biomass_panel_600dpi.png"
    ),
    file.path(
      repo_root,
      "biomass",
      "Plant_Biomass_ANOVA_results",
      "10_plant_biomass_panel_600dpi.png"
    )
  ),
  "Figure 1A plant biomass panel"
)

tax_rds_file <- find_first_existing(
  c(
    file.path(repo_root, "taxonomy", "results", "Figure1_NatCom", "Fig1B_Taxonomy_PCoA_NatCom_panel.rds")
  ),
  "Figure 1B taxonomy PCoA RDS"
)

path_rds_file <- find_first_existing(
  c(
    file.path(repo_root, "function", "results", "merged", "Fig1C_Pathway_PCoA_NatCom_panel.rds")
  ),
  "Figure 1C pathway PCoA RDS"
)

message("Using Figure 1A plant biomass panel: ", biomass_file)
message("Using Figure 1B taxonomy RDS: ", tax_rds_file)
message("Using Figure 1C pathway RDS: ", path_rds_file)

# Load Panels B and C as editable ggplot objects. Their source pipelines use
# matching dimensions and typography and save them without embedded legends.
pB_source <- readRDS(tax_rds_file)
pC_source <- readRDS(path_rds_file)

if (!inherits(pB_source, "ggplot")) {
  stop("Figure 1B RDS does not contain a ggplot object: ", tax_rds_file, call. = FALSE)
}
if (!inherits(pC_source, "ggplot")) {
  stop("Figure 1C RDS does not contain a ggplot object: ", path_rds_file, call. = FALSE)
}

pB_source <- pB_source + ggplot2::theme(legend.position = "none")
pC_source <- pC_source + ggplot2::theme(legend.position = "none")

# Remove any panel tags stored inside older RDS objects. The merger draws the
# final uppercase labels itself so only one A, B or C appears.
pB_source <- pB_source + ggplot2::theme(plot.tag = ggplot2::element_blank())
pC_source <- pC_source + ggplot2::theme(plot.tag = ggplot2::element_blank())

# Copy all source panels into manuscript_figures/Main_Figures.
copy_to_main_figures(biomass_file)
copy_to_main_figures(c(tax_rds_file, path_rds_file))

# Reserve a narrow strip above each image for the panel letter. This keeps the
# letters at the top-left of their own panels and prevents collisions with the
# top y-axis tick labels. B and C use identical image boxes and scales.
make_image_panel <- function(
  image_file,
  label,
  label_x,
  image_height,
  label_y = 0.995
) {
  cowplot::ggdraw() +
    cowplot::draw_image(
      image_file,
      x = 0,
      y = 0,
      width = 1,
      height = image_height,
      scale = 1
    ) +
    cowplot::draw_label(
      toupper(label),
      x = label_x,
      y = label_y,
      hjust = 0,
      vjust = 1,
      size = 16,
      fontface = "bold"
    )
}

pA <- make_image_panel(
  biomass_file,
  label = "A",
  label_x = 0.075,
  image_height = 0.95,
  label_y = 0.93
)
make_plot_panel <- function(plot_object, label, label_x, plot_height = 0.93) {
  cowplot::ggdraw() +
    cowplot::draw_plot(
      plot_object,
      x = 0,
      y = 0,
      width = 1,
      height = plot_height
    ) +
    cowplot::draw_label(
      toupper(label),
      x = label_x,
      y = 0.995,
      hjust = 0,
      vjust = 1,
      size = 16,
      fontface = "bold"
    )
}

pB <- make_plot_panel(pB_source, label = "B", label_x = 0.22)
pC <- make_plot_panel(pC_source, label = "C", label_x = 0.22)

# Create a six-environment vector shared legend matching both source panels.
environment_levels <- c(
  "Eyrewell_Forest",
  "Kowhai (Irrigated)",
  "Kowhai (Rainfed)",
  "LU_H8",
  "Rolleston",
  "West_Coast"
)
environment_colours <- c(
  Eyrewell_Forest = "#D55E00",
  "Kowhai (Irrigated)" = "#8C510A",
  "Kowhai (Rainfed)" = "#E69F00",
  LU_H8 = "#009E73",
  Rolleston = "#0072B2",
  West_Coast = "#CC79A7"
)

# The older RDS plots contain all six environments in their data, but their
# point layers mapped colour to the five-level Location variable. Remap those
# stored ggplot objects; the underlying PCoA coordinates are unchanged.
use_six_environment_colours <- function(plot, panel) {
  if (!all(c("Environment", "Compartment") %in% names(plot$data))) {
    stop("Figure 1", panel, " RDS lacks Environment/Compartment data. ",
         "Regenerate this RDS with the updated source pipeline.", call. = FALSE)
  }
  environment <- trimws(as.character(plot$data$Environment))
  environment[environment %in% c("Kowhai Irrigated", "Kowhai_Irrigated")] <-
    "Kowhai (Irrigated)"
  environment[environment %in% c("Kowhai Rainfed", "Kowhai_Rainfed")] <-
    "Kowhai (Rainfed)"
  missing <- setdiff(environment_levels, unique(environment))
  unexpected <- setdiff(unique(environment), environment_levels)
  if (length(missing) || length(unexpected)) {
    stop("Figure 1", panel, " field-environment data mismatch. Missing: ",
         paste(missing, collapse = ", "), "; unexpected: ",
         paste(unexpected, collapse = ", "), ". Regenerate this RDS with ",
         "the updated source pipeline.", call. = FALSE)
  }
  plot$data$Environment <- factor(environment, levels = environment_levels)
  point_layers <- which(vapply(
    plot$layers, function(layer) inherits(layer$geom, "GeomPoint"), logical(1)
  ))
  if (length(point_layers) != 1L) {
    stop("Figure 1", panel, " must have exactly one point layer for recolouring.",
         call. = FALSE)
  }
  plot$layers[[point_layers]]$mapping$colour <-
    ggplot2::aes(colour = Environment)$colour
  # Replacing the old five-location scale may emit a harmless ggplot notice.
  plot <- suppressMessages(plot + ggplot2::scale_colour_manual(
    values = environment_colours,
    breaks = environment_levels,
    drop = FALSE
  ))
  point_colours <- unique(ggplot2::ggplot_build(plot)$data[[point_layers]]$colour)
  if (!setequal(point_colours, unname(environment_colours))) {
    stop("Figure 1", panel, " six-colour validation failed. Observed: ",
         paste(point_colours, collapse = ", "), call. = FALSE)
  }
  message("Figure 1", panel, ": mapped six field environments to six colours.")
  plot
}
pB_source <- use_six_environment_colours(pB_source, "B")
pC_source <- use_six_environment_colours(pC_source, "C")

# Refresh wrappers after the source ggplots have been corrected.
pB <- make_plot_panel(pB_source, label = "B", label_x = 0.22)
pC <- make_plot_panel(pC_source, label = "C", label_x = 0.22)

legend_data <- expand.grid(
  Environment = factor(environment_levels, levels = environment_levels),
  Compartment = factor(
    c("Rhizosphere", "Root"),
    levels = c("Rhizosphere", "Root")
  )
)

legend_plot <- ggplot2::ggplot(
  legend_data,
  ggplot2::aes(
    x = 1,
    y = 1,
    colour = Environment,
    shape = Compartment
  )
) +
  ggplot2::geom_point(size = 3.0) +
  ggplot2::scale_colour_manual(
    values = environment_colours,
    breaks = environment_levels,
    labels = c(
      "Eyrewell Forest",
      "Kowhai (Irrigated)",
      "Kowhai (Rainfed)",
      "LU_H8",
      "Rolleston",
      "West Coast"
    ),
    drop = FALSE
  ) +
  ggplot2::scale_shape_manual(
    values = c(Rhizosphere = 16, Root = 17),
    drop = FALSE
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      title = "Field environment",
      order = 1,
      override.aes = list(shape = 16, size = 3.0)
    ),
    shape = ggplot2::guide_legend(
      title = "Compartment",
      order = 2,
      override.aes = list(colour = "black", size = 3.0)
    )
  ) +
  ggplot2::theme_void(base_size = 11) +
  ggplot2::theme(
    legend.position = "right",
    legend.justification = "center",
    legend.title = ggplot2::element_text(size = 11, face = "bold"),
    legend.text = ggplot2::element_text(size = 9.5, face = "bold"),
    legend.key.height = grid::unit(4.8, "mm"),
    legend.key.width = grid::unit(4.5, "mm"),
    legend.spacing.y = grid::unit(1.5, "mm"),
    legend.margin = ggplot2::margin(0, 3, 0, 0, unit = "mm")
  )

p_shared_legend <- cowplot::get_legend(legend_plot)

lower_row <- cowplot::plot_grid(
  pB,
  pC,
  p_shared_legend,
  nrow = 1,
  # Panels B and C use equal layout widths and retain equal height.
  rel_widths = c(1.10, 1.10, 0.72),
  align = "h",
  axis = "tb"
)

fig1 <- cowplot::plot_grid(
  pA,
  lower_row,
  ncol = 1,
  rel_heights = c(1.55, 1.0)
)

output_stem <- "Figure1_NatCom_merged_v24"
png_output <- file.path(out_dir, paste0(output_stem, ".png"))
pdf_output <- file.path(out_dir, paste0(output_stem, ".pdf"))
rds_output <- file.path(out_dir, paste0(output_stem, ".rds"))
components_rds_output <- file.path(
  out_dir,
  paste0(output_stem, "_components.rds")
)

# Save both the complete layout and its components for later modification.
# The B and C source pipelines also save their original ggplot objects as RDS,
# so axes, labels, scales and themes can be edited without rebuilding the data.
saveRDS(fig1, rds_output)
saveRDS(
  list(
    figure = fig1,
    panel_A = pA,
    panel_B = pB,
    panel_C = pC,
    panel_B_source = pB_source,
    panel_C_source = pC_source,
    shared_legend = p_shared_legend,
    source_files = list(
      biomass_script = biomass_script,
      biomass = biomass_file,
      taxonomy_rds = tax_rds_file,
      pathway_rds = path_rds_file
    ),
    layout = list(
      lower_row_relative_widths = c(1.10, 1.10, 0.72),
      figure_relative_heights = c(1.55, 1.0),
      export_width_mm = 250,
      export_height_mm = 225
    )
  ),
  components_rds_output
)

ggplot2::ggsave(
  png_output,
  fig1,
  width = 250,
  height = 225,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  pdf_output,
  fig1,
  width = 250,
  height = 225,
  units = "mm",
  useDingbats = FALSE,
  bg = "white"
)

copy_to_main_figures(c(pdf_output, png_output, rds_output, components_rds_output))

message("Saved merged Figure 1 PNG: ", png_output)
message("Saved merged Figure 1 PDF: ", pdf_output)
message("Saved editable merged Figure 1 RDS: ", rds_output)
message("Saved editable Figure 1 components RDS: ", components_rds_output)
message("Copied Figure 1 outputs to: ", main_fig_dir)
message("Done.")
