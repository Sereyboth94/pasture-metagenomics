# ============================================================
# Figure 1 standalone pipeline
# Study sites map (Figure 1A) + merged Figure 1 panels
# Nature Communications style
#
# Repository location:
#   location/run_location_figure1_pipeline.R
#
# Run from repository root in RStudio:
#   source("location/run_location_figure1_pipeline.R")
#
# Or from command line:
#   Rscript location/run_location_figure1_pipeline.R
#
# Expected repository structure:
#   Github/
#   ├── location/
#   │   ├── run_location_figure1_pipeline.R
#   │   └── results/
#   ├── taxonomy/
#   │   └── results/
#   └── function/
#       └── results/
#
# Outputs:
#   location/results/Fig1A_Study_Sites_Map_with_Canterbury_Inset.pdf
#   location/results/Fig1A_Study_Sites_Map_with_Canterbury_Inset.png
#   location/results/Fig1A_Study_Sites_Map_with_Canterbury_Inset.svg
#   location/results/Figure1_NatCom_merged.pdf
#   location/results/Figure1_NatCom_merged.png
#
# Manuscript figure copies:
#   manuscript_figures/Main_Figures/Fig1A_Study_Sites_Map_with_Canterbury_Inset.pdf
#   manuscript_figures/Main_Figures/Fig1A_Study_Sites_Map_with_Canterbury_Inset.png
#   manuscript_figures/Main_Figures/Fig1A_Study_Sites_Map_with_Canterbury_Inset.svg
#   manuscript_figures/Main_Figures/Fig1B_Taxonomy_PCoA_NatCom_panel.png
#   manuscript_figures/Main_Figures/Fig1C_Pathway_PCoA_NatCom_panel.png
#   manuscript_figures/Main_Figures/Figure1_NatCom_merged.pdf
#   manuscript_figures/Main_Figures/Figure1_NatCom_merged.png
# ============================================================

# ----------------------------
# 0) Packages
# ----------------------------
required_pkgs <- c(
  "sf",
  "ggplot2",
  "rnaturalearth",
  "rnaturalearthdata",
  "ggrepel",
  "ggspatial",
  "cowplot",
  "magick",
  "grid"
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
  library(sf)
  library(ggplot2)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggrepel)
  library(ggspatial)
  library(cowplot)
  library(magick)
  library(grid)
})

# ----------------------------
# 1) Repository-relative paths
# ----------------------------
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = TRUE))
  }

  # Fallback for interactive RStudio use when source() does not expose script path
  return(normalizePath(file.path(getwd(), "location", "run_location_figure1_pipeline.R"),
                       winslash = "/", mustWork = FALSE))
}

script_path <- get_script_path()
location_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(location_dir, ".."), winslash = "/", mustWork = FALSE)

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
  # Copy matching PDF/PNG/SVG versions when they exist.
  base <- tools::file_path_sans_ext(primary_file)
  sidecars <- paste0(base, c(".pdf", ".png", ".svg"))
  copy_to_main_figures(sidecars)
}

message("Repository root: ", repo_root)
message("Location output directory: ", out_dir)
message("Main manuscript figure directory: ", main_fig_dir)

# ----------------------------
# 2) Study site coordinates
# ----------------------------
sites <- data.frame(
  Location = c("Eyrewell_Forest", "LU_H8", "Kowhai", "Rolleston", "West_Coast"),
  Label = c("Eyrewell Forest", "LU_H8", "Kowhai", "Rolleston", "West Coast"),
  Latitude = c(-43.45203, -43.64875, -43.63625, -43.58647, -42.29572),
  Longitude = c(172.33886, 172.45306, 172.46994, 172.35739, 171.53972)
)

loc_cols <- c(
  "Eyrewell_Forest" = "#D55E00",
  "Kowhai"          = "#E69F00",
  "LU_H8"           = "#009E73",
  "Rolleston"       = "#0072B2",
  "West_Coast"      = "#CC79A7"
)

# ----------------------------
# 3) Generate Figure 1A map
# ----------------------------
nz <- rnaturalearth::ne_countries(
  country = "New Zealand",
  scale = "medium",
  returnclass = "sf"
)

sites_sf <- sf::st_as_sf(
  sites,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

nz <- sf::st_transform(nz, 2193)
sites_sf <- sf::st_transform(sites_sf, 2193)

coords <- sf::st_coordinates(sites_sf)
sites_xy <- cbind(sf::st_drop_geometry(sites_sf), coords)

canterbury_box <- sf::st_as_sfc(
  sf::st_bbox(
    c(
      xmin = 172.26,
      xmax = 172.52,
      ymin = -43.69,
      ymax = -43.43
    ),
    crs = 4326
  )
)

canterbury_box <- sf::st_transform(sf::st_sf(geometry = canterbury_box), 2193)

theme_map <- function(base_size = 7) {
  ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.3),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

p_main <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = nz, fill = "grey96", colour = "grey50", linewidth = 0.25) +
  ggplot2::geom_sf(
    data = sites_sf[sites_sf$Location == "West_Coast", ],
    ggplot2::aes(colour = Location),
    size = 2.8
  ) +
  ggplot2::geom_text(
    data = subset(sites_xy, Location == "West_Coast"),
    ggplot2::aes(x = X - 25000, y = Y + 30000, label = Label, colour = Location),
    size = 1.8,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::geom_sf(
    data = canterbury_box,
    fill = NA,
    colour = "black",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  ggplot2::scale_colour_manual(values = loc_cols, guide = "none") +
  ggplot2::coord_sf(
    xlim = c(1240000, 1700000),
    ylim = c(5040000, 5490000),
    expand = FALSE
  ) +
  ggspatial::annotation_scale(
    location = "bl",
    width_hint = 0.28,
    text_cex = 0.55,
    line_width = 0.25,
    unit_category = "metric"
  ) +
  ggspatial::annotation_north_arrow(
    location = "br",
    which_north = "true",
    height = grid::unit(0.55, "cm"),
    width = grid::unit(0.55, "cm"),
    pad_x = grid::unit(0.15, "cm"),
    pad_y = grid::unit(0.15, "cm"),
    style = ggspatial::north_arrow_fancy_orienteering(
      line_width = 0.25,
      text_size = 5
    )
  ) +
  theme_map()

p_inset <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = nz, fill = "grey96", colour = "grey50", linewidth = 0.22) +
  ggplot2::geom_sf(
    data = sites_sf[sites_sf$Location != "West_Coast", ],
    ggplot2::aes(colour = Location),
    size = 2.2
  ) +
  ggrepel::geom_text_repel(
    data = subset(sites_xy, Location != "West_Coast"),
    ggplot2::aes(x = X, y = Y, label = Label, colour = Location),
    size = 1.5,
    fontface = "bold",
    segment.size = 0.12,
    box.padding = 0.10,
    point.padding = 0.08,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = loc_cols, guide = "none") +
  ggplot2::coord_sf(
    xlim = sf::st_bbox(canterbury_box)[c("xmin", "xmax")],
    ylim = sf::st_bbox(canterbury_box)[c("ymin", "ymax")],
    expand = FALSE
  ) +
  theme_map(base_size = 6)

final_map <- cowplot::plot_grid(
  p_main,
  p_inset,
  nrow = 1,
  rel_widths = c(3.6, 1.0),
  align = "h"
)

ggplot2::ggsave(
  file.path(out_dir, "Fig1A_Study_Sites_Map_with_Canterbury_Inset.pdf"),
  final_map,
  width = 150,
  height = 75,
  units = "mm",
  useDingbats = FALSE
)

ggplot2::ggsave(
  file.path(out_dir, "Fig1A_Study_Sites_Map_with_Canterbury_Inset.png"),
  final_map,
  width = 150,
  height = 75,
  units = "mm",
  dpi = 600
)

ggplot2::ggsave(
  file.path(out_dir, "Fig1A_Study_Sites_Map_with_Canterbury_Inset.svg"),
  final_map,
  width = 150,
  height = 75,
  units = "mm"
)

copy_to_main_figures(file.path(
  out_dir,
  c(
    "Fig1A_Study_Sites_Map_with_Canterbury_Inset.pdf",
    "Fig1A_Study_Sites_Map_with_Canterbury_Inset.png",
    "Fig1A_Study_Sites_Map_with_Canterbury_Inset.svg"
  )
))

message("Saved Figure 1A map outputs to: ", out_dir)
message("Copied Figure 1A map outputs to: ", main_fig_dir)

# ----------------------------
# 4) Merge Figure 1A, 1B, and 1C
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

map_file <- file.path(out_dir, "Fig1A_Study_Sites_Map_with_Canterbury_Inset.png")

tax_file <- find_first_existing(
  c(
    file.path(out_dir, "Fig1B_Taxonomy_PCoA_NatCom_panel.png"),
    file.path(repo_root, "location", "data", "Fig1B_Taxonomy_PCoA_NatCom_panel.png"),
    file.path(repo_root, "taxonomy", "results", "Fig1B_Taxonomy_PCoA_NatCom_panel.png"),
    file.path(repo_root, "taxonomy", "results", "Figure1_NatCom", "Fig1B_Taxonomy_PCoA_NatCom_panel.png"),
    file.path(repo_root, "taxonomy", "results", "outputs_taxonomy_both_sqrtBC", "Figure1_NatCom", "Fig1B_Taxonomy_PCoA_NatCom_panel.png")
  ),
  "Figure 1B taxonomy PCoA panel"
)

path_file <- find_first_existing(
  c(
    file.path(out_dir, "Fig1C_Pathway_PCoA_NatCom_panel.png"),
    file.path(repo_root, "location", "data", "Fig1C_Pathway_PCoA_NatCom_panel.png"),
    file.path(repo_root, "function", "results", "Fig1C_Pathway_PCoA_NatCom_panel.png"),
    file.path(repo_root, "function", "results", "Figure1_NatCom_Pathway_PCoA", "Fig1C_Pathway_PCoA_NatCom_panel.png")
  ),
  "Figure 1C pathway PCoA panel"
)

message("Using Figure 1A map: ", map_file)
message("Using Figure 1B taxonomy panel: ", tax_file)
message("Using Figure 1C pathway panel: ", path_file)

# Copy Figure 1B and 1C source panels into manuscript_figures/Main_Figures
copy_matching_sidecars(tax_file)
copy_matching_sidecars(path_file)

map_img  <- magick::image_read(map_file)
tax_img  <- magick::image_read(tax_file)
path_img <- magick::image_read(path_file)

map_grob  <- grid::rasterGrob(as.raster(map_img),  interpolate = TRUE)
tax_grob  <- grid::rasterGrob(as.raster(tax_img),  interpolate = TRUE)
path_grob <- grid::rasterGrob(as.raster(path_img), interpolate = TRUE)

pA <- cowplot::ggdraw() + cowplot::draw_grob(map_grob)
pB <- cowplot::ggdraw() + cowplot::draw_grob(tax_grob)
pC <- cowplot::ggdraw() + cowplot::draw_grob(path_grob)

fig1 <- cowplot::plot_grid(
  pA,
  cowplot::plot_grid(
    pB,
    pC,
    nrow = 1,
    labels = c("B", "C"),
    label_size = 14
  ),
  ncol = 1,
  labels = c("A", ""),
  label_size = 14,
  rel_heights = c(0.8, 1.2)
)

ggplot2::ggsave(
  file.path(out_dir, "Figure1_NatCom_merged.png"),
  fig1,
  width = 180,
  height = 180,
  units = "mm",
  dpi = 600
)

ggplot2::ggsave(
  file.path(out_dir, "Figure1_NatCom_merged.pdf"),
  fig1,
  width = 180,
  height = 180,
  units = "mm",
  useDingbats = FALSE
)

copy_to_main_figures(file.path(
  out_dir,
  c(
    "Figure1_NatCom_merged.pdf",
    "Figure1_NatCom_merged.png"
  )
))

message("Saved merged Figure 1 outputs to: ", out_dir)
message("Copied Figure 1 outputs to: ", main_fig_dir)
message("Done.")
