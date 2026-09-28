# Plant biomass analysis: one-factor and two-factor ANOVA
# Designed for Plant_Biomass(2).xlsx
#
# Experimental structure in the supplied workbook:
#   - between-block factor: Location (6 levels)
#   - within-block factor: Treatment (Control and Trichoderma)
#   - block/pair identifier: Replicate nested within Location (8 per location)
#   - response: Weight
#
# The primary analysis is a two-factor linear mixed model. It preserves the
# Control/Trichoderma pairing by including Location:Replicate as a random intercept.
# Site-specific one-factor tests are paired/repeated-measures ANOVAs. With two
# treatment levels, their F statistic equals the squared paired-t statistic.

# Run once if required:
# install.packages(c(
#   "readxl", "dplyr", "tidyr", "ggplot2", "lme4", "lmerTest",
#   "emmeans", "car"
# ))

required_packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "lme4", "lmerTest",
  "emmeans", "car"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(car)
})

# -----------------------------------------------------------------------------
# 1. User settings
# -----------------------------------------------------------------------------

input_directory <- paste0(
  "C:/Users/soths/OneDrive - Lincoln University/Writing/NatCom/Github/",
  "Biomass"
)

# The first name matches the newly supplied workbook. The second preserves
# compatibility if the file is renamed to remove the Windows download suffix.
input_candidates <- file.path(
  input_directory,
  c("Plant_Biomass(2).xlsx", "Plant_Biomass.xlsx")
)
input_file <- input_candidates[file.exists(input_candidates)][1]
sheet_name <- "Sheet1"
# Plotmath superscript avoids Windows graphics-device encoding warnings for ⁻.
y_axis_label <- expression("Plant biomass (kg ha"^{-1}*")")

if (length(input_file) == 0 || is.na(input_file)) {
  stop(
    "Input workbook not found. Expected one of: ",
    paste(input_candidates, collapse = "; ")
  )
}

output_dir <- file.path(dirname(input_file), "Plant_Biomass_ANOVA_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Sum-to-zero contrasts make Type III tests interpretable in the presence of
# the Location x Treatment interaction.
old_contrasts <- getOption("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))

# -----------------------------------------------------------------------------
# 2. Import and validate the experimental design
# -----------------------------------------------------------------------------

dat <- read_excel(input_file, sheet = sheet_name) |>
  as.data.frame()

required_columns <- c("Location", "Treatment", "Replicate", "Weight")
missing_columns <- setdiff(required_columns, names(dat))

if (length(missing_columns) > 0) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

if (anyNA(dat[required_columns])) {
  stop("Missing values were found in one or more required columns.")
}

if (!is.numeric(dat$Weight)) {
  stop("Weight must be numeric.")
}

dat <- dat |>
  mutate(
    Location = trimws(as.character(Location)),
    Treatment = ifelse(
      tolower(trimws(as.character(Treatment))) %in% c("panch", "trichoderma"),
      "Trichoderma", trimws(as.character(Treatment))
    )
  )

location_order <- c(
  "Eyrewell_Forest", "Kowhai_Irrigated", "Kowhai_Rainfed",
  "LU_H8", "Rolleston", "West_Coast"
)

unexpected_locations <- setdiff(unique(dat$Location), location_order)
if (length(unexpected_locations) > 0) {
  location_order <- c(location_order, sort(unexpected_locations))
}
location_order <- location_order[location_order %in% unique(dat$Location)]

dat <- dat |>
  mutate(
    Location = factor(Location, levels = location_order),
    Treatment = factor(Treatment, levels = c("Control", "Trichoderma")),
    Replicate = factor(Replicate),
    Block = interaction(Location, Replicate, drop = TRUE),
    Cell = interaction(Location, Treatment, drop = TRUE)
  )

if (anyNA(dat$Treatment)) {
  stop("Treatment must contain only 'Control' and 'Trichoderma'.")
}

duplicate_keys <- dat |>
  count(Location, Treatment, Replicate, name = "n") |>
  filter(n != 1)

if (nrow(duplicate_keys) > 0) {
  print(duplicate_keys)
  stop("Each Location x Treatment x Replicate combination must occur once.")
}

pair_check <- dat |>
  count(Block, Treatment) |>
  tidyr::complete(Block, Treatment, fill = list(n = 0))

if (any(pair_check$n != 1)) {
  print(pair_check |> filter(n != 1))
  stop("Every block must contain exactly one observation for each treatment.")
}

write.csv(dat, file.path(output_dir, "01_clean_data.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Descriptive statistics
# -----------------------------------------------------------------------------

descriptive <- dat |>
  group_by(Location, Treatment) |>
  summarise(
    n = n(),
    mean = mean(Weight),
    sd = sd(Weight),
    se = sd / sqrt(n),
    ci95_lower = mean - qt(0.975, df = n - 1) * se,
    ci95_upper = mean + qt(0.975, df = n - 1) * se,
    median = median(Weight),
    .groups = "drop"
  )

write.csv(
  descriptive,
  file.path(output_dir, "02_descriptive_statistics.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. Primary analysis: two-factor mixed-model ANOVA
# -----------------------------------------------------------------------------

# Location and Treatment are fixed effects. Block is Replicate nested within
# Location and is random, preserving the paired Control/Trichoderma observations.
fit_mixed <- lmer(
  Weight ~ Location * Treatment + (1 | Block),
  data = dat,
  REML = TRUE
)

anova_raw <- anova(
  fit_mixed,
  type = 3,
  ddf = "Satterthwaite"
)

anova_table <- data.frame(
  Term = rownames(anova_raw),
  anova_raw,
  row.names = NULL,
  check.names = FALSE
)

# Approximate partial eta-squared calculated from F and its degrees of freedom.
anova_table$partial_eta_squared <- with(
  anova_table,
  (`F value` * NumDF) / ((`F value` * NumDF) + DenDF)
)

write.csv(
  anova_table,
  file.path(output_dir, "03_two_way_mixed_ANOVA.csv"),
  row.names = FALSE
)

capture.output(
  summary(fit_mixed),
  file = file.path(output_dir, "03a_mixed_model_full_summary.txt")
)

writeLines(
  paste("Singular fit:", lme4::isSingular(fit_mixed, tol = 1e-5)),
  con = file.path(output_dir, "03b_singularity_check.txt")
)

# Classical repeated-measures ANOVA as a transparent cross-check.
# Location is between blocks; Treatment is within blocks.
fit_rm_aov <- aov(
  Weight ~ Location * Treatment + Error(Block / Treatment),
  data = dat
)

capture.output(
  summary(fit_rm_aov),
  file = file.path(output_dir, "04_two_way_repeated_measures_ANOVA.txt")
)

# -----------------------------------------------------------------------------
# 5. Planned post-hoc comparisons
# -----------------------------------------------------------------------------

# Primary simple effects: Trichoderma versus Control within each location.
# Holm correction controls family-wise error across all locations.
treatment_by_location <- emmeans(fit_mixed, ~ Treatment | Location)
treatment_contrasts <- pairs(
  treatment_by_location,
  reverse = TRUE,
  adjust = "holm"
) |>
  summary(infer = c(TRUE, TRUE)) |>
  as.data.frame()

write.csv(
  treatment_contrasts,
  file.path(output_dir, "05_treatment_within_location_Holm.csv"),
  row.names = FALSE
)

# Secondary comparisons among locations within each treatment.
# Tukey correction is appropriate for all pairwise location comparisons.
location_by_treatment <- emmeans(fit_mixed, ~ Location | Treatment)
location_contrasts <- pairs(
  location_by_treatment,
  adjust = "tukey"
) |>
  summary(infer = c(TRUE, TRUE)) |>
  as.data.frame()

write.csv(
  location_contrasts,
  file.path(output_dir, "06_location_within_treatment_Tukey.csv"),
  row.names = FALSE
)

estimated_means <- emmeans(fit_mixed, ~ Location * Treatment) |>
  summary(infer = c(TRUE, TRUE)) |>
  as.data.frame()

write.csv(
  estimated_means,
  file.path(output_dir, "07_estimated_marginal_means.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 6. Site-specific one-factor repeated-measures ANOVAs
# -----------------------------------------------------------------------------

# With exactly two treatment levels, a paired t test and a one-factor
# repeated-measures ANOVA are identical: F(1, n - 1) = t^2.
paired_wide <- dat |>
  select(Location, Replicate, Treatment, Weight) |>
  pivot_wider(names_from = Treatment, values_from = Weight)

one_way_results <- paired_wide |>
  group_by(Location) |>
  group_modify(function(.x, .y) {
    test <- t.test(.x$Trichoderma, .x$Control, paired = TRUE)
    difference <- .x$Trichoderma - .x$Control

    data.frame(
      n_pairs = length(difference),
      control_mean = mean(.x$Control),
      control_sd = sd(.x$Control),
      trichoderma_mean = mean(.x$Trichoderma),
      trichoderma_sd = sd(.x$Trichoderma),
      mean_difference_Trichoderma_minus_Control = mean(difference),
      # Percent increase between treatment means; the displayed SD values
      # describe biomass in kg ha^-1, not percentages.
      percent_change_from_Control = 100 * mean(difference) / mean(.x$Control),
      ci95_lower = unname(test$conf.int[1]),
      ci95_upper = unname(test$conf.int[2]),
      F_value = unname(test$statistic)^2,
      df1 = 1,
      df2 = unname(test$parameter),
      p_value = test$p.value,
      Cohen_dz = mean(difference) / sd(difference)
    )
  }) |>
  ungroup() |>
  mutate(p_value_Holm = p.adjust(p_value, method = "holm"))

write.csv(
  one_way_results,
  file.path(output_dir, "08_one_way_ANOVA_within_each_location.csv"),
  row.names = FALSE
)

# Save the conventional aov output for each location as an additional record.
one_way_models <- lapply(levels(dat$Location), function(site) {
  site_data <- droplevels(filter(dat, Location == site))
  aov(Weight ~ Treatment + Error(Replicate), data = site_data)
})
names(one_way_models) <- levels(dat$Location)

capture.output(
  lapply(one_way_models, summary),
  file = file.path(output_dir, "08a_one_way_ANOVA_full_output.txt")
)

# -----------------------------------------------------------------------------
# 7. Model diagnostics
# -----------------------------------------------------------------------------

diagnostic_data <- dat |>
  mutate(
    fitted_value = fitted(fit_mixed),
    residual = resid(fit_mixed),
    standardised_residual = as.numeric(scale(residual))
  )

shapiro_result <- shapiro.test(diagnostic_data$residual)
levene_result <- car::leveneTest(
  residual ~ Cell,
  data = diagnostic_data,
  center = median
)

capture.output(
  list(
    residual_Shapiro_Wilk = shapiro_result,
    residual_Levene_test_across_cells = levene_result
  ),
  file = file.path(output_dir, "09_diagnostic_tests.txt")
)

p_residual <- ggplot(
  diagnostic_data,
  aes(x = fitted_value, y = standardised_residual, colour = Location)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey45") +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = "Fitted value", y = "Standardised residual", colour = "Location") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

p_qq <- ggplot(diagnostic_data, aes(sample = residual)) +
  stat_qq(size = 2, alpha = 0.8) +
  stat_qq_line(linewidth = 0.6, colour = "#0072B2") +
  labs(x = "Theoretical quantile", y = "Residual quantile") +
  theme_classic(base_size = 11)

ggsave(
  file.path(output_dir, "09a_residuals_vs_fitted.pdf"),
  p_residual,
  width = 7.2,
  height = 5.0,
  units = "in"
)

ggsave(
  file.path(output_dir, "09b_residual_QQ.pdf"),
  p_qq,
  width = 5.5,
  height = 5.0,
  units = "in"
)

# -----------------------------------------------------------------------------
# 8. Publication-quality paired-data figure
# -----------------------------------------------------------------------------

treatment_colours <- c(Control = "#4D4D4D", Trichoderma = "#0072B2")

location_labels <- c(
  Kowhai_Rainfed = "Kowhai Rainfed",
  Kowhai_Irrigated = "Kowhai Irrigated",
  Rolleston = "Rolleston",
  Eyrewell_Forest = "Eyrewell Forest",
  West_Coast = "West Coast",
  LU_H8 = "LU_H8"
)

# Two lines per field: treatment mean ± sample SD (kg ha^-1), with the
# Trichoderma percentage calculated from the two treatment means.
biomass_annotations <- one_way_results |>
  transmute(
    Location,
    x = 1.5,
    y = max(dat$Weight) * 1.12,
    label = sprintf(
      "atop(bold('Control: %s ± %s'), bolditalic(Trichoderma)~bold(': %s ± %s (%+.1f%%)'))",
      formatC(control_mean, format = "f", digits = 0, big.mark = ","),
      formatC(control_sd, format = "f", digits = 0, big.mark = ","),
      formatC(trichoderma_mean, format = "f", digits = 0, big.mark = ","),
      formatC(trichoderma_sd, format = "f", digits = 0, big.mark = ","),
      percent_change_from_Control
    )
  )

p_biomass <- ggplot(
  dat,
  aes(x = Treatment, y = Weight, group = Block)
) +
  geom_line(colour = "grey65", linewidth = 0.70, alpha = 0.85) +
  geom_point(
    aes(colour = Treatment),
    size = 2.7,
    alpha = 0.92
  ) +
  geom_point(
    data = descriptive,
    aes(x = Treatment, y = mean),
    inherit.aes = FALSE,
    shape = 23,
    size = 4.2,
    stroke = 1.0,
    colour = "black",
    fill = "white"
  ) +
  geom_text(
    data = biomass_annotations,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.7,
    parse = TRUE
  ) +
  facet_wrap(
    ~ Location,
    ncol = 3,
    nrow = 2,
    labeller = as_labeller(location_labels)
  ) +
  scale_colour_manual(values = treatment_colours) +
  scale_x_discrete(
    labels = function(x) parse(text = ifelse(
      x == "Trichoderma", "italic(Trichoderma)", "'Control'"
    ))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.22))) +
  labs(
    x = NULL,
    y = y_axis_label,
    colour = "Treatment"
  ) +
  theme_classic(base_size = 13, base_family = "sans") +
  theme(
  strip.background = element_blank(),
  strip.text = element_text(face = "bold", size = 14),
  axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
  axis.text.y = element_text(size = 11),
  axis.title.y = element_text(size = 13, face = "bold"),
  panel.spacing = grid::unit(1.0, "lines"),
  plot.margin = margin(10, 10, 10, 10),
  panel.background = element_rect(fill = "white", colour = NA),
  plot.background = element_rect(fill = "white", colour = NA),
  legend.position = "none"
)

ggsave(
  file.path(output_dir, "10_plant_biomass_paired_plot.pdf"),
  p_biomass,
  width = 10,
  height = 5.4,
  units = "in",
  bg = "white"
)

ggsave(
  file.path(output_dir, "10_plant_biomass_paired_plot_600dpi.tiff"),
  p_biomass,
  width = 10,
  height = 5.4,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# High-resolution white-background PNG for combining with other figure panels.
ggsave(
  file.path(output_dir, "10_plant_biomass_panel_600dpi.png"),
  p_biomass,
  width = 10,
  height = 5.4,
  units = "in",
  dpi = 600,
  bg = "white"
)

# -----------------------------------------------------------------------------
# 9. Reproducibility record and concise console output
# -----------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "11_session_info.txt")
)

options(contrasts = old_contrasts)

cat("\nTwo-way mixed-model ANOVA\n")
print(anova_table)

cat("\nTrichoderma versus Control within each location (Holm-adjusted)\n")
print(treatment_contrasts)

cat("\nSite-specific one-factor repeated-measures ANOVAs\n")
print(one_way_results)

cat("\nResults saved to:", normalizePath(output_dir), "\n")
