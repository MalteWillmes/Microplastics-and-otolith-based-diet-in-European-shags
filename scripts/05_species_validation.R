# Compare model predictions with expert species identifications
# Comparison is based on numbers of otoliths per site and category

rm(list = ls())

library(tidyverse)
library(readxl)
library(here)

validation_seed <- 20261002
set.seed(validation_seed)

# Settings ----------------------------------------------------------------

selected_sites <- tibble(
  site_key = c(
    "F2R1",
    "F3R1",
    "F3R6",
    "F3R8",
    "F4R8",
    "F4R10"
  ),
  site = c(
    "F2-R1",
    "F3-R1",
    "F3-R6",
    "F3-R8",
    "F4-R8",
    "F4-R10"
  )
)

category_key_levels <- c(
  "cod",
  "saithe",
  "haddock",
  "poor_cod",
  "norway_pout",
  "gadoid_unknown",
  "lesser_sandeel",
  "wrasse",
  "other",
  "unknown"
)

category_labels <- c(
  "cod" = "Cod",
  "saithe" = "Saithe",
  "haddock" = "Haddock",
  "poor_cod" = "Poor Cod",
  "norway_pout" = "Norway Pout",
  "gadoid_unknown" = "Unknown gadoid",
  "lesser_sandeel" = "Lesser Sandeel",
  "wrasse" = "Wrasse sp.",
  "other" = "Other species",
  "unknown" = "Unknown"
)

category_levels <- unname(
  category_labels[
    category_key_levels
  ]
)

model_category_lookup <- c(
  "cod" = "cod",
  "saithe" = "saithe",
  "haddock" = "haddock",
  "poor_cod" = "poor_cod",
  "norway_pout" = "norway_pout",
  "gadoids_unknown" = "gadoid_unknown",
  "lesser_sand_eel" = "lesser_sandeel",
  "wrasse" = "wrasse"
)

validation_colors <- c(
  "Expert" = "#0072B2",
  "Model" = "#E69F00"
)

category_colors <- c(
  "Cod" = "#0072B2",
  "Saithe" = "#E69F00",
  "Haddock" = "#009E73",
  "Poor Cod" = "#CC79A7",
  "Norway Pout" = "#56B4E9",
  "Unknown gadoid" = "#999999",
  "Lesser Sandeel" = "#F0E442",
  "Wrasse sp." = "#D55E00",
  "Other species" = "#8C510A",
  "Unknown" = "#222222"
)

# Standardize site names --------------------------------------------------

clean_site <- function(x) {
  
  site_key <- x %>%
    as.character() %>%
    str_to_upper() %>%
    str_replace_all(
      "[^A-Z0-9]",
      ""
    )
  
  case_when(
    site_key == "R6F3" ~ "F3R6",
    site_key == "R8F4" ~ "F4R8",
    TRUE ~ site_key
  )
}

# Integer axis breaks -----------------------------------------------------

integer_breaks <- function(x) {
  breaks <- pretty(x)
  breaks[
    breaks >= 0 &
      abs(breaks - round(breaks)) < 1e-8
  ]
}

# Safe correlation --------------------------------------------------------

safe_cor <- function(
    x,
    y,
    method = "spearman"
) {
  
  keep <- is.finite(x) &
    is.finite(y)
  
  x <- x[keep]
  y <- y[keep]
  
  if (
    length(x) < 3 ||
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    return(
      NA_real_
    )
  }
  
  cor(
    x,
    y,
    method = method
  )
}

# Read data ---------------------------------------------------------------

predicted_data <- read_csv(
  here(
    "outputs",
    "predicted_data.csv"
  ),
  show_col_types = FALSE
)

expert_file <- here(
  "data",
  "expert_reader_compare.xlsx"
)

expert_data <- read_excel(
  expert_file,
  sheet = "data"
)

# Model predictions -------------------------------------------------------

model_data <- predicted_data %>%
  filter(
    prediction_valid,
    !is.na(predicted_species),
    !is.na(sample_group)
  ) %>%
  mutate(
    site_key = clean_site(
      sample_group
    ),
    predicted_species = as.character(
      predicted_species
    ),
    category_key = unname(
      model_category_lookup[
        predicted_species
      ]
    )
  )

unexpected_model_species <- model_data %>%
  filter(
    is.na(category_key)
  ) %>%
  distinct(
    predicted_species
  ) %>%
  pull(
    predicted_species
  )

if (length(unexpected_model_species) > 0) {
  stop(
    "Unexpected model species found: ",
    paste(
      unexpected_model_species,
      collapse = ", "
    )
  )
}

model_data <- model_data %>%
  inner_join(
    selected_sites,
    by = "site_key"
  ) %>%
  mutate(
    category = factor(
      unname(
        category_labels[
          category_key
        ]
      ),
      levels = category_levels
    )
  )

cat("\nModel sites included:\n")

print(
  model_data %>%
    count(
      site,
      name = "n_model"
    ),
  n = Inf
)

cat("\nModel species to validation categories:\n")

model_category_mapping <- model_data %>%
  count(
    predicted_species,
    category_key,
    category,
    name = "n"
  ) %>%
  arrange(
    category,
    predicted_species
  )

print(
  model_category_mapping,
  n = Inf
)

write_csv(
  model_category_mapping,
  here(
    "outputs",
    "model_validation_category_mapping.csv"
  )
)

model_counts <- model_data %>%
  count(
    site,
    category,
    name = "model_n"
  )

# Expert identifications --------------------------------------------------
# The organized data sheet contains the authoritative Category column.
# Species is retained only as the original expert description for auditing.

required_expert_columns <- c(
  "ID",
  "Number of ol",
  "Species",
  "Category"
)

missing_expert_columns <- setdiff(
  required_expert_columns,
  names(expert_data)
)

if (length(missing_expert_columns) > 0) {
  stop(
    "Missing required columns in expert data sheet: ",
    paste(
      missing_expert_columns,
      collapse = ", "
    )
  )
}

expert_validation <- expert_data %>%
  transmute(
    site_key = clean_site(
      ID
    ),
    expert_species_raw = str_squish(
      as.character(
        Species
      )
    ),
    category_key = str_to_lower(
      str_squish(
        as.character(
          Category
        )
      )
    ),
    n_raw = as.character(
      `Number of ol`
    ),
    n = suppressWarnings(
      as.numeric(
        as.character(
          `Number of ol`
        )
      )
    )
  )

unexpected_expert_categories <- expert_validation %>%
  filter(
    !is.na(category_key),
    category_key != "",
    !category_key %in% category_key_levels
  ) %>%
  distinct(
    category_key
  ) %>%
  pull(
    category_key
  )

if (length(unexpected_expert_categories) > 0) {
  stop(
    "Unexpected Category values in expert data sheet: ",
    paste(
      unexpected_expert_categories,
      collapse = ", "
    )
  )
}

expert_validation <- expert_validation %>%
  inner_join(
    selected_sites,
    by = "site_key"
  ) %>%
  mutate(
    category = factor(
      unname(
        category_labels[
          category_key
        ]
      ),
      levels = category_levels
    )
  )

cat("\nExpert sites included:\n")

print(
  expert_validation %>%
    group_by(
      site
    ) %>%
    summarise(
      n_rows = n(),
      n_otoliths = sum(
        n,
        na.rm = TRUE
      ),
      .groups = "drop"
    ),
  n = Inf
)

cat("\nExpert categories from organized data sheet:\n")

expert_category_mapping <- expert_validation %>%
  group_by(
    expert_species_raw,
    category_key,
    category
  ) %>%
  summarise(
    n_rows = n(),
    n_otoliths = sum(
      n,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    category,
    expert_species_raw
  )

print(
  expert_category_mapping,
  n = Inf
)

write_csv(
  expert_category_mapping,
  here(
    "outputs",
    "expert_validation_category_mapping.csv"
  )
)

# Check expert counts -----------------------------------------------------

expert_unquantified <- expert_validation %>%
  filter(
    is.na(n)
  ) %>%
  select(
    site,
    expert_species_raw,
    category_key,
    n_raw,
    category
  )

write_csv(
  expert_unquantified,
  here(
    "outputs",
    "expert_validation_unquantified.csv"
  )
)

if (nrow(expert_unquantified) > 0) {
  print(
    expert_unquantified,
    n = Inf
  )
  stop(
    "One or more expert counts are not numeric. ",
    "Resolve these entries before running the agreement analysis."
  )
}

expert_missing_category <- expert_validation %>%
  filter(
    is.na(category_key) |
      category_key == "" |
      is.na(category)
  ) %>%
  select(
    site,
    expert_species_raw,
    category_key,
    n_raw
  )

if (nrow(expert_missing_category) > 0) {
  print(
    expert_missing_category,
    n = Inf
  )
  stop(
    "One or more expert rows have a missing Category value."
  )
}

# Expert counts -----------------------------------------------------------

expert_counts <- expert_validation %>%
  group_by(
    site,
    category
  ) %>%
  summarise(
    expert_n = sum(n),
    .groups = "drop"
  )

# Complete comparison -----------------------------------------------------

comparison <- expand_grid(
  site = selected_sites$site,
  category = factor(
    category_levels,
    levels = category_levels
  )
) %>%
  left_join(
    model_counts,
    by = c(
      "site",
      "category"
    )
  ) %>%
  left_join(
    expert_counts,
    by = c(
      "site",
      "category"
    )
  ) %>%
  mutate(
    model_n = replace_na(
      model_n,
      0L
    ),
    
    expert_n = replace_na(
      expert_n,
      0
    ),
    
    site = factor(
      site,
      levels = selected_sites$site
    ),
    
    category = factor(
      category,
      levels = category_levels
    )
  ) %>%
  arrange(
    site,
    category
  )

write_csv(
  comparison,
  here(
    "outputs",
    "expert_model_species_validation.csv"
  )
)

# Counts, proportions and differences ------------------------------------

agreement_data <- comparison %>%
  group_by(
    site
  ) %>%
  mutate(
    model_total = sum(
      model_n
    ),
    
    expert_total = sum(
      expert_n
    ),
    
    model_prop = if_else(
      model_total > 0,
      model_n / model_total,
      NA_real_
    ),
    
    expert_prop = if_else(
      expert_total > 0,
      expert_n / expert_total,
      NA_real_
    ),
    
    count_difference =
      model_n -
      expert_n,
    
    absolute_count_difference =
      abs(
        count_difference
      ),
    
    proportion_difference =
      model_prop -
      expert_prop,
    
    difference_pp =
      100 *
      proportion_difference,
    
    absolute_difference_pp =
      abs(
        difference_pp
      )
  ) %>%
  ungroup()

write_csv(
  agreement_data,
  here(
    "outputs",
    "expert_model_agreement_data.csv"
  )
)

# Site-level agreement ----------------------------------------------------

agreement_by_site <- agreement_data %>%
  group_by(
    site
  ) %>%
  summarise(
    expert_total = first(
      expert_total
    ),
    
    model_total = first(
      model_total
    ),
    
    total_count_difference =
      model_total -
      expert_total,
    
    mean_bias_count = mean(
      count_difference
    ),
    
    MAE_count = mean(
      absolute_count_difference
    ),
    
    RMSE_count = sqrt(
      mean(
        count_difference^2
      )
    ),
    
    MAE_pp = mean(
      absolute_difference_pp
    ),
    
    RMSE_pp = sqrt(
      mean(
        difference_pp^2
      )
    ),
    
    max_absolute_difference_pp = max(
      absolute_difference_pp
    ),
    
    exact_count_matches = sum(
      count_difference == 0
    ),
    
    bray_curtis_counts =
      sum(
        abs(
          model_n -
            expert_n
        )
      ) /
      sum(
        model_n +
          expert_n
      ),
    
    bray_curtis_composition =
      sum(
        abs(
          model_prop -
            expert_prop
        ),
        na.rm = TRUE
      ) /
      sum(
        model_prop +
          expert_prop,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

print(
  agreement_by_site,
  n = Inf
)

write_csv(
  agreement_by_site,
  here(
    "outputs",
    "expert_model_agreement_by_site.csv"
  )
)

# Category-level agreement ------------------------------------------------

agreement_by_category <- agreement_data %>%
  group_by(
    category
  ) %>%
  summarise(
    expert_total = sum(
      expert_n
    ),
    
    model_total = sum(
      model_n
    ),
    
    total_count_difference =
      model_total -
      expert_total,
    
    mean_bias_count = mean(
      count_difference
    ),
    
    MAE_count = mean(
      absolute_count_difference
    ),
    
    RMSE_count = sqrt(
      mean(
        count_difference^2
      )
    ),
    
    mean_bias_pp = mean(
      difference_pp
    ),
    
    MAE_pp = mean(
      absolute_difference_pp
    ),
    
    RMSE_pp = sqrt(
      mean(
        difference_pp^2
      )
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    category
  )

print(
  agreement_by_category,
  n = Inf
)

write_csv(
  agreement_by_category,
  here(
    "outputs",
    "expert_model_agreement_by_category.csv"
  )
)

# Overall agreement -------------------------------------------------------

agreement_overall <- agreement_data %>%
  summarise(
    n_sites = n_distinct(
      site
    ),
    
    n_site_category_comparisons =
      n(),
    
    expert_total = sum(
      expert_n
    ),
    
    model_total = sum(
      model_n
    ),
    
    total_count_difference =
      model_total -
      expert_total,
    
    mean_bias_count = mean(
      count_difference
    ),
    
    MAE_count = mean(
      absolute_count_difference
    ),
    
    RMSE_count = sqrt(
      mean(
        count_difference^2
      )
    ),
    
    exact_count_match_proportion = mean(
      count_difference == 0
    ),
    
    MAE_pp = mean(
      absolute_difference_pp
    ),
    
    RMSE_pp = sqrt(
      mean(
        difference_pp^2
      )
    ),
    
    max_absolute_difference_pp = max(
      absolute_difference_pp
    ),
    
    spearman_count = safe_cor(
      expert_n,
      model_n
    ),
    
    spearman_proportion = safe_cor(
      expert_prop,
      model_prop
    )
  ) %>%
  mutate(
    mean_site_bray_curtis_counts =
      mean(
        agreement_by_site$
          bray_curtis_counts,
        na.rm = TRUE
      ),
    
    mean_site_bray_curtis_composition =
      mean(
        agreement_by_site$
          bray_curtis_composition,
        na.rm = TRUE
      )
  )

print(
  agreement_overall
)

write_csv(
  agreement_overall,
  here(
    "outputs",
    "expert_model_agreement_overall.csv"
  )
)

# Site totals -------------------------------------------------------------

site_totals <- agreement_by_site %>%
  select(
    site,
    expert_total,
    model_total,
    total_count_difference
  )

write_csv(
  site_totals,
  here(
    "outputs",
    "expert_model_species_validation_site_totals.csv"
  )
)

# Side-by-side count plot -------------------------------------------------

plot_data <- comparison %>%
  select(
    site,
    category,
    Model = model_n,
    Expert = expert_n
  ) %>%
  pivot_longer(
    cols = c(
      Model,
      Expert
    ),
    names_to = "method",
    values_to = "n"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "Expert",
        "Model"
      )
    ),
    
    category = factor(
      category,
      levels = category_levels
    )
  )

p_validation <- ggplot(
  plot_data,
  aes(
    x = category,
    y = n,
    fill = method
  )
) +
  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  facet_wrap(
    ~site,
    nrow = 2,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = validation_colors
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = integer_breaks,
    expand = expansion(
      mult = c(
        0,
        0.05
      )
    )
  ) +
  labs(
    x = NULL,
    y = "Number of otoliths",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 12
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "top"
  )

p_validation

ggsave(
  here(
    "outputs",
    "expert_model_species_validation.png"
  ),
  p_validation,
  width = 13,
  height = 7
)

# 1:1 count agreement -----------------------------------------------------

p_count_agreement <- ggplot(
  agreement_data,
  aes(
    x = expert_n,
    y = model_n
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point(
    size = 2.5
  ) +
  geom_text(
    aes(
      label = site
    ),
    vjust = -0.7,
    size = 2.5,
    check_overlap = TRUE
  ) +
  facet_wrap(
    ~category,
    scales = "free"
  ) +
  scale_x_continuous(
    breaks = integer_breaks
  ) +
  scale_y_continuous(
    breaks = integer_breaks
  ) +
  labs(
    x = "Expert count",
    y = "Model count"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )

p_count_agreement

ggsave(
  here(
    "outputs",
    "expert_model_count_agreement.png"
  ),
  p_count_agreement,
  width = 12,
  height = 8
)

# 1:1 proportional agreement ---------------------------------------------

p_proportion_agreement <- ggplot(
  agreement_data,
  aes(
    x = expert_prop,
    y = model_prop
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point(
    size = 2.5
  ) +
  geom_text(
    aes(
      label = site
    ),
    vjust = -0.7,
    size = 2.5,
    check_overlap = TRUE
  ) +
  facet_wrap(
    ~category
  ) +
  scale_x_continuous(
    labels = scales::percent,
    limits = c(
      0,
      1
    )
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(
      0,
      1
    )
  ) +
  labs(
    x = "Expert proportion",
    y = "Model proportion"
  ) +
  coord_equal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )

p_proportion_agreement

ggsave(
  here(
    "outputs",
    "expert_model_proportion_agreement.png"
  ),
  p_proportion_agreement,
  width = 12,
  height = 8
)

# Paired species-composition bars ----------------------------------------

composition_plot_data <- agreement_data %>%
  select(
    site,
    category,
    Expert = expert_prop,
    Model = model_prop
  ) %>%
  pivot_longer(
    cols = c(
      Expert,
      Model
    ),
    names_to = "method",
    values_to = "proportion"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "Expert",
        "Model"
      )
    ),
    
    category = factor(
      category,
      levels = category_levels
    )
  )

write_csv(
  composition_plot_data,
  here(
    "outputs",
    "expert_model_composition_by_site.csv"
  )
)

p_composition <- ggplot(
  composition_plot_data,
  aes(
    x = method,
    y = proportion,
    fill = category
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.25,
    width = 0.72
  ) +
  facet_wrap(
    ~site,
    nrow = 2
  ) +
  scale_fill_manual(
    values = category_colors,
    breaks = category_levels,
    drop = FALSE,
    name = NULL
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(
      0,
      1
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  labs(
    x = NULL,
    y = "Species composition"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 12
    ),
    legend.position = "right"
  )

p_composition

ggsave(
  here(
    "outputs",
    "expert_model_composition_by_site.png"
  ),
  p_composition,
  width = 11,
  height = 7
)

# Difference heatmap ------------------------------------------------------

heatmap_data <- agreement_data %>%
  mutate(
    category = factor(
      category,
      levels = rev(
        category_levels
      )
    ),
    
    difference_label = sprintf(
      "%+.1f",
      difference_pp
    ),
    
    text_color = if_else(
      absolute_difference_pp >= 20,
      "white",
      "black"
    )
  )

p_difference_heatmap <- ggplot(
  heatmap_data,
  aes(
    x = site,
    y = category,
    fill = difference_pp
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = difference_label,
      color = text_color
    ),
    size = 3
  ) +
  scale_color_identity() +
  scale_fill_gradient2(
    low = "#5E3C99",
    mid = "white",
    high = "#E66101",
    midpoint = 0,
    name = "Model - expert\n(percentage points)"
  ) +
  labs(
    x = "Site",
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_difference_heatmap

ggsave(
  here(
    "outputs",
    "expert_model_difference_heatmap.png"
  ),
  p_difference_heatmap,
  width = 9,
  height = 6
)