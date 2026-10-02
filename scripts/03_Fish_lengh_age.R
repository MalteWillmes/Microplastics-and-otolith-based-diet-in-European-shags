# Reconstruct fish sizes and age classes from predicted otolith species

rm(list = ls())

library(tidyverse)
library(here)

# Settings ----------------------------------------------------------------

analysis_seed <- 20261002
set.seed(analysis_seed)

analysis_species_levels <- c(
  "cod",
  "saithe",
  "haddock",
  "poor_cod",
  "norway_pout",
  "gadoids_unknown",
  "lesser_sand_eel",
  "wrasse",
  "other_species",
  "unknown"
)

species_labels <- c(
  "cod" = "Cod",
  "saithe" = "Saithe",
  "haddock" = "Haddock",
  "poor_cod" = "Poor Cod",
  "norway_pout" = "Norway Pout",
  "gadoids_unknown" = "Unknown gadoid",
  "lesser_sand_eel" = "Lesser Sandeel",
  "wrasse" = "Wrasse sp.",
  "other_species" = "Other species",
  "unknown" = "Unknown"
)

species_palette <- c(
  "cod" = "#0072B2",
  "saithe" = "#E69F00",
  "haddock" = "#009E73",
  "poor_cod" = "#CC79A7",
  "norway_pout" = "#56B4E9",
  "gadoids_unknown" = "#999999",
  "lesser_sand_eel" = "#F0E442",
  "wrasse" = "#D55E00",
  "other_species" = "#8C510A",
  "unknown" = "#222222"
)

plot_seed <- analysis_seed

# Explicit species-age order ---------------------------------------------

species_age_levels <- c(
  "cod_0",
  "cod_1",
  "cod_2",
  "cod_3",
  "cod",
  "saithe_0",
  "saithe_1",
  "saithe_2",
  "saithe",
  "haddock",
  "poor_cod_0",
  "poor_cod_1",
  "poor_cod",
  "norway_pout",
  "gadoids_unknown",
  "lesser_sand_eel",
  "wrasse",
  "other_species",
  "unknown"
)

# Read data ---------------------------------------------------------------

fish_data_readin <- read_csv(
  here(
    "outputs",
    "predicted_data.csv"
  ),
  show_col_types = FALSE
)

# Check species -----------------------------------------------------------

unexpected_species <- fish_data_readin %>%
  filter(
    prediction_valid,
    !is.na(predicted_species)
  ) %>%
  distinct(
    predicted_species
  ) %>%
  filter(
    !predicted_species %in%
      analysis_species_levels
  ) %>%
  pull(
    predicted_species
  )

if (length(unexpected_species) > 0) {
  stop(
    "Unexpected predicted species found: ",
    paste(
      unexpected_species,
      collapse = ", "
    )
  )
}

# Reconstruct fish length and age ----------------------------------------

fish_data <- fish_data_readin %>%
  select(
    any_of(
      c(
        "Fish_ID",
        "folder",
        "picname",
        "key",
        "cal"
      )
    ),
    dataset,
    category,
    sample_group,
    OL,
    OW,
    OA,
    OP,
    prediction_valid,
    any_of(
      c(
        "model_predicted_species",
        "uncertain_gadoid",
        "gadoid_probability"
      )
    ),
    predicted_species,
    top_probability,
    second_species,
    second_probability,
    probability_margin
  ) %>%
  mutate(
    species = as.character(
      predicted_species
    ),
    
    FL_mm = case_when(
      species == "cod" ~
        0.41 + OL * 22.44,
      
      species == "poor_cod" &
        OL > 0.5 ~
        -49.9 + 28.091 * OL,
      
      species == "haddock" ~
        8.785 * OL^1.38,
      
      species == "norway_pout" ~
        -42.6 + 29.522 * OL,
      
      species == "saithe" &
        OL < 5.5 ~
        -4.24 + 23.5 * OL,
      
      species == "saithe" &
        OL >= 5.5 ~
        8.97297 * OL^1.53,
      
      species == "lesser_sand_eel" ~
        (
          (18.76 + 45.75 * OL) +
            (-4.024 + 56.84 * OL)
        ) / 2,
      
      TRUE ~
        NA_real_
    ),
    
    age = case_when(
      species == "cod" &
        FL_mm < 150 ~
        0L,
      
      species == "cod" &
        FL_mm < 250 ~
        1L,
      
      species == "cod" &
        FL_mm < 300 ~
        2L,
      
      species == "cod" &
        FL_mm >= 300 ~
        3L,
      
      species == "saithe" &
        FL_mm < 120 ~
        0L,
      
      species == "saithe" &
        FL_mm < 250 ~
        1L,
      
      species == "saithe" &
        FL_mm >= 250 ~
        2L,
      
      species == "poor_cod" &
        FL_mm < 80 ~
        0L,
      
      species == "poor_cod" &
        FL_mm >= 80 ~
        1L,
      
      TRUE ~
        NA_integer_
    ),
    
    species = factor(
      species,
      levels = analysis_species_levels
    ),
    
    age_f = factor(
      age,
      levels = 0:3
    ),
    
    species_age = case_when(
      is.na(species) ~
        NA_character_,
      
      species == "gadoids_unknown" ~
        "gadoids_unknown",
      
      is.na(age) ~
        as.character(species),
      
      TRUE ~
        paste0(
          species,
          "_",
          age
        )
    ),
    
    species_age = factor(
      species_age,
      levels = species_age_levels
    )
  ) %>%
  select(
    -predicted_species
  )

# Summary -----------------------------------------------------------------

length_summary <- fish_data %>%
  group_by(
    species,
    .drop = FALSE
  ) %>%
  summarise(
    n = n(),
    
    n_length_estimated = sum(
      !is.na(FL_mm)
    ),
    
    n_age_estimated = sum(
      !is.na(age)
    ),
    
    mean_FL_mm = mean(
      FL_mm,
      na.rm = TRUE
    ),
    
    min_FL_mm = min(
      FL_mm,
      na.rm = TRUE
    ),
    
    max_FL_mm = max(
      FL_mm,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    across(
      c(
        mean_FL_mm,
        min_FL_mm,
        max_FL_mm
      ),
      ~ replace(
        .x,
        !is.finite(.x),
        NA_real_
      )
    )
  ) %>%
  arrange(
    species
  )

print(
  length_summary
)

write_csv(
  length_summary,
  here(
    "outputs",
    "fish_length_age_summary.csv"
  )
)

# Unknown gadoid summary --------------------------------------------------

gadoid_unknown_summary <- fish_data %>%
  filter(
    prediction_valid
  ) %>%
  summarise(
    n_total = n(),
    
    n_gadoids_unknown = sum(
      species == "gadoids_unknown",
      na.rm = TRUE
    ),
    
    proportion_gadoids_unknown =
      n_gadoids_unknown /
      n_total
  )

cat(
  "\nUnknown gadoid classifications:\n"
)

print(
  gadoid_unknown_summary
)

write_csv(
  gadoid_unknown_summary,
  here(
    "outputs",
    "gadoids_unknown_summary.csv"
  )
)

# Unknown gadoids by category and sample ---------------------------------

gadoid_unknown_by_sample <- fish_data %>%
  filter(
    prediction_valid,
    species == "gadoids_unknown"
  ) %>%
  count(
    category,
    sample_group,
    name = "n_gadoids_unknown"
  ) %>%
  arrange(
    category,
    sample_group
  )

write_csv(
  gadoid_unknown_by_sample,
  here(
    "outputs",
    "gadoids_unknown_by_sample.csv"
  )
)

# Species-age colors and labels ------------------------------------------

species_age_lookup <- tibble(
  species_age = species_age_levels
) %>%
  mutate(
    species = case_when(
      str_detect(
        species_age,
        "^cod"
      ) ~
        "cod",
      
      str_detect(
        species_age,
        "^saithe"
      ) ~
        "saithe",
      
      str_detect(
        species_age,
        "^haddock"
      ) ~
        "haddock",
      
      str_detect(
        species_age,
        "^poor_cod"
      ) ~
        "poor_cod",
      
      str_detect(
        species_age,
        "^norway_pout"
      ) ~
        "norway_pout",
      
      species_age ==
        "gadoids_unknown" ~
        "gadoids_unknown",
      
      str_detect(
        species_age,
        "^lesser_sand_eel"
      ) ~
        "lesser_sand_eel",
      
      str_detect(
        species_age,
        "^wrasse"
      ) ~
        "wrasse",
      
      species_age ==
        "other_species" ~
        "other_species",
      
      species_age ==
        "unknown" ~
        "unknown"
    ),
    
    age = case_when(
      str_detect(
        species_age,
        "_[0-9]+$"
      ) ~
        as.integer(
          str_extract(
            species_age,
            "[0-9]+$"
          )
        ),
      
      TRUE ~
        NA_integer_
    ),
    
    color = unname(
      species_palette[
        species
      ]
    ),
    
    species_label = recode(
      species,
      !!!species_labels
    ),
    
    age_label = case_when(
      species == "poor_cod" &
        age == 1 ~
        "1+",
      
      species == "saithe" &
        age == 2 ~
        "2+",
      
      species == "cod" &
        age == 3 ~
        "3+",
      
      is.na(age) ~
        NA_character_,
      
      TRUE ~
        as.character(age)
    ),
    
    plot_label = if_else(
      is.na(age),
      species_label,
      paste(
        species_label,
        age_label
      )
    )
  )

species_age_palette <- set_names(
  species_age_lookup$color,
  species_age_lookup$species_age
)

species_age_label_vec <- set_names(
  species_age_lookup$plot_label,
  species_age_lookup$species_age
)

# Species-age count summary -----------------------------------------------

species_age_count_summary <- fish_data %>%
  filter(
    prediction_valid,
    !is.na(species_age)
  ) %>%
  count(
    species_age,
    name = "n_otoliths",
    .drop = FALSE
  ) %>%
  mutate(
    species_age = as.character(
      species_age
    )
  ) %>%
  left_join(
    species_age_lookup %>%
      select(
        species_age,
        species,
        species_label,
        age,
        age_label,
        plot_label
      ),
    by = "species_age"
  ) %>%
  filter(
    n_otoliths > 0
  ) %>%
  mutate(
    species_age = factor(
      species_age,
      levels = species_age_levels
    ),
    age_class = if_else(
      is.na(age_label),
      "Not assigned",
      age_label
    )
  ) %>%
  arrange(
    species_age
  ) %>%
  transmute(
    species = species_label,
    age_class,
    species_age = as.character(
      species_age
    ),
    n_otoliths
  )

cat(
  "\nSpecies-age counts:\n"
)

print(
  species_age_count_summary,
  n = Inf
)

write_csv(
  species_age_count_summary,
  here(
    "outputs",
    "species_age_count_summary.csv"
  )
)

# Length and age histogram ------------------------------------------------

p_age_size_histo <- fish_data %>%
  filter(
    !is.na(FL_mm)
  ) %>%
  ggplot(
    aes(
      x = FL_mm,
      fill = species
    )
  ) +
  geom_histogram(
    binwidth = 10,
    color = "black"
  ) +
  facet_wrap(
    ~species,
    scales = "free_y",
    labeller = as_labeller(
      species_labels
    )
  ) +
  scale_fill_manual(
    values = species_palette,
    breaks = analysis_species_levels,
    labels = species_labels,
    drop = FALSE
  ) +
  labs(
    x = "Estimated fish length (mm)",
    y = "Fish count"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.background = element_blank()
  )

p_age_size_histo

ggsave(
  here(
    "outputs",
    "age_size_histo.png"
  ),
  p_age_size_histo,
  width = 8,
  height = 4
)

# Pooled species-age proportions and Wilson CI ---------------------------

z <- qnorm(
  0.975
)

species_age_summary <- fish_data %>%
  filter(
    prediction_valid,
    !is.na(species_age),
    !is.na(category)
  ) %>%
  count(
    category,
    species_age,
    name = "n",
    .drop = FALSE
  ) %>%
  complete(
    category,
    species_age = factor(
      species_age_levels,
      levels = species_age_levels
    ),
    fill = list(
      n = 0
    )
  ) %>%
  group_by(
    category
  ) %>%
  mutate(
    N = sum(n),
    
    proportion =
      n / N,
    
    denominator =
      1 +
      z^2 / N,
    
    center =
      (
        proportion +
          z^2 / (2 * N)
      ) /
      denominator,
    
    half_width =
      z *
      sqrt(
        proportion *
          (1 - proportion) /
          N +
          z^2 /
          (4 * N^2)
      ) /
      denominator,
    
    CI_lower = pmax(
      0,
      center -
        half_width
    ),
    
    CI_upper = pmin(
      1,
      center +
        half_width
    )
  ) %>%
  ungroup() %>%
  mutate(
    category = recode(
      category,
      "No_MP" =
        "No Microplastics",
      "MP" =
        "Microplastics"
    ),
    
    category = factor(
      category,
      levels = c(
        "No Microplastics",
        "Microplastics"
      )
    ),
    
    species_age = factor(
      species_age,
      levels = species_age_levels
    )
  ) %>%
  arrange(
    category,
    species_age
  )

write_csv(
  species_age_summary,
  here(
    "outputs",
    "species_age_pooled_summary.csv"
  )
)

# Drop globally empty species-age classes from final figure ---------------

species_age_present <- species_age_summary %>%
  group_by(
    species_age
  ) %>%
  summarise(
    total_n = sum(n),
    .groups = "drop"
  ) %>%
  filter(
    total_n > 0
  ) %>%
  arrange(
    species_age
  ) %>%
  pull(
    species_age
  ) %>%
  as.character()

species_age_plot_data <- species_age_summary %>%
  filter(
    as.character(
      species_age
    ) %in%
      species_age_present
  ) %>%
  mutate(
    species_age = factor(
      as.character(
        species_age
      ),
      levels =
        species_age_present
    )
  )

# Final species-age figure ------------------------------------------------

p_species_age <- ggplot(
  species_age_plot_data,
  aes(
    x = species_age,
    y = proportion,
    fill = species_age
  )
) +
  geom_col(
    width = 0.75,
    color = "black",
    linewidth = 0.4
  ) +
  geom_errorbar(
    aes(
      ymin = CI_lower,
      ymax = CI_upper
    ),
    width = 0.2,
    linewidth = 0.5
  ) +
  facet_grid(
    . ~ category
  ) +
  scale_fill_manual(
    values =
      species_age_palette[
        species_age_present
      ],
    breaks =
      species_age_present,
    labels =
      species_age_label_vec[
        species_age_present
      ],
    drop = TRUE
  ) +
  scale_x_discrete(
    limits =
      species_age_present,
    labels =
      species_age_label_vec[
        species_age_present
      ],
    drop = TRUE
  ) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(
      mult = c(
        0,
        0.05
      )
    )
  ) +
  labs(
    x = "Species and age-class",
    y = "Proportion of otoliths"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 14
    ),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

p_species_age

ggsave(
  here(
    "outputs",
    "species_age_pooled_proportion.png"
  ),
  p_species_age,
  width = 10,
  height = 6
)

# Species-age composition by sample --------------------------------------

sample_group_age_data <- fish_data %>%
  filter(
    prediction_valid,
    !is.na(species_age)
  ) %>%
  count(
    category,
    sample_group,
    species_age,
    name = "n",
    .drop = FALSE
  ) %>%
  group_by(
    category,
    sample_group
  ) %>%
  mutate(
    proportion =
      n / sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    species_age = factor(
      species_age,
      levels = species_age_levels
    )
  )

p_predict_bar_group_age <- ggplot(
  sample_group_age_data,
  aes(
    x = sample_group,
    y = proportion,
    fill = species_age
  )
) +
  geom_col(
    color = "black",
    position = position_stack(
      reverse = TRUE
    )
  ) +
  facet_grid(
    cols = vars(
      category
    ),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = species_age_palette,
    breaks = species_age_levels,
    labels = species_age_label_vec,
    name = "Predicted species / age",
    drop = TRUE
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    x = "Sample group",
    y = "Proportion of otoliths"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_predict_bar_group_age

ggsave(
  here(
    "outputs",
    "p_predict_bar_group_age.png"
  ),
  p_predict_bar_group_age,
  width = 9,
  height = 4
)

# Length by predicted species --------------------------------------------

p_length_species <- fish_data %>%
  filter(
    !is.na(FL_mm)
  ) %>%
  ggplot(
    aes(
      x = species,
      y = FL_mm,
      fill = species
    )
  ) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    position = position_jitter(
      width = 0.15,
      seed = plot_seed
    ),
    alpha = 0.4,
    size = 1
  ) +
  scale_fill_manual(
    values = species_palette,
    breaks = analysis_species_levels,
    labels = species_labels,
    drop = FALSE
  ) +
  scale_x_discrete(
    limits = analysis_species_levels,
    labels = species_labels,
    drop = TRUE
  ) +
  labs(
    x = "Predicted species",
    y = "Estimated fish length (mm)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_length_species

ggsave(
  here(
    "outputs",
    "fish_length_by_species.png"
  ),
  p_length_species,
  width = 8,
  height = 5
)

# Export ------------------------------------------------------------------

write_csv(
  fish_data,
  here(
    "outputs",
    "fish_data.csv"
  )
)