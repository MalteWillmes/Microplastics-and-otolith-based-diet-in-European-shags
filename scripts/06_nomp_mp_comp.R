# Test species and age-class composition between MP and No-MP sites

rm(list = ls())

library(tidyverse)
library(here)

# Settings ----------------------------------------------------------------

analysis_seed <- 20261002
set.seed(analysis_seed)

species_levels <- c(
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

# Ordered analysis groups -------------------------------------------------

analysis_lookup <- tribble(
  ~analysis_id,          ~species,              ~analysis_type,   ~label,
  "cod_total",           "cod",                 "Species total",  "Cod",
  "cod_0",               "cod",                 "Age class",      "Cod 0",
  "cod_1",               "cod",                 "Age class",      "Cod 1",
  "cod_2",               "cod",                 "Age class",      "Cod 2",
  "cod_3",               "cod",                 "Age class",      "Cod 3+",
  
  "saithe_total",        "saithe",              "Species total",  "Saithe",
  "saithe_0",            "saithe",              "Age class",      "Saithe 0",
  "saithe_1",            "saithe",              "Age class",      "Saithe 1",
  "saithe_2",            "saithe",              "Age class",      "Saithe 2+",
  
  "haddock",             "haddock",             "Species total",  "Haddock",
  
  "poor_cod_total",      "poor_cod",            "Species total",  "Poor Cod",
  "poor_cod_0",          "poor_cod",            "Age class",      "Poor Cod 0",
  "poor_cod_1",          "poor_cod",            "Age class",      "Poor Cod 1+",
  
  "norway_pout",         "norway_pout",         "Species total",  "Norway Pout",
  "gadoids_unknown",     "gadoids_unknown",     "Species total",  "Unknown gadoid",
  "lesser_sand_eel",     "lesser_sand_eel",     "Species total",  "Lesser Sandeel",
  "wrasse",              "wrasse",              "Species total",  "Wrasse sp.",
  "other_species",       "other_species",       "Species total",  "Other species",
  "unknown",             "unknown",             "Species total",  "Unknown"
) %>%
  mutate(
    color = unname(
      species_palette[
        species
      ]
    )
  )

# Read data ---------------------------------------------------------------

fish_data <- read_csv(
  here(
    "outputs",
    "fish_data.csv"
  ),
  show_col_types = FALSE
) %>%
  filter(
    prediction_valid,
    !is.na(species),
    !is.na(category),
    !is.na(sample_group)
  ) %>%
  mutate(
    species = as.character(
      species
    ),
    
    species_age = as.character(
      species_age
    ),
    
    category = factor(
      category,
      levels = c(
        "No_MP",
        "MP"
      )
    )
  )

if (
  !all(
    c(
      "No_MP",
      "MP"
    ) %in%
    fish_data$category
  )
) {
  stop(
    "Both No_MP and MP samples are required."
  )
}

# Check species -----------------------------------------------------------

unexpected_species <- setdiff(
  unique(
    fish_data$species
  ),
  species_levels
)

if (
  length(
    unexpected_species
  ) > 0
) {
  stop(
    "Unexpected species found: ",
    paste(
      unexpected_species,
      collapse = ", "
    )
  )
}

# Site totals -------------------------------------------------------------

site_totals <- fish_data %>%
  count(
    category,
    sample_group,
    name = "N"
  )

cat(
  "\nSites by treatment:\n"
)

print(
  site_totals %>%
    count(
      category,
      name = "n_sites"
    )
)

# Create species-total observations --------------------------------------

species_total_data <- fish_data %>%
  mutate(
    analysis_id = case_when(
      species == "cod" ~
        "cod_total",
      
      species == "poor_cod" ~
        "poor_cod_total",
      
      species == "saithe" ~
        "saithe_total",
      
      TRUE ~
        species
    )
  ) %>%
  select(
    category,
    sample_group,
    analysis_id
  )

# Create age-specific observations ---------------------------------------

age_class_data <- fish_data %>%
  filter(
    species %in%
      c(
        "cod",
        "poor_cod",
        "saithe"
      ),
    species_age %in%
      analysis_lookup$analysis_id
  ) %>%
  transmute(
    category,
    sample_group,
    analysis_id =
      species_age
  )

# Combine species totals and age classes ---------------------------------

analysis_data <- bind_rows(
  species_total_data,
  age_class_data
)

# Drop analysis groups absent from the entire dataset ---------------------

analysis_groups_present <- analysis_data %>%
  distinct(
    analysis_id
  ) %>%
  pull(
    analysis_id
  )

analysis_lookup <- analysis_lookup %>%
  filter(
    analysis_id %in%
      analysis_groups_present
  )

analysis_levels <- analysis_lookup$analysis_id

analysis_data <- analysis_data %>%
  mutate(
    analysis_id = factor(
      analysis_id,
      levels = analysis_levels
    )
  )

# Counts at each site -----------------------------------------------------

site_counts <- analysis_data %>%
  count(
    category,
    sample_group,
    analysis_id,
    name = "n",
    .drop = FALSE
  ) %>%
  complete(
    nesting(
      category,
      sample_group
    ),
    analysis_id = factor(
      analysis_levels,
      levels = analysis_levels
    ),
    fill = list(
      n = 0
    )
  ) %>%
  left_join(
    site_totals,
    by = c(
      "category",
      "sample_group"
    )
  ) %>%
  left_join(
    analysis_lookup %>%
      select(
        analysis_id,
        species,
        analysis_type,
        label
      ),
    by = "analysis_id"
  ) %>%
  arrange(
    category,
    sample_group,
    analysis_id
  )

write_csv(
  site_counts,
  here(
    "outputs",
    "species_composition_site_counts.csv"
  )
)

# Fit one quasibinomial GLM per analysis group ----------------------------

fit_composition <- function(dat) {
  
  if (sum(dat$n) == 0) {
    
    return(
      tibble(
        n_sites = nrow(dat),
        n_sites_present = 0,
        estimate = NA_real_,
        SE = NA_real_,
        odds_ratio = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  fit <- glm(
    cbind(
      n,
      N - n
    ) ~ category,
    family = quasibinomial(),
    data = dat
  )
  
  coef_table <- summary(
    fit
  )$coefficients
  
  if (
    !"categoryMP" %in%
    rownames(
      coef_table
    )
  ) {
    
    return(
      tibble(
        n_sites = nrow(dat),
        n_sites_present = sum(
          dat$n > 0
        ),
        estimate = NA_real_,
        SE = NA_real_,
        odds_ratio = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  estimate <- coef_table[
    "categoryMP",
    "Estimate"
  ]
  
  SE <- coef_table[
    "categoryMP",
    "Std. Error"
  ]
  
  p_value <- coef_table[
    "categoryMP",
    "Pr(>|t|)"
  ]
  
  critical_value <- qt(
    0.975,
    df = df.residual(
      fit
    )
  )
  
  tibble(
    n_sites = nrow(
      dat
    ),
    
    n_sites_present = sum(
      dat$n > 0
    ),
    
    estimate =
      estimate,
    
    SE =
      SE,
    
    odds_ratio = exp(
      estimate
    ),
    
    CI_lower = exp(
      estimate -
        critical_value *
        SE
    ),
    
    CI_upper = exp(
      estimate +
        critical_value *
        SE
    ),
    
    p_value =
      p_value
  )
}

glm_results <- site_counts %>%
  select(
    category,
    sample_group,
    analysis_id,
    n,
    N
  ) %>%
  group_by(
    analysis_id
  ) %>%
  group_modify(
    ~ fit_composition(
      .x
    )
  ) %>%
  ungroup()

# Observed pooled proportions --------------------------------------------

pooled_proportions <- site_counts %>%
  group_by(
    analysis_id,
    category
  ) %>%
  summarise(
    n = sum(
      n
    ),
    
    N = sum(
      N
    ),
    
    proportion =
      n / N,
    
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from =
      category,
    
    values_from = c(
      n,
      N,
      proportion
    ),
    
    names_glue =
      "{.value}_{category}"
  )

# Final results -----------------------------------------------------------

glm_results <- glm_results %>%
  left_join(
    pooled_proportions,
    by = "analysis_id"
  ) %>%
  left_join(
    analysis_lookup,
    by = "analysis_id"
  ) %>%
  mutate(
    difference =
      proportion_MP -
      proportion_No_MP,
    
    difference_pp =
      100 *
      difference,
    
    p_adjusted = p.adjust(
      p_value,
      method = "BH"
    ),
    
    significant =
      p_adjusted < 0.05,
    
    direction = case_when(
      difference > 0 ~
        "Higher in MP",
      
      difference < 0 ~
        "Higher in No-MP",
      
      TRUE ~
        "No difference"
    ),
    
    analysis_id = factor(
      analysis_id,
      levels = analysis_levels
    )
  ) %>%
  select(
    analysis_id,
    label,
    species,
    analysis_type,
    n_sites,
    n_sites_present,
    n_No_MP,
    N_No_MP,
    proportion_No_MP,
    n_MP,
    N_MP,
    proportion_MP,
    difference,
    difference_pp,
    odds_ratio,
    CI_lower,
    CI_upper,
    estimate,
    SE,
    p_value,
    p_adjusted,
    significant,
    direction
  ) %>%
  arrange(
    analysis_id
  )

cat(
  "\nGLM results:\n"
)

print(
  glm_results,
  n = Inf
)

write_csv(
  glm_results,
  here(
    "outputs",
    "species_composition_glm_results.csv"
  )
)

# Also retain previous output filename ------------------------------------

write_csv(
  glm_results,
  here(
    "outputs",
    "species_age_glm_results.csv"
  )
)

# Species-total results ---------------------------------------------------

species_total_results <- glm_results %>%
  filter(
    analysis_type ==
      "Species total"
  )

cat(
  "\nSpecies-total results:\n"
)

print(
  species_total_results,
  n = Inf
)

write_csv(
  species_total_results,
  here(
    "outputs",
    "species_total_glm_results.csv"
  )
)

# Age-class results -------------------------------------------------------

age_class_results <- glm_results %>%
  filter(
    analysis_type ==
      "Age class"
  )

cat(
  "\nAge-class results:\n"
)

print(
  age_class_results,
  n = Inf
)

write_csv(
  age_class_results,
  here(
    "outputs",
    "species_age_class_glm_results.csv"
  )
)

# Significant results -----------------------------------------------------

cat(
  "\nSpecies and age-class groups differing between MP and No-MP sites",
  "\nafter BH correction across all tests:\n\n"
)

significant_results <- glm_results %>%
  filter(
    significant
  )

if (
  nrow(
    significant_results
  ) == 0
) {
  
  cat(
    "None at adjusted P < 0.05.\n"
  )
  
} else {
  
  print(
    significant_results,
    n = Inf
  )
}

# Effect-size plot --------------------------------------------------------

plot_lookup <- analysis_lookup %>%
  filter(
    analysis_id %in%
      glm_results$analysis_id
  )

plot_data <- glm_results %>%
  filter(
    is.finite(
      odds_ratio
    ),
    is.finite(
      CI_lower
    ),
    is.finite(
      CI_upper
    ),
    odds_ratio > 0,
    CI_lower > 0,
    CI_upper > 0
  ) %>%
  mutate(
    label = factor(
      label,
      levels = rev(
        plot_lookup$label
      )
    ),
    
    species = factor(
      species,
      levels = species_levels
    )
  )

p_glm <- ggplot(
  plot_data,
  aes(
    x = odds_ratio,
    y = label,
    color = species
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = 2,
    color = "grey50"
  ) +
  geom_errorbar(
    aes(
      xmin = CI_lower,
      xmax = CI_upper
    ),
    width = 0.2,
    orientation = "y",
    linewidth = 0.6
  ) +
  geom_point(
    aes(
      shape = analysis_type
    ),
    size = 3
  ) +
  scale_color_manual(
    values =
      species_palette,
    breaks =
      species_levels,
    labels =
      species_labels[species_levels],
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "Species total" = 16,
      "Age class" = 17
    )
  ) +
  scale_x_log10() +
  labs(
    x = "Odds ratio (MP vs No-MP)",
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid.minor =
      element_blank(),
    
    panel.grid.major.y =
      element_blank(),
    
    legend.position =
      "top"
  )

p_glm

ggsave(
  here(
    "outputs",
    "species_composition_glm_effects.png"
  ),
  p_glm,
  width = 7,
  height = 7
)