# Predict unknown otoliths to species
# Uses final random forest model from reference dataset

rm(list = ls())

library(tidyverse)
library(shapeR)
library(caret)

# Settings ----------------------------------------------------------------

prediction_seed <- 20261002
set.seed(prediction_seed)

model_file <- "outputs/model_eval/final_species_model.rds"
output_file <- "outputs/predicted_data.csv"

mp_file <- "shape_extraction/MP_otoliths/shapedata_mp_samples.RData"
no_mp_file <- "shape_extraction/noMP_otoliths/shapedata_no_mp_samples.RData"
no_mp_file2 <- "shape_extraction/noMP_otoliths_newsample/shapedata_no_mp_samples2.RData"

gadoid_threshold <- 0.50

display_species_order <- c(
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

display_species_labels <- c(
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

display_species_colors <- c(
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

model_species_order <- c(
  "cod",
  "saithe",
  "haddock",
  "poor_cod",
  "norway_pout",
  "lesser_sand_eel",
  "wrasse"
)

gadoid_species <- c(
  "cod",
  "saithe",
  "haddock",
  "poor_cod",
  "norway_pout"
)

prediction_species_levels <- display_species_order[
  display_species_order %in%
    c(
      model_species_order,
      "gadoids_unknown"
    )
]

prediction_species_labels <- display_species_labels[
  prediction_species_levels
]

prediction_species_colors <- display_species_colors[
  prediction_species_levels
]

# Read final model --------------------------------------------------------

model_bundle <- readRDS(
  model_file
)

final_model <- model_bundle$model
final_predictors <- model_bundle$predictors
model_name <- model_bundle$model_name
fitted_species_levels <- model_bundle$species_levels

unexpected_species <- setdiff(
  fitted_species_levels,
  model_species_order
)

missing_species <- setdiff(
  model_species_order,
  fitted_species_levels
)

if (length(unexpected_species) > 0) {
  stop(
    "Unexpected species found in the fitted model: ",
    paste(
      unexpected_species,
      collapse = ", "
    )
  )
}

if (length(missing_species) > 0) {
  stop(
    "Expected species missing from the fitted model: ",
    paste(
      missing_species,
      collapse = ", "
    )
  )
}

# Use the canonical order regardless of the order stored in the model.
species_levels <- model_species_order
expected_model_species <- model_species_order

gadoid_species <- intersect(
  gadoid_species,
  species_levels
)

cat(
  "\nModel:",
  model_name,
  "\n"
)

cat(
  "Predictor set:",
  model_bundle$predictor_set,
  "\n"
)

cat(
  "Class balancing:",
  model_bundle$balance,
  "\n"
)

cat(
  "Number of predictors:",
  length(final_predictors),
  "\n"
)

cat(
  "Species:",
  paste(
    species_levels,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Gadoids:",
  paste(
    gadoid_species,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Unknown-gadoid threshold:",
  gadoid_threshold,
  "\n\n"
)

print(
  model_bundle$best_tune
)

# Model information -------------------------------------------------------

model_info <- tibble(
  model = model_name,
  predictor_set = model_bundle$predictor_set,
  balance = model_bundle$balance,
  prediction_seed = prediction_seed,
  n_predictors = length(final_predictors),
  gadoid_threshold = gadoid_threshold,
  gadoid_species = paste(
    gadoid_species,
    collapse = "; "
  ),
  species = paste(
    species_levels,
    collapse = "; "
  )
)

write_csv(
  model_info,
  "outputs/prediction_model_info.csv"
)

# Read shapeR RData object ------------------------------------------------

read_shape_object <- function(path) {
  
  if (!file.exists(path)) {
    stop(
      "File not found: ",
      path
    )
  }
  
  env <- new.env()
  
  objects_loaded <- load(
    path,
    envir = env
  )
  
  if ("shape" %in% objects_loaded) {
    return(
      env$shape
    )
  }
  
  if (length(objects_loaded) == 1) {
    return(
      env[[objects_loaded]]
    )
  }
  
  shape_candidates <- objects_loaded[
    map_lgl(
      objects_loaded,
      ~ {
        x <- env[[.x]]
        
        isS4(x) &&
          "shape.coef.raw" %in%
          methods::slotNames(x)
      }
    )
  ]
  
  if (length(shape_candidates) == 1) {
    return(
      env[[shape_candidates]]
    )
  }
  
  stop(
    "Could not identify the shapeR object in: ",
    path,
    "\nObjects found: ",
    paste(
      objects_loaded,
      collapse = ", "
    )
  )
}

# Extract shapeR samples --------------------------------------------------

read_shape_samples <- function(
    path,
    category,
    sample_column,
    dataset
) {
  
  shape_object <- read_shape_object(
    path
  )
  
  shape_object <- setFilter(
    shape_object
  )
  
  master <- getMasterlist(
    shape_object
  ) %>%
    as_tibble()
  
  if (!sample_column %in% names(master)) {
    stop(
      "Column '",
      sample_column,
      "' not found in: ",
      path
    )
  }
  
  if (!"key" %in% names(master)) {
    
    if (
      !all(
        c(
          "folder",
          "picname"
        ) %in% names(master)
      )
    ) {
      stop(
        "Cannot construct sample key in: ",
        path
      )
    }
    
    master <- master %>%
      mutate(
        key = paste(
          folder,
          picname,
          sep = ";"
        )
      )
  }
  
  if (!"cal" %in% names(master)) {
    master <- master %>%
      mutate(
        cal = NA_real_
      )
  }
  
  # Shape coefficients ----------------------------------------------------
  
  wavelet_raw <- getWavelet(
    shape_object
  )
  
  fourier_raw <- getFourier(
    shape_object
  )
  
  stopifnot(
    nrow(master) == nrow(wavelet_raw),
    nrow(master) == nrow(fourier_raw)
  )
  
  wavelet_data <- wavelet_raw %>%
    as.data.frame() %>%
    as_tibble()
  
  names(wavelet_data) <- paste0(
    "wavelet_",
    make.names(
      names(wavelet_data),
      unique = TRUE
    )
  )
  
  fourier_data <- fourier_raw %>%
    as.data.frame() %>%
    as_tibble()
  
  names(fourier_data) <- paste0(
    "fourier_",
    make.names(
      names(fourier_data),
      unique = TRUE
    )
  )
  
  # Raw morphometrics -----------------------------------------------------
  
  raw_data <- shape_object@shape.coef.raw %>%
    as.data.frame()
  
  required_raw <- c(
    "otolith.area",
    "otolith.length",
    "otolith.width",
    "otolith.perimeter"
  )
  
  if (
    !all(
      required_raw %in%
      names(raw_data)
    )
  ) {
    stop(
      "Raw morphometric columns are missing from: ",
      path
    )
  }
  
  raw_data <- raw_data %>%
    select(
      all_of(
        required_raw
      )
    )
  
  # Match raw measurements to master list --------------------------------
  
  if (
    !is.null(
      rownames(raw_data)
    ) &&
    all(
      master$key %in%
      rownames(raw_data)
    )
  ) {
    
    raw_measurements <- raw_data %>%
      rownames_to_column(
        "key"
      ) %>%
      as_tibble() %>%
      rename(
        raw_OA = otolith.area,
        raw_OL = otolith.length,
        raw_OW = otolith.width,
        raw_OP = otolith.perimeter
      )
    
    sample_info <- master %>%
      left_join(
        raw_measurements,
        by = "key"
      )
    
  } else {
    
    if (
      nrow(raw_data) !=
      nrow(master)
    ) {
      stop(
        "Could not align raw measurements with master list in: ",
        path
      )
    }
    
    raw_measurements <- raw_data %>%
      as_tibble() %>%
      rename(
        raw_OA = otolith.area,
        raw_OL = otolith.length,
        raw_OW = otolith.width,
        raw_OP = otolith.perimeter
      )
    
    sample_info <- bind_cols(
      master,
      raw_measurements
    )
  }
  
  # Sample information ----------------------------------------------------
  
  sample_info <- sample_info %>%
    mutate(
      category = .env$category,
      
      sample_group = str_squish(
        as.character(
          .data[[sample_column]]
        )
      ),
      
      dataset = .env$dataset,
      
      OL = otolith.length,
      OW = otolith.width,
      OA = otolith.area,
      OP = otolith.perimeter,
      
      AR = raw_OL / raw_OW,
      
      C = raw_OP^2 / raw_OA,
      
      E = (
        raw_OL -
          raw_OW
      ) /
        (
          raw_OL +
            raw_OW
        ),
      
      FF = 4 *
        pi *
        raw_OA /
        raw_OP^2,
      
      RE = raw_OA /
        (
          raw_OL *
            raw_OW
        ),
      
      RO = 4 *
        raw_OA /
        (
          pi *
            raw_OL^2
        )
    ) %>%
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
      AR,
      C,
      E,
      FF,
      RE,
      RO
    )
  
  bind_cols(
    sample_info,
    wavelet_data,
    fourier_data
  )
}

# Read unknown samples ----------------------------------------------------

MP_readin <- read_shape_samples(
  path = mp_file,
  category = "MP",
  sample_column = "MP_sample",
  dataset = "MP"
)

no_MP_readin <- read_shape_samples(
  path = no_mp_file,
  category = "No_MP",
  sample_column = "noMP_sample",
  dataset = "No_MP_1"
)

no_MP_readin2 <- read_shape_samples(
  path = no_mp_file2,
  category = "No_MP",
  sample_column = "noMP_sample",
  dataset = "No_MP_2"
)

# Combine samples ---------------------------------------------------------

sample_data <- bind_rows(
  MP_readin,
  no_MP_readin,
  no_MP_readin2
)

site_summary <- sample_data %>%
  filter(
    !is.na(sample_group)
  ) %>%
  count(
    category,
    sample_group,
    name = "n_otoliths"
  ) %>%
  arrange(
    category,
    sample_group
  )

cat(
  "\nSites loaded:\n"
)

print(
  site_summary,
  n = Inf
)

write_csv(
  site_summary,
  "outputs/prediction_site_summary.csv"
)

write_csv(
  sample_data,
  "outputs/sample_data.csv"
)

cat(
  "\nTotal unknown otoliths:",
  nrow(sample_data),
  "\n"
)

# Calibration summary -----------------------------------------------------

calibration_summary <- sample_data %>%
  summarise(
    n = n(),
    
    n_with_calibration = sum(
      !is.na(cal)
    ),
    
    n_without_calibration = sum(
      is.na(cal)
    ),
    
    n_with_size = sum(
      if_all(
        all_of(
          c(
            "OL",
            "OW",
            "OA",
            "OP"
          )
        ),
        ~ !is.na(.x) &
          is.finite(.x)
      )
    ),
    
    n_with_shape_indices = sum(
      if_all(
        all_of(
          c(
            "AR",
            "FF",
            "RE",
            "RO"
          )
        ),
        ~ !is.na(.x) &
          is.finite(.x)
      )
    )
  )

cat(
  "\nCalibration and morphometric data:\n"
)

print(
  calibration_summary
)

# Check model predictors --------------------------------------------------

missing_predictors <- setdiff(
  final_predictors,
  names(sample_data)
)

if (length(missing_predictors) > 0) {
  stop(
    "The following predictors required by the model are missing:\n",
    paste(
      missing_predictors,
      collapse = "\n"
    )
  )
}

# Prepare prediction data -------------------------------------------------

prediction_data <- sample_data %>%
  select(
    all_of(
      final_predictors
    )
  ) %>%
  as.data.frame()

valid_prediction <- prediction_data %>%
  as_tibble() %>%
  transmute(
    valid = if_all(
      everything(),
      ~ !is.na(.x) &
        is.finite(.x)
    )
  ) %>%
  pull(
    valid
  )

cat(
  "\nSamples available for prediction:",
  sum(valid_prediction),
  "of",
  nrow(sample_data),
  "\n"
)

if (any(!valid_prediction)) {
  
  cat(
    "\nSamples excluded because of missing/non-finite predictors:\n"
  )
  
  print(
    sample_data %>%
      filter(
        !valid_prediction
      ) %>%
      select(
        any_of(
          c(
            "Fish_ID",
            "folder",
            "picname",
            "dataset",
            "category",
            "sample_group",
            "cal"
          )
        )
      ),
    n = Inf
  )
}

# Prediction probabilities -----------------------------------------------
# A fixed seed is set immediately before prediction. Final species classes
# are then derived from the probability matrix with deterministic tie
# breaking, so identical model/data inputs always give identical classes.

probability_matrix <- matrix(
  NA_real_,
  nrow = nrow(sample_data),
  ncol = length(species_levels),
  dimnames = list(
    NULL,
    species_levels
  )
)

if (any(valid_prediction)) {
  
  set.seed(prediction_seed)
  
  probability_valid <- predict(
    final_model,
    newdata = prediction_data[
      valid_prediction,
      ,
      drop = FALSE
    ],
    type = "prob"
  )
  
  probability_matrix[
    valid_prediction,
    colnames(probability_valid)
  ] <- as.matrix(
    probability_valid
  )
}

# Deterministic model class from probabilities ---------------------------

model_predicted_species <- rep(
  NA_character_,
  nrow(sample_data)
)

if (any(valid_prediction)) {
  
  probability_for_class <- probability_matrix[
    valid_prediction,
    model_species_order,
    drop = FALSE
  ]
  
  winning_class <- max.col(
    probability_for_class,
    ties.method = "first"
  )
  
  model_predicted_species[
    valid_prediction
  ] <- model_species_order[
    winning_class
  ]
}

probability_data <- probability_matrix %>%
  as_tibble()

names(probability_data) <- paste0(
  "prob_",
  names(probability_data)
)

# Prediction confidence ---------------------------------------------------

get_prediction_confidence <- function(x) {
  
  valid <- which(
    !is.na(x) &
      is.finite(x)
  )
  
  if (length(valid) == 0) {
    
    return(
      tibble(
        top_probability = NA_real_,
        second_species = NA_character_,
        second_probability = NA_real_,
        probability_margin = NA_real_
      )
    )
  }
  
  ordered <- valid[
    order(
      -x[valid],
      valid
    )
  ]
  
  top_probability <- x[
    ordered[1]
  ]
  
  if (length(ordered) >= 2) {
    
    second_species <- names(x)[
      ordered[2]
    ]
    
    second_probability <- x[
      ordered[2]
    ]
    
  } else {
    
    second_species <- NA_character_
    second_probability <- NA_real_
  }
  
  tibble(
    top_probability = top_probability,
    second_species = second_species,
    second_probability = second_probability,
    probability_margin =
      top_probability -
      second_probability
  )
}

prediction_confidence <- map_dfr(
  seq_len(
    nrow(probability_matrix)
  ),
  ~ get_prediction_confidence(
    probability_matrix[
      .x,
    ]
  )
)

# Total probability assigned to gadoids ----------------------------------

gadoid_probability <- rep(
  NA_real_,
  nrow(sample_data)
)

if (length(gadoid_species) > 0) {
  
  gadoid_probability[
    valid_prediction
  ] <- rowSums(
    probability_matrix[
      valid_prediction,
      gadoid_species,
      drop = FALSE
    ]
  )
}

# Apply unknown-gadoid rule -----------------------------------------------

uncertain_gadoid <- (
  valid_prediction &
    model_predicted_species %in%
    gadoid_species &
    prediction_confidence$
    top_probability <
    gadoid_threshold
)

predicted_species <- case_when(
  !valid_prediction ~
    NA_character_,
  
  uncertain_gadoid ~
    "gadoids_unknown",
  
  TRUE ~
    model_predicted_species
)

# Check predictions -------------------------------------------------------

unexpected_predictions <- predicted_species[
  !is.na(predicted_species) &
    !predicted_species %in%
    prediction_species_levels
] %>%
  unique()

if (length(unexpected_predictions) > 0) {
  stop(
    "Unexpected predicted species: ",
    paste(
      unexpected_predictions,
      collapse = ", "
    )
  )
}

# Set consistent species order -------------------------------------------

predicted_species <- factor(
  predicted_species,
  levels = prediction_species_levels
)

model_predicted_species <- factor(
  model_predicted_species,
  levels = expected_model_species
)

# Gadoid reclassification summary ----------------------------------------

gadoid_reclass_summary <- tibble(
  model_predicted_species =
    model_predicted_species,
  
  predicted_species =
    predicted_species,
  
  top_probability =
    prediction_confidence$
    top_probability,
  
  gadoid_probability =
    gadoid_probability,
  
  uncertain_gadoid =
    uncertain_gadoid
) %>%
  filter(
    uncertain_gadoid
  ) %>%
  count(
    model_predicted_species,
    name = "n_reclassified"
  )

cat(
  "\nGadoids reclassified as unknown gadoid:\n"
)

print(
  gadoid_reclass_summary,
  n = Inf
)

write_csv(
  gadoid_reclass_summary,
  "outputs/gadoid_reclassification_summary.csv"
)

# Final prediction dataset ------------------------------------------------

predicted_data <- bind_cols(
  
  sample_data %>%
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
      AR,
      C,
      E,
      FF,
      RE,
      RO
    ),
  
  tibble(
    prediction_valid =
      valid_prediction,
    
    model_predicted_species =
      model_predicted_species,
    
    predicted_species =
      predicted_species,
    
    uncertain_gadoid =
      uncertain_gadoid,
    
    gadoid_probability =
      gadoid_probability
  ),
  
  prediction_confidence,
  
  probability_data
)

# Export predictions ------------------------------------------------------

write_csv(
  predicted_data,
  output_file
)

cat(
  "\nPrediction file successfully written to:\n",
  normalizePath(
    output_file,
    winslash = "/"
  ),
  "\n"
)

# Gadoid diagnostics ------------------------------------------------------

gadoid_diagnostics <- predicted_data %>%
  filter(
    prediction_valid,
    model_predicted_species %in%
      gadoid_species
  ) %>%
  select(
    any_of(
      c(
        "Fish_ID",
        "picname"
      )
    ),
    category,
    sample_group,
    model_predicted_species,
    predicted_species,
    top_probability,
    second_species,
    second_probability,
    probability_margin,
    gadoid_probability,
    uncertain_gadoid
  ) %>%
  arrange(
    top_probability
  )

write_csv(
  gadoid_diagnostics,
  "outputs/gadoid_prediction_diagnostics.csv"
)

# Overall prediction summary ---------------------------------------------

predicted_data_overview <- predicted_data %>%
  filter(
    prediction_valid
  ) %>%
  count(
    predicted_species,
    name = "n",
    .drop = FALSE
  ) %>%
  mutate(
    proportion =
      n /
      sum(n)
  ) %>%
  arrange(
    predicted_species
  )

print(
  predicted_data_overview,
  n = Inf
)

write_csv(
  predicted_data_overview,
  "outputs/predicted_data_overview.csv"
)

# Prediction summary by MP category --------------------------------------

prediction_category_summary <- predicted_data %>%
  filter(
    prediction_valid
  ) %>%
  count(
    category,
    predicted_species,
    name = "n",
    .drop = FALSE
  ) %>%
  group_by(
    category
  ) %>%
  mutate(
    proportion =
      n /
      sum(n)
  ) %>%
  ungroup()

write_csv(
  prediction_category_summary,
  "outputs/prediction_category_summary.csv"
)

# Prediction summary by sample group -------------------------------------

prediction_group_summary <- predicted_data %>%
  filter(
    prediction_valid
  ) %>%
  count(
    category,
    sample_group,
    predicted_species,
    name = "n",
    .drop = FALSE
  ) %>%
  group_by(
    category,
    sample_group
  ) %>%
  mutate(
    proportion =
      n /
      sum(n)
  ) %>%
  ungroup()

write_csv(
  prediction_group_summary,
  "outputs/prediction_group_summary.csv"
)

# Confidence summary ------------------------------------------------------

prediction_confidence_summary <- predicted_data %>%
  filter(
    prediction_valid
  ) %>%
  group_by(
    predicted_species,
    .drop = FALSE
  ) %>%
  summarise(
    n = n(),
    
    mean_probability = mean(
      top_probability,
      na.rm = TRUE
    ),
    
    median_probability = median(
      top_probability,
      na.rm = TRUE
    ),
    
    mean_margin = mean(
      probability_margin,
      na.rm = TRUE
    ),
    
    median_margin = median(
      probability_margin,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

write_csv(
  prediction_confidence_summary,
  "outputs/prediction_confidence_summary.csv"
)

# Prediction count plot ---------------------------------------------------

p_predict_bar <- predicted_data %>%
  filter(
    prediction_valid
  ) %>%
  ggplot(
    aes(
      x = predicted_species,
      fill = predicted_species
    )
  ) +
  geom_bar(
    color = "black"
  ) +
  facet_grid(
    cols = vars(
      category
    )
  ) +
  scale_fill_manual(
    values =
      prediction_species_colors,
    breaks =
      prediction_species_levels,
    labels =
      prediction_species_labels,
    drop = FALSE
  ) +
  scale_x_discrete(
    breaks =
      prediction_species_levels,
    labels =
      prediction_species_labels,
    drop = FALSE
  ) +
  labs(
    x = "Predicted species",
    y = "Number of otoliths"
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

p_predict_bar

ggsave(
  "outputs/p_predict_bar.png",
  p_predict_bar,
  width = 11,
  height = 4
)

# Proportional prediction plot -------------------------------------------

prediction_category_plot <- predicted_data %>%
  filter(
    prediction_valid
  ) %>%
  count(
    category,
    predicted_species,
    name = "n",
    .drop = FALSE
  ) %>%
  group_by(
    category
  ) %>%
  mutate(
    proportion =
      n /
      sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    category = factor(
      category,
      levels = c(
        "No_MP",
        "MP"
      )
    )
  )

p_predict_bar_prop <- ggplot(
  prediction_category_plot,
  aes(
    x = category,
    y = proportion,
    fill = predicted_species
  )
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values =
      prediction_species_colors,
    breaks =
      prediction_species_levels,
    labels =
      prediction_species_labels,
    name =
      "Predicted species",
    drop = FALSE
  ) +
  scale_y_continuous(
    labels =
      scales::percent,
    limits = c(
      0,
      1
    ),
    expand =
      expansion(
        mult = c(
          0,
          0.02
        )
      )
  ) +
  labs(
    x = NULL,
    y = "Proportion of otoliths"
  ) +
  theme_bw() +
  theme(
    panel.grid =
      element_blank()
  )

p_predict_bar_prop

ggsave(
  "outputs/p_predict_bar_prop.png",
  p_predict_bar_prop,
  width = 7,
  height = 5
)

# Predictions by sample group --------------------------------------------

p_predict_bar_group <- predicted_data %>%
  filter(
    prediction_valid
  ) %>%
  ggplot(
    aes(
      x = sample_group,
      fill = predicted_species
    )
  ) +
  geom_bar(
    color = "black"
  ) +
  facet_grid(
    cols = vars(
      category
    ),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values =
      prediction_species_colors,
    breaks =
      prediction_species_levels,
    labels =
      prediction_species_labels,
    name =
      "Predicted species",
    drop = FALSE
  ) +
  labs(
    x = "Sample group",
    y = "Number of otoliths"
  ) +
  theme_bw() +
  theme(
    panel.grid =
      element_blank(),
    strip.background =
      element_blank(),
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

p_predict_bar_group

ggsave(
  "outputs/p_predict_bar_group.png",
  p_predict_bar_group,
  width = 11,
  height = 5
)

# Prediction confidence plot ---------------------------------------------

p_confidence <- predicted_data %>%
  filter(
    prediction_valid
  ) %>%
  ggplot(
    aes(
      x = predicted_species,
      y = top_probability,
      fill = predicted_species
    )
  ) +
  geom_hline(
    yintercept =
      gadoid_threshold,
    linetype =
      "dashed"
  ) +
  geom_boxplot(
    outlier.shape =
      NA
  ) +
  geom_jitter(
    position = position_jitter(
      width = 0.15,
      height = 0,
      seed = prediction_seed
    ),
    alpha = 0.3,
    size = 1
  ) +
  scale_fill_manual(
    values =
      prediction_species_colors,
    breaks =
      prediction_species_levels,
    labels =
      prediction_species_labels,
    drop = FALSE
  ) +
  scale_x_discrete(
    breaks =
      prediction_species_levels,
    labels =
      prediction_species_labels,
    drop = FALSE
  ) +
  labs(
    x = "Predicted species",
    y = "Prediction probability"
  ) +
  theme_bw() +
  theme(
    legend.position =
      "none",
    panel.grid =
      element_blank(),
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

p_confidence

ggsave(
  "outputs/prediction_confidence.png",
  p_confidence,
  width = 11,
  height = 5
)