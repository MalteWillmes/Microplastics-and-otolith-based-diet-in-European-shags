# Otolith shape species classification
# Random forest models using outline coefficients and morphometrics

rm(list = ls())

library(tidyverse)
library(shapeR)
library(vegan)
library(caret)
library(ranger)

# Settings ----------------------------------------------------------------

shape_file <- "shape_extraction/reference_otoliths/shapedata_refoto.rds"
out_dir <- "outputs/model_eval"

# Set TRUE to restore the saved workspace and skip the long model run.
# Set FALSE whenever data, predictors, model settings, or model code change.
reload_saved_environment <- TRUE
environment_file <- file.path(
  out_dir,
  "model_building_environment.RData"
)

# Canonical species order, labels, and colors used throughout the project.
# Only the seven species in model_species_order are reference-model classes.
# Unknown gadoid, Other species, and Unknown are downstream categories.
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

model_species_order <- c(
  "cod",
  "saithe",
  "haddock",
  "poor_cod",
  "norway_pout",
  "lesser_sand_eel",
  "wrasse"
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

species_palette_all <- c(
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

if (reload_saved_environment) {
  
  if (!file.exists(environment_file)) {
    stop(
      "Saved model environment not found: ",
      environment_file
    )
  }
  
  load(
    environment_file,
    envir = .GlobalEnv
  )
  
  cat(
    "\nSaved model-building environment loaded from:\n",
    environment_file,
    "\n"
  )
  
} else {
  
  # split_seed <- 20261001
  # cv_seed <- 20261001
  split_seed <- 20261002
  cv_seed <- 20261002
  
  n_repeats <- 10
  n_trees <- 2000
  
  excluded_species <- "vahls_eelpout"
  
  
  # Read shapeR object ------------------------------------------------------
  
  shape <- readRDS(shape_file)
  
  # Reset any stored shapeR filter
  shape <- setFilter(shape)
  
  master <- getMasterlist(shape) %>%
    as_tibble()
  
  if (!"key" %in% names(master)) {
    master <- master %>%
      mutate(
        key = paste(
          folder,
          picname,
          sep = ";"
        )
      )
  }
  
  wavelet_raw <- getWavelet(shape)
  fourier_raw <- getFourier(shape)
  
  stopifnot(
    nrow(master) == nrow(wavelet_raw),
    nrow(master) == nrow(fourier_raw)
  )
  
  # Remove coefficient columns containing non-finite values -----------------
  
  wavelet_raw <- wavelet_raw[
    ,
    apply(
      wavelet_raw,
      2,
      function(x) all(is.finite(x))
    ),
    drop = FALSE
  ]
  
  fourier_raw <- fourier_raw[
    ,
    apply(
      fourier_raw,
      2,
      function(x) all(is.finite(x))
    ),
    drop = FALSE
  ]
  
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
  
  wavelet_names <- names(wavelet_data)
  fourier_names <- names(fourier_data)
  
  # Outline reconstruction diagnostic --------------------------------------
  
  outline_reconstruction <- estimate.outline.reconstruction(
    shape
  )
  
  png(
    file.path(
      out_dir,
      "outline_reconstruction.png"
    ),
    width = 10,
    height = 5,
    units = "in",
    res = 300
  )
  
  outline.reconstruction.plot(
    outline_reconstruction,
    ref.w.level = 5,
    ref.f.harmonics = 12,
    max.num.harmonics = 32
  )
  
  dev.off()
  
  # Raw morphometrics -------------------------------------------------------
  
  raw_data <- shape@shape.coef.raw %>%
    as.data.frame()
  
  required_raw <- c(
    "otolith.area",
    "otolith.length",
    "otolith.width",
    "otolith.perimeter"
  )
  
  if (!all(required_raw %in% names(raw_data))) {
    stop(
      "Required raw morphometric columns are missing from shape.coef.raw."
    )
  }
  
  raw_data <- raw_data %>%
    select(
      all_of(required_raw)
    )
  
  if (
    !is.null(rownames(raw_data)) &&
    all(master$key %in% rownames(raw_data))
  ) {
    
    raw_measurements <- raw_data %>%
      rownames_to_column("key") %>%
      as_tibble() %>%
      rename(
        raw_OA = otolith.area,
        raw_OL = otolith.length,
        raw_OW = otolith.width,
        raw_OP = otolith.perimeter
      )
    
    sample_data <- master %>%
      left_join(
        raw_measurements,
        by = "key"
      )
    
  } else {
    
    if (nrow(raw_data) != nrow(master)) {
      stop(
        "Could not align raw morphometrics with the shapeR master list."
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
    
    sample_data <- bind_cols(
      master,
      raw_measurements
    )
  }
  
  # Morphometrics -----------------------------------------------------------
  
  sample_data <- sample_data %>%
    mutate(
      
      # Calibrated absolute size
      OL = otolith.length,
      OW = otolith.width,
      OA = otolith.area,
      OP = otolith.perimeter,
      
      # Scale-independent indices from raw measurements
      AR = raw_OL / raw_OW,
      C = raw_OP^2 / raw_OA,
      E = (raw_OL - raw_OW) /
        (raw_OL + raw_OW),
      FF = 4 * pi * raw_OA /
        raw_OP^2,
      RE = raw_OA /
        (raw_OL * raw_OW),
      RO = 4 * raw_OA /
        (pi * raw_OL^2)
    )
  
  write_csv(
    sample_data,
    file.path(
      out_dir,
      "reference_morphometrics.csv"
    )
  )
  
  # Calibration summary ----------------------------------------------------
  
  calibration_summary <- sample_data %>%
    filter(
      !species %in% excluded_species
    ) %>%
    mutate(
      species = factor(
        species,
        levels = model_species_order
      )
    ) %>%
    group_by(species, .drop = FALSE) %>%
    summarise(
      n = n(),
      n_calibrated = sum(
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
      .groups = "drop"
    ) %>%
    filter(
      n > 0
    ) %>%
    arrange(
      species
    )
  
  print(
    calibration_summary,
    n = Inf
  )
  
  write_csv(
    calibration_summary,
    file.path(
      out_dir,
      "reference_calibration_summary.csv"
    )
  )
  
  # Species palette and order -----------------------------------------------
  
  observed_species <- sample_data %>%
    filter(
      !species %in% excluded_species
    ) %>%
    distinct(species) %>%
    pull(species) %>%
    as.character()
  
  unexpected_species <- setdiff(
    observed_species,
    model_species_order
  )
  
  if (length(unexpected_species) > 0) {
    stop(
      "Unexpected reference species found: ",
      paste(
        unexpected_species,
        collapse = ", "
      )
    )
  }
  
  species_levels_all <- model_species_order[
    model_species_order %in% observed_species
  ]
  
  missing_reference_species <- setdiff(
    model_species_order,
    observed_species
  )
  
  if (length(missing_reference_species) > 0) {
    cat(
      "\nReference species in the model order but absent from the reference data:",
      paste(
        missing_reference_species,
        collapse = ", "
      ),
      "\n"
    )
  }
  
  species_palette <- species_palette_all[
    species_levels_all
  ]
  
  # Morphometric figure -----------------------------------------------------
  
  oto_indi <- sample_data %>%
    filter(
      !species %in% excluded_species
    ) %>%
    mutate(
      species = factor(
        species,
        levels = species_levels_all
      )
    ) %>%
    select(
      species,
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
    ) %>%
    pivot_longer(
      cols = -species,
      names_to = "parameter",
      values_to = "value"
    )
  
  p_oto_indi <- ggplot(
    oto_indi,
    aes(
      x = species,
      y = value,
      fill = species
    )
  ) +
    geom_boxplot(
      outlier.alpha = 0.3
    ) +
    facet_wrap(
      ~parameter,
      scales = "free_y",
      ncol = 2
    ) +
    scale_fill_manual(
      values = species_palette,
      breaks = species_levels_all,
      labels = species_labels[species_levels_all],
      drop = FALSE
    ) +
    scale_x_discrete(
      limits = species_levels_all,
      labels = species_labels[species_levels_all],
      drop = FALSE
    ) +
    labs(
      x = NULL,
      y = "Value"
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
  
  p_oto_indi
  
  ggsave(
    file.path(
      out_dir,
      "otolith_parameters_boxplots.png"
    ),
    p_oto_indi,
    width = 10,
    height = 14
  )
  
  # CAP analysis ------------------------------------------------------------
  
  cap_keep <- !master$species %in%
    excluded_species
  
  cap_wavelet <- wavelet_raw[
    cap_keep,
    ,
    drop = FALSE
  ]
  
  cap_species <- factor(
    master$species[
      cap_keep
    ],
    levels = species_levels_all
  )
  
  cap_res <- capscale(
    cap_wavelet ~ cap_species,
    distance = "euclidean"
  )
  
  cap_test <- anova(
    cap_res,
    permutations = 999
  )
  
  cap_axis_test <- anova(
    cap_res,
    by = "axis",
    permutations = 999
  )
  
  print(cap_test)
  print(cap_axis_test)
  
  capture.output(
    cap_test,
    file = file.path(
      out_dir,
      "CAP_permutation_test.txt"
    )
  )
  
  capture.output(
    cap_axis_test,
    file = file.path(
      out_dir,
      "CAP_axis_tests.txt"
    )
  )
  
  cap_eigen <- eigenvals(
    cap_res,
    model = "constrained"
  )
  
  cap_eigen_ratio <-
    cap_eigen /
    sum(cap_eigen)
  
  cap_scores <- scores(
    cap_res,
    display = "sites",
    choices = 1:2
  ) %>%
    as.data.frame() %>%
    as_tibble()
  
  names(cap_scores)[1:2] <- c(
    "CAP1",
    "CAP2"
  )
  
  cap_scores <- cap_scores %>%
    mutate(
      species = cap_species
    )
  
  cap_centroids <- cap_scores %>%
    group_by(species) %>%
    summarise(
      CAP1 = mean(CAP1),
      CAP2 = mean(CAP2),
      .groups = "drop"
    ) %>%
    mutate(
      species_label = species_labels[
        as.character(species)
      ]
    )
  
  CAP1_lab <- paste0(
    "CAP1 (",
    round(
      cap_eigen_ratio[1] * 100,
      1
    ),
    "%)"
  )
  
  CAP2_lab <- paste0(
    "CAP2 (",
    round(
      cap_eigen_ratio[2] * 100,
      1
    ),
    "%)"
  )
  
  p_cap <- ggplot(
    cap_scores,
    aes(
      CAP1,
      CAP2,
      color = species
    )
  ) +
    geom_point(
      alpha = 0.2
    ) +
    geom_point(
      data = cap_centroids,
      size = 3
    ) +
    geom_text(
      data = cap_centroids,
      aes(
        label = species_label
      ),
      nudge_y = 0.02,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = species_palette,
      breaks = species_levels_all,
      labels = species_labels[species_levels_all],
      drop = FALSE
    ) +
    labs(
      x = CAP1_lab,
      y = CAP2_lab,
      color = "Species"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    )
  
  p_cap
  
  ggsave(
    file.path(
      out_dir,
      "CAP_species.png"
    ),
    p_cap,
    width = 8,
    height = 6
  )
  
  # Classification dataset --------------------------------------------------
  
  classification_base <- sample_data %>%
    select(
      any_of(
        c(
          "Fish_ID",
          "key",
          "cal"
        )
      ),
      species,
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
  
  classification_data_raw <- bind_cols(
    classification_base,
    wavelet_data,
    fourier_data
  ) %>%
    filter(
      !species %in% excluded_species
    ) %>%
    mutate(
      species = factor(
        species,
        levels = species_levels_all
      )
    )
  
  # Predictors --------------------------------------------------------------
  
  size_predictors <- c(
    "OL",
    "OW",
    "OA",
    "OP"
  )
  
  shape_indices <- c(
    "AR",
    "FF",
    "RE",
    "RO"
  )
  
  morphometric_predictors <- c(
    size_predictors,
    shape_indices
  )
  
  # C and E remain in descriptive outputs but are not predictors because:
  # C is the inverse form of FF and E is determined by AR.
  
  predictor_sets <- list(
    
    Wavelet =
      wavelet_names,
    
    Fourier =
      fourier_names,
    
    Wavelet_Fourier =
      c(
        wavelet_names,
        fourier_names
      ),
    
    Morphometrics =
      morphometric_predictors,
    
    Wavelet_Morphometrics =
      c(
        wavelet_names,
        morphometric_predictors
      ),
    
    Fourier_Morphometrics =
      c(
        fourier_names,
        morphometric_predictors
      ),
    
    Wavelet_Fourier_Morphometrics =
      c(
        wavelet_names,
        fourier_names,
        morphometric_predictors
      )
  )
  
  # Use same specimens for every candidate model ---------------------------
  
  all_candidate_predictors <- predictor_sets %>%
    unlist() %>%
    unique()
  
  excluded_incomplete <- classification_data_raw %>%
    filter(
      !if_all(
        all_of(
          all_candidate_predictors
        ),
        ~ !is.na(.x) &
          is.finite(.x)
      )
    )
  
  if (nrow(excluded_incomplete) > 0) {
    
    cat(
      "\nReference samples excluded because calibration or model predictors are missing:\n"
    )
    
    print(
      excluded_incomplete %>%
        count(
          species,
          name = "n_excluded"
        ),
      n = Inf
    )
  }
  
  classification_data <- classification_data_raw %>%
    filter(
      if_all(
        all_of(
          all_candidate_predictors
        ),
        ~ !is.na(.x) &
          is.finite(.x)
      )
    ) %>%
    mutate(
      species = droplevels(
        species
      )
    )
  
  cat(
    "\nSamples retained for classification:\n"
  )
  
  print(
    table(
      classification_data$species
    )
  )
  
  if (
    any(
      table(
        classification_data$species
      ) < 3
    )
  ) {
    stop(
      "At least one species has fewer than three usable reference samples."
    )
  }
  
  # Remove near-zero variance predictors -----------------------------------
  
  remove_nzv <- function(
    data,
    predictors
  ) {
    
    nzv <- nearZeroVar(
      data[
        ,
        predictors,
        drop = FALSE
      ]
    )
    
    if (length(nzv) > 0) {
      predictors <- predictors[
        -nzv
      ]
    }
    
    predictors
  }
  
  # Train/test split --------------------------------------------------------
  
  set.seed(
    split_seed
  )
  
  train_index <- createDataPartition(
    classification_data$species,
    p = 0.70,
    list = FALSE
  )
  
  classdata_train <- classification_data[
    train_index,
    ,
    drop = FALSE
  ]
  
  classdata_test <- classification_data[
    -train_index,
    ,
    drop = FALSE
  ]
  
  cat(
    "\nTraining samples:\n"
  )
  
  print(
    table(
      classdata_train$species
    )
  )
  
  cat(
    "\nTest samples:\n"
  )
  
  print(
    table(
      classdata_test$species
    )
  )
  
  predictor_sets <- map(
    predictor_sets,
    ~ remove_nzv(
      classdata_train,
      .x
    )
  )
  
  # Export candidate predictor sets ----------------------------------------
  
  predictor_set_summary <- imap_dfr(
    predictor_sets,
    ~ tibble(
      predictor_set = .y,
      n_predictors = length(.x),
      predictors = paste(
        .x,
        collapse = "; "
      )
    )
  )
  
  write_csv(
    predictor_set_summary,
    file.path(
      out_dir,
      "candidate_predictor_sets.csv"
    )
  )
  
  # Cross-validation --------------------------------------------------------
  
  min_class_n <- min(
    table(
      classdata_train$species
    )
  )
  
  n_folds <- min(
    5,
    min_class_n
  )
  
  if (n_folds < 2) {
    stop(
      "At least one species has fewer than two training samples."
    )
  }
  
  set.seed(
    cv_seed
  )
  
  cv_index <- createMultiFolds(
    classdata_train$species,
    k = n_folds,
    times = n_repeats
  )
  
  ctrl_down <- trainControl(
    method = "repeatedcv",
    number = n_folds,
    repeats = n_repeats,
    index = cv_index,
    classProbs = TRUE,
    savePredictions = "final",
    sampling = "down",
    allowParallel = FALSE
  )
  
  ctrl_weight <- trainControl(
    method = "repeatedcv",
    number = n_folds,
    repeats = n_repeats,
    index = cv_index,
    classProbs = TRUE,
    savePredictions = "final",
    allowParallel = FALSE
  )
  
  # Tuning ------------------------------------------------------------------
  
  make_grid <- function(
    n_predictors
  ) {
    
    mtry_values <-
      sqrt(n_predictors) *
      c(
        0.5,
        1,
        2
      ) %>%
      round() %>%
      pmax(1) %>%
      pmin(n_predictors) %>%
      unique()
    
    expand.grid(
      mtry = mtry_values,
      splitrule = c(
        "gini",
        "extratrees"
      ),
      min.node.size = c(
        1,
        5
      ),
      stringsAsFactors = FALSE
    )
  }
  
  make_case_weights <- function(y) {
    
    class_n <- table(y)
    
    weights <-
      1 /
      class_n[
        as.character(y)
      ]
    
    as.numeric(
      weights /
        mean(weights)
    )
  }
  
  # Ranger fitting function -------------------------------------------------
  
  fit_ranger <- function(
    data,
    predictors,
    balance
  ) {
    
    model_data <- data %>%
      select(
        species,
        all_of(
          predictors
        )
      ) %>%
      as.data.frame()
    
    tune_grid <- make_grid(
      length(
        predictors
      )
    )
    
    set.seed(
      cv_seed
    )
    
    if (balance == "Downsampling") {
      
      train(
        species ~ .,
        data = model_data,
        method = "ranger",
        metric = "Kappa",
        tuneGrid = tune_grid,
        importance = "permutation",
        num.trees = n_trees,
        num.threads = 1,
        trControl = ctrl_down
      )
      
    } else {
      
      train(
        species ~ .,
        data = model_data,
        method = "ranger",
        metric = "Kappa",
        tuneGrid = tune_grid,
        importance = "permutation",
        num.trees = n_trees,
        num.threads = 1,
        weights = make_case_weights(
          model_data$species
        ),
        trControl = ctrl_weight
      )
    }
  }
  
  # Fit candidate models ----------------------------------------------------
  
  models <- list()
  
  model_metadata <- tibble(
    model = character(),
    predictor_set = character(),
    balance = character()
  )
  
  for (
    predictor_set_name in
    names(predictor_sets)
  ) {
    
    for (
      balance_method in
      c(
        "Downsampling",
        "Weighting"
      )
    ) {
      
      model_name_current <- paste(
        predictor_set_name,
        balance_method,
        sep = "_"
      )
      
      cat(
        "\nFitting:",
        model_name_current,
        "\n"
      )
      
      models[[
        model_name_current
      ]] <- fit_ranger(
        data = classdata_train,
        predictors =
          predictor_sets[[
            predictor_set_name
          ]],
        balance =
          balance_method
      )
      
      model_metadata <- bind_rows(
        model_metadata,
        tibble(
          model =
            model_name_current,
          predictor_set =
            predictor_set_name,
          balance =
            balance_method
        )
      )
    }
  }
  
  # Cross-validation performance -------------------------------------------
  
  cv_performance <- imap_dfr(
    models,
    function(
    fit,
    model_name_current
    ) {
      
      best_result <- fit$results %>%
        as_tibble() %>%
        inner_join(
          fit$bestTune %>%
            as_tibble(),
          by = names(
            fit$bestTune
          )
        ) %>%
        slice(1)
      
      tibble(
        model =
          model_name_current,
        Accuracy =
          best_result$Accuracy,
        AccuracySD =
          best_result$AccuracySD,
        Kappa =
          best_result$Kappa,
        KappaSD =
          best_result$KappaSD
      )
    }
  ) %>%
    left_join(
      model_metadata,
      by = "model"
    ) %>%
    arrange(
      desc(Kappa)
    )
  
  print(
    cv_performance,
    n = Inf
  )
  
  write_csv(
    cv_performance,
    file.path(
      out_dir,
      "model_CV_performance.csv"
    )
  )
  
  # CV comparison plot ------------------------------------------------------
  
  model_resamples <- resamples(
    models
  )
  
  capture.output(
    summary(
      model_resamples
    ),
    file = file.path(
      out_dir,
      "model_CV_summary.txt"
    )
  )
  
  png(
    file.path(
      out_dir,
      "model_CV_comparison.png"
    ),
    width = 10,
    height = 7,
    units = "in",
    res = 300
  )
  
  print(
    bwplot(
      model_resamples,
      metric = "Kappa"
    )
  )
  
  dev.off()
  
  # Select final model ------------------------------------------------------
  
  best_model_name <- cv_performance %>%
    slice_max(
      Kappa,
      n = 1,
      with_ties = FALSE
    ) %>%
    pull(model)
  
  best_metadata <- model_metadata %>%
    filter(
      model ==
        best_model_name
    )
  
  best_predictor_set <-
    best_metadata$predictor_set
  
  best_balance <-
    best_metadata$balance
  
  final_model <-
    models[[
      best_model_name
    ]]
  
  final_predictors <-
    predictor_sets[[
      best_predictor_set
    ]]
  
  cat(
    "\nSelected model:",
    best_model_name,
    "\n"
  )
  
  cat(
    "Predictor set:",
    best_predictor_set,
    "\n"
  )
  
  cat(
    "Balance method:",
    best_balance,
    "\n"
  )
  
  cat(
    "Predictors:",
    length(final_predictors),
    "\n"
  )
  
  print(
    final_model$bestTune
  )
  
  # Withheld test evaluation ------------------------------------------------
  
  test_data <- classdata_test %>%
    select(
      all_of(
        final_predictors
      )
    ) %>%
    as.data.frame()
  
  predict_test <- predict(
    final_model,
    newdata = test_data,
    type = "raw"
  )
  
  Class_summary_test <- confusionMatrix(
    data = predict_test,
    reference =
      classdata_test$species,
    mode = "everything"
  )
  
  print(
    Class_summary_test
  )
  
  write.csv(
    Class_summary_test$table,
    file.path(
      out_dir,
      "Class_summary_test.csv"
    ),
    row.names = TRUE
  )
  
  # Overall test performance ------------------------------------------------
  
  test_metrics <- tibble(
    model =
      best_model_name,
    
    predictor_set =
      best_predictor_set,
    
    balance =
      best_balance,
    
    n_test =
      nrow(
        classdata_test
      ),
    
    accuracy = unname(
      Class_summary_test$overall[
        "Accuracy"
      ]
    ),
    
    accuracy_ci_lower = unname(
      Class_summary_test$overall[
        "AccuracyLower"
      ]
    ),
    
    accuracy_ci_upper = unname(
      Class_summary_test$overall[
        "AccuracyUpper"
      ]
    ),
    
    kappa = unname(
      Class_summary_test$overall[
        "Kappa"
      ]
    ),
    
    macro_f1 = mean(
      Class_summary_test$byClass[
        ,
        "F1"
      ],
      na.rm = TRUE
    ),
    
    macro_balanced_accuracy = mean(
      Class_summary_test$byClass[
        ,
        "Balanced Accuracy"
      ],
      na.rm = TRUE
    )
  )
  
  print(
    test_metrics
  )
  
  write_csv(
    test_metrics,
    file.path(
      out_dir,
      "test_metrics.csv"
    )
  )
  
  # Species-specific performance -------------------------------------------
  
  species_support <- classdata_test %>%
    count(
      species,
      name = "N_test"
    ) %>%
    mutate(
      species =
        as.character(
          species
        )
    )
  
  species_performance <-
    Class_summary_test$byClass %>%
    as.data.frame() %>%
    rownames_to_column(
      "species"
    ) %>%
    as_tibble() %>%
    mutate(
      species = str_remove(
        species,
        "^Class: "
      )
    ) %>%
    left_join(
      species_support,
      by = "species"
    ) %>%
    mutate(
      species = factor(
        species,
        levels = species_levels_all
      )
    ) %>%
    arrange(
      species
    )
  
  write_csv(
    species_performance,
    file.path(
      out_dir,
      "species_classification_performance.csv"
    )
  )
  
  # Confusion matrix --------------------------------------------------------
  
  class_plot_data <-
    Class_summary_test$table %>%
    as.data.frame() %>%
    mutate(
      Reference = factor(
        as.character(Reference),
        levels = species_levels_all
      ),
      Prediction = factor(
        as.character(Prediction),
        levels = species_levels_all
      )
    ) %>%
    group_by(
      Reference
    ) %>%
    mutate(
      percent =
        100 *
        Freq /
        sum(Freq)
    ) %>%
    ungroup()
  
  p_confmat <- ggplot(
    class_plot_data,
    aes(
      x = Reference,
      y = Prediction,
      fill = percent
    )
  ) +
    geom_tile() +
    geom_text(
      aes(
        label = sprintf(
          "%.1f",
          percent
        )
      )
    ) +
    scale_fill_gradient(
      low = "white",
      high = "#3575b5",
      limits = c(
        0,
        100
      )
    ) +
    scale_x_discrete(
      limits = species_levels_all,
      labels = species_labels[species_levels_all],
      drop = FALSE
    ) +
    scale_y_discrete(
      limits = species_levels_all,
      labels = species_labels[species_levels_all],
      drop = FALSE
    ) +
    labs(
      x = "Reference species",
      y = "Predicted species",
      fill = "Percent"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  p_confmat
  
  ggsave(
    file.path(
      out_dir,
      "confusion_matrix.png"
    ),
    p_confmat,
    width = 8,
    height = 7
  )
  
  # Variable importance -----------------------------------------------------
  
  varimp_data <- varImp(
    final_model,
    scale = TRUE
  )$importance %>%
    as.data.frame() %>%
    rownames_to_column(
      "variable"
    ) %>%
    as_tibble()
  
  if (
    !"Overall" %in%
    names(varimp_data)
  ) {
    
    varimp_data <- varimp_data %>%
      mutate(
        Overall = rowMeans(
          across(
            where(
              is.numeric
            )
          ),
          na.rm = TRUE
        )
      )
  }
  
  varimp_data <- varimp_data %>%
    arrange(
      desc(
        Overall
      )
    )
  
  write_csv(
    varimp_data,
    file.path(
      out_dir,
      "feature_importance.csv"
    )
  )
  
  p_varimp <- varimp_data %>%
    slice_head(
      n = 20
    ) %>%
    mutate(
      variable = fct_reorder(
        variable,
        Overall
      )
    ) %>%
    ggplot(
      aes(
        x = Overall,
        y = variable
      )
    ) +
    geom_col() +
    labs(
      x = "Permutation importance",
      y = NULL
    ) +
    theme_classic()
  
  p_varimp
  
  ggsave(
    file.path(
      out_dir,
      "feature_importance.png"
    ),
    p_varimp,
    width = 7,
    height = 7
  )
  
  # Manuscript model table --------------------------------------------------
  
  selected_cv <- cv_performance %>%
    filter(
      model ==
        best_model_name
    )
  
  manuscript_model_table <- tibble(
    
    Model =
      best_model_name,
    
    Predictor_set =
      best_predictor_set,
    
    Balance =
      best_balance,
    
    N_total =
      nrow(
        classification_data
      ),
    
    N_training =
      nrow(
        classdata_train
      ),
    
    N_test =
      nrow(
        classdata_test
      ),
    
    CV_accuracy =
      selected_cv$Accuracy,
    
    CV_accuracy_SD =
      selected_cv$AccuracySD,
    
    CV_kappa =
      selected_cv$Kappa,
    
    CV_kappa_SD =
      selected_cv$KappaSD,
    
    Test_accuracy =
      test_metrics$accuracy,
    
    Test_accuracy_CI =
      sprintf(
        "%.3f–%.3f",
        test_metrics$accuracy_ci_lower,
        test_metrics$accuracy_ci_upper
      ),
    
    Test_kappa =
      test_metrics$kappa,
    
    Macro_F1 =
      test_metrics$macro_f1,
    
    Macro_balanced_accuracy =
      test_metrics$
      macro_balanced_accuracy
  ) %>%
    mutate(
      across(
        where(
          is.numeric
        ),
        ~ round(
          .x,
          3
        )
      )
    )
  
  print(
    manuscript_model_table
  )
  
  write_csv(
    manuscript_model_table,
    file.path(
      out_dir,
      "manuscript_model_performance.csv"
    )
  )
  
  # Manuscript species table ------------------------------------------------
  
  manuscript_species_table <-
    species_performance %>%
    select(
      species,
      N_test,
      Sensitivity,
      Specificity,
      `Pos Pred Value`,
      F1,
      `Balanced Accuracy`
    ) %>%
    rename(
      Recall =
        Sensitivity,
      Precision =
        `Pos Pred Value`,
      Balanced_accuracy =
        `Balanced Accuracy`
    ) %>%
    mutate(
      across(
        c(
          Recall,
          Specificity,
          Precision,
          F1,
          Balanced_accuracy
        ),
        ~ round(
          .x,
          3
        )
      )
    ) %>%
    mutate(
      species = factor(
        species,
        levels = species_levels_all
      )
    ) %>%
    arrange(
      species
    )
  
  print(
    manuscript_species_table
  )
  
  write_csv(
    manuscript_species_table,
    file.path(
      out_dir,
      "manuscript_species_performance.csv"
    )
  )
  
  # Save final model --------------------------------------------------------
  
  model_bundle <- list(
    
    model =
      final_model,
    
    model_name =
      best_model_name,
    
    representation =
      best_predictor_set,
    
    predictor_set =
      best_predictor_set,
    
    balance =
      best_balance,
    
    pca =
      FALSE,
    
    pca_threshold =
      NA_real_,
    
    predictors =
      final_predictors,
    
    size_predictors =
      size_predictors,
    
    shape_indices =
      shape_indices,
    
    species_levels =
      levels(
        classification_data$species
      ),
    
    species_palette =
      species_palette[
        levels(
          classification_data$species
        )
      ],
    
    species_labels =
      species_labels[
        levels(
          classification_data$species
        )
      ],
    
    display_species_order =
      display_species_order,
    
    display_species_labels =
      species_labels[
        display_species_order
      ],
    
    display_species_palette =
      species_palette_all[
        display_species_order
      ],
    
    best_tune =
      final_model$bestTune,
    
    excluded_species =
      excluded_species,
    
    split_seed =
      split_seed,
    
    cv_seed =
      cv_seed,
    
    n_folds =
      n_folds,
    
    n_repeats =
      n_repeats,
    
    n_trees =
      n_trees
  )
  
  saveRDS(
    model_bundle,
    file.path(
      out_dir,
      "final_species_model.rds"
    )
  )
  
  cat(
    "\nFinal model saved to:\n",
    file.path(
      out_dir,
      "final_species_model.rds"
    ),
    "\n"
  )
  
  # Save model-building environment -----------------------------------------
  
  objects_to_save <- setdiff(
    ls(
      envir = .GlobalEnv
    ),
    c(
      "reload_saved_environment",
      "objects_to_save"
    )
  )
  
  save(
    list = objects_to_save,
    file = environment_file,
    envir = .GlobalEnv
  )
  
  rm(objects_to_save)
  
  cat(
    "\nModel-building environment saved to:\n",
    environment_file,
    "\nSet reload_saved_environment <- TRUE at the top of the script to restore it without rerunning the models.\n"
  )
  
}

# Reapply current display settings after either a fresh fit or environment load.
# This allows spelling/order changes to propagate to plots without refitting models.
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

species_palette_all <- c(
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

# Keep the model-only order aligned to the current project order.
model_species_order <- display_species_order[
  display_species_order %in% c(
    "cod",
    "saithe",
    "haddock",
    "poor_cod",
    "norway_pout",
    "lesser_sand_eel",
    "wrasse"
  )
]

# If a saved environment was loaded, re-order the display objects without changing
# the fitted model itself.
if (exists("classification_data")) {
  species_levels_all <- model_species_order[
    model_species_order %in% levels(classification_data$species)
  ]
}

species_palette <- species_palette_all[
  species_levels_all
]

# Regenerate order-sensitive figures --------------------------------------
# These plots are rebuilt after either a fresh model fit or a saved-
# environment load so changes to species order, labels, or colors propagate
# without rerunning the random-forest models.

# Otolith metric boxplots -------------------------------------------------

oto_indi <- sample_data %>%
  filter(
    !species %in% excluded_species
  ) %>%
  mutate(
    species = factor(
      as.character(species),
      levels = species_levels_all
    )
  ) %>%
  select(
    species,
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
  ) %>%
  pivot_longer(
    cols = -species,
    names_to = "parameter",
    values_to = "value"
  )

p_oto_indi <- ggplot(
  oto_indi,
  aes(
    x = species,
    y = value,
    fill = species
  )
) +
  geom_boxplot(
    outlier.alpha = 0.3
  ) +
  facet_wrap(
    ~parameter,
    scales = "free_y",
    ncol = 2
  ) +
  scale_fill_manual(
    values = species_palette,
    breaks = species_levels_all,
    labels = species_labels[species_levels_all],
    drop = FALSE
  ) +
  scale_x_discrete(
    limits = species_levels_all,
    labels = species_labels[species_levels_all],
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Value"
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

p_oto_indi

ggsave(
  file.path(
    out_dir,
    "otolith_parameters_boxplots.png"
  ),
  p_oto_indi,
  width = 10,
  height = 14
)

# CAP plot ----------------------------------------------------------------

cap_scores <- cap_scores %>%
  mutate(
    species = factor(
      as.character(species),
      levels = species_levels_all
    )
  )

cap_centroids <- cap_scores %>%
  group_by(species) %>%
  summarise(
    CAP1 = mean(CAP1),
    CAP2 = mean(CAP2),
    .groups = "drop"
  ) %>%
  mutate(
    species_label = species_labels[
      as.character(species)
    ]
  )

p_cap <- ggplot(
  cap_scores,
  aes(
    CAP1,
    CAP2,
    color = species
  )
) +
  geom_point(
    alpha = 0.2
  ) +
  geom_point(
    data = cap_centroids,
    size = 3
  ) +
  geom_text(
    data = cap_centroids,
    aes(
      label = species_label
    ),
    nudge_y = 0.02,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = species_palette,
    breaks = species_levels_all,
    labels = species_labels[species_levels_all],
    drop = FALSE
  ) +
  labs(
    x = CAP1_lab,
    y = CAP2_lab,
    color = "Species"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )

p_cap

ggsave(
  file.path(
    out_dir,
    "CAP_species.png"
  ),
  p_cap,
  width = 8,
  height = 6
)

# Mean reconstructed otolith shapes ---------------------------------------
# This block runs after either a fresh model fit or a saved-environment load,
# so the figure can be regenerated without refitting the models.
#
# Custom plotting is used here instead of plotWaveletShape() and
# plotFourierShape() so that plotting limits are calculated across all species.
# This prevents larger reconstructed outlines and angle labels from being
# clipped when species differ substantially in mean otolith size.

shape_plot <- shape

shape_plot@master.list$species_plot <- factor(
  species_labels[
    as.character(
      shape_plot@master.list$species
    )
  ],
  levels = unname(
    species_labels[
      species_levels_all
    ]
  )
)

shape_plot <- setFilter(
  shape_plot,
  as.character(
    shape_plot@master.list$species
  ) %in%
    species_levels_all
)

shape_plot_colors <- unname(
  species_palette[
    species_levels_all
  ]
)

inverse_wavelet <- getFromNamespace(
  ".shapeR.inverse.wavelet",
  "shapeR"
)

inverse_fourier <- getFromNamespace(
  ".shapeR.iefourier",
  "shapeR"
)

plot_reconstructed_shapes <- function(
    object,
    class_name,
    reconstruction = c("wavelet", "fourier"),
    colors,
    lwd = 2,
    lty = 1
) {
  reconstruction <- match.arg(reconstruction)
  
  classes <- object@master.list[[class_name]]
  keep <- object@filter & !is.na(classes)
  
  class_levels <- levels(
    droplevels(
      factor(
        classes[keep],
        levels = levels(classes)
      )
    )
  )
  
  if (length(class_levels) == 0) {
    stop("No reference species available for shape reconstruction.")
  }
  
  outlines <- vector(
    "list",
    length(class_levels)
  )
  
  names(outlines) <- class_levels
  
  if (reconstruction == "wavelet") {
    for (i in seq_along(class_levels)) {
      class_i <- class_levels[i]
      ind <- keep & classes == class_i
      
      mean_coef <- apply(
        object@wavelet.coef[ind, , drop = FALSE],
        2,
        mean,
        na.rm = TRUE
      )
      
      mean_radius <- mean(
        object@master.list$mean.radii[ind],
        na.rm = TRUE
      )
      
      reconstructed <- inverse_wavelet(
        mean_coef,
        mean_radius
      )
      
      outlines[[i]] <- tibble(
        x = reconstructed$X,
        y = reconstructed$Y
      )
    }
  } else {
    fourier_coef <- cbind(
      -1,
      0,
      0,
      object@fourier.coef
    )
    
    fourier_seq <- seq(
      1,
      12 * 4,
      by = 4
    )
    
    for (i in seq_along(class_levels)) {
      class_i <- class_levels[i]
      ind <- keep & classes == class_i
      
      mean_coef <- apply(
        fourier_coef[ind, , drop = FALSE],
        2,
        mean,
        na.rm = TRUE
      )
      
      reconstructed <- inverse_fourier(
        mean_coef[fourier_seq],
        mean_coef[fourier_seq + 1],
        mean_coef[fourier_seq + 2],
        mean_coef[fourier_seq + 3],
        12,
        64 * 4
      )
      
      outlines[[i]] <- tibble(
        x = reconstructed$x,
        y = reconstructed$y
      )
    }
  }
  
  all_x <- unlist(
    map(
      outlines,
      "x"
    ),
    use.names = FALSE
  )
  
  all_y <- unlist(
    map(
      outlines,
      "y"
    ),
    use.names = FALSE
  )
  
  x_range <- range(
    all_x,
    finite = TRUE
  )
  
  y_range <- range(
    all_y,
    finite = TRUE
  )
  
  x_span <- diff(x_range)
  y_span <- diff(y_range)
  
  if (!is.finite(x_span) || x_span == 0) {
    x_span <- 1
  }
  
  if (!is.finite(y_span) || y_span == 0) {
    y_span <- 1
  }
  
  xlim <- x_range +
    c(-1, 1) *
    x_span *
    0.18
  
  ylim <- y_range +
    c(-1, 1) *
    y_span *
    0.22
  
  center_x <- mean(
    map_dbl(
      outlines,
      ~ mean(.x$x, na.rm = TRUE)
    )
  )
  
  center_y <- mean(
    map_dbl(
      outlines,
      ~ mean(.x$y, na.rm = TRUE)
    )
  )
  
  plot(
    NA,
    NA,
    type = "n",
    xlim = xlim,
    ylim = ylim,
    xlab = "",
    ylab = "",
    axes = FALSE,
    frame.plot = FALSE,
    asp = 1,
    xaxs = "i",
    yaxs = "i"
  )
  
  abline(
    h = center_y,
    v = center_x,
    lty = 2,
    col = "grey"
  )
  
  for (i in seq_along(outlines)) {
    lines(
      outlines[[i]]$x,
      outlines[[i]]$y,
      col = colors[i],
      lwd = lwd,
      lty = lty
    )
  }
  
  text(
    center_x,
    ylim[2] - 0.05 * diff(ylim),
    "90\u00b0",
    cex = 1.05
  )
  
  text(
    center_x,
    ylim[1] + 0.05 * diff(ylim),
    "270\u00b0",
    cex = 1.05
  )
  
  text(
    xlim[2] - 0.05 * diff(xlim),
    center_y,
    "0\u00b0",
    cex = 1.05
  )
  
  text(
    xlim[1] + 0.06 * diff(xlim),
    center_y,
    "180\u00b0",
    cex = 1.05
  )
  
  legend(
    "bottomleft",
    legend = class_levels,
    col = colors[seq_along(class_levels)],
    lty = lty,
    lwd = lwd,
    cex = 0.85,
    bty = "n",
    inset = 0.01
  )
}

png(
  file.path(
    out_dir,
    "otolith_shape_reconstructions.png"
  ),
  width = 8,
  height = 11,
  units = "in",
  res = 300
)

par(
  mfrow = c(2, 1),
  mar = c(2.5, 2.5, 3, 2.5),
  mgp = c(2.2, 0.7, 0)
)

plot_reconstructed_shapes(
  shape_plot,
  "species_plot",
  reconstruction = "wavelet",
  colors = shape_plot_colors,
  lwd = 2,
  lty = 1
)

mtext(
  "A",
  side = 3,
  adj = 0,
  line = 0.4,
  font = 2,
  cex = 1.4
)

mtext(
  "Wavelet reconstruction",
  side = 3,
  adj = 0.5,
  line = 0.4,
  cex = 1.1
)

plot_reconstructed_shapes(
  shape_plot,
  "species_plot",
  reconstruction = "fourier",
  colors = shape_plot_colors,
  lwd = 2,
  lty = 1
)

mtext(
  "B",
  side = 3,
  adj = 0,
  line = 0.4,
  font = 2,
  cex = 1.4
)

mtext(
  "Fourier reconstruction",
  side = 3,
  adj = 0.5,
  line = 0.4,
  cex = 1.1
)

dev.off()
