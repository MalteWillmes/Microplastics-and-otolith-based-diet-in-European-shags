# Export all predicted samples to Excel

rm(list = ls())

library(tidyverse)
library(here)
library(openxlsx)

predicted_data <- read_csv(
  here("outputs", "predicted_data.csv"),
  show_col_types = FALSE
)

all_predicted_samples <- predicted_data %>%
  filter(
    prediction_valid,
    !is.na(predicted_species)
  ) %>%
  transmute(
    category,
    site = sample_group,
    Fish_ID,
    picname,
    predicted_species,
    top_probability
  ) %>%
  arrange(
    category,
    site,
    picname
  )

write.xlsx(
  all_predicted_samples,
  here("outputs", "all predicted samples.xlsx"),
  sheetName = "all predicted samples",
  overwrite = TRUE
)