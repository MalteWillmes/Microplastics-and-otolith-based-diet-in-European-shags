# Otolith shape analysis - No MP samples
# No fish size adjustment

library(tidyverse)
library(shapeR)

# Project -----------------------------------------------------------------
rm(list = ls())
project_dir <- here::here()

shape <- shapeR(
  project_dir,
  "FISH_noMP.csv"
)

# Check imported metadata -------------------------------------------------

shape@master.list.org[["Fish_ID"]]

# Outline extraction ------------------------------------------------------

shape <- detect.outline(
  shape,
  threshold = 0.2,
  write.outline.w.org = FALSE,
  mouse.click = FALSE
)

# Inspect individual outline if needed
# show.original.with.outline(
#   shape,
#   "f3r1",
#   "f3r1-7"
# )

# Remove incorrect outlines if needed
shape <- remove.outline(
  shape,
  "f3r1",
  "f3r1-7"
)

shape <- remove.outline(
  shape,
  "f4r7",
  "f4r7-10"
)

shape <- remove.outline(
  shape,
  "f4r7",
  "f4r7-12"
)

shape <- remove.outline(
  shape,
  "f4r7",
  "f4r7-18"
)
# Contour smoothing -------------------------------------------------------

shape <- smoothout(
  shape,
  n = 100
)

# Shape coefficients ------------------------------------------------------

shape <- generateShapeCoefficients(shape)

# Reconstruction quality -------------------------------------------------

reconstruction <- estimate.outline.reconstruction(shape)

outline.reconstruction.plot(
  reconstruction,
  max.num.harmonics = 15
)

# Select coefficient resolution ------------------------------------------

n_wavelet_levels <- 5
n_fourier_harmonics <- 12

# Link metadata and calibration ------------------------------------------

shape <- enrich.master.list(
  shape,
  folder_name = "folder",
  pic_name = "picname",
  calibration = "cal",
  n.wavelet.levels = n_wavelet_levels,
  n.fourier.freq = n_fourier_harmonics
)

# Extract analysis data ---------------------------------------------------

master <- getMasterlist(shape)
measurements <- getMeasurements(shape)
wavelet <- getWavelet(shape)
fourier <- getFourier(shape)

# Check matched data ------------------------------------------------------

master

# Check sample groups if available
if ("noMP_sample" %in% names(master)) {
  print(table(master$noMP_sample))
}

# Shape measurements ------------------------------------------------------

measurements

# Unstandardized wavelet coefficients ------------------------------------

plotWavelet(
  shape,
  level = n_wavelet_levels,
  class.name = NULL,
  useStdcoef = FALSE
)

# Unstandardized Fourier coefficients ------------------------------------

plotFourier(
  shape,
  class.name = NULL,
  useStdcoef = FALSE
)

# Save --------------------------------------------------------------------

saveRDS(
  shape,
  file = "shapedata_no_mp_samples.rds"
)

