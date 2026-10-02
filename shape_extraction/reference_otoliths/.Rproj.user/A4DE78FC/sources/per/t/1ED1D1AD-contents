# Otolith shape analysis
# No fish size adjustment

library(tidyverse)
library(shapeR)

# Project -----------------------------------------------------------------
rm(list = ls())
project_dir <- here::here()

shape <- shapeR(
  project_dir,
  "FISH.csv"
)

# Check imported metadata -------------------------------------------------

shape@master.list.org[["Fish_ID"]]

# Outline extraction ------------------------------------------------------

shape <- detect.outline(
  shape,
  threshold = 0.2,
  write.outline.w.org = F,
  mouse.click = F
)

# Inspect individual outline if needed
# show.original.with.outline(
#   shape,
#   "sypike",
#   "ref-syp4"
# )

# Remove incorrect outline if needed
 shape <- remove.outline(shape, "labridae", "Labridae ts-sk 24-37_3")
 shape <- remove.outline(shape, "labridae", "Labridae ts-sk 24-37_5")
 shape <- remove.outline(shape, "labridae", "Labridae ts-sk 24-37_8")
 shape <- remove.outline(shape, "labridae", "Labridae ts-sk 24-92_2")
 shape <- remove.outline(shape, "labridae", "Labridae ts-sk 24-92_3")
 shape <- remove.outline(shape, "labridae", "Labridae ts-sk 24-92_4")
 shape <- remove.outline(shape, "torsk", "Cod_ts_sk_24-88_2")
 shape <- remove.outline(shape, "sei", "refsei-4") 
 
 
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

table(master$species)

# Shape measurements ------------------------------------------------------

measurements

tapply(
  measurements$otolith.area,
  master$species,
  mean,
  na.rm = TRUE
)

# Mean shapes by species --------------------------------------------------

plotWaveletShape(
  shape,
  "species",
  show.angle = TRUE,
  lwd = 2,
  lty = 1
)

plotFourierShape(
  shape,
  "species",
  show.angle = TRUE,
  lwd = 2,
  lty = 1
)

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
  file = "shapedata_refoto.rds"
)