library(magick)
library(tidyverse)

input_dir  <- "C:/shape_analysis/original_images"
output_dir <- "C:/shape_analysis/bw_images"

threshold_value <- "35%"

files <- list.files(
  path = input_dir,
  pattern = "\\.(tif|tiff)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

walk(files, function(file) {
  
  relative_path <- str_remove(
    normalizePath(file, winslash = "/"),
    fixed(paste0(normalizePath(input_dir, winslash = "/"), "/"))
  )
  
  output_file <- file.path(
    output_dir,
    str_replace(
      relative_path,
      regex("\\.(tif|tiff)$", ignore_case = TRUE),
      ".jpg"
    )
  )
  
  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Read only the first TIFF frame/page
  img <- image_read(paste0(file, "[0]"))
  
  # Convert to grayscale
  gray <- image_convert(
    img,
    colorspace = "gray"
  )
  
  image_destroy(img)
  
  # Pixels above threshold -> white
  white <- image_threshold(
    gray,
    type = "white",
    threshold = threshold_value
  )
  
  image_destroy(gray)
  
  # Pixels below threshold -> black
  bw <- image_threshold(
    white,
    type = "black",
    threshold = threshold_value
  )
  
  image_destroy(white)
  
  # White otolith on black background
  image_write(
    bw,
    path = output_file,
    format = "jpg",
    quality = 100
  )
  
  image_destroy(bw)
  
  gc()
})

cat(
  "Converted",
  length(files),
  "images using threshold",
  threshold_value,
  "\n"
)