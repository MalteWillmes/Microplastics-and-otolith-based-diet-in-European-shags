library(tidyverse)
library(exiftoolr)

input_dir <- "C:/shape_analysis/original_images"
output_file <- "C:/shape_analysis/shapeR_calibration.csv"

tif_files <- list.files(
  path = input_dir,
  pattern = "\\.(tif|tiff)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(tif_files) == 0) {
  stop("No TIFF files found in: ", input_dir)
}

metadata <- exif_read(
  path = tif_files,
  tags = c(
    "FileName",
    "Directory",
    "ImageWidth",
    "ImageHeight",
    "ImageDescription"
  ),
  quiet = TRUE
) %>%
  as_tibble()

extract_ome_value <- function(x, field) {
  str_match(
    x,
    paste0(field, '="([^"]+)"')
  )[, 2]
}

convert_to_um <- function(value, unit) {
  
  value <- as.numeric(value)
  
  unit <- unit %>%
    str_to_lower() %>%
    str_replace_all("μ", "µ") %>%
    str_trim()
  
  case_when(
    unit %in% c(
      "µm", "um", "micron", "microns",
      "micrometer", "micrometers",
      "micrometre", "micrometres"
    ) ~ value,
    
    unit %in% c(
      "nm", "nanometer", "nanometers",
      "nanometre", "nanometres"
    ) ~ value / 1000,
    
    unit %in% c(
      "mm", "millimeter", "millimeters",
      "millimetre", "millimetres"
    ) ~ value * 1000,
    
    unit %in% c(
      "cm", "centimeter", "centimeters",
      "centimetre", "centimetres"
    ) ~ value * 10000,
    
    unit %in% c(
      "m", "meter", "meters",
      "metre", "metres"
    ) ~ value * 1e6,
    
    TRUE ~ NA_real_
  )
}

calibration <- metadata %>%
  mutate(
    source_file = normalizePath(
      SourceFile,
      winslash = "/",
      mustWork = FALSE
    ),
    file_name = basename(source_file),
    picname = tools::file_path_sans_ext(file_name),
    folder = basename(dirname(source_file)),
    
    physical_size_x = extract_ome_value(
      ImageDescription,
      "PhysicalSizeX"
    ),
    
    physical_size_y = extract_ome_value(
      ImageDescription,
      "PhysicalSizeY"
    ),
    
    physical_size_x_unit = extract_ome_value(
      ImageDescription,
      "PhysicalSizeXUnit"
    ),
    
    physical_size_y_unit = extract_ome_value(
      ImageDescription,
      "PhysicalSizeYUnit"
    ),
    
    pixel_size_x_um = convert_to_um(
      physical_size_x,
      physical_size_x_unit
    ),
    
    pixel_size_y_um = convert_to_um(
      physical_size_y,
      physical_size_y_unit
    ),
    
    cal_x = 1000 / pixel_size_x_um,
    cal_y = 1000 / pixel_size_y_um,
    
    xy_difference_percent =
      abs(cal_x - cal_y) /
      ((cal_x + cal_y) / 2) * 100,
    
    cal = case_when(
      !is.na(cal_x) &
        !is.na(cal_y) &
        xy_difference_percent <= 1 ~
        (cal_x + cal_y) / 2,
      
      TRUE ~ NA_real_
    ),
    
    status = case_when(
      is.na(ImageDescription) ~
        "No ImageDescription",
      
      is.na(physical_size_x) &
        is.na(physical_size_y) ~
        "No OME physical pixel size",
      
      is.na(physical_size_x_unit) |
        is.na(physical_size_y_unit) ~
        "Pixel-size unit missing",
      
      is.na(pixel_size_x_um) |
        is.na(pixel_size_y_um) ~
        "Unknown pixel-size unit",
      
      xy_difference_percent > 1 ~
        "X/Y pixel size mismatch",
      
      is.na(cal) ~
        "Calibration could not be calculated",
      
      TRUE ~
        "OK"
    )
  ) %>%
  select(
    folder,
    picname,
    file_name,
    source_file,
    ImageWidth,
    ImageHeight,
    physical_size_x,
    physical_size_x_unit,
    physical_size_y,
    physical_size_y_unit,
    pixel_size_x_um,
    pixel_size_y_um,
    cal_x,
    cal_y,
    cal,
    status
  )

write_csv(
  calibration,
  output_file
)

status_summary <- calibration %>%
  count(status, sort = TRUE)

calibration_summary <- calibration %>%
  filter(status == "OK") %>%
  count(
    pixel_size_x_um,
    cal,
    sort = TRUE
  ) %>%
  arrange(pixel_size_x_um)

problem_files <- calibration %>%
  filter(status != "OK") %>%
  select(
    folder,
    picname,
    physical_size_x,
    physical_size_x_unit,
    physical_size_y,
    physical_size_y_unit,
    status
  )

print(status_summary)
print(calibration_summary)
print(problem_files)