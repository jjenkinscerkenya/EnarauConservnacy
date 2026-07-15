# Input manifest for Objective 1's manually-downloaded Dynamic World raster exports
# (DW_INPUT_RASTER_DIR, see 00_config.R). Filenames here match what's actually present on disk
# after manual download: the user stripped each export cell's Drive filename prefix
# (seasonal_categorical_/conversion_pressure_/connectivity_inputs_/transitions_) so every local
# file starts directly at "DW_" -- e.g. Drive's
# connectivity_inputs_DW_connectivity_inputs_current_2022_2025_project.tif is locally
# DW_connectivity_inputs_current_2022_2025_project.tif. Confirmed 2026-07-15 against the actual
# 15-file first-batch download. The underlying band/class scheme is unchanged from
# scripts/python/notebooks/historical_change_detection.ipynb (export cells c34, c38, c42).

expand_seasonal_manifest <- function() {
  years <- 2016:2025
  seasons <- c("wet", "dry")
  grid <- expand.grid(year = years, season = seasons, stringsAsFactors = FALSE)
  grid$habitat_class_file <- file.path(
    DW_INPUT_RASTER_DIR,
    sprintf("DW_class_%s_%d_project.tif", grid$season, grid$year)
  )
  grid$pressure_file <- file.path(
    DW_INPUT_RASTER_DIR,
    sprintf("DW_pressure_%s_%d_project.tif", grid$season, grid$year)
  )
  grid
}

period_manifest <- data.frame(
  token = DW_PERIOD_TOKENS,
  stack_file = file.path(
    DW_INPUT_RASTER_DIR,
    sprintf("DW_connectivity_inputs_%s_project.tif", DW_PERIOD_TOKENS)
  ),
  class_file = file.path(
    DW_INPUT_RASTER_DIR,
    sprintf("DW_class_%s_project.tif", DW_PERIOD_TOKENS)
  ),
  pressure_file = file.path(
    DW_INPUT_RASTER_DIR,
    sprintf("DW_pressure_%s_project.tif", DW_PERIOD_TOKENS)
  ),
  stringsAsFactors = FALSE
)

# 9-band stack order, matching config.py's DW_DERIVED_BANDS exactly
DW_STACK_BAND_NAMES <- c(
  "natural_prob", "woody_prob", "grass_prob", "conversion_pressure_prob",
  "hard_conversion_prob", "bare_degradation_prob", "water_wetland_prob",
  "top1_prob", "valid_obs_count"
)

TRANSITION_FILES <- c(
  baseline_to_current = file.path(DW_INPUT_RASTER_DIR, "DW_transition_baseline_2016_2018_to_current_2022_2025_project.tif"),
  pre_to_current       = file.path(DW_INPUT_RASTER_DIR, "DW_transition_pre_2019_2021_to_current_2022_2025_project.tif")
)

MASK_FILES <- c(
  analysis_mask           = file.path(DW_INPUT_RASTER_DIR, "DW_connectivity_analysis_mask_project.tif"),
  natural_habitat_mask     = file.path(DW_INPUT_RASTER_DIR, "DW_natural_habitat_mask_current_2022_2025_project.tif"),
  high_quality_source_mask = file.path(DW_INPUT_RASTER_DIR, "DW_high_quality_source_mask_current_2022_2025_project.tif")
)

#' Check which expected input rasters are present locally.
#' @param strict If TRUE, stop() on any missing file. If FALSE, warn() and return the missing
#'   list -- use FALSE for the partial-download smoke test (Verification step 3 of the plan).
check_inputs_present <- function(strict = TRUE) {
  seasonal <- expand_seasonal_manifest()
  files <- c(
    period_manifest$stack_file, period_manifest$class_file, period_manifest$pressure_file,
    seasonal$habitat_class_file, seasonal$pressure_file,
    TRANSITION_FILES, MASK_FILES
  )
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    msg <- sprintf(
      "%d/%d expected input rasters missing from %s:\n%s",
      length(missing), length(files), DW_INPUT_RASTER_DIR,
      paste(missing, collapse = "\n")
    )
    if (strict) stop(msg) else warning(msg, call. = FALSE)
  } else {
    message(sprintf("All %d expected input rasters present in %s.", length(files), DW_INPUT_RASTER_DIR))
  }
  invisible(missing)
}

#' Read a habitat_class raster, verify it matches Objective 1's actual (empirically confirmed)
#' scheme, and normalize its NoData representation to real NA.
#'
#' Empirically confirmed 2026-07-15 against the real downloaded rasters (not just the notebook
#' source, which implied a `-9999` sentinel that turns out not to appear in practice -- likely
#' clipped to 0 by the export's small unsigned integer pixel type): habitat_class rasters carry
#' NO registered GDAL NoData tag (`terra::NAflag()` reads NaN) and contain literal value `0` for
#' "no classification / outside AOI / fails the QA gate" (0's pixel count on the current-period
#' raster matches DW_connectivity_analysis_mask_project.tif's own "invalid" pixel count exactly).
#' `-9999` itself never appears as literal data in the files actually downloaded. This function
#' converts 0 -> NA immediately so every downstream consumer (recode, QA coverage, metrics) can
#' rely on ordinary is.na() semantics without re-deriving this -- do not read habitat_class
#' rasters via bare terra::rast() elsewhere in this project.
read_habitat_raster <- function(path) {
  if (!file.exists(path)) stop("Habitat raster not found: ", path)
  r <- terra::rast(path)
  crs_code <- tryCatch(terra::crs(r, describe = TRUE)$code, error = function(e) NA)
  if (is.na(crs_code) || crs_code != "32736") {
    warning(path, ": CRS did not resolve to EPSG:32736 (got ", crs_code, ") -- verify before trusting downstream results")
  }
  vals <- terra::values(r, na.rm = TRUE)
  if (length(vals) > 0) {
    if (any(vals == NODATA_SENTINEL)) {
      warning(path, ": literal -9999 present -- unexpected given empirical findings, investigate before trusting this file")
    }
    if (any(!vals %in% 0:8)) {
      stop(path, ": value outside the 0-8 habitat class scheme found (0 = no classification/outside AOI, 1-8 = real classes)")
    }
  }
  r <- terra::subst(r, 0, NA)
  names(r) <- "habitat_class"
  r
}

#' Read a pressure_class raster (values 0/1/2 -- 0 = "Low" is a REAL class here, unlike
#' habitat_class's 0, so this function does NOT substitute 0 -> NA. Pressure rasters' own NoData
#' representation has not yet been empirically verified (no DW_pressure_*.tif files were in the
#' first download batch) -- verify against the same approach used for read_habitat_raster() above
#' once they're available, rather than assuming NoData looks the same across both raster types.
read_pressure_raster <- function(path) {
  if (!file.exists(path)) stop("Pressure raster not found: ", path)
  r <- terra::rast(path)
  vals <- terra::values(r, na.rm = TRUE)
  if (length(vals) > 0 && any(vals == NODATA_SENTINEL)) {
    warning(path, ": literal -9999 present -- unexpected given habitat_class's empirical findings, investigate")
  }
  r
}

#' Read the 9-band period probability/derived stack and apply DW_STACK_BAND_NAMES.
read_period_stack <- function(path) {
  if (!file.exists(path)) stop("Period stack raster not found: ", path)
  r <- terra::rast(path)
  if (terra::nlyr(r) != length(DW_STACK_BAND_NAMES)) {
    stop(path, ": expected ", length(DW_STACK_BAND_NAMES), " bands, found ", terra::nlyr(r))
  }
  names(r) <- DW_STACK_BAND_NAMES
  r
}

#' Read + reproject a site boundary GeoJSON to PROJECT_CRS.
#' Source GeoJSONs are WGS84 lon/lat (confirmed from their embedded coordinates) -- always
#' reproject before any raster crop/mask/area operation.
read_site_boundary <- function(site_id) {
  path <- AOI_PATHS[[site_id]]
  if (is.null(path)) stop("Unknown site_id: ", site_id)
  sf::st_read(path, quiet = TRUE) |> sf::st_transform(crs = PROJECT_CRS)
}
