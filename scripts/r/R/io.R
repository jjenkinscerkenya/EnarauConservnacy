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
#' Corrected 2026-07-18: every exported raster DOES carry a correctly-registered `-9999` GDAL
#' NoData tag (confirmed via direct gdalinfo-style inspection) -- an earlier
#' `terra::NAflag()`-based check that suggested otherwise was itself misleading (NAflag() reports
#' a user-set override, not whether the file's own embedded tag is present/honored; terra
#' silently applies the embedded tag on read, which is why literal `-9999` never shows up in
#' `terra::values()`). The real, separate defect (root-caused to an `eetools.io.
#' export_image_to_drive()` bug, since fixed upstream for future exports -- see the write-back
#' note `2026-07-17-gee-utils-eetools-unmask-samefootprint-bug.md`) is that habitat_class/
#' pressure_class rasters exported directly after a terminal `.clip(project_geom)` also contain a
#' literal `0` fill OUTSIDE the true (irregular, buffered) AOI polygon but inside the export's
#' rectangular bounding box -- a hard-coded Earth Engine exporter fallback that bypasses the
#' NoData tag entirely. For habitat_class this is harmless to detect: `0` is never a real class
#' (classes are 1-8), so it's unambiguous regardless of position. This function converts `0 ->
#' NA` immediately so every downstream consumer (recode, QA coverage, metrics) can rely on
#' ordinary is.na() semantics without re-deriving this -- do not read habitat_class rasters via
#' bare terra::rast() elsewhere in this project. (Contrast with `read_pressure_raster()`, where
#' `0` is a real class value and this trick does NOT work -- that function needs an actual
#' geometric clip instead.)
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
#' habitat_class's 0). Resolved 2026-07-18: pressure_class rasters have the same literal-0
#' exterior-fringe defect as habitat_class (see read_habitat_raster()'s docs for the root cause),
#' but the 0->NA substitution trick that works for habitat_class does NOT work here, since 0 is a
#' real "Low" pressure class both inside the true AOI and in the bogus fringe -- indistinguishable
#' by value alone. This function instead masks to the reconstructed project_geom polygon
#' (build_project_geom()/project_geom_vect(), below), converting everything outside the true AOI
#' to NA regardless of value; real interior 0/1/2 survive untouched. The upstream eetools bug
#' this works around has since been fixed (see the write-back note
#' 2026-07-17-gee-utils-eetools-unmask-samefootprint-bug.md) -- this clip is only needed for
#' rasters exported before that fix landed; a fresh re-export would not need it.
read_pressure_raster <- function(path) {
  if (!file.exists(path)) stop("Pressure raster not found: ", path)
  r <- terra::rast(path)
  r <- terra::mask(r, project_geom_vect())
  vals <- terra::values(r, na.rm = TRUE)
  if (length(vals) > 0) {
    if (any(vals == NODATA_SENTINEL)) {
      warning(path, ": literal -9999 present -- unexpected, investigate")
    }
    if (any(!vals %in% 0:2)) {
      stop(path, ": value outside the 0-2 pressure class scheme found (0=Low, 1=Moderate, 2=High)")
    }
  }
  names(r) <- "pressure_class"
  r
}

#' Read the 9-band period probability/derived stack and apply DW_STACK_BAND_NAMES.
#' Confirmed 2026-07-18: does NOT have the exterior-fringe defect (unlike habitat_class/
#' pressure_class) -- these rasters are built via .select() on an .addBands()-augmented
#' composite, and every one examined comes out correctly masked with no code change needed here.
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

# Reconstruction of the notebook's `project_geom` (union of the 4 SITES polygons, buffered by
# STUDY_AREA_BUFFER_M), used to work around the pressure_class exterior-fringe defect -- see
# read_pressure_raster()'s docs. This is a PLANAR buffer (sf::st_buffer() in EPSG:32736), while
# GEE's own `.buffer(150, maxError=10)` is a GEODESIC buffer with up to 10m tolerance -- expect a
# sub-pixel-to-few-pixel discrepancy at the boundary, immaterial at 10m resolution for masking
# purposes, not a correctness concern for this use.
.project_geom_cache <- new.env(parent = emptyenv())

#' Build (or return the cached) reconstructed project_geom polygon, as an sf object in
#' PROJECT_CRS. Memoized per R session since every numbered script reads many rasters that each
#' need it -- call with force_refresh = TRUE only if SITES/AOI_PATHS somehow change mid-session
#' (never happens in normal script execution).
build_project_geom <- function(force_refresh = FALSE) {
  if (!force_refresh && !is.null(.project_geom_cache$geom)) return(.project_geom_cache$geom)
  geoms <- do.call(c, lapply(SITES$site_id, function(sid) sf::st_geometry(read_site_boundary(sid))))
  buffered <- sf::st_buffer(sf::st_union(sf::st_make_valid(geoms)), dist = STUDY_AREA_BUFFER_M)
  .project_geom_cache$geom <- sf::st_sf(geometry = buffered)
  .project_geom_cache$geom
}

#' terra::vect() version of build_project_geom(), for direct use in terra::mask().
project_geom_vect <- function(force_refresh = FALSE) {
  terra::vect(build_project_geom(force_refresh))
}

#' Read a transition_code raster (from_class*10 + to_class, valid codes are both digits 1-8,
#' i.e. 11-88) and validate it. Confirmed 2026-07-18: these rasters do NOT have the
#' exterior-fringe defect (built via arithmetic on two already-clipped habitat_class images,
#' which for reasons not fully pinned down avoids the bug -- see the write-back note) -- no clip
#' needed here, just the same validation discipline read_habitat_raster() already applies.
read_transition_raster <- function(path) {
  if (!file.exists(path)) stop("Transition raster not found: ", path)
  r <- terra::rast(path)
  vals <- terra::values(r, na.rm = TRUE)
  if (length(vals) > 0) {
    valid_codes <- as.vector(outer(1:8, 1:8, function(a, b) a * 10 + b))
    if (any(!vals %in% valid_codes)) {
      stop(path, ": transition code outside the valid from*10+to (1-8 x 1-8) scheme found -- ",
           "investigate before trusting this file")
    }
  }
  names(r) <- "transition_code"
  r
}
