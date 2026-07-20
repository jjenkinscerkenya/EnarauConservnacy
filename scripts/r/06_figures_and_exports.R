# Objective 3, step 6: linkage priority score synthesis, Objective 2 cross-check.
#
# RUN AS: cd scripts/r && Rscript 06_figures_and_exports.R
#
# Depends on 02 (transition tables), 03 (site metrics), 04 (moving-window rasters -- at the 500m
# radius specifically, see LINKAGE_SCORE_RADIUS_M below), and 05 (patch importance/graph). Run
# those first; this script skips any component gracefully (with a message) if its upstream
# output isn't present yet, but the final linkage_priority_score.tif requires all four.
#
# Figures are NOT generated here -- all Objective 3 plotting lives in
# scripts/python/plotting/landscape_metrics.ipynb (pandas/matplotlib/seaborn), reading this
# script's and 03's table outputs from outputs/tables/, for consistency with Objectives 1-2's
# plotting convention (kept in Python, not duplicated per-language).

source("00_config.R")
source("R/io.R")
source("R/recode.R")
source("R/scoring.R")

message("=== 06_figures_and_exports ===")

LINKAGE_SCORE_RADIUS_M <- 500  # which 04 radius feeds local_connectivity_score below

current_class_path <- period_manifest$class_file[period_manifest$token == "current_2022_2025"]
current_stack_path <- period_manifest$stack_file[period_manifest$token == "current_2022_2025"]
transition_raster_path <- TRANSITION_FILES["baseline_to_current"]

have_current_class <- file.exists(current_class_path)
have_current_stack <- file.exists(current_stack_path)
have_transition_raster <- file.exists(transition_raster_path)
have_local_connectivity <- file.exists(file.path(LANDSCAPE_RASTER_DIR, sprintf("local_connectivity_change_score_w%dm.tif", LINKAGE_SCORE_RADIUS_M)))
have_patch_importance <- file.exists(file.path(VECTORS_DIR, "natural_patches_current.gpkg"))

if (have_current_class && have_transition_raster && have_local_connectivity && have_patch_importance && have_current_stack) {

  current_r <- read_habitat_raster(current_class_path)
  current_bin <- make_natural_binary(current_r)
  grid_template <- current_bin

  # ---- 1. Persistent or recovered natural habitat (from the transition raster directly) ----
  transition_r <- read_transition_raster(transition_raster_path)
  stable_natural_codes <- c(11, 12, 13, 21, 22, 23, 31, 32, 33)
  recovered_codes <- c(41, 42, 43, 51, 52, 53, 61, 62, 63)
  persistent_or_recovered <- terra::classify(
    transition_r,
    rcl = matrix(c(
      cbind(stable_natural_codes, 1),
      cbind(recovered_codes, 1)
    ), ncol = 2)
  )
  persistent_or_recovered <- terra::ifel(is.na(persistent_or_recovered), 0, persistent_or_recovered)
  persistent_or_recovered <- terra::resample(persistent_or_recovered, grid_template, method = "near")

  # ---- 2. Local connectivity change score (from 04, resampled to the current grid) ----
  local_connectivity_r <- terra::rast(file.path(LANDSCAPE_RASTER_DIR, sprintf("local_connectivity_change_score_w%dm.tif", LINKAGE_SCORE_RADIUS_M)))
  local_connectivity_r <- terra::resample(local_connectivity_r, grid_template, method = "bilinear")
  local_connectivity_score <- normalize01(terra::values(local_connectivity_r))
  local_connectivity_r <- terra::setValues(grid_template, local_connectivity_score)

  # ---- 3. Patch importance (rasterized from 05's polygons) ----
  patches_sf <- sf::st_read(file.path(VECTORS_DIR, "natural_patches_current.gpkg"), quiet = TRUE)
  patch_importance_r <- terra::rasterize(terra::vect(patches_sf), grid_template, field = "patch_importance_score", background = 0)

  # ---- 4. Low conversion pressure (from the current-period 9-band stack) ----
  current_stack <- read_period_stack(current_stack_path)
  low_pressure_r <- 1 - terra::resample(current_stack[["conversion_pressure_prob"]], grid_template, method = "bilinear")

  # ---- 5. Bottleneck relevance (rasterized betweenness + articulation-point boost, mid threshold) ----
  bottleneck_r <- terra::rasterize(
    terra::vect(patches_sf), grid_template,
    field = "betweenness", background = 0
  )
  bottleneck_r <- normalize01(terra::values(bottleneck_r))
  bottleneck_r <- terra::setValues(grid_template, bottleneck_r)

  # ---- 6. Corridor proximity (distance decay from corridor_p1 U corridor_p2) ----
  corridor_sf <- do.call(rbind, lapply(c("corridor_p1", "corridor_p2"), function(sid) {
    b <- read_site_boundary(sid)
    sf::st_sf(site_id = sid, geometry = sf::st_geometry(sf::st_union(b)))
  }))
  corridor_proximity_r <- distance_decay_score(grid_template, corridor_sf)

  # ---- Combine ----
  linkage_priority_score <- compute_linkage_priority_score(
    persistent_or_recovered_natural = persistent_or_recovered,
    local_connectivity = local_connectivity_r,
    patch_importance = normalize01(terra::values(patch_importance_r)) |> (\(v) terra::setValues(grid_template, v))(),
    low_conversion_pressure = low_pressure_r,
    bottleneck_relevance = bottleneck_r,
    corridor_proximity = corridor_proximity_r
  )
  names(linkage_priority_score) <- "linkage_priority_score"
  terra::writeRaster(linkage_priority_score, file.path(LANDSCAPE_RASTER_DIR, "linkage_priority_score.tif"), overwrite = TRUE)

  # Candidate linkage areas: top quartile of the score, dissolved to polygons for reporting.
  score_threshold <- stats::quantile(terra::values(linkage_priority_score, na.rm = TRUE), 0.75, na.rm = TRUE)
  candidate_mask <- terra::ifel(linkage_priority_score >= score_threshold, 1, NA)
  candidate_poly <- terra::as.polygons(candidate_mask, dissolve = TRUE, na.rm = TRUE)
  sf::st_write(sf::st_as_sf(candidate_poly), file.path(VECTORS_DIR, "candidate_linkage_areas.gpkg"), delete_dsn = TRUE, quiet = TRUE)

  message("Wrote linkage_priority_score.tif and candidate_linkage_areas.gpkg.")
} else {
  message("Skipping linkage_priority_score synthesis -- missing one or more upstream inputs:")
  message("  current class raster: ", have_current_class)
  message("  current 9-band stack: ", have_current_stack)
  message("  baseline->current transition raster: ", have_transition_raster)
  message("  04's local_connectivity_change_score_w", LINKAGE_SCORE_RADIUS_M, "m.tif: ", have_local_connectivity)
  message("  05's natural_patches_current.gpkg: ", have_patch_importance)
  message("Run 02-05 first (04 needs a full, non-smoke-test run at radius ", LINKAGE_SCORE_RADIUS_M, "m).")
}

# ---- Objective 2 cross-check (table-level join, not per-pixel raster math -- 30m vs 10m
# resolution mismatch; see the plan's Objective 2->3 handoff notes) ----
message("=== Objective 2 cross-check ===")

net_natural_change_path <- file.path(TABLES_DIR, "landscape_net_natural_change_by_site.csv")
nbr_event_path <- file.path(TABLES_DIR, "landtrendr_nbrseg_event_summary_by_site.csv")
msavi2_event_path <- file.path(TABLES_DIR, "landtrendr_msavi2seg_event_summary_by_site.csv")

if (file.exists(net_natural_change_path) && file.exists(nbr_event_path) && file.exists(msavi2_event_path)) {
  net_natural_change <- readr::read_csv(net_natural_change_path, show_col_types = FALSE)
  nbr_events <- readr::read_csv(nbr_event_path, show_col_types = FALSE)
  msavi2_events <- readr::read_csv(msavi2_event_path, show_col_types = FALSE)

  # area_ha_sum in these CSVs is raw hectares, NOT pre-normalized by site area -- normalize here.
  site_areas_ha <- net_natural_change$site_area_ha
  names(site_areas_ha) <- net_natural_change$site_id

  normalize_landtrendr <- function(df, run_label) {
    if (!"area_ha_sum" %in% names(df)) return(NULL)
    df$pct_of_site_area <- 100 * df$area_ha_sum / site_areas_ha[df$site_id]
    df$run <- run_label
    df
  }
  nbr_norm <- normalize_landtrendr(nbr_events, "nbr_dry")
  msavi2_norm <- normalize_landtrendr(msavi2_events, "msavi2_wet")
  landtrendr_normalized <- rbind(nbr_norm, msavi2_norm)

  crosscheck <- net_natural_change |>
    dplyr::select(site_id, site_name, net_change_baseline_to_current_pct_of_site, net_change_pre_to_current_pct_of_site) |>
    dplyr::left_join(
      landtrendr_normalized |> dplyr::select(site_id, run, change_type, pct_of_site_area) |>
        tidyr::pivot_wider(names_from = c(run, change_type), values_from = pct_of_site_area, names_prefix = "landtrendr_"),
      by = "site_id"
    )

  # Load-bearing caveat check: does Enarau's natural-habitat-expansion signal also appear at
  # Mbokishi (the untouched reference site)? If so, it's more likely regional/rainfall-driven
  # than Enarau-specific conservation activity -- a site *diverging* from Mbokishi is the
  # defensible localized-effect claim (per the Objective 2->3 handoff notes).
  mbokishi_change <- crosscheck$net_change_baseline_to_current_pct_of_site[crosscheck$site_id == "mbokishi"]
  crosscheck$diverges_from_mbokishi_pct_points <- crosscheck$net_change_baseline_to_current_pct_of_site - mbokishi_change
  crosscheck$caveat <- ifelse(
    crosscheck$site_id != "mbokishi" & sign(crosscheck$net_change_baseline_to_current_pct_of_site) == sign(mbokishi_change) &
      abs(crosscheck$diverges_from_mbokishi_pct_points) < 2,
    "CAUTION: change direction/magnitude matches Mbokishi (untouched reference site) within 2 pct points -- likely regional/rainfall-driven, not necessarily Enarau-specific conservation effect.",
    "Diverges from Mbokishi's own trajectory -- more defensible as a localized effect."
  )

  readr::write_csv(crosscheck, file.path(TABLES_DIR, "landscape_vs_objective2_crosscheck_by_site.csv"))
  message("Wrote landscape_vs_objective2_crosscheck_by_site.csv. Caveat flags:")
  print(crosscheck[, c("site_id", "net_change_baseline_to_current_pct_of_site", "caveat")])
} else {
  message("Skipping Objective 2 cross-check -- run 02_area_transition_summaries.R first, and confirm ",
          "landtrendr_nbrseg_event_summary_by_site.csv / landtrendr_msavi2seg_event_summary_by_site.csv ",
          "are present in outputs/tables/ (from Objective 2).")
}

message("=== 06_figures_and_exports complete ===")
message("Figures: run scripts/python/plotting/landscape_metrics.ipynb to (re)generate plots from ",
        "this script's and 03's table outputs.")
