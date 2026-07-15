# Objective 3, step 4: moving-window local connectivity change maps (Level 3a).
#
# RUN AS: cd scripts/r && Rscript 04_moving_window_connectivity.R
#
# PROJECT-WIDE FIRST, clip only for reporting/figures afterward -- a real corridor pinch-point
# can straddle a site boundary; masking to site before window_lsm() would truncate it. See the
# masking-order rule in the plan.
#
# window_lsm() is the single most expensive step in this pipeline (PD/CLUMPY require per-window
# patch delineation, far more costly than PLAND/ED's focal-sum approach). DO NOT run this
# project-wide/full-radius on a whim -- run the smoke test below first (Verification step 5 of
# the plan) and read off the printed timing estimate before committing to a full run.

source("00_config.R")
source("R/io.R")
source("R/recode.R")

# ---- Smoke-test controls: set SMOKE_TEST_SITE to a site_id to crop to that site + a 150 m
# buffer before timing a single radius, instead of running the full project extent. Leave NULL
# for the real production run only after the timing smoke test has been reviewed. ----
SMOKE_TEST_SITE <- "corridor_p1"   # set to NULL for the full production run
SMOKE_TEST_RADII_M <- 500          # single radius to time during the smoke test

message("=== 04_moving_window_connectivity ===")

current_class_path <- period_manifest$class_file[period_manifest$token == "current_2022_2025"]
baseline_class_path <- period_manifest$class_file[period_manifest$token == "baseline_2016_2018"]

if (!file.exists(current_class_path) || !file.exists(baseline_class_path)) {
  stop("Current and/or baseline period class rasters not yet downloaded -- run 01_prepare_inputs.R first and check the manifest.")
}

current_r <- read_habitat_raster(current_class_path)
baseline_r <- read_habitat_raster(baseline_class_path)

current_bin <- make_natural_binary(current_r)
baseline_bin <- make_natural_binary(baseline_r)

if (!is.null(SMOKE_TEST_SITE)) {
  message("SMOKE TEST MODE: cropping to '", SMOKE_TEST_SITE, "' + 150 m buffer, radius = ", SMOKE_TEST_RADII_M, " m only.")
  site_boundary <- read_site_boundary(SMOKE_TEST_SITE)
  site_buffered <- sf::st_buffer(site_boundary, dist = 150)
  current_bin <- terra::crop(current_bin, terra::vect(site_buffered))
  baseline_bin <- terra::crop(baseline_bin, terra::vect(site_buffered))
  radii_to_run <- SMOKE_TEST_RADII_M
} else {
  radii_to_run <- MOVING_WINDOW_RADII_M
}

window_metrics <- c("lsm_l_pland", "lsm_l_ed", "lsm_l_clumpy", "lsm_l_pd")

run_window_lsm_for_period <- function(r_bin, radius_m) {
  win_size <- radius_m %/% 10 + 1  # matches the vault plan doc's own 500m->51x51 reference example
  win <- matrix(1, nrow = win_size, ncol = win_size)
  t0 <- Sys.time()
  result <- landscapemetrics::window_lsm(
    landscape = r_bin, window = win, what = window_metrics, directions = 8, progress = TRUE
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("  window_lsm() at radius %dm (%dx%d cells) took %.1f sec.", radius_m, win_size, win_size, elapsed))
  result
}

for (radius_m in radii_to_run) {
  message("=== Radius ", radius_m, "m ===")
  current_windows <- run_window_lsm_for_period(current_bin, radius_m)
  baseline_windows <- run_window_lsm_for_period(baseline_bin, radius_m)

  if (!is.null(SMOKE_TEST_SITE)) {
    message("Smoke test complete for radius ", radius_m, "m -- review timing above before running the full pipeline. ",
            "No output rasters written in smoke-test mode.")
    next
  }

  # window_lsm() returns one raster layer per metric, named after the metric.
  pland_change <- current_windows[["lsm_l_pland"]] - baseline_windows[["lsm_l_pland"]]
  ed_change    <- current_windows[["lsm_l_ed"]] - baseline_windows[["lsm_l_ed"]]
  pd_change    <- current_windows[["lsm_l_pd"]] - baseline_windows[["lsm_l_pd"]]
  clumpy_change <- current_windows[["lsm_l_clumpy"]] - baseline_windows[["lsm_l_clumpy"]]

  scale_raster <- function(r) (r - terra::global(r, "mean", na.rm = TRUE)[1, 1]) / terra::global(r, "sd", na.rm = TRUE)[1, 1]
  local_connectivity_change_score <- scale_raster(pland_change) + scale_raster(clumpy_change) -
    scale_raster(ed_change) - scale_raster(pd_change)
  names(local_connectivity_change_score) <- "local_connectivity_change_score"

  suffix <- sprintf("_w%dm", radius_m)
  terra::writeRaster(current_windows[["lsm_l_pland"]], file.path(LANDSCAPE_RASTER_DIR, paste0("local_natural_prop_current", suffix, ".tif")), overwrite = TRUE)
  terra::writeRaster(baseline_windows[["lsm_l_pland"]], file.path(LANDSCAPE_RASTER_DIR, paste0("local_natural_prop_baseline", suffix, ".tif")), overwrite = TRUE)
  terra::writeRaster(pland_change, file.path(LANDSCAPE_RASTER_DIR, paste0("local_natural_prop_change_baseline_to_current", suffix, ".tif")), overwrite = TRUE)
  terra::writeRaster(ed_change, file.path(LANDSCAPE_RASTER_DIR, paste0("local_edge_density_change_baseline_to_current", suffix, ".tif")), overwrite = TRUE)
  terra::writeRaster(pd_change, file.path(LANDSCAPE_RASTER_DIR, paste0("local_patch_density_change_baseline_to_current", suffix, ".tif")), overwrite = TRUE)
  terra::writeRaster(local_connectivity_change_score, file.path(LANDSCAPE_RASTER_DIR, paste0("local_connectivity_change_score", suffix, ".tif")), overwrite = TRUE)

  message("Wrote radius-", radius_m, "m rasters to ", LANDSCAPE_RASTER_DIR)
}

message("=== 04_moving_window_connectivity complete ===")
if (!is.null(SMOKE_TEST_SITE)) {
  message("This was a SMOKE TEST run (SMOKE_TEST_SITE = '", SMOKE_TEST_SITE, "'). ",
          "Set SMOKE_TEST_SITE <- NULL at the top of this script for the full production run, ",
          "after confirming the timing above extrapolates to an acceptable full-extent runtime.")
}
