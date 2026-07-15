# Objective 3, step 3: site-level fragmentation/connectivity metrics (Level 2).
#
# RUN AS: cd scripts/r && Rscript 03_landscape_metrics.R
#
# MASKS TO SITE FIRST (crop + mask each raster to each site polygon before calculate_lsm()) --
# this is Level 2, answering "what does fragmentation look like within this site." Level 3
# (04/05) must NOT follow this pattern -- see the masking-order rule in the plan.
#
# Processes only the seasonal/period rasters actually present locally (missing files are skipped
# with a message, not a hard failure) -- run 01_prepare_inputs.R first to see what's missing.

source("00_config.R")
source("R/io.R")
source("R/recode.R")
source("R/qa.R")
source("R/metrics.R")

message("=== 03_landscape_metrics ===")

coverage_year_season_path <- file.path(TABLES_DIR, "landscape_valid_pixel_coverage_by_site_year_season.csv")
coverage_period_path <- file.path(TABLES_DIR, "landscape_valid_pixel_coverage_by_site_period.csv")
coverage_year_season <- if (file.exists(coverage_year_season_path)) readr::read_csv(coverage_year_season_path, show_col_types = FALSE) else NULL
coverage_period <- if (file.exists(coverage_period_path)) readr::read_csv(coverage_period_path, show_col_types = FALSE) else NULL

site_boundaries <- lapply(SITES$site_id, read_site_boundary)
names(site_boundaries) <- SITES$site_id

#' Run class + binary + entropy metrics for one habitat_class raster path, across all 4 sites,
#' excluding any site below the coverage threshold (if a coverage table is supplied).
run_metrics_for_raster <- function(habitat_class_path, id_values, coverage_table = NULL) {
  r <- read_habitat_raster(habitat_class_path)
  r_full <- recode_full_habitat(r)
  r_bin  <- make_natural_binary(r)

  rows <- lapply(SITES$site_id, function(sid) {
    if (!is.null(coverage_table)) {
      match_rows <- coverage_table[coverage_table$site_id == sid, , drop = FALSE]
      for (nm in names(id_values)) match_rows <- match_rows[match_rows[[nm]] == id_values[[nm]], , drop = FALSE]
      if (nrow(match_rows) > 0 && isTRUE(match_rows$below_threshold[1])) {
        message("  Skipping ", sid, " (", paste(id_values, collapse = "/"), ") -- below valid-pixel coverage threshold.")
        return(NULL)
      }
    }
    site_vect <- terra::vect(site_boundaries[[sid]])
    r_full_site <- mask_to_site(r_full, site_vect)
    r_bin_site  <- mask_to_site(r_bin, site_vect)

    class_m <- calculate_class_metrics(r_full_site)
    class_m <- cbind(site_id = sid, class_m, as.data.frame(id_values))

    binary_m <- calculate_binary_metrics(r_bin_site)
    binary_m <- cbind(site_id = sid, binary_m, as.data.frame(id_values))

    entropy_full <- calculate_entropy_pilot(r_full_site)
    entropy_full <- cbind(site_id = sid, landscape_type = "full_class", entropy_full, as.data.frame(id_values))
    entropy_bin <- calculate_entropy_pilot(r_bin_site)
    entropy_bin <- cbind(site_id = sid, landscape_type = "binary_natural", entropy_bin, as.data.frame(id_values))

    list(class = class_m, binary = binary_m, entropy = dplyr::bind_rows(entropy_full, entropy_bin))
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  list(
    class = dplyr::bind_rows(lapply(rows, `[[`, "class")),
    binary = dplyr::bind_rows(lapply(rows, `[[`, "binary")),
    entropy = dplyr::bind_rows(lapply(rows, `[[`, "entropy"))
  )
}

seasonal_manifest <- expand_seasonal_manifest()
seasonal_available <- seasonal_manifest[file.exists(seasonal_manifest$habitat_class_file), ]
period_available <- period_manifest[file.exists(period_manifest$class_file), ]

all_class <- list(); all_binary <- list(); all_entropy <- list()

if (nrow(seasonal_available) > 0) {
  message("Processing ", nrow(seasonal_available), " seasonal per-year rasters...")
  for (i in seq_len(nrow(seasonal_available))) {
    res <- run_metrics_for_raster(
      seasonal_available$habitat_class_file[i],
      id_values = list(year = seasonal_available$year[i], season = seasonal_available$season[i]),
      coverage_table = coverage_year_season
    )
    all_class[[length(all_class) + 1]] <- res$class
    all_binary[[length(all_binary) + 1]] <- res$binary
    all_entropy[[length(all_entropy) + 1]] <- res$entropy
  }
} else {
  message("No seasonal rasters available yet -- skipping seasonal metrics.")
}

if (nrow(period_available) > 0) {
  message("Processing ", nrow(period_available), " period composite rasters...")
  for (i in seq_len(nrow(period_available))) {
    res <- run_metrics_for_raster(
      period_available$class_file[i],
      id_values = list(period = period_available$token[i]),
      coverage_table = coverage_period
    )
    all_class[[length(all_class) + 1]] <- res$class
    all_binary[[length(all_binary) + 1]] <- res$binary
    all_entropy[[length(all_entropy) + 1]] <- res$entropy
  }
} else {
  message("No period composite rasters available yet -- skipping period metrics.")
}

# bind_rows (not rbind) -- seasonal rows carry year/season columns, period rows carry a period
# column instead; rbind() can't reconcile the mismatched schemas but bind_rows() aligns them and
# fills the other set's columns with NA, which is exactly what downstream code expects (e.g. the
# metric-change summary below filters on !is.na(period), which only works if seasonal rows
# contribute a real (NA) period column rather than the column being entirely absent).
class_metrics <- dplyr::bind_rows(all_class)
binary_metrics <- dplyr::bind_rows(all_binary)
entropy_metrics <- dplyr::bind_rows(all_entropy)

if (!is.null(class_metrics) && nrow(class_metrics) > 0) {
  readr::write_csv(class_metrics, file.path(TABLES_DIR, "landscape_fragmentation_metrics_full_class_by_site_year_season.csv"))
}
if (!is.null(binary_metrics) && nrow(binary_metrics) > 0) {
  readr::write_csv(binary_metrics, file.path(TABLES_DIR, "landscape_connectivity_metrics_binary_natural_by_site_year_season.csv"))
}
if (!is.null(entropy_metrics) && nrow(entropy_metrics) > 0) {
  readr::write_csv(entropy_metrics, file.path(TABLES_DIR, "landscape_entropy_pilot_by_site_year_season.csv"))
}

# ---- Correlation screen (binary-natural metric set + entropy pilot, pooled across all series) ----
if (!is.null(binary_metrics) && nrow(binary_metrics) > 0) {
  message("=== Correlation screen ===")
  id_cols <- intersect(c("site_id", "year", "season", "period"), names(binary_metrics))
  screen <- screen_metric_correlation(binary_metrics, id_cols = id_cols)
  readr::write_csv(
    as.data.frame(as.table(screen$correlation_matrix)),
    file.path(TABLES_DIR, "landscape_metric_correlation_matrix.csv")
  )
  if (nrow(screen$flagged_pairs) > 0) {
    message(nrow(screen$flagged_pairs), " metric pair(s) flagged with |r| > ", CORRELATION_FLAG_THRESHOLD, " -- review before finalizing the reporting metric set:")
    print(screen$flagged_pairs)
  } else {
    message("No metric pairs exceeded the |r| > ", CORRELATION_FLAG_THRESHOLD, " correlation flag threshold.")
  }
}

# ---- Metric change summary (period-to-period deltas) ----
if (!is.null(binary_metrics) && "period" %in% names(binary_metrics) && nrow(binary_metrics) > 0) {
  period_binary <- binary_metrics[!is.na(binary_metrics$period), ]
  wide <- tidyr::pivot_wider(
    period_binary[, c("site_id", "period", "metric", "value")],
    names_from = period, values_from = value
  )
  if (all(c("baseline_2016_2018", "current_2022_2025") %in% names(wide))) {
    wide$change_baseline_to_current <- wide$current_2022_2025 - wide$baseline_2016_2018
    readr::write_csv(
      wide[, c("site_id", "metric", "baseline_2016_2018", "current_2022_2025", "change_baseline_to_current")],
      file.path(TABLES_DIR, "landscape_metric_change_baseline_to_current.csv")
    )
  }
  if (all(c("pre_2019_2021", "current_2022_2025") %in% names(wide))) {
    wide$change_pre_to_current <- wide$current_2022_2025 - wide$pre_2019_2021
    readr::write_csv(
      wide[, c("site_id", "metric", "pre_2019_2021", "current_2022_2025", "change_pre_to_current")],
      file.path(TABLES_DIR, "landscape_metric_change_pre_to_current.csv")
    )
  }
}

message("=== 03_landscape_metrics complete ===")
