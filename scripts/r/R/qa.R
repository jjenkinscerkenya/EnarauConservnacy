# Valid-pixel-coverage QA: a site-year with only ~55% valid pixel
# coverage produced artificial patch breaks that inflated NP/PD/ED. No standalone valid_obs_count
# raster exists per-year/seasonal (only inside the 9-band period stacks) -- masking already
# happened upstream at classification time (Objective 1's min-observation/top1_prob gates), so the
# NA fraction on the habitat_class raster itself IS the QA signal; no separate QA band needed.

#' Compute valid-pixel coverage for one habitat_class raster, restricted to one site's polygon.
#' Uses read_habitat_raster() (not bare terra::rast()) so literal 0 ("no classification /
#' outside AOI", see that function's docs) is correctly counted as invalid, not valid -- reading
#' the raw file directly here would silently overcount coverage.
#' @param habitat_class_path Path to a habitat_class GeoTIFF.
#' @param site_vect A terra SpatVector (or sf object) of the site boundary, already in PROJECT_CRS.
compute_valid_coverage <- function(habitat_class_path, site_vect) {
  r <- read_habitat_raster(habitat_class_path)
  if (inherits(site_vect, "sf")) site_vect <- terra::vect(site_vect)
  r_crop <- terra::crop(r, site_vect)
  vals <- terra::extract(r_crop, site_vect, ID = FALSE)[[1]]
  total <- length(vals)
  valid <- sum(!is.na(vals))
  data.frame(
    valid_pixel_count = valid,
    total_pixel_count = total,
    coverage_pct = if (total > 0) valid / total else NA_real_,
    below_threshold = if (total > 0) (valid / total) < VALID_PIXEL_COVERAGE_MIN else NA
  )
}

#' Run compute_valid_coverage() across every (site x year x season) row of a manifest data.frame
#' (expects a `habitat_class_file` column plus identifying columns to carry through).
compute_valid_coverage_table <- function(manifest, id_cols) {
  rows <- lapply(seq_len(nrow(manifest)), function(i) {
    site_ids <- SITES$site_id
    do.call(rbind, lapply(site_ids, function(sid) {
      site_vect <- terra::vect(read_site_boundary(sid))
      cov <- compute_valid_coverage(manifest$habitat_class_file[i], site_vect)
      cbind(
        site_id = sid,
        manifest[i, id_cols, drop = FALSE],
        cov
      )
    }))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Print the excluded (below-threshold) rows loudly -- never a silent drop.
report_excluded_rows <- function(coverage_table, context_label) {
  excluded <- coverage_table[coverage_table$below_threshold %in% TRUE, ]
  if (nrow(excluded) > 0) {
    message(sprintf(
      "[QA] %d row(s) excluded from %s for <%.0f%% valid-pixel coverage:",
      nrow(excluded), context_label, 100 * VALID_PIXEL_COVERAGE_MIN
    ))
    print(excluded)
  } else {
    message(sprintf("[QA] No rows excluded from %s -- all site/period-season combinations meet the %.0f%% coverage threshold.",
                     context_label, 100 * VALID_PIXEL_COVERAGE_MIN))
  }
  invisible(excluded)
}
