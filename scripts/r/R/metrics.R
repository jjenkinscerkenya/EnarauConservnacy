# calculate_lsm() wrappers, correlation screen, and the entropy-metric pilot check
# (Nowosad & Stepinski 2019 -- lsm_l_ent/lsm_l_mutinf as a cheap, weakly-correlated complexity
# axis, per the vault plan doc's open question).
#
# IMPORTANT interpretation note (confirmed empirically 2026-07-15): PLAND/ED/AI/etc. are
# percentages of the VALID (classified, non-NA) landscape extent within each site crop, not of
# the site's full nominal polygon area. Roughly 15% of Enarau's own polygon has no valid
# classification at all (cross-checked against Objective 1's own dw_area_by_class_by_site_period
# .csv: summing habitat_area_ha_1..8 gives ~982 ha against a 1161 ha polygon). This means a
# PLAND of e.g. 89% for Enarau means "89% of the CLASSIFIED area is natural," not "89% of Enarau."
# 01_prepare_inputs.R's valid-pixel-coverage table gives the classified/polygon ratio needed to
# reconcile the two -- always report both when quoting a PLAND-style percentage in the final
# deliverable.

#' Mask a project-wide raster to one site's polygon (crop + mask). Used by 03 (Level 2), which
#' masks to site BEFORE computing metrics -- do NOT reuse this for Level 3 (04/05), which must
#' compute on the full project extent first and clip only for reporting (see plan's masking-order
#' rule: a real patch/pinch-point can straddle a site boundary).
mask_to_site <- function(r, site_vect) {
  if (inherits(site_vect, "sf")) site_vect <- terra::vect(site_vect)
  r_crop <- terra::crop(r, site_vect)
  terra::mask(r_crop, site_vect)
}

#' Class-level fragmentation metrics on the full habitat-class raster (masked to site).
calculate_class_metrics <- function(r_site) {
  landscapemetrics::calculate_lsm(
    landscape = r_site,
    what = SELECTED_CLASS_METRICS,
    directions = 8,
    edge_depth = EDGE_DEPTH_CELLS
  )
}

#' Class-level connectivity metrics on the binary natural-habitat raster (masked to site),
#' filtered to class == 1 (natural).
calculate_binary_metrics <- function(r_site_binary) {
  m <- landscapemetrics::calculate_lsm(
    landscape = r_site_binary,
    what = SELECTED_BINARY_METRICS,
    directions = 8,
    edge_depth = EDGE_DEPTH_CELLS
  )
  m[m$class == 1, ]
}

#' Landscape-level entropy pilot (cheap -- no patch delineation required).
calculate_entropy_pilot <- function(r_site) {
  landscapemetrics::calculate_lsm(
    landscape = r_site,
    what = c("lsm_l_ent", "lsm_l_mutinf"),
    directions = 8
  )
}

#' Pivot a long metrics table (columns: site_id, ..., metric, value) to wide (one column per
#' metric) and flag |r| > CORRELATION_FLAG_THRESHOLD pairs. Does NOT auto-drop anything -- a
#' human picks which metric to keep per "tells a similar ecological story" (source doc's own
#' criterion, not purely mechanical).
screen_metric_correlation <- function(metrics_long, id_cols) {
  wide <- tidyr::pivot_wider(
    metrics_long[, c(id_cols, "metric", "value")],
    names_from = metric, values_from = value
  )
  numeric_cols <- wide[, setdiff(names(wide), id_cols), drop = FALSE]
  numeric_cols <- numeric_cols[, vapply(numeric_cols, is.numeric, logical(1)), drop = FALSE]
  corr <- suppressWarnings(stats::cor(numeric_cols, use = "pairwise.complete.obs"))
  corr_long <- as.data.frame(as.table(corr))
  names(corr_long) <- c("metric_a", "metric_b", "r")
  corr_long <- corr_long[corr_long$metric_a != corr_long$metric_b, ]
  corr_long$abs_r <- abs(corr_long$r)
  flagged <- corr_long[!is.na(corr_long$abs_r) & corr_long$abs_r > CORRELATION_FLAG_THRESHOLD, ]
  flagged <- flagged[order(-flagged$abs_r), ]
  list(correlation_matrix = corr, flagged_pairs = flagged)
}
