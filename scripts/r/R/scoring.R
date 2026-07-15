# Patch importance score (05) and linkage priority score (06). Both use GLOBAL normalization
# across all 4 sites' patches/pixels pooled together, not per-site -- appropriate since the goal
# is cross-site corridor prioritization (which patches matter most for the whole network), not a
# within-site ranking. Per-site normalization would make each site's single largest patch always
# score 1.0 locally regardless of true landscape-scale importance.

#' Min-max normalize to [0, 1]. NA-safe (ignores NA in range calc; propagates NA through).
normalize01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0 || !is.finite(diff(rng))) {
    return(rep(if (length(x) > 0) 0.5 else numeric(0), length(x)))
  }
  (x - rng[1]) / diff(rng)
}

#' Combine patch-level inputs into the weighted patch_importance_score.
#' @param area_ha,core_area_ha,enn_m,betweenness,mean_pressure Numeric vectors, one per patch,
#'   all the same length/order.
compute_patch_importance_score <- function(area_ha, core_area_ha, enn_m, betweenness, mean_pressure) {
  betweenness <- ifelse(is.na(betweenness), 0, betweenness)  # isolated/degree-0 patches score 0, not NA
  w <- PATCH_IMPORTANCE_WEIGHTS
  w["area"] * normalize01(area_ha) +
    w["core"] * normalize01(core_area_ha) +
    w["isolation"] * (1 - normalize01(enn_m)) +
    w["betweenness"] * normalize01(betweenness) +
    w["pressure"] * (1 - normalize01(mean_pressure))
}

#' Linear distance-decay score from a reference geometry: 1 inside/at the boundary, decaying
#' linearly to 0 at `decay_m`.
distance_decay_score <- function(r_or_points, reference_sf, decay_m = CORRIDOR_PROXIMITY_DECAY_M) {
  ref_union <- sf::st_union(reference_sf)
  if (inherits(r_or_points, "SpatRaster")) {
    dist_r <- terra::distance(r_or_points, terra::vect(ref_union))
    score <- 1 - terra::clamp(dist_r / decay_m, lower = 0, upper = 1)
    names(score) <- "corridor_proximity_score"
    return(score)
  }
  dist_v <- as.numeric(sf::st_distance(r_or_points, ref_union))
  pmax(0, 1 - dist_v / decay_m)
}

#' Combine the six components into the linkage_priority_score, per LINKAGE_SCORE_WEIGHTS.
#' All six inputs must already be normalized to a comparable [0, 1]-ish scale (rasters aligned
#' to the same grid) before calling this.
compute_linkage_priority_score <- function(persistent_or_recovered_natural, local_connectivity,
                                            patch_importance, low_conversion_pressure,
                                            bottleneck_relevance, corridor_proximity) {
  w <- LINKAGE_SCORE_WEIGHTS
  w["persistent_or_recovered_natural"] * persistent_or_recovered_natural +
    w["local_connectivity"] * local_connectivity +
    w["patch_importance"] * patch_importance +
    w["low_conversion_pressure"] * low_conversion_pressure +
    w["bottleneck_relevance"] * bottleneck_relevance +
    w["corridor_proximity"] * corridor_proximity
}
