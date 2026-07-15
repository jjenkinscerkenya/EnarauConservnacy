# Patch delineation, patch-level metrics, and the igraph Euclidean patch-graph connectivity
# screen. Runs on the full project-extent binary natural-habitat raster (see masking-order rule
# in 00_config.R/plan) -- a corridor patch spanning two sites must not be truncated.

#' Delineate patches from a binary natural-habitat raster (1 = natural, 0 = converted, NA =
#' excluded), filter to >= MIN_PATCH_AREA_HA, and return both the cleaned patch-ID raster and its
#' dissolved polygons.
#'
#' zeroAsNA = TRUE is required here: without it, terra::patches() also delineates "patches" out
#' of the 0 (converted/non-natural) background, which would corrupt both the area filter and the
#' patch_id numbering. With zeroAsNA = TRUE, only class-1 (natural) cells receive patch IDs --
#' the same connected-component labeling landscapemetrics::calculate_lsm(level="patch") uses
#' internally for this class, on the same directions=8 setting, so the resulting patch_id values
#' are expected to align with calculate_patch_metrics()'s own `id` column. This is verified in
#' 05_patch_importance_graph.R (patch count cross-check) rather than assumed silently.
delineate_patches <- function(r_bin_natural) {
  patch_id <- terra::patches(r_bin_natural, directions = 8, zeroAsNA = TRUE)
  pixel_area_ha <- prod(terra::res(patch_id)) / 10000
  patch_freq <- as.data.frame(terra::freq(patch_id, useNA = "no"))
  patch_freq$area_ha <- patch_freq$count * pixel_area_ha
  keep_ids <- patch_freq$value[patch_freq$area_ha >= MIN_PATCH_AREA_HA]

  patch_id_clean <- terra::ifel(patch_id %in% keep_ids, patch_id, NA)
  patch_poly <- terra::as.polygons(patch_id_clean, dissolve = TRUE, na.rm = TRUE)
  names(patch_poly) <- "patch_id"

  list(patch_id_raster = patch_id_clean, patch_polygons = patch_poly, area_table = patch_freq)
}

#' Patch-level landscapemetrics (AREA/CORE/SHAPE/ENN) on the binary raster.
calculate_patch_metrics <- function(r_bin_natural) {
  landscapemetrics::calculate_lsm(
    landscape = r_bin_natural,
    what = c("lsm_p_area", "lsm_p_core", "lsm_p_shape", "lsm_p_enn"),
    directions = 8,
    edge_depth = EDGE_DEPTH_CELLS
  )
}

#' Assign each patch polygon to its primary site by largest-area overlap; flag patches spanning
#' more than one site (expected/meaningful for corridor-crossing patches, not an error).
attribute_patches_to_sites <- function(patch_poly, sites_sf) {
  patch_sf <- sf::st_as_sf(patch_poly)
  inter <- suppressWarnings(sf::st_intersection(patch_sf, sites_sf))
  inter$overlap_area_m2 <- as.numeric(sf::st_area(inter))

  by_patch <- split(inter, inter$patch_id)
  summary_rows <- lapply(by_patch, function(rows) {
    primary <- rows[which.max(rows$overlap_area_m2), ]
    data.frame(
      patch_id = primary$patch_id[1],
      primary_site_id = primary$site_id[1],
      spans_multiple_sites = length(unique(rows$site_id)) > 1,
      site_ids = paste(sort(unique(rows$site_id)), collapse = ";")
    )
  })
  do.call(rbind, summary_rows)
}

#' Mean conversion pressure within a buffer around each patch polygon.
mean_pressure_around_patches <- function(patch_poly, pressure_raster, buffer_m = PATCH_PRESSURE_BUFFER_M) {
  patch_sf <- sf::st_as_sf(patch_poly)
  buffered <- sf::st_buffer(patch_sf, dist = buffer_m)
  vals <- terra::extract(pressure_raster, terra::vect(buffered), fun = mean, na.rm = TRUE, ID = FALSE)
  data.frame(patch_id = patch_sf$patch_id, mean_pressure = vals[[1]])
}

#' Build an igraph patch graph at one distance threshold. Edges are polygon-to-polygon Euclidean
#' distance (sf::st_distance), not centroid distance -- a large, irregular patch is closer to its
#' neighbor at its nearest edge than its centroid-to-centroid distance would suggest.
build_patch_graph <- function(patch_poly, distance_threshold_m) {
  patch_sf <- sf::st_as_sf(patch_poly)
  n <- nrow(patch_sf)
  dist_mat <- sf::st_distance(patch_sf)
  dist_mat <- units::drop_units(dist_mat)

  edges <- which(dist_mat <= distance_threshold_m & upper.tri(dist_mat), arr.ind = TRUE)
  edge_df <- data.frame(
    from = patch_sf$patch_id[edges[, 1]],
    to   = patch_sf$patch_id[edges[, 2]],
    distance_m = dist_mat[edges]
  )

  g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = data.frame(name = patch_sf$patch_id))

  data.frame(
    patch_id = patch_sf$patch_id,
    degree = igraph::degree(g),
    betweenness = igraph::betweenness(g),
    component = igraph::components(g)$membership,
    is_articulation_point = patch_sf$patch_id %in% names(igraph::articulation_points(g))
  ) -> node_metrics

  list(graph = g, node_metrics = node_metrics, edges = edge_df)
}
