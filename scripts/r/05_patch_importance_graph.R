# Objective 3, step 5: patch delineation, patch metrics, and the igraph patch-graph
# connectivity screen (Level 3b/c).
#
# RUN AS: cd scripts/r && Rscript 05_patch_importance_graph.R
#
# PROJECT-WIDE FIRST (same masking-order rule as 04) -- a corridor patch spanning Enarau/Corridor
# Phase 1/2/Mbokishi must not be truncated at a site boundary before delineation. Site attribution
# happens AFTER delineation, by largest-area overlap.
#
# Does not depend on 04 -- patch_importance_score uses area/core/isolation/betweenness/pressure
# only, no moving-window input.

source("00_config.R")
source("R/io.R")
source("R/recode.R")
source("R/patch_graph.R")
source("R/scoring.R")

message("=== 05_patch_importance_graph ===")

current_class_path <- period_manifest$class_file[period_manifest$token == "current_2022_2025"]
current_pressure_path <- period_manifest$pressure_file[period_manifest$token == "current_2022_2025"]

if (!file.exists(current_class_path) || !file.exists(current_pressure_path)) {
  stop("Current-period class/pressure rasters not yet downloaded -- run 01_prepare_inputs.R first and check the manifest.")
}

current_r <- read_habitat_raster(current_class_path)
current_bin <- make_natural_binary(current_r)
pressure_r <- read_pressure_raster(current_pressure_path)

message("Delineating patches (>= ", MIN_PATCH_AREA_HA, " ha, 8-connectivity)...")
patches_result <- delineate_patches(current_bin)
patch_poly <- patches_result$patch_polygons
message(nrow(patch_poly), " patches retained.")

message("Computing patch-level metrics (edge_depth = ", EDGE_DEPTH_CELLS, " cell)...")
patch_metrics <- calculate_patch_metrics(current_bin)
readr::write_csv(patch_metrics, file.path(TABLES_DIR, "landscape_patch_metrics_current.csv"))

# Cross-check: terra::patches() (used for patch_poly above) and calculate_lsm(level="patch")
# (used for patch_metrics) are independent calls that must agree on patch identity -- see the
# comment in R/patch_graph.R's delineate_patches(). Verify rather than assume.
n_poly_patches <- nrow(patch_poly)
n_metric_patches <- length(unique(patch_metrics$id[patch_metrics$class == 1]))
if (n_poly_patches != n_metric_patches) {
  warning(sprintf(
    "Patch count mismatch: %d polygons from terra::patches() vs %d patch IDs from calculate_lsm() -- patch_id join below may be unreliable.",
    n_poly_patches, n_metric_patches
  ))
} else {
  message("Patch count cross-check OK: ", n_poly_patches, " patches agree between terra::patches() and calculate_lsm().")
}

message("Attributing patches to primary sites...")
sites_sf <- do.call(rbind, lapply(SITES$site_id, function(sid) {
  b <- read_site_boundary(sid)
  sf::st_sf(site_id = sid, geometry = sf::st_geometry(sf::st_union(b)))
}))
site_attribution <- attribute_patches_to_sites(patch_poly, sites_sf)

message("Computing mean conversion pressure within ", PATCH_PRESSURE_BUFFER_M, "m of each patch...")
patch_pressure <- mean_pressure_around_patches(patch_poly, pressure_r, buffer_m = PATCH_PRESSURE_BUFFER_M)

# ---- Patch graph at each distance threshold ----
graph_results <- lapply(PATCH_GRAPH_DISTANCE_THRESHOLDS_M, function(thresh) {
  message("Building patch graph at ", thresh, "m distance threshold...")
  build_patch_graph(patch_poly, thresh)
})
names(graph_results) <- paste0("w", PATCH_GRAPH_DISTANCE_THRESHOLDS_M, "m")

graph_metrics_all <- do.call(rbind, lapply(names(graph_results), function(nm) {
  cbind(threshold = nm, graph_results[[nm]]$node_metrics)
}))
readr::write_csv(graph_metrics_all, file.path(TABLES_DIR, "landscape_patch_graph_metrics_current.csv"))

# ---- Patch importance score (uses the 500m threshold's betweenness as the "stepping-stone role"
# term -- the source doc's own recommended default starting radius for local patch structure) ----
mid_threshold_name <- paste0("w", MOVING_WINDOW_RADII_M[MOVING_WINDOW_RADII_M == 500], "m")
mid_graph <- graph_results[[mid_threshold_name]]$node_metrics

# patch_metrics already has all four (area/core/shape/enn) pivoted long by `metric`; reshape once.
patch_metrics_wide <- tidyr::pivot_wider(
  patch_metrics[patch_metrics$class == 1, c("id", "metric", "value")],
  names_from = metric, values_from = value
)
names(patch_metrics_wide)[names(patch_metrics_wide) == "id"] <- "patch_id"

importance_input <- patch_metrics_wide |>
  dplyr::left_join(mid_graph[, c("patch_id", "betweenness")], by = "patch_id") |>
  dplyr::left_join(patch_pressure, by = "patch_id") |>
  dplyr::left_join(site_attribution, by = "patch_id")

importance_input$patch_importance_score <- compute_patch_importance_score(
  area_ha = importance_input$area,
  core_area_ha = importance_input$core,
  enn_m = importance_input$enn,
  betweenness = importance_input$betweenness,
  mean_pressure = importance_input$mean_pressure
)

readr::write_csv(importance_input, file.path(TABLES_DIR, "landscape_patch_importance_scores_current.csv"))

# ---- Vector outputs ----
patch_poly_sf <- sf::st_as_sf(patch_poly) |>
  dplyr::left_join(importance_input, by = "patch_id")

sf::st_write(patch_poly_sf, file.path(VECTORS_DIR, "natural_patches_current.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# Multi-layer patch graph gpkg: one "nodes" layer + one "edges_<threshold>" layer per threshold.
nodes_sf <- sf::st_as_sf(patch_poly) |> dplyr::left_join(importance_input, by = "patch_id")
sf::st_write(nodes_sf, file.path(VECTORS_DIR, "patch_graph_current.gpkg"), layer = "nodes", delete_dsn = TRUE, quiet = TRUE)
for (nm in names(graph_results)) {
  edges <- graph_results[[nm]]$edges
  if (nrow(edges) > 0) {
    patch_centroids <- sf::st_centroid(sf::st_as_sf(patch_poly))
    edge_lines <- lapply(seq_len(nrow(edges)), function(i) {
      from_geom <- patch_centroids[patch_centroids$patch_id == edges$from[i], ]
      to_geom <- patch_centroids[patch_centroids$patch_id == edges$to[i], ]
      sf::st_sfc(sf::st_linestring(rbind(sf::st_coordinates(from_geom), sf::st_coordinates(to_geom))))
    })
    edges_sf <- sf::st_sf(edges, geometry = do.call(c, edge_lines), crs = sf::st_crs(patch_poly))
    sf::st_write(edges_sf, file.path(VECTORS_DIR, "patch_graph_current.gpkg"), layer = paste0("edges_", nm), append = TRUE, quiet = TRUE)
  }
}

# High-priority stepping-stone patches: top-quartile patch_importance_score AND degree >= 2 at
# the mid threshold (i.e. genuinely a stepping stone between >=2 other patches, not just large).
importance_threshold <- stats::quantile(importance_input$patch_importance_score, 0.75, na.rm = TRUE)
stepping_stones <- patch_poly_sf[
  patch_poly_sf$patch_importance_score >= importance_threshold & patch_poly_sf$betweenness > 0,
]
sf::st_write(stepping_stones, file.path(VECTORS_DIR, "high_priority_stepping_stone_patches.gpkg"), delete_dsn = TRUE, quiet = TRUE)

message("=== 05_patch_importance_graph complete ===")
message(nrow(stepping_stones), " high-priority stepping-stone patches identified.")
