# Objective 3, step 2: area/transition summaries.
#
# RUN AS: cd scripts/r && Rscript 02_area_transition_summaries.R
#
# Reads EXISTING zonal-statistics CSVs (already computed server-side in GEE,
# independent of which rasters were exported) rather than recomputing area-by-class from raw
# rasters -- these CSVs already cover annual/seasonal AND period-level area-by-class, plus
# transition areas. No raster reads in this script.

source("00_config.R")

message("=== 02_area_transition_summaries ===")

area_year_season <- readr::read_csv(file.path(TABLES_DIR, "dw_area_by_class_by_site_year_season.csv"), show_col_types = FALSE)
area_period      <- readr::read_csv(file.path(TABLES_DIR, "dw_area_by_class_by_site_period.csv"), show_col_types = FALSE)
transition_area  <- readr::read_csv(file.path(TABLES_DIR, "dw_transition_area_by_site.csv"), show_col_types = FALSE)

habitat_cols <- paste0("habitat_area_ha_", 1:8)
pressure_cols <- paste0("pressure_area_ha_", 0:2)

# ---- Tidy (long) area-by-class tables ----
tidy_area <- function(df, id_cols) {
  df |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(habitat_cols),
      names_to = "habitat_class", names_prefix = "habitat_area_ha_",
      values_to = "area_ha"
    ) |>
    dplyr::mutate(
      habitat_class = as.integer(habitat_class),
      habitat_label = DW_HABITAT_CLASS_LABELS[as.character(habitat_class)],
      habitat_group = dplyr::case_when(
        habitat_class %in% NATURAL_CLASSES ~ "natural",
        habitat_class %in% CONVERSION_CLASSES ~ "conversion",
        TRUE ~ "excluded"
      )
    ) |>
    dplyr::select(dplyr::all_of(c(id_cols, "site_id", "site_name")), habitat_class, habitat_label, habitat_group, area_ha)
}

area_year_season_tidy <- tidy_area(area_year_season, id_cols = c("year", "season"))
area_period_tidy <- tidy_area(area_period, id_cols = "period")

readr::write_csv(area_year_season_tidy, file.path(TABLES_DIR, "landscape_area_by_class_tidy_site_year_season.csv"))
readr::write_csv(area_period_tidy, file.path(TABLES_DIR, "landscape_area_by_class_tidy_site_period.csv"))

# ---- Net natural habitat change per site (current - baseline, current - pre) ----
natural_area_by_period <- area_period_tidy |>
  dplyr::filter(habitat_group == "natural") |>
  dplyr::group_by(site_id, site_name, period) |>
  dplyr::summarise(natural_area_ha = sum(area_ha), .groups = "drop")

# site_area_ha needs the actual reprojected boundary areas (source geojsons are WGS84 lon/lat,
# so this must reproject, not read raw coordinate extents).
site_areas_ha <- vapply(SITES$site_id, function(sid) {
  as.numeric(sum(sf::st_area(sf::st_transform(sf::st_read(SITES$path[SITES$site_id == sid], quiet = TRUE), PROJECT_CRS)))) / 10000
}, numeric(1))

net_natural_change <- natural_area_by_period |>
  tidyr::pivot_wider(names_from = period, values_from = natural_area_ha, names_prefix = "natural_ha_") |>
  dplyr::mutate(
    site_area_ha = site_areas_ha[site_id],
    net_change_baseline_to_current_ha = natural_ha_current - natural_ha_baseline,
    net_change_pre_to_current_ha = natural_ha_current - natural_ha_pre,
    net_change_baseline_to_current_pct_of_site = 100 * net_change_baseline_to_current_ha / site_area_ha,
    net_change_pre_to_current_pct_of_site = 100 * net_change_pre_to_current_ha / site_area_ha
  )

readr::write_csv(net_natural_change, file.path(TABLES_DIR, "landscape_net_natural_change_by_site.csv"))

# ---- Transition grouping (source doc Sec.9.2) ----
transition_groups <- list(
  stable_natural       = c(11, 12, 13, 21, 22, 23, 31, 32, 33),
  natural_to_cropland  = c(14, 24, 34),
  natural_to_built     = c(15, 25, 35),
  natural_to_bare      = c(16, 26, 36),
  pressure_to_natural  = c(41, 42, 43, 51, 52, 53, 61, 62, 63),
  stable_pressure      = c(44, 45, 46, 54, 55, 56, 64, 65, 66)
)
# Structural-change sub-flags, not mutually exclusive with stable_natural above  
# used for woodland/grassland shifts.
structural_change_groups <- list(
  woody_to_grass_or_mixed = c(12, 13),
  grass_to_woody_or_mixed = c(21, 23)
)

transition_cols <- grep("^transition_area_ha_", names(transition_area), value = TRUE)

transition_long <- transition_area |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(transition_cols),
    names_to = "transition_code", names_prefix = "transition_area_ha_",
    values_to = "area_ha"
  ) |>
  dplyr::mutate(transition_code = as.integer(transition_code))

group_lookup <- utils::stack(transition_groups) |>
  dplyr::rename(transition_code = values, category = ind) |>
  dplyr::mutate(transition_code = as.integer(as.character(transition_code)))

structural_lookup <- utils::stack(structural_change_groups) |>
  dplyr::rename(transition_code = values, structural_category = ind) |>
  dplyr::mutate(transition_code = as.integer(as.character(transition_code)))

transition_summary <- transition_long |>
  dplyr::left_join(group_lookup, by = "transition_code") |>
  dplyr::left_join(structural_lookup, by = "transition_code") |>
  dplyr::mutate(category = dplyr::coalesce(as.character(category), "other")) |>
  dplyr::group_by(site_id, site_name, comparison, category) |>
  dplyr::summarise(area_ha = sum(area_ha), .groups = "drop") |>
  dplyr::mutate(
    site_area_ha = site_areas_ha[site_id],
    pct_of_site_area = 100 * area_ha / site_area_ha
  )

readr::write_csv(transition_summary, file.path(TABLES_DIR, "landscape_transition_summary_by_site_period.csv"))

message("=== 02_area_transition_summaries complete ===")
message("Wrote: landscape_area_by_class_tidy_site_year_season.csv, landscape_area_by_class_tidy_site_period.csv, ",
        "landscape_net_natural_change_by_site.csv, landscape_transition_summary_by_site_period.csv")
