# Objective 3, step 1: pre-flight checks + derived-raster prep.
#
# RUN AS: cd scripts/r && Rscript 01_prepare_inputs.R
# (NOT `Rscript scripts/r/01_prepare_inputs.R` from the repo root -- renv only activates when
# the working directory is scripts/r/ itself. See 00_config.R's header comment.)
#
# Requires the Objective 1 rasters listed in the plan to be manually downloaded from Google
# Drive folder CERK_Enarau_DW_HistoricalChange into outputs/rasters/dynamic_world/ first. If only
# a partial download is present, this script still runs (STRICT_INPUT_CHECK = FALSE below) and
# reports what's missing rather than hard-failing -- see Verification step 3 of the plan.

source("00_config.R")
source("R/io.R")
source("R/recode.R")
source("R/qa.R")

STRICT_INPUT_CHECK <- FALSE  # set TRUE once all 60 files are downloaded

message("=== 01_prepare_inputs: checking input manifest ===")
missing_files <- check_inputs_present(strict = STRICT_INPUT_CHECK)

# ---- Site boundaries: reproject WGS84 lon/lat -> PROJECT_CRS ----
message("=== Reprojecting site boundaries to ", PROJECT_CRS, " ===")
site_boundaries <- lapply(SITES$site_id, read_site_boundary)
names(site_boundaries) <- SITES$site_id
site_areas_ha <- vapply(site_boundaries, function(s) as.numeric(sum(sf::st_area(s))) / 10000, numeric(1))
message("Site areas (ha): ", paste(sprintf("%s=%.1f", names(site_areas_ha), site_areas_ha), collapse = ", "))

# ---- NoData / CRS sanity + landscapemetrics::check_landscape() on the current-period class raster ----
current_class_path <- period_manifest$class_file[period_manifest$token == "current_2022_2025"]
if (file.exists(current_class_path)) {
  message("=== NoData/CRS sanity check on current-period habitat_class raster ===")
  current_class_r <- read_habitat_raster(current_class_path)
  check_result <- landscapemetrics::check_landscape(recode_full_habitat(current_class_r))
  print(check_result)
  # check_landscape()'s OK column is a display glyph (character U+2714 heavy check mark),
  # not logical -- compare by codepoint rather than treating it as TRUE/FALSE or risking a
  # source-encoding mismatch on the literal glyph.
  if (!all(vapply(check_result$OK, utf8ToInt, integer(1)) == 10004L)) {
    warning("landscapemetrics::check_landscape() flagged issues -- review before trusting downstream metrics.")
  }

  # ---- Mask-reuse agreement check (current period only) ----
  if (file.exists(MASK_FILES["natural_habitat_mask"])) {
    message("=== Mask-reuse agreement check (current period) ===")
    natural_mask_r <- terra::rast(MASK_FILES["natural_habitat_mask"])
    plain_binary <- make_natural_binary(current_class_r)
    reconciled_binary <- build_current_binary_natural(current_class_r, natural_mask_r)
    check_mask_reuse_agreement(plain_binary, reconciled_binary)
  } else {
    message("Skipping mask-reuse check -- natural_habitat_mask file not yet downloaded.")
  }
} else {
  message("Skipping NoData/CRS sanity + mask-reuse checks -- current-period habitat_class raster not yet downloaded.")
}

# ---- Valid-pixel-coverage QA (TGBS_Kwale 80% rule) ----
seasonal_manifest <- expand_seasonal_manifest()
seasonal_available <- seasonal_manifest[file.exists(seasonal_manifest$habitat_class_file), ]
period_available <- period_manifest[file.exists(period_manifest$class_file), ]

if (nrow(seasonal_available) > 0) {
  message("=== Valid-pixel coverage QA: seasonal per-year rasters (", nrow(seasonal_available), "/", nrow(seasonal_manifest), " available) ===")
  coverage_year_season <- compute_valid_coverage_table(seasonal_available, id_cols = c("year", "season"))
  readr::write_csv(coverage_year_season, file.path(TABLES_DIR, "landscape_valid_pixel_coverage_by_site_year_season.csv"))
  report_excluded_rows(coverage_year_season, "seasonal per-year rasters")
} else {
  message("No seasonal rasters downloaded yet -- skipping seasonal coverage QA.")
}

if (nrow(period_available) > 0) {
  message("=== Valid-pixel coverage QA: period composites (", nrow(period_available), "/", nrow(period_manifest), " available) ===")
  period_available_renamed <- period_available
  names(period_available_renamed)[names(period_available_renamed) == "class_file"] <- "habitat_class_file"
  # id_cols = "period" (not "token") so this table's ID column matches what 03/04/06 look up by --
  # values are still the 5 DW_PERIOD_TOKENS-style tokens (baseline_2016_2018 etc.), a different
  # value domain than Objective 1's own 3-valued "period" column in dw_area_by_class_by_site_period.csv
  # (consumed separately in 02) -- don't join those two tables on this column name without checking.
  names(period_available_renamed)[names(period_available_renamed) == "token"] <- "period"
  coverage_period <- compute_valid_coverage_table(period_available_renamed, id_cols = "period")
  readr::write_csv(coverage_period, file.path(TABLES_DIR, "landscape_valid_pixel_coverage_by_site_period.csv"))
  report_excluded_rows(coverage_period, "period composites")
} else {
  message("No period composite rasters downloaded yet -- skipping period coverage QA.")
}

message("=== 01_prepare_inputs complete ===")
if (length(missing_files) > 0) {
  message(length(missing_files), " input file(s) still missing -- re-run after downloading the rest from Drive.")
}
