# Mirrors C:/Users/harre/Repos/CERK/EnarauConservnacy/config.py -- config.py is the Python
# source of truth for this repo. R cannot `import config.py` (no reticulate in use anywhere in
# this repo, by design -- see CLAUDE.md's "mixed-language stack, not pure Python"). If config.py's
# PROJECT_CRS, DW_PERIODS, DW_HABITAT_CLASS_LABELS, SITES, or STUDY_AREA_BUFFER_M change, this
# file must be updated to match by hand -- there is no automated drift check.
#
# renv note: this project's renv (scripts/r/renv.lock) only activates when the R process's
# working directory is scripts/r/ itself. Always run scripts as:
#   cd scripts/r && Rscript 01_prepare_inputs.R
# NOT `Rscript scripts/r/01_prepare_inputs.R` from the repo root -- the latter leaves the working
# directory at the repo root (no .Rprofile there), so renv never activates and whatever is on the
# machine's global R library silently runs instead.

find_repo_root <- function() {
  candidates <- c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "..", ".."))
  for (d in candidates) {
    if (file.exists(file.path(d, "config.py"))) return(normalizePath(d))
  }
  stop(
    "Could not locate EnarauConservancy repo root (config.py not found) from ", getwd(),
    ". Run this script with working directory = scripts/r/ (see renv note above)."
  )
}

REPO_ROOT <- find_repo_root()

#################### FILE PATH HANDLING ####################
DATA_DIR             <- file.path(REPO_ROOT, "data")
OUTPUTS_DIR           <- file.path(REPO_ROOT, "outputs")
RASTER_DIR            <- file.path(OUTPUTS_DIR, "rasters")
PLOTS_DIR             <- file.path(OUTPUTS_DIR, "plots")
TABLES_DIR            <- file.path(OUTPUTS_DIR, "tables")
LANDSCAPE_RASTER_DIR  <- file.path(RASTER_DIR, "landscape_metrics")  # Objective 3's OWN raster outputs
DW_INPUT_RASTER_DIR   <- file.path(RASTER_DIR, "dynamic_world")      # manually-downloaded Objective 1 inputs
VECTORS_DIR           <- file.path(OUTPUTS_DIR, "vectors")

for (d in c(LANDSCAPE_RASTER_DIR, VECTORS_DIR, TABLES_DIR, PLOTS_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

#################### PROJECT-WIDE SETTINGS ####################
# Confirmed 2026-07-06, WGS 84 / UTM zone 36S -- same correction as config.py's own comment
# (the Objective 1 plan document's EPSG:32637 is the wrong hemisphere for this AOI).
PROJECT_CRS <- "EPSG:32736"
STUDY_AREA_BUFFER_M <- 150  # documentation parity only; already baked into downloaded raster extents

#################### AOI BOUNDARIES / SITE METADATA ####################
AOI_PATHS <- list(
  enarau      = file.path(DATA_DIR, "enarau_conservancy.geojson"),
  mbokishi    = file.path(DATA_DIR, "mbokishi_conservancy.geojson"),
  corridor_p1 = file.path(DATA_DIR, "phase_1_corridor.geojson"),
  corridor_p2 = file.path(DATA_DIR, "phase_2_corridor.geojson")
)

SITES <- data.frame(
  site_id   = c("enarau", "mbokishi", "corridor_p1", "corridor_p2"),
  site_name = c("Enarau Conservancy", "Mbokishi Conservancy", "Corridor Phase 1", "Corridor Phase 2"),
  path      = unlist(AOI_PATHS, use.names = FALSE),
  stringsAsFactors = FALSE
)

#################### DYNAMIC WORLD CLASS SCHEME (Objective 1) ####################
# There is no class 0 -- classify_habitat() never emits it. NoData is the raster's own -9999
# sentinel (see NODATA_SENTINEL below), not a literal class value.
DW_HABITAT_CLASS_LABELS <- c(
  `1` = "Woody", `2` = "Grassland", `3` = "Mixed natural", `4` = "Cropland",
  `5` = "Built", `6` = "Bare/degraded", `7` = "Water/flooded veg", `8` = "Uncertain"
)
NATURAL_CLASSES    <- c(1, 2, 3)
CONVERSION_CLASSES <- c(4, 5, 6)
EXCLUDED_CLASSES   <- c(7, 8)  # + raster NA

DW_PRESSURE_CLASS_LABELS <- c(`0` = "Low", `1` = "Moderate", `2` = "High")

#################### PERIODS ####################
# Matches config.py's DW_PERIODS (Objective 1); Objective 2's own PERIODS uses different year
# ranges (Landsat baseline 1984-2000) and is NOT the same dict -- don't conflate the two.
DW_PERIODS <- list(
  baseline = c(2016, 2018),
  pre      = c(2019, 2021),
  current  = c(2022, 2025)
)
# Literal filename tokens used by Objective 1's raster exports (see scripts/r/R/io.R)
DW_PERIOD_TOKENS <- c(
  "baseline_2016_2018", "pre_2019_2021", "current_2022_2025",
  "current_wet_2022_2025", "current_dry_2022_2025"
)

#################### RASTER / QA CONSTANTS ####################
NODATA_SENTINEL <- -9999  # eetools.io.export_image_to_drive's unmask() sentinel; verified in 01_prepare_inputs.R
                          # that terra reads this as NA on load, not literal data.
VALID_PIXEL_COVERAGE_MIN <- 0.80  # TGBS_Kwale precedent: site-years below this show artificial
                                  # patch breaks that inflate NP/PD/ED -- exclude from Level 2/3.

#################### LANDSCAPE METRICS PARAMETERS ####################
# Package default (1 cell = 10 m at this resolution). Provisional, same treatment as config.py's
# DW_HABITAT_THRESHOLDS -- not yet literature-justified, must stay fixed across the whole project
# once chosen so CORE/patch-importance numbers remain comparable across periods and sites.
EDGE_DEPTH_CELLS <- 1

MOVING_WINDOW_RADII_M <- c(250, 500, 1000)
PATCH_GRAPH_DISTANCE_THRESHOLDS_M <- c(250, 500, 1000)
MIN_PATCH_AREA_HA <- 1
PATCH_PRESSURE_BUFFER_M <- 100        # "mean conversion pressure around patch" buffer distance
CORRIDOR_PROXIMITY_DECAY_M <- 2000    # linear decay-to-zero distance from corridor_p1 U corridor_p2

#################### CORRELATION / METRIC SETS ####################
SELECTED_CLASS_METRICS <- c(
  "lsm_c_ca", "lsm_c_pland", "lsm_c_pd", "lsm_c_lpi", "lsm_c_ed",
  "lsm_c_ai", "lsm_c_clumpy", "lsm_c_cohesion", "lsm_c_enn_mn", "lsm_c_mesh"
)
SELECTED_BINARY_METRICS <- setdiff(SELECTED_CLASS_METRICS, "lsm_c_ca")  # CA not meaningful for a 1-class binary landscape
CORRELATION_FLAG_THRESHOLD <- 0.85

#################### LINKAGE PRIORITY SCORE WEIGHTS (06) ####################
LINKAGE_SCORE_WEIGHTS <- c(
  persistent_or_recovered_natural = 0.25,
  local_connectivity              = 0.20,
  patch_importance                = 0.20,
  low_conversion_pressure         = 0.15,
  bottleneck_relevance             = 0.10,
  corridor_proximity              = 0.10
)
PATCH_IMPORTANCE_WEIGHTS <- c(
  area        = 0.30,
  core        = 0.20,
  isolation   = 0.20,  # applied as inverse-normalized ENN
  betweenness = 0.15,
  pressure    = 0.15   # applied as inverse-normalized pressure
)
