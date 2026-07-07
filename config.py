from pathlib import Path

#################### FILE PATH HANDLING #######################

REPO_ROOT = Path(__file__).resolve().parents[0]
DATA_DIR = REPO_ROOT / "data"
OUTPUTS_DIR = REPO_ROOT / "outputs"
RASTER_DIR = OUTPUTS_DIR / "rasters"
PLOTS_DIR = OUTPUTS_DIR / "plots"
TABLES_DIR = OUTPUTS_DIR / "tables"
LANDSCAPE_RASTER_DIR = RASTER_DIR / "landscape_metrics"

# Output directory names
OUTPUT_DIRS = {
    "plots": PLOTS_DIR,
    "landscape_metrics": LANDSCAPE_RASTER_DIR,
    "tables": TABLES_DIR,
    "rasters": RASTER_DIR,
}

#################### AOI BOUNDARIES ##########################
AOI_PATHS = {
    "enarau": DATA_DIR / "enarau_conservancy.geojson",
    "mbokishi": DATA_DIR / "mbokishi_conservancy.geojson",
    "phase_1_corridor": DATA_DIR / "phase_1_corridor.geojson",
    "phase_2_corridor": DATA_DIR / "phase_2_corridor.geojson",
}

#################### SITE METADATA ##########################
SITES = [
    {
        "site_id": "enarau",
        "site_name": "Enarau Conservancy",
        "path": AOI_PATHS["enarau"],
    },
    {
        "site_id": "mbokishi",
        "site_name": "Mbokishi Conservancy",
        "path": AOI_PATHS["mbokishi"],
    },
    {
        "site_id": "corridor_p1",
        "site_name": "Corridor Phase 1",
        "path": AOI_PATHS["phase_1_corridor"],
    },
    {
        "site_id": "corridor_p2",
        "site_name": "Corridor Phase 2",
        "path": AOI_PATHS["phase_2_corridor"],
    },
]
STUDY_AREA_BUFFER_M = (
    150  # buffer around the union of all sites for project-wide exports
)

#################### DYNAMIC WORLD HISTORICAL CHANGE (Objective 1) ##########################
DW_CRS = "EPSG:32736"  # WGS 84 / UTM zone 36S -- the plan's EPSG:32637 (UTM 37N) is the wrong hemisphere
DW_EXPORT_FOLDER = "CERK_Enarau_DW_HistoricalChange"
DW_EXCLUDED_PROB_BANDS = [
    "snow_and_ice"
]  # filtered out of eetools.constants.DW_PROBABILITY_BANDS

DW_DERIVED_BANDS = [
    "natural_prob",
    "woody_prob",
    "grass_prob",
    "conversion_pressure_prob",
    "hard_conversion_prob",
    "bare_degradation_prob",
    "water_wetland_prob",
    "top1_prob",
    "valid_obs_count",
]

DW_YEAR_START = 2016
DW_YEAR_END = 2025
DW_SEASONS = ["wet", "dry", "annual"]
DW_PERIODS = {
    "baseline": (2016, 2018),
    "pre": (2019, 2021),
    "current": (2022, 2025),
}
# Recalibrated 2026-07-06 (round 1: was 6/3; round 2, same day: 4/2 -> 3/1) -- visual QA after
# round 1 still showed a moderate amount of masked (no-class) pixels. valid_obs_count is a
# per-pixel count of cloud-free Sentinel-2 acquisitions in the window; the wet season in
# particular (Mar-May, the rainiest months) frequently doesn't clear more than 1-2 cloud-free
# passes at this AOI. The plan itself (objective-1 doc, §6) anticipates this: "Adjust these
# thresholds if cloud masking or image availability is too restrictive." At min_obs=1 the
# "composite" is effectively whatever single image was available -- still transparently a
# median, just over a window of 1, and still not a final calibration; revisit against
# high-resolution imagery per the plan's open questions.
DW_MIN_OBS_ANNUAL = 3  # was 6 -> 4 -> 3
DW_MIN_OBS_SEASONAL = 1  # was 3 -> 2 -> 1
DW_COVERAGE_WARNING_PCT = 60  # was 80 -> 70 -> 60 -- QA/reporting threshold only (see coverage_flag); does not mask any pixels

# Habitat classification thresholds -- round 2 recalibration (2026-07-06), following a code fix
# to classify_habitat's rule precedence (see that function's docstring in the notebook): rules
# were being applied crops-first/natural-last, and since each `.where()` overwrites matches, the
# LAST rule silently won -- so cropland pixels with moderate natural_prob were being overwritten
# to "mixed natural"/grassland/uncertain regardless of threshold values. That is now fixed by
# reordering (crops has highest effective precedence) and by dropping a redundant
# dominance-margin check on the "mixed natural" catch-all that left skewed-but-real vegetation
# signal with no valid class. The threshold nudges below are a secondary, smaller pass on top of
# that fix, for the residual moderate masking/uncertainty visual QA still showed afterward.
# Original (2026-07 plan) values noted per key; still starting values pending visual calibration
# against high-resolution imagery; do not treat as final.
DW_HABITAT_THRESHOLDS = {
    "crops_min": 0.28,  # was 0.45 -> 0.40 -> 0.35 -> 0.25 -- crops is a single raw DW band,
    # while natural_prob/woody_prob are SUMS of 2-3 bands; since all 9 class probabilities sum
    # to 1 per pixel, an aggregated band starts from a structurally higher ceiling than any
    # single band, so a similar absolute floor systematically favors natural/woody/grass over
    # crops even when crops is the clear plurality (e.g. crops=0.28 with the rest thinly spread
    # across trees/shrub/grass/bare). RISK: crops has the highest effective precedence (see
    # classify_habitat), so this floor alone decides cropland with no competing-signal check --
    # watch for false-positive cropland at Mbokishi (the reference/intact-habitat site) after
    # this change; walk back toward ~0.28-0.30 if it appears there.
    "built_min": 0.38,  # was 0.35 -> 0.30 -- built-up pixels are often small/mixed, per the plan's own note
    "water_wetland_min": 0.35,  # was 0.45 -> 0.40
    "bare_min": 0.35,  # was 0.40 -> 0.35
    "woody_min": 0.40,  # was 0.45 -> 0.35 -> 0.30
    "grass_min": 0.35,  # was 0.45 -> 0.35 -> 0.30
    "woody_grass_margin": 0.06,  # was 0.10 -> 0.08 -> 0.06 -- smaller dominance gap needed to call woody vs. grass
    "natural_min": 0.30,  # was 0.55 -> 0.40 -> 0.35 -> 0.30 -- the "mixed natural habitat"
    # catch-all; a true fallback (no margin restriction, see notebook classify_habitat) so this
    # is the main lever for how much moderate-confidence vegetation signal counts as classifiable
    "top1_min": 0.25,  # was 0.35 -> 0.30 -> 0.25 -- DW confidence proxy; sand/bare/arid surfaces
    # are known to depress top1_prob by scoring across multiple classes at once (see
    # wiki/tools/dynamic-world.md)
}
DW_HABITAT_CLASS_CODES = list(
    range(1, 9)
)  # classify_habitat never emits 0 (NoData/outside AOI)

# Conversion-pressure thresholds
DW_PRESSURE_THRESHOLDS = {"moderate_min": 0.35, "high_min": 0.55}
DW_PRESSURE_CLASS_CODES = [0, 1, 2]

# All from-class x to-class combinations of the 8 habitat classes.
DW_TRANSITION_CODES = [
    f * 10 + t for f in DW_HABITAT_CLASS_CODES for t in DW_HABITAT_CLASS_CODES
]

# Objective 4 connectivity-mask thresholds (connectivity-input handoff amendment)
DW_CONNECTIVITY_THRESHOLDS = {
    "natural_prob_min": 0.60,
    "conversion_pressure_max": 0.30,
    "top1_prob_min": 0.35,
}

# Visualization presets
DW_HABITAT_CLASS_VIS = {
    "min": 0,
    "max": 8,
    "palette": [
        "#cccccc",  # 0 nodata
        "#397d49",  # 1 woody
        "#88b053",  # 2 grassland
        "#a3c586",  # 3 mixed natural
        "#e49635",  # 4 cropland
        "#c4281b",  # 5 built
        "#a59b8f",  # 6 bare/degraded
        "#419bdf",  # 7 water/flooded vegetation
        "#b39fe1",  # 8 uncertain
    ],
}
DW_PRESSURE_VIS = {
    "min": 0,
    "max": 2,
    "palette": ["#1a9850", "#fee08b", "#d73027"],
}
# transition_code = from_class * 10 + to_class, where both from_class and to_class are
# DW_HABITAT_CLASS_CODES (1-8, see DW_HABITAT_CLASS_VIS above for what each code means). "min"
# and "max" are the lowest/highest codes that can occur (11 = woody->woody i.e. persistent
# woody; 88 = uncertain->uncertain i.e. persistent uncertain) -- NOT a severity scale. EE/geemap
# linearly interpolates the 6 palette colors across that 11-88 numeric range, so color is driven
# mainly by the FROM class (the tens digit) and only subtly by the TO class (the ones digit,
# worth only ~1/77th of the range) -- the ramp is a visual differentiator, not a "loss=red,
# gain=green" encoding. Approximate stop values under linear interpolation (step = 77/5 = 15.4):
DW_TRANSITION_VIS = {
    "min": 11,
    "max": 88,
    "palette": [
        "#800026",  # ~11.0 -- from_class 1 (woody) origin, e.g. persistent woody (11)
        "#f03b20",  # ~26.4 -- from_class 2 (grassland) origin, e.g. grass->cropland (24)
        "#fd8d3c",  # ~41.8 -- from_class 4 (cropland) origin, e.g. cropland->woody recovery (41)
        "#ffffb2",  # ~57.2 -- from_class 5-6 (built/bare) origin
        "#c7e9b4",  # ~72.6 -- from_class 7 (water/flooded) origin
        "#41b6c4",  # ~88.0 -- from_class 8 (uncertain) origin, e.g. persistent uncertain (88)
    ],
}
