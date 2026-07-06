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
    500  # buffer around the union of all sites for project-wide exports
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
DW_MIN_OBS_ANNUAL = 6
DW_MIN_OBS_SEASONAL = 3
DW_COVERAGE_WARNING_PCT = (
    80  # flag site-years below this valid-pixel coverage (TGBS_Kwale precedent)
)

# Habitat classification thresholds (plan §7) -- starting values pending visual calibration
# against high-resolution imagery; do not treat as final.
DW_HABITAT_THRESHOLDS = {
    "crops_min": 0.45,
    "built_min": 0.35,
    "water_wetland_min": 0.45,
    "bare_min": 0.40,
    "woody_min": 0.45,
    "grass_min": 0.45,
    "woody_grass_margin": 0.10,
    "natural_min": 0.55,
    "top1_min": 0.35,
}
DW_HABITAT_CLASS_CODES = list(
    range(1, 9)
)  # classify_habitat never emits 0 (NoData/outside AOI)

# Conversion-pressure thresholds (plan §8)
DW_PRESSURE_THRESHOLDS = {"moderate_min": 0.35, "high_min": 0.55}
DW_PRESSURE_CLASS_CODES = [0, 1, 2]

# All from-class x to-class combinations of the 8 habitat classes (plan §8 transition coding).
DW_TRANSITION_CODES = [
    f * 10 + t for f in DW_HABITAT_CLASS_CODES for t in DW_HABITAT_CLASS_CODES
]

# Objective 4 connectivity-mask thresholds (connectivity-input handoff amendment)
DW_CONNECTIVITY_THRESHOLDS = {
    "natural_prob_min": 0.60,
    "conversion_pressure_max": 0.30,
    "top1_prob_min": 0.35,
}

# Visualization presets (geemap / ee visParams)
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
DW_TRANSITION_VIS = {
    "min": 11,
    "max": 88,
    "palette": [
        "#800026",
        "#f03b20",
        "#fd8d3c",
        "#ffffb2",
        "#c7e9b4",
        "#41b6c4",
    ],
}
