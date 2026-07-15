# Derived-raster recodes from a habitat_class raster. Expect values 1-8 with NA elsewhere if the
# raster was loaded via read_habitat_raster() (which already substitutes literal 0 -> NA -- see
# that function's docs for why). The explicit "0 -> NA" mappings below are defense-in-depth for
# any caller that passes a raw, not-yet-normalized raster directly.

#' Full habitat-class raster for composition/fragmentation metrics: masks out water/flooded
#' vegetation (7) and uncertain (8), keeping classes 1-6.
recode_full_habitat <- function(r) {
  # terra::classify()'s default (others = NULL) leaves any value not in `rcl` unchanged -- exactly
  # what's wanted here for classes 1-6 (not listed below), so no `others=` argument is needed.
  out <- terra::classify(r, rcl = matrix(c(0, NA, 7, NA, 8, NA), ncol = 2, byrow = TRUE))
  names(out) <- "habitat_class"
  out
}

#' Binary natural-habitat raster: 1 = natural (woody/grassland/mixed natural), 0 = conversion
#' pressure classes (cropland/built/bare-degraded), NA = no classification/water/flooded veg/uncertain.
make_natural_binary <- function(r) {
  out <- terra::classify(
    r,
    rcl = matrix(c(0, NA, 1, 1, 2, 1, 3, 1, 4, 0, 5, 0, 6, 0, 7, NA, 8, NA), ncol = 2, byrow = TRUE)
  )
  names(out) <- "natural_binary"
  out
}

#' Current-period binary natural-habitat raster, reconciled against Objective 1's own
#' natural_habitat_mask export.
#'
#' Objective 1's natural_habitat_mask file is 1/NA-only (built via .selfMask()) -- it cannot by
#' itself represent "0 = converted but valid", so it can't be used as a drop-in replacement for
#' make_natural_binary(). Instead: build the 0-class from a plain reclass of the current
#' habitat_class raster (as make_natural_binary() does for baseline/pre, which have no equivalent
#' Objective 1 export), then DEMOTE any pixel the plain reclass called "natural" (1) to NA
#' wherever Objective 1's stricter mask (top1_prob/valid_obs_count QA gate, see
#' DW_CONNECTIVITY_THRESHOLDS in config.py) does NOT also confirm it as natural. Never promote a
#' pixel the plain reclass called 0 or NA.
#'
#' Net effect: identical to plain make_natural_binary(current_class_r) except a small number of
#' natural-class edge pixels that fail the stricter QA gate get demoted to NA rather than left as
#' 1. Verify the two rasters agree on >95% of valid pixels (see 01_prepare_inputs.R) -- a bigger
#' gap means the QA gate is dropping more than expected and current-period connectivity metrics
#' need investigation before being trusted.
build_current_binary_natural <- function(current_class_r, natural_habitat_mask_r) {
  plain <- make_natural_binary(current_class_r)
  # natural_habitat_mask_r is 1/NA; align it to plain's grid before comparing.
  mask_aligned <- terra::resample(natural_habitat_mask_r, plain, method = "near")
  downgrade <- is.na(mask_aligned) & (plain == 1)
  out <- terra::ifel(downgrade, NA, plain)
  names(out) <- "natural_binary_current_reconciled"
  out
}

#' % of valid (non-NA) pixels where `plain` and `reconciled` agree, restricted to pixels valid in
#' at least one of the two rasters. Called from 01_prepare_inputs.R; prints a loud PASS/investigate
#' message rather than failing silently.
check_mask_reuse_agreement <- function(plain, reconciled) {
  p <- terra::values(plain, na.rm = FALSE)
  r <- terra::values(reconciled, na.rm = FALSE)
  either_valid <- !is.na(p) | !is.na(r)
  agree <- (is.na(p) & is.na(r)) | (!is.na(p) & !is.na(r) & p == r)
  pct_agree <- if (sum(either_valid) > 0) sum(agree & either_valid) / sum(either_valid) else NA_real_
  status <- if (!is.na(pct_agree) && pct_agree > 0.95) "PASS" else "INVESTIGATE"
  message(sprintf(
    "[%s] Current-period binary-natural mask-reuse agreement: %.2f%% of valid pixels (expect >95%%).",
    status, 100 * pct_agree
  ))
  invisible(pct_agree)
}
