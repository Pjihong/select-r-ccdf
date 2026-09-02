# selectrgev 0.2.1 (2026-09-02)

Updates to the **RGEV11** (nonstationary) code path.
Source files received as `rsel_rgev11_all_13Aug26.R` (13 Aug 2026) and
`rgev11_fit_park_28Feb26.R` (28 Feb 2026); merged into the undated files
in `R/` so that the package keeps a single definition per function.

## Bug fixes

* `rsel.rgev11()`: the MLE was fitted on the **full** `xdat` matrix at every
  step instead of the first `sth` columns. All calls now pass
  `as.matrix(xdat[, 1:sth])` to `rgev11.fit.park()` and `ns.trsf.rlos.Gum()`.
* `ns.trsf.rlos.Gum()`: `maxr` was read before being defined
  (`if (is.null(maxr)) maxr = ncol(xdat)`); it is now set unconditionally
  from `ncol(xdat)`.

## Changes in behaviour

* `rsel.rgev11()`, `method = "ed"` and `method = "spacing"`, `r = 1`:
  the p-value was hard-coded to 1. It is now an actual Cramer-von Mises
  test of the Gumbel-transformed first order statistic against `punif`,
  so selection can now stop at `r = 0`.
* `rsel.rgev11()`: the penalised refit triggered by `kappa < -0.5` now uses
  `reltol = 1e-5` with `num_inits` unchanged (previously `num_inits * 2`).
* `rsel.rgev11()`: default `num_inits` 15 -> 10.
* `rsel.rgev11()`: the returned object gains `$sig` (the significance level).
* `rsel.rgev11()`, `qqplot = TRUE`, `method = "spacing"`: an `r = 1` panel is
  now drawn from the GEV-CDF transform.
* `rgev11.fit.park()`: default `maxit` 1000 -> 200.

## Known limitation

* The `seq.cut = FALSE` branch of `rsel.rgev11()` is commented out in this
  version, so only the sequential-stopping path (`seq.cut = TRUE`, the
  default) is functional. Calling with `seq.cut = FALSE` returns an empty
  `$result`.

# selectrgev 0.2.0 (2026-08-31)

## New features

* New main program `rsel.rgev.SeqStop()` (`R/rsel.rgev.SeqStop.R`):
  sequential selection of the optimal r with stopping rules — raw p-values
  first, then **ForwardStop** and **StrongStop** adjusted for sequential
  multiple testing.
* New GOF helper `rccdf.gof.new()` (`R/rccdf.gof.new.R`):
  conditional complementary-CDF transform for the r-th order statistic.
* New GOF helper `spacing.gof()` (`R/spacing.gof.rgev.R`):
  spacing-based transform (Tawn approach).

## Updates

* `rgev.fit.park()`: tuned default `ntry` values
  (30 → 10 in `rgevmle.park()`, 20 → 10 in `gev.max.consT()`).
* `rk3d.fit.park()`: tuned defaults (`num_inits` 50 → 20, `maxit` 2000 → 1000,
  `reltol` 1e-8 → 1e-5) and revised parameter bounds for `xi` and `h`.
* `NAMESPACE`: export `rsel.rgev.SeqStop()`, `rccdf.gof.new()`, `spacing.gof()`.
* `DESCRIPTION`: declare the `evmr`, `ismev`, and `extRemes` dependencies.

## Housekeeping

* Removed date-suffixed duplicate files from `R/`
  (old versions remain available in the git history).
* Removed `R/sancheong.RDS`; the dataset lives in `data/`.

# selectrgev 0.1.0 (2026-02-28)

* Initial release: `rsel.rgev()` (stationary RGEV) and `rsel.rgev11()`
  (nonstationary RGEV11) with `"ed"`, `"ccdf"`, and `"spacing"` GOF methods.
* Included dataset: `sancheong` — r-largest annual daily rainfall,
  Sancheong station (KMA ASOS 285), 1972–2022.
