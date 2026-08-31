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
