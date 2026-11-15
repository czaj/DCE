# Changelog

## Unreleased

- Report the quality of the covariance matrix automatically. MNL, MXL, LC,
  LCMXL, HMNL, HMXL, GMXL, HLC, MIMIC, LML, MDCEV and MMDCEV now print, next to
  the existing `Hessian` row and into the `.xlsx` output, the median
  `SE(BHHH)/SE(H)`, the smallest eigenvalue of the observed Hessian,
  `max |gradient|`, the largest Newton step in standard errors and the number
  of free parameters without a finite standard error; the same numbers are in
  `Results.diagnostics`. A median ratio outside `[0.8, 1.25]` raises a warning:
  the information matrix equality fails, so the point is most likely a maximum
  of one particular set of draws rather than of the likelihood, and the
  standard errors are not fit for a table. Where the comparison cannot be made
  - the reported standard errors are the BHHH ones, the gradient is numerical,
  or the outer product is singular - the output says so instead of showing a
  healthy ratio. The implementation is `Tools/hessianDiagnostics.m`, with the
  bookkeeping the modules share in `Tools/dceHessDiagnostics.m`; it only reads
  matrices the modules have already computed, costs at most one extra gradient
  evaluation (`HessEstFix` 3 and 4, where no jacobian is at hand), and leaves
  estimates and standard errors unchanged. The count data models in
  [CDM](https://github.com/czaj/CDM) report the same diagnostics. Test:
  `DCE/tests/test_hessian_diagnostics.m`.

- Changed behaviour: the Hessian settings are no longer overridden silently.
  Weights together with the analytical Hessian raise an error instead of
  resetting every weight to 1. Each compatibility override in MXL (WTP space,
  covariates of means or scale, non-normal distributions, non-linear
  transformations, missing data, numerical gradient) names the setting it
  changes and both values, and keeps the ex-post estimator the user asked for:
  it turns off `EstimOpt.ApproxHess`, which only concerns what the optimizer is
  given, and changes `EstimOpt.HessEstFix` only when it is 4 - the analytical
  Hessian that those specifications do not have - which becomes 5 with an
  analytical gradient and 0 without one. Standard errors requested with
  `EstimOpt.ApproxHess = 0` and `HessEstFix` 3 or 5 in those specifications
  therefore change: they now come from the requested estimator instead of
  silently from the optimization Hessian. `HessEstFix = 5` works in LC and
  LCMXL, where it used to fail, and `HessEstFix = 4` falls back there to 3 with
  an explicit warning. In HLC the numerical Hessian of `HessEstFix = 3` was
  assigned to `Results.jacobian`, which left `Results.hess` unset and failed
  after the whole estimation; it goes into `Results.hess` now.

- Use scrambled Sobol defaults with 1024 draws, skip 0 and leap 0. A shared
  Tools generator preserves each respondent's draw prefix across draw counts
  and supports the explicit historical global layout. Keep user-specified
  counts and seeds; use power-of-two demand defaults. Update DCE and Tools
  together. See [draw settings and replication](DRAWS.md).

- Add opt-in active-coordinate optimization, reuse of positive definite starting
  curvature and fit-only screening for HMXL, HMNL and MIMIC. Preserve nonzero
  fixed anchors in equality-constrained hybrid estimation. Correct the HMXL/HMNL
  latent-interaction mask with mean covariates/full covariance, ordered-cutoff
  delta standard errors and HLC structural-block indexing in native output,
  and modern fmincon Hessian
  logging. Add a start-rebasing helper for changed draws and OLS/ordered items.
  Defaults and likelihood kernels are unchanged. See the
  [estimation guide](Hybrid/ESTIMATION.md) and
  [validation report](Hybrid/OPTIMIZATION_VALIDATION.md).

- Enable automatic double-precision GPU selection in MXL for supported
  normal/lognormal models, including missing data and varying Xm/Xs.
  Detect MATLAB/device support and memory, compare warmed CPU/GPU calls,
  and retain CPU when GPU is unavailable, slower or fails. Optional
  `EstimOpt.GPU` values are `auto` (default), `cpu` and `gpu`.
  See the [production GPU report](MXL/GPU_AUTO_REPORT.md).

- Add a separate double-precision MXL GPU experiment and correctness/timing
  checks, including missing choices, varying Xm/Xs and a CH optimizer check.
  On D7/R2026b the resident-data GPU medians are 0.058/0.780 seconds for
  CH/pooled, versus 0.321/4.301 seconds with three CPU workers. These initial
  prototype measurements are historical; production integration is described
  above. See the [GPU test report](MXL/GPU_TEST_REPORT.md).

- Extend the compact normal/lognormal choice kernel to missing alternatives and
  tasks, task/alternative-specific mean and scale covariates, HMXL with latent
  scale and diagonal/full coefficient covariance, and all LCMXL classes.
  Preserve public inputs/options and fallback distributions/Hessian workflows.
  See [the extended report](MXL/EXTENDED_MEMORY_REPORT.md) for current scope and
  validation: all 184 extended comparisons passed; final MXL CH/pooled medians
  are 0.162/4.34 seconds and active-worker fault rates meet the reduction target.
  The initial complete-data measurements remain in the historical
  [memory report](MXL/MEMORY_REPORT.md). That CPU-only change did not add GPU.

## 6bf56b8 (2026-10-02)

- Reduce MXL likelihood/gradient allocations for complete-data normal and
  lognormal models in preference and WTP space, with diagonal/full covariance
  and respondent-specific mean covariates. Matrix contractions replace expanded
  gradient arrays; immutable inputs are cached once per existing parallel worker.
  Serial execution and all existing function/option interfaces are retained.
- Add reproducible HEAD comparisons and per-process allocation benchmarks.
  See [the historical memory report](MXL/MEMORY_REPORT.md) for measured results,
  initial scope and its then-existing CT-specific mean/scale limitation, and
  [the GPU assessment](MXL/GPU_ASSESSMENT.md)
  for hardware feasibility. No GPU backend is introduced.
