# Changelog

## Unreleased

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
