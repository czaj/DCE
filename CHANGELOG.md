# Changelog

## Unreleased

- Reduce MXL likelihood/gradient allocations for complete-data normal and
  lognormal models in preference and WTP space, with diagonal/full covariance
  and respondent-specific mean covariates. Matrix contractions replace expanded
  gradient arrays; immutable inputs are cached once per existing parallel worker.
  Serial execution and all existing function/option interfaces are retained.
- Add reproducible HEAD comparisons and per-process allocation benchmarks.
  See [the memory report](MXL/MEMORY_REPORT.md) for measured results, scope and
  the existing CT-specific mean/scale limitation, and [the GPU assessment](MXL/GPU_ASSESSMENT.md)
  for hardware feasibility. No GPU backend is introduced.
