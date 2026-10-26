# Hybrid optimization validation — 27 October 2026

The change builds on DCE `8c49be4`, with Tools `feda9f5`. It does not replace the blocked CPU/GPU
likelihood introduced earlier. New optimization options are opt-in. No private
respondent data, model workspaces, correspondence or machine paths are included
in the repository.

## Correctness checks

MATLAB R2026b, double precision:

- `test_hybrid_optimization`: exact nonzero fixed parameters; active-coordinate
  and whitened BFGS; objective/gradient/Hessian chain rules; trust-region with a
  supplied Hessian; numerical-gradient optimization; full-vector callbacks and
  user stop; rejection of indefinite initial curvature; cutoff delta method,
  clipping and latent-link indices for FullCov 0/1/2.
- `test_hybrid_public_options`: six final-output paths through public HMXL,
  HMNL and MIMIC (active optimization and constrained SQP), plus three screening
  paths using curvature. The checks cover observed mean interactions, latent
  masks, fixed nonzero coefficients, covariance offsets in two ordinal items
  and absence of inferential fields from screening results.
- `test_hlc_cutoff_output`: the public HLC path with two ordinal items and
  nonzero fixed class/loadings checks the corrected cumulative thresholds and
  their covariance slices. This is an output regression, not validation of a
  scientifically estimated HLC model.
- `test_hybrid_rebase`: 14 cases spanning HMXL FullCov 0/1/2, HMNL and MIMIC,
  two LVs, mean interactions, and OLS/ordered measurements. At identical raw
  latent values, measurement probabilities and choice indices are preserved to
  less than 1e-12 after a draw-count change.
- Existing CPU regression: 72 likelihood/score cases covering seven measurement
  families, missing items/choices, singleton dimensions and FullCov 2 passed.
  Maximum relative central-difference score error was 6.50e-11; block-gradient
  discrepancy was below 5e-16. Some old legacy-kernel comparisons retain known
  numerical defects; passing does not mean that every legacy value is correct.
- Existing `test_hlc_memory`, `test_hybrid_auto`, and
  `test_hybrid_estimators` also passed. These include CPU fallback and the
  previous public-estimator draw-layout checks.
- In a real three-LV HMXL with 713 respondents and 163 parameters, the 64-draw
  check reproduced the reference objective/scores to 8.01e-16/6.45e-15 relative
  error. All-coordinate finite-difference error was 1.63e-6. No respondent
  probability was floored in that check.

The tiny public-estimator fixtures intentionally stop after one iteration and
are not identified research models. Their standard errors are checked only as
algebraic output transformations, never presented as valid scientific inference.

## Reproduce the public checks

Add current DCE and Tools (including subfolders) to the MATLAB path, then run in
a dedicated session:

```matlab
test_hybrid_optimization;
test_hybrid_rebase;
test_hybrid_public_options(fullfile(tempdir,'dce_hybrid_public'));
test_hlc_cutoff_output(fullfile(tempdir,'dce_hlc_cutoffs'));
test_hybrid_memory(fullfile(tempdir,'dce_hybrid_memory'),0);
test_hlc_memory;
test_hybrid_auto;
test_hybrid_estimators(fullfile(tempdir,'dce_hybrid_estimators'));
```

The research-data comparison requires authorized access to its input MAT and is
not reproducible from public synthetic fixtures alone. Its aggregate results
are a stress-case benchmark, not a broad claim about all HCMs.

