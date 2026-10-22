# Blocked hybrid likelihoods

Update **both DCE and Tools**. Add both repositories (including subfolders) to
the MATLAB path. Tools supplies `hybridLatentMoments`.

HMXL, HMNL, MIMIC and HLC now evaluate respondent blocks. The default path
keeps the existing likelihood arguments and respondent-level output order.
It supports normal, lognormal and fixed choice coefficients, preference/WTP
space, varying mean/scale covariates, latent scale effects, missing choice
tasks/alternatives and missing measurement items. HMXL supports FullCov 0/1/2,
including analytical derivatives with multiple correlated latent errors.
HLC continues to expose a value-only likelihood.

The shared measurement kernel handles OLS, multinomial logit, ordered probit,
Poisson, negative binomial, ZIP and ZINB. Each estimator's existing input
validation still applies. Other choice distributions retain the legacy path.

## Options

- `EstimOpt.GPU = 'auto'` (default): compare warmed CPU and GPU evaluations,
  require numerical agreement, and use GPU only if at least 10% faster.
  Value-only and gradient decisions are separate.
- `EstimOpt.GPU = 'cpu'` or `0`: CPU only, using an existing process pool if
  present. No pool is created by the likelihood.
- `EstimOpt.GPU = 'gpu'` or `1`: request GPU after the same numerical check.
  An unavailable/unsupported device or insufficient GPU memory falls back to
  CPU. Calculations remain double precision.
- `EstimOpt.HybridBlockSize`: optional positive integer overriding automatic
  block sizing. `0` selects the old implementation for diagnostic comparisons,
  including its known numerical defects; it is not the recommended default.
- `EstimOpt.HessEstFix = 5`: HMXL, HMNL and MIMIC now support central differences
  of their analytical gradients through Tools' `hessianFromGradient`, with
  observation weights and inactive-parameter handling. The usual robust
  covariance option remains applicable.

`clear hybrid_likelihood` releases the current cache and resets GPU selection.
Changed data, options, model or process pool replace the cache; changed
parameters reuse it. A GPU device change resets its resident copy and decisions.

## Memory and normalization

Latent variables retain the **global** mean and sample standard deviation over
all respondents and draws, with denominator `NP*NRep-1`. They are never
normalized independently within blocks. Tools computes the centered moments
of structural covariates and draws without expanding the covariates. Those
small moments give the normalization derivatives for each current parameter
vector. Structural and latent covariance scores are contracted after weighting
by the joint choice/measurement posterior, avoiding global derivative tensors.

Measurement categories also come from the whole sample. A block containing
only one category, a singleton respondent or missing measurements cannot
change the parameter layout. Draw-level likelihoods and measurement scores
exist only for the current block. The page-based choice kernel is shared with
the existing GPU MXL implementation.

The original draws remain resident, with one copy per process worker; this is
not shared physical memory. The change concerns likelihood and score
evaluation, not posterior-score/R-squared reporting or all legacy distribution
paths. GPU results are gathered before returning to the optimizer.

## Correctness fixes

Missing measurement items contribute a neutral probability and zero score
for every measurement type. Ordered-probit upper tails avoid subtraction of
two CDFs close to one. Derivatives respect clipped thresholds/scales and the
`RealMin >= 2` probability floor: a locally constant floored objective has
zero derivative. A mining check found two floored respondents for whom the
old code returned nonzero scores, while all other scores agreed.

MIMIC's null-model initialization now calls the existing `LL_hmnl0` function;
the previous name did not exist. Its numerical Hessian is stored in
`Results.hess` rather than `Results.jacobian`. Pseudo-random and Latin
hypercube draws now retain exactly the latent dimensions, as its existing
quasi-random path already did.

See [measured results and reproduction](VALIDATION.md).
