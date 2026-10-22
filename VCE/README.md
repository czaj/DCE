# MDCEV and mixed MDCEV: blocked likelihoods

Update **both DCE and Tools**, including their subfolders on the MATLAB path.
These are the volume-choice estimators currently provided by DCE. They share
one double-precision likelihood/score implementation; their public arguments
and output ordering are unchanged. MDCEV returns one negative log likelihood
per task, and MMDCEV returns one per respondent after joint panel integration.

## What changed

- Likelihoods run in bounded blocks. Random coefficients and mean shifters
  are no longer expanded over all respondents and draws together. Attribute
  gradients use matrix products instead of an alternatives-by-tasks-by-draws-
  by-parameters tensor. Profile scores are integrated before expansion over
  profile covariates.
- MDCEV now supplies analytical gradients. Both estimators default to
  `NumGrad=0`; explicitly requested numerical gradients remain available.
- Both estimators support `HessEstFix=5`, central differences of weighted
  analytical gradients using Tools' `hessianFromGradient`. Ordinary BHHH
  (`HessEstFix=1`) also uses analytical scores when `NumGrad=0`.
- The sign of MMDCEV's objective-based numerical Hessian (`HessEstFix=3`)
  is corrected. Observation weights are no longer silently replaced with
  ones when `ApproxHess=0`. Panel weights must be constant within a respondent.
- The Jacobian sum is stabilized in log space. Alpha-profile arithmetic
  preserves `1-alpha=exp(-fit)` instead of subtracting a rounded alpha from
  one, including when the profile index is large.

Normal, lognormal and fixed mixed coefficients, diagonal/full covariance,
alpha/gamma/mixed profiles, grouped profile parameters and covariates are
covered. Logical indicator covariates are converted to doubles once in the
input cache. Existing restrictions on missing tasks/alternatives remain:
quantities must be finite and nonnegative, prices positive, and each task
must contain at least one consumed good.

## Options and execution

```matlab
EstimOpt.NumGrad = 0;
EstimOpt.GPU = 'auto'; % default: measure and compare CPU/GPU
EstimOpt.HessEstFix = 5; % optional Hessian after optimization
```

`GPU='auto'` uses a GPU only after warmed measurements show at least a 10%
speed improvement and matching values/gradients. `GPU='cpu'` (or `0`)
disables GPU evaluation; `GPU='gpu'` (or `1`) requests it after the same
correctness check. Unsupported/unavailable devices and insufficient memory
fall back to CPU. Value-only and gradient decisions are calibrated separately.

`VCEBlockSize` optionally sets a positive integer number of respondents
(MMDCEV) or tasks (MDCEV) per block. Automatic sizing targets approximately
64 MiB of CPU temporaries and available GPU memory. An existing process pool
is used; the likelihood never starts a pool. Each worker keeps its own cached
input copy. Changed inputs/options/pool replace the cache, while changed
parameters reuse it. `clear vce_likelihood` releases the cache.

Raw draws and returned respondent scores still occupy memory proportional
to sample size. This is not shared physical memory between workers. Input
draw-generation methods and their random sequences have not been changed.
First-call setup/calibration costs are separate from warmed evaluation times.

## Demand prediction

`MDCEV_demand` and `MMDCEV_demand` now generate utility errors and solve budget
equations in bounded blocks, using Tools' shared `mdcevDemand`. Bisection
works on all task/draw pairs in a block together. MMDCEV selects `NSim` draws
before constructing random coefficients, instead of constructing all `NRep`
coefficient draws and then discarding most of them. Sobol scrambling and draw
ordering are preserved for the same random state.

These prediction functions currently run on CPU; `GPU` controls likelihood
evaluation. `VCEBlockSize` also bounds prediction blocks. The original random
draw matrix may still be regenerated in full when it is absent from saved
MMDCEV results. Invalid profile parameters and infeasible numeraire requests
raise an error rather than silently returning invalid demand.

See [validation and measured performance](VALIDATION.md). Existing inference
options retain their meanings: `HessEstFix=5` does not itself request clustered
or sandwich covariance, and no new analytical Hessian (`HessEstFix=4`) is added.
