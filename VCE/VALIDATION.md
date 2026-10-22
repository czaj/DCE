# VCE validation and performance, 23 October 2026

The changes cover the volume-choice estimators in this toolbox: MDCEV and
MMDCEV, their likelihood/gradient/Hessian entry points, and demand prediction.
Update both DCE and Tools. See [options and limitations](README.md).

## Benchmark method

Baseline: DCE `e930e8b79f91b0f7f25e01096e1be5b641ecc2f6` and Tools
`fc5cd4bec212661b4bde6d70d41f965105a0f407`. Each version ran in a fresh
MATLAB R2026b process, with two CPU threads, no process pool, two warm-up
evaluations, and at least five timed evaluations over at least 20 seconds.
CPU pairs ran sequentially on the same idle node. Times below are medians
for **one likelihood plus gradient evaluation**, not complete estimation.

Synthetic fixtures have three goods, six tasks per respondent, three utility
attributes, full covariance for MMDCEV, normal/lognormal coefficients, mean
covariates, and a grouped mixed alpha/gamma profile with a profile covariate.
MMDCEV has 17 parameters; MDCEV has 11. The old MDCEV path uses numerical
gradients, as its former public estimator required; the new path is analytical.
Both MMDCEV paths use analytical gradients.

Memory is sampled process private commit from Windows `GetProcessMemoryInfo`
every 200 ms during timed evaluations. It includes MATLAB and cached inputs;
it excludes warm-up peaks. MiB means 2^20 bytes. Page faults include soft
faults and must not be interpreted as disk paging.

## Measured CPU results

| Model and workload | Node | Before (s) | After (s) | Speedup | Before peak (MiB) | After peak (MiB) |
|---|---|---:|---:|---:|---:|---:|
| MDCEV, 10,000 respondents / 60,000 tasks | D4 | 0.43338 | 0.04134 | 10.48x | 795.72 | 831.51 |
| MMDCEV, 1,000 respondents, 2,000 draws | D3 | 4.79555 | 1.68893 | 2.84x | 3,782.44 | 859.36 |
| MMDCEV, 1,000 respondents, 10,000 draws | D3 | 24.35066 | 8.27092 | 2.94x | 15,903.39 | 1,006.23 |

D3: Intel Core i7-7820X, 120 GB RAM. D4: Intel Core i9-7900X, 128 GB RAM
(cluster inventory).
MMDCEV sampled peak private memory fell by 77.3% and 93.7%, respectively.
MDCEV's small test used **35.79 MiB more** (+4.5%) because of cached analytical
inputs. Its speed gain is not a claim of reduced total process memory.

Page faults per evaluation fell from 2,145,073 to 682,901 for MMDCEV at 2,000
draws, from 10,724,336 to 2,469,218 at 10,000 draws, and from 204,183 to 18,108
for MDCEV. These counters describe allocation/memory activity, not hard faults.

Full per-task/respondent outputs were compared, rather than only summed
objectives:

| Pair | Maximum absolute likelihood difference | Maximum absolute score difference |
|---|---:|---:|
| MMDCEV, 2,000 draws: old/new CPU | 7.11e-15 | 2.31e-14 |
| MMDCEV, 10,000 draws: old/new CPU | 3.55e-15 | 2.29e-13 |
| MDCEV: old numerical/new analytical | 1.78e-15 | 1.93e-10 |
| MMDCEV, 2,000 draws: new CPU/GPU | 3.55e-15 | 1.95e-14 |

## GPU validation

On D7, NVIDIA RTX 5060 (driver 617.42), the 1,000-respondent/2,000-draw
MMDCEV fixture took **0.13575 s** per warmed likelihood/gradient evaluation.
The test asserted that the GPU actually executed the calculation. This is
a standalone measurement: D7 differs from D3 and had an existing CPU research
job, so it is **not a controlled CPU/GPU speedup comparison**. Host process
memory counters do not measure device VRAM.

All 47 correctness fixtures also passed on the GPU with block sizes 1 and 3
(94 comparisons), for both value-only and gradient calls.
Maximum relative CPU/GPU errors were 2.67e-16 for values and 9.20e-16 for
gradients. Routing tests cover
CPU, automatic selection, forced GPU, repeated calls, changed coefficients,
changed data, invalid options, and CPU fallback on D4's unavailable device.

## Correctness and compatibility

- 47 likelihood fixtures passed central finite differences for every
  parameter, block-size invariance, weighted objective/score/BHHH checks,
  inactive parameters, and serial versus three-process-worker comparisons.
  Maximum relative finite-difference error was 4.21e-10. Coverage includes
  alpha/gamma/mixed profiles, grouped parameters, profile/mean covariates,
  diagonal/full covariance, fixed/normal/lognormal coefficients, logical
  covariates, one respondent/draw/attribute, and a 200-task panel.
- 45 ordinary fixtures agreed with old likelihood values within 2.37e-16
  relative error. Two additional fixtures use alpha-profile index 50:
  the old subtraction `1-(1-exp(-50))` loses the positive value and yields
  invalid likelihoods, so these are checked against finite differences and
  CPU/GPU agreement instead of requiring the old result.
- Six public-estimator smoke runs checked `HessEstFix=1,3,5` for both models,
  with unequal respondent weights and an inactive parameter. Objective-based
  and gradient-based Hessians agreed within 1.82e-10 relative error. These
  one-iteration, small-sample checks do not establish convergence or valid
  inferential results for a substantive model.
- Twelve demand cases preserved the original scrambled Sobol sequence and
  agreed with old predictions within 7.11e-15 absolute error. Maximum budget
  discrepancy was 4.27e-14. Zero-income and feasible/infeasible numeraire
  cases passed. These tiny checks are not a demand performance benchmark.
- `test_mdcev_mmdcev` passed. Explicit numerical-gradient and value-only
  optimizer paths, and rejection of inconsistent within-panel weights,
  are covered by `test_vce_routing`.
- Local checks on two existing saved milk smoke-test fits (11 respondents)
  passed: MDCEV has 44 parameters, MMDCEV 59 parameters and 64 draws. Maximum
  relative finite-difference error was 4.96e-9, and old/new likelihood error
  was below 2.42e-16. The saved data stayed local and are not in either repo.
  This was not re-estimation of the full milk/flour studies.

## Reproduction

Put both toolboxes and their subfolders on the MATLAB path. Export the old
`probs_mdcev.m`, `probs_mmdcev.m`, `MDCEV_demand.m`, and `MMDCEV_demand.m`
from the DCE baseline commit above into a separate `baseline` folder. Prefix
each filename and its first-line function name with `baseline_`. The old
functions are self-contained; do not add the entire old toolbox to the path.

```matlab
maxNumCompThreads(2);
test_vce_memory('output/checks','baseline',3);
test_vce_gpu('output/checks/vce.mat','output/gpu'); % supported NVIDIA device
test_vce_estimators('output/public');
test_vce_demand('baseline','output/demand');
test_vce_routing;
test_mdcev_mmdcev;
```

Run each benchmark call in a fresh process, using a new output directory:

```matlab
bench_vce_memory('output/old2k','baseline','MMDCEV','legacy',1000,2000);
bench_vce_memory('output/new2k','baseline','MMDCEV','cpu',1000,2000);
% Repeat with 10000 draws; for MDCEV use NP=10000 and R=1.
```

`bench_vce_memory` uses the existing Windows native-memory sampler
`tests/measure_mxl_memory.ps1`. Raw benchmark JSON/MAT files, memory samples,
and correctness CSVs for this run are retained outside the repositories in
`_vce_validation_20261023/output` under the workspace.

Raw random draws and returned scores still scale with sample size, and each
CPU worker owns an input copy. Prediction is CPU-only. Full-data estimation,
bootstrap/inference accuracy, and demand memory/speed at production scale
were not benchmarked here.
