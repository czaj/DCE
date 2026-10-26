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
  their covariance slices. A structural-equation variant reproduced the old
  output failure when class and LV counts differed, then passed after fixing
  the structural block offset and its LV headers. This is an output regression,
  not validation of a
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

## Matched optimization experiment

A three-LV, full-covariance HMXL stress case has 713 respondents, 163 total and
157 free parameters, and 128 scrambled Sobol draws per respondent. The objective,
weights, guards and draws are identical within each start comparison. All methods
use analytical respondent scores. The two starts are a construction from MIMIC
and a smaller HCM (negative LL 9876.529913348), and a later unfinished checkpoint
(9724.828351813). The latter's observed Hessian was already available.

The first comparison ran on a Xeon E5-2687W v4 with three process workers and two
computational threads per worker/client. Other cluster jobs were running; these
are single-run elapsed times, not repeated idle-node medians. Pool creation,
initial cache warm-up and final inference are outside the timed optimization.
Limits were 1000 iterations, 20000 evaluations and a 360-second callback cap;
a callback can stop only after an iteration. All use optimality tolerance 1e-5,
function/step tolerances 1e-8, interpreted in each solver's coordinates.

### Continuing the unfinished checkpoint

| Method | Final negative LL | Iterations | Evaluations | Seconds | Raw max gradient | Exit |
|---|---:|---:|---:|---:|---:|---:|
| Full-vector BFGS | 9724.810126679 | 1000 | 1004 | 227.77 | 0.002828 | 0 |
| Active, centered BFGS | 9724.810126679 | 1000 | 1004 | 227.74 | 0.002828 | 0 |
| Active L-BFGS | 9724.811097296 | 1000 | 1084 | 338.75 | 0.2221 | 0 |
| Equality-constrained SQP | 9724.822171963 | 42 | 1054 | 364.18 | 2.523 | -1 |
| Trust-region / BHHH | 9724.820616616 | 12 | 13 | 4.73 | 3.633 | 3 |
| BFGS / previous observed-H metric | 9724.810123293 | 3 | 4 | 1.69 | 0.000336 | 1 |

The full-vector and active BFGS calls reached virtually the same objective but
hit the iteration cap. Removing fixed dimensions alone did not accelerate this
case. The available observed-H metric reached a slightly better objective in
three iterations/four evaluations. This is a continuation benefit, not a claim
of a 100-fold complete estimation speed-up from scratch. The cost of obtaining
the starting Hessian must be counted when it is not already available.

Trust-region/BHHH's positive exit code did not establish stationarity: its raw
gradient remained large and the objective was worse. L-BFGS and SQP did not
improve the result within the common budgets. Do not select an algorithm on its
exit flag or speed alone.

### Starting from the MIMIC/smaller-HCM construction

| Method | Final negative LL | Iterations | Evaluations | Seconds | Raw max gradient | Exit |
|---|---:|---:|---:|---:|---:|---:|
| Full-vector BFGS | 9724.889732608 | 1000 | 1027 | 234.67 | 0.251 | 0 |
| Active, centered BFGS | 9724.895865888 | 984 | 1014 | 360.13 | 0.3771 | -1 |
| Active L-BFGS | 9727.572865291 | 1000 | 1078 | 323.17 | 2.113 | 0 |
| Equality-constrained SQP | 9811.454024549 | 47 | 1086 | 366.02 | 156.9 | -1 |
| Trust-region / BHHH | 9741.373439967 | 1000 | 1001 | 283.20 | 11.05 | 0 |

Centering active coordinates changes the initial line-search/stopping scale, so
its path need not match full-vector BFGS even though the statistical objective
is unchanged. These limited runs are comparisons of search behavior; an
iteration/time cap is not a valid final estimate.

### A cheap BHHH metric for BFGS

In a separate experiment on an i9-7920X, using the same draws and starts, the
initial respondent-score outer product supplied a fixed coordinate metric to
BFGS. This differs from supplying BHHH at every trust-region iteration. Hardware
and concurrent load differ from the table above: compare likelihood/evaluation
counts directly, but do not form a cross-machine wall-time speed ratio.

| Start | Final negative LL | Iterations | Evaluations | Optimization seconds | Metric seconds | Raw max gradient | Exit |
|---|---:|---:|---:|---:|---:|---:|---:|
| cold | 9731.043213714 | 105 | 107 | 39.38 | 1.927 | 0.006305 | 1 |
| partial | 9724.810123294 | 54 | 56 | 19.50 | 0.396 | 0.001319 | 1 |

The cheap metric reached essentially the same objective from the checkpoint,
but stopped at a worse objective from the more distant construction. A metric
can accelerate convergence without selecting the best basin. Keep other
start/method combinations rather than treating one fast run as sufficient.

### Independent endpoint checks

Fresh observed Hessians were calculated by central differences of analytical
scores, with one score row per respondent and the corresponding weighted
sandwich. No ridge or pseudoinverse was used.

| Candidate | Observed H minimum eigenvalue | rcond(H) | rcond(sandwich) | Maximum Newton correction / model SE | Local numeric check |
|---|---:|---:|---:|---:|---|
| Observed-H metric, checkpoint | 0.237841 | 8.48e-08 | 3.02e-09 | 4.64e-06 | PASS |
| BHHH metric, cold | 0.0149445 | 3.48e-09 | 6.69e-12 | 7.81e-05 | PASS |
| BHHH metric, partial | 0.237841 | 8.48e-08 | 3.02e-09 | 1.35e-05 | PASS |

For the observed-H continuation, the independent final Hessian took 138.71 s
on a different processor (i9-7900X); this was a validation/inference cost, not
part of the 1.69 s optimization timing. Positive definiteness and a small
Newton correction check local numerical stationarity, not global optimality,
construct validity or adequate integration precision. The raw conditioning is
poor enough that units and scaling still matter.

The original batch harness attempted to create an already existing process
pool after completing its first method. Those completed results were retained;
a separate batch helper closed pools between remaining methods, skipped
completed outputs, and preserved the original source/manifest. No fit was
restarted or duplicated to obtain a preferred result.

## Scope and limits

The tests support the implementation and the staged workflow in the
[estimation guide](ESTIMATION.md). They do not establish a universal optimizer
ranking. There was no broad high-draw or many-data-set search, no new GPU timing
of the optimizer adapter, and no validation on earlier MATLAB releases.
Existing GPU kernels were unchanged. HLC retains its value-only likelihood and
legacy optimization; its shared equality/cutoff/logging fixes were tested, but
the new optimization options are for HMXL/HMNL/MIMIC.

All low-draw benchmark outcomes are retained, including capped and inferior
ones. Final scientific models still require multistart, integration stability,
measurement identification and valid respondent-level inference.
