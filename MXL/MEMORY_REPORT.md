# MXL memory allocation report

Date: 2026-10-02. Baseline: `dd6f704`. Branch: `mxl-memory`.

## Scope and status

This change reduces repeated allocation and worker transfer in `LL_mxl`
without changing its arguments, outputs, `MXL`, or `EstimOpt` options.
The frozen toolbox and data under
`C:\Users\miq\Documents\lasy\replication_package` were read-only throughout.
Work and measurements started after both NEE multi-start queues finished and
their 19 expected result files were received and validated. Benchmark batches
were run sequentially on an otherwise computation-idle computer.

The numerical comparisons and the CH re-estimation passed. Both CH and pooled
performance series meet the worker page-fault reduction target with shorter
evaluation times. The before/after measurements use R2026a consistently;
the user's current default is R2026b. Its separate numerical and estimation
validation also passed and is reported below, without mixing releases in the
R2026a performance ratios.

## Implementation

The new path handles complete data, one or two outputs, diagonal or full
covariance, preference or WTP space, and distributions `Dist = -1, 0, 1`
(fixed, normal, lognormal). It requires `mCT = 0`, `NVarS = 0`, no nonlinear
transformation, and no Johnson transformation. Respondent-specific mean
covariates are supported. WTP mappings must refer to the configured cost
parameters, including models with several mapped costs.

Other distributions, missing data, choice-task/alternative-specific mean
covariates, scale covariates, and analytic Hessians keep the existing paths.
The existing special all-normal, diagonal WTP path for `mCT > 0` is also kept.
These fallback paths were not rewritten as part of the memory optimization.

For each respondent, let `X` be the alternative-by-attribute design matrix,
`P` the alternative-by-draw probability matrix, `E` the standard-normal
draws, and `panel` the product of chosen probabilities across choice tasks.
The former attribute-by-task-by-draw expansion is replaced by the contraction

```matlab
D = sum(X(chosen,:),1)' - X'*P;
```

Distribution and WTP chain rules are applied to the resulting
`NVarA x NRep` array. In the full-covariance gradient, the weighted derivative
is contracted with `E'` before selecting the Cholesky entries. This replaces
the `NVarCholesky x NRep` temporary with a small `NVarA x NVarA` product.
The diagonal-covariance branch uses an elementwise reduction of the same
compact arrays. A subnormal-probability guard retains the legacy multiplication
order when `RealMin = 0` and the simulated probability lies below `realmin`.

Coefficients are generated for each respondent on the worker. The client no
longer constructs the full `NVarA x (NRep*NP)` coefficient matrix in this path.
For the pooled forest case, one such double matrix would contain
1,215,840,000 bytes (1.216 GB); the corresponding draw matrix has the same
size. Static `YY`, `XXa`, `XXm`, and `err` are stored in a
`parallel.pool.Constant`, so each worker receives one full copy per dataset
rather than new coefficient/draw slices at every evaluation. The trade-off is
resident static data per worker, not shared memory between workers.

The cache checks both the pool identity and exact input contents using
`isequaln`, not only array dimensions. It refreshes when any cached input
changes and retains at most one dataset. A fallback evaluation may leave the
previous cached dataset resident. Replacing the dataset or evaluating the new
path without a pool releases the previous constant; `clear LL_mxl` also
releases the cached state. The function uses the current pool if present and
does not create one. Without a pool, `parfor (n = 1:NP,0)` executes locally.

Respondent blocking was not added: the matrix contractions remove the large
per-respondent buffers while retaining the existing respondent-level work
distribution. The measurements below support this smaller implementation.
No GPU path was implemented or benchmarked; see the separate
[GPU assessment](GPU_ASSESSMENT.md).

## Numerical validation

The comparison harness extracts `dd6f704:MXL/LL_mxl.m` outside the frozen
replication directory and renames only its entry function. Old and new
functions receive identical parameters, data, options, weights, and Sobol
draws. Both per-person outputs and the entire weighted gradient are compared,
not just selected gradient entries. The reported LL is
`-sum(W.*f)` because `LL_mxl` returns negative log likelihood by person.

Relative errors use `abs(delta)/max(1,abs(reference))` for LL and
`max(abs(delta))/max(1,max(abs(reference)))` for vectors or arrays. Assertions
require LL/value errors at most `1e-8` and gradient errors at most `1e-6`.
Tests also check finite outputs and agreement of value-only and gradient
branches. Rows labelled serial have no pool; parallel rows use three process
workers.

### Core cases

The requested raw `NEWFOREX_DCE_demo.mat` is used with its long-format `Y`
choice indicator and `SKIP` missing indicator, following the demo design.
The forest design and sample filters come from `load_forest_data.m`; CH is
`estim_sample & country_model == "CH"`, and pooled uses `Weq`. Parameters
come from the published CH and pooled MXL results. The published MAT contains
estimates but no saved `INPUT`, `EstimOpt`, or optimizer options, so the
fixtures reconstruct those from the replication specification and toolbox
defaults. Forest tests use 17 random coefficients and 1,000 Sobol draws.

| Case | NP | LL, old and new | Max relative gradient error, serial/3 workers |
| --- | ---: | ---: | ---: |
| Demo preference, diagonal covariance | 789 | -6928.16606386632 | 2.22e-14 |
| Demo preference, full covariance | 789 | -6842.47599520024 | 4.62e-14 |
| Demo WTP, diagonal covariance | 789 | -7051.30135955028 | 2.60e-15 |
| Demo WTP, full covariance | 789 | -7024.18210528086 | 8.48e-16 |
| Forest CH, WTP/full covariance | 644 | -4793.81636376126 | 1.58e-12 |
| Forest pooled, WTP/full covariance, Weq | 8940 | -80435.0119931802 | 2.91e-15 |

All 12 core comparisons passed. The weighted LL difference was zero in these
runs; CH rounds to the required `-4793.8164`. Across all core rows, the maximum
per-person value difference was `3.55e-15`, the maximum absolute per-person
gradient difference was `1.07e-12`, and the maximum relative per-person
gradient difference was `7.72e-14`. The largest relative value-only difference
was `2.13e-16`.

The following are single warmed evaluation timings recorded during correctness
checks, not independent performance-series medians. All times are seconds.

| Case | Gradient serial HEAD/new | Gradient 3 workers HEAD/new | Value serial HEAD/new | Value 3 workers HEAD/new |
| --- | ---: | ---: | ---: | ---: |
| Demo preference, diagonal | 0.664 / 0.282 | 0.772 / 0.203 | 0.395 / 0.184 | 0.184 / 0.166 |
| Demo preference, full | 0.759 / 0.296 | 0.737 / 0.192 | 0.215 / 0.181 | 0.167 / 0.163 |
| Demo WTP, diagonal | 0.966 / 0.350 | 0.939 / 0.188 | 0.246 / 0.200 | 0.184 / 0.159 |
| Demo WTP, full | 0.902 / 0.322 | 0.932 / 0.185 | 0.252 / 0.219 | 0.182 / 0.167 |
| Forest CH | 2.507 / 0.417 | 2.460 / 0.178 | 0.324 / 0.218 | 0.221 / 0.146 |
| Forest pooled | 33.951 / 5.335 | 31.005 / 2.047 | 3.870 / 3.050 | 2.919 / 1.710 |

### Supplemental and fallback checks

All 16 supplemental comparisons passed: two mapped WTP costs, respondent mean
covariates, same-shape changes to `XXa`, `err`, `YY`, and `XXm`, diagonal
covariance, a fixed coefficient, and `RealMin = 1` underflow handling. The
maximum relative gradient error was `2.35e-16`; values were identical.

The covariate matrix contains 102 checks across serial and three-worker modes:
preference/WTP, diagonal/full covariance, no/respondent/task/alternative mean
covariates, and no/respondent/task scale covariates. Scale values are constant
across alternatives within a task, matching the existing scale-gradient
representation. It also checks the existing all-normal diagonal `mCT` WTP
path and subnormal probabilities with clipping disabled.

There were 70 successful numerical gradient comparisons (maximum relative
error `2.66e-16`) and 32 unchanged legacy gradient exceptions. The latter
combine task/alternative-specific mean covariates with scale covariates:
both HEAD and the candidate raise the same
`MATLAB:getReshapeDims:notSameNumel` error with the same message. They count
as preserved behavior, not successful gradient support. Value-only evaluations
passed numerically for all 102 cases. The two subnormal cases reproduce the
expected complete gradients `[1 0 2 0]` and `[1 0 2 0 0]`.

The existing repository test `test_mdcev_mmdcev` also passed.

### CH re-estimation

The full MXL estimation entry was run with three workers, starting at the
published CH parameters. It follows the trust-region/user-Hessian stage and
quasi-Newton stage, using reconstructed replication settings. This is a
published-optimum reproduction check, not a new multi-start search or evidence
that the published optimum is globally best.

| Workers | Optimizer time (s) | Published LL | Re-estimated LL | Relative LL error | Relative parameter error |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 31.7373 | -4793.81636376125 | -4793.81636367476 | 1.80e-11 | 3.96e-7 |

Both assertions passed. The absolute LL difference was `8.65e-8`.

### R2026b validation

The same validation was repeated on the current default MATLAB R2026b,
version `26.2.0.3386108`, in a dedicated session. The run exited successfully
and `validate_R2026b.log` ends with `R2026b VALIDATION PASSED`. All result CSV
and MAT files were present and nonempty before this section was added.

| Check | R2026b result |
| --- | --- |
| 12 core comparisons, serial and 3 workers | Passed; weighted LL difference 0; maximum relative weighted-gradient error 1.58e-12 |
| Per-person core outputs | Maximum absolute value difference 3.55e-15; maximum absolute gradient difference 1.07e-12; maximum relative gradient error 7.72e-14 |
| 16 supplemental comparisons | Passed; values identical; maximum relative gradient error 2.35e-16 |
| 102 covariate/fallback checks | 70 numerical gradient comparisons passed, maximum relative error 2.66e-16; 32 identical existing reshape exceptions; all 102 value-only comparisons passed |
| Full CH estimation, 3 workers | 32.1302 s; LL -4793.81636367476; relative LL error 1.80e-11; relative parameter error 3.96e-7 |
| Existing `test_mdcev_mmdcev` | Passed |

The complete saved R2026a benchmark outputs were also checked in this R2026b
validation session. `performance/equivalence.csv` compares every per-person
value and gradient entry for HEAD+parfor and the new implementation against
actual HEAD, for both CH and pooled. All four comparisons passed. Maximum
absolute gradient differences for HEAD+parfor were `8.88e-15` (CH) and
`2.00e-15` (pooled); those for the new implementation were `1.07e-12` and
`1.91e-13`. This is a numerical check of the recorded benchmark results, not
a second performance benchmark on R2026b.

R2026b `checkcode` completed without syntax errors for the changed MATLAB
source and harnesses. It reports advisory warnings for existing fallback
code, broadcast handles/options, the estimation harness's existing global
backup interface, small growing test-result arrays, and benchmark PID lookup.
The native Windows counter self-check also passed.

## Performance method

Hardware: Windows 10, i9-13900KS (24 cores/32 threads), 192 GB DDR5-4200.
MATLAB: R2026a Update 3 (`26.1.0.3276743`) with three process workers.
Each implementation runs in a fresh batch session on identical fixture inputs.
`pooled_HEAD_initial` generated and saved the pooled fixture before warm-up;
the diagnostic and new runs loaded that same saved `C`. Thus their timed
inputs, including parameters and draws, are identical. Data preparation, the
`DataCleanDCE` separation-check LP, draw generation, pool startup, JIT, and
initial constant transfer are outside the evaluation timings. The separation
LP is a setup diagnostic and does not change the likelihood definition; its
potentially long preparation time is not attributed to likelihood evaluation
or Sobol generation. The test fixture now sets `CheckSeparation = 0` to skip
that pre-estimation LP on future fixture construction, without changing the
numerical likelihood inputs or definition.

Two warm-up evaluations precede the measured series. The series contains at
least five value-plus-gradient evaluations and at least 20 seconds of calls.
The reported time is the median of individual calls. A PID handshake identifies
the client and all three workers before the 1 Hz sampler starts.

Windows `GetProcessMemoryInfo` supplies cumulative page faults, working set,
private bytes, and native lifetime memory peaks. Fault rates are counter deltas
divided by elapsed sample time; DWORD wrap is handled. These are **all page
fault types**, not an isolation of demand-zero faults. The measurements show
allocation-related fault reduction but do not measure the exact demand-zero
fraction or bytes zeroed by the kernel.

Sampled memory peaks exclude warm-up and can miss sub-second spikes. Native
per-process lifetime peaks include startup and warm-up and are reported
separately. Total sampled peaks are maxima of simultaneous process sums, not
sums of independent peaks. Summed working sets may double-count shared pages;
summed private bytes are the more useful process-memory comparison here.

### Baseline distinction

The actual `dd6f704` value-plus-gradient respondent loops use `for`, so its
three pool workers are largely idle in these benchmarks. The frozen NEE copy
was inspected read-only and has `parfor` gradient loops. A separate diagnostic
`HEAD+parfor` copy changes only the two gradient respondent loops of HEAD to
`parfor`. It is useful for the worker allocation target, but is **not** the
unmodified HEAD baseline or a direct benchmark of the frozen toolbox.

### Evaluation time and sampled total memory

MiB means `2^20` bytes. The CH HEAD row uses the completed cached-fixture run;
fixture caching affects preparation, not the timed likelihood calls.

| Case | Implementation | Calls | Median evaluation (s) | Sampled peak private (MiB) | Sampled peak working set (MiB) |
| --- | --- | ---: | ---: | ---: | ---: |
| CH, NP=644 | HEAD | 8 | 2.5066 | 4072.3 | 3952.9 |
| CH, NP=644 | HEAD+parfor diagnostic | 11 | 1.8615 | 6735.4 | 6587.7 |
| CH, NP=644 | New | 129 | 0.1546 | 4241.7 | 4101.1 |
| Pooled, NP=8940 | HEAD | 5 | 36.8380 | 9738.7 | 9589.3 |
| Pooled, NP=8940 | HEAD+parfor diagnostic | 5 | 24.2897 | 44618.7 | 42759.4 |
| Pooled, NP=8940 | New | 10 | 2.0437 | 8819.6 | 8679.1 |

For CH the new median is 16.22 times shorter than actual HEAD and 12.04 times
shorter than HEAD+parfor. Sampled total private memory is 37.0% lower than
HEAD+parfor, but 4.2% higher than HEAD, whose gradient work stays on the client.
This is consistent with retaining static data on the newly active workers.

For pooled, the new median is 18.03 times shorter than actual HEAD and 11.89
times shorter than HEAD+parfor. Sampled total private memory is 9.4% lower
than HEAD and 80.2% lower than HEAD+parfor. All three performance-series
implementations returned the same pooled LL, `-80435.01199318023`.

### CH page faults and per-process memory

Fault rates below are per-process means over the sampled steady-state window;
the peak column is the largest sampled rate, not an instantaneous peak.
Memory columns are private commit in MiB.

| Implementation | Process | Mean faults/s | Peak sampled faults/s | Sampled peak private | Native lifetime peak private |
| --- | --- | ---: | ---: | ---: | ---: |
| HEAD | Client | 837440 | 869667 | 1437.1 | 1501.6 |
| HEAD | Worker 1 | 8.4 | 163.4 | 891.8 | 902.8 |
| HEAD | Worker 2 | 7.5 | 136.7 | 905.0 | 913.4 |
| HEAD | Worker 3 | 8.1 | 157.5 | 838.4 | 847.3 |
| HEAD+parfor | Client | 326228 | 624329 | 2522.4 | 2730.4 |
| HEAD+parfor | Worker 1 | 489941 | 625854 | 1416.8 | 1557.6 |
| HEAD+parfor | Worker 2 | 485020 | 611609 | 1410.4 | 1557.9 |
| HEAD+parfor | Worker 3 | 485818 | 618171 | 1414.8 | 1559.7 |
| New | Client | 131.9 | 886.2 | 1189.4 | 1333.7 |
| New | Worker 1 | 351.9 | 7130.5 | 1026.7 | 1170.5 |
| New | Worker 2 | 306.2 | 6194.6 | 1016.9 | 1164.2 |
| New | Worker 3 | 375.4 | 6588.0 | 1017.9 | 1170.4 |

Compared with the parallel diagnostic baseline, the mean worker rates fall by
1392, 1584, and 1294 times, respectively. CH therefore exceeds the requested
fivefold reduction on every active worker without an evaluation-time increase.
The idle HEAD workers are not used to claim a worker reduction.

### Pooled page faults and per-process memory

The column definitions and units are the same as in the CH table.

| Implementation | Process | Mean faults/s | Peak sampled faults/s | Sampled peak private (MiB) | Native lifetime peak private (MiB) |
| --- | --- | ---: | ---: | ---: | ---: |
| HEAD | Client | 789874 | 1265213 | 7049.3 | 7049.4 |
| HEAD | Worker 1 | 12.8 | 2167.9 | 898.9 | 905.0 |
| HEAD | Worker 2 | 12.3 | 2094.4 | 897.1 | 902.8 |
| HEAD | Worker 3 | 15.2 | 2179.8 | 893.5 | 893.5 |
| HEAD+parfor | Client | 354563 | 1281173 | 21012.8 | 21203.2 |
| HEAD+parfor | Worker 1 | 527420 | 1424904 | 9749.2 | 9749.3 |
| HEAD+parfor | Worker 2 | 527401 | 1067428 | 8074.7 | 8245.9 |
| HEAD+parfor | Worker 3 | 527289 | 1211891 | 9750.2 | 9750.4 |
| New | Client | 12179.1 | 20076.6 | 2398.5 | 4763.7 |
| New | Worker 1 | 4068.9 | 6917.4 | 2139.2 | 4531.5 |
| New | Worker 2 | 4067.3 | 6918.2 | 2146.5 | 4538.7 |
| New | Worker 3 | 4069.9 | 6923.8 | 2136.9 | 4528.8 |

Compared with HEAD+parfor, the mean pooled worker rates fall by 129.6, 129.7,
and 129.6 times. Every active worker exceeds the fivefold target. The native
lifetime peaks are higher than the sampled new steady-state peaks because
they include pool startup and initial data transfer, which the steady-state
sampler deliberately excludes.

### Faults normalized by evaluation count

As the new implementation completes more calls per second, rates are also
normalized by the number of evaluations. These are total worker faults over
the sampled series divided by completed calls, not direct allocation counts
and not demand-zero-only counts. The sampled window can include a short
handshake/completion interval around the timed calls.

| Case | Implementation | Worker 1 faults/evaluation | Worker 2 faults/evaluation | Worker 3 faults/evaluation |
| --- | --- | ---: | ---: | ---: |
| CH | HEAD, workers largely idle | 21.1 | 18.9 | 20.4 |
| CH | HEAD+parfor diagnostic | 942073 | 932490 | 934025 |
| CH | New | 57.8 | 50.3 | 61.6 |
| Pooled | HEAD, workers largely idle | 473.6 | 457.8 | 563.0 |
| Pooled | HEAD+parfor diagnostic | 12863244 | 12862291 | 12859563 |
| Pooled | New | 8623.9 | 8619.4 | 8625.0 |

## Raw results and reproduction

All raw outputs are outside the frozen replication package, under
`C:\Users\miq\Documents\_dce_memory`:

- `correctness_v1/correctness.csv` and `.mat`: 12 core comparisons, complete
  old/new per-person values and gradients, parameters, options, and weights.
- `correctness_v1/supplemental_0.csv`, `supplemental_3.csv`, and their MAT files:
  16 supplemental comparisons.
- `covariates_v2/covariates.csv` and `.mat`: 102 fallback/covariate checks,
  including explicit baseline and candidate exception records.
- `CH_estimation_v1/CH_estimation.csv` and `.mat`: published and re-estimated
  CH results, including both optimizer stages.
- `correctness_R2026b`, `covariates_R2026b`, and `CH_estimation_R2026b`:
  repeated R2026b checks with the same CSV/MAT layout, including both
  supplemental files. `validate_R2026b.log` records the release and final
  validation marker.
- `performance/equivalence.csv`: four complete-array comparisons of the
  recorded CH/pooled diagnostic and new benchmark outputs against HEAD.
- `performance/CH_HEAD_cached`, `CH_HEAD_parfor`, and `CH_new`: `summary.json`,
  `evaluation.json`, `process_samples.tsv`, `process_summary.tsv`, MATLAB log,
  and result MAT files. JSON contains every process PID and unrounded metrics.
- `performance/pooled_HEAD_initial`, `pooled_HEAD_parfor`, and `pooled_new`:
  completed pooled series with the same output layout as the CH series.
  `performance/input_pooled.mat` is the shared fixture generated before the
  initial HEAD warm-up and subsequently loaded unchanged.

Runnable entry points, each in a dedicated MATLAB session with no existing
pool:

```matlab
addpath('C:\Users\miq\Documents\MATLAB\GitHub\DCE\tests');
out = 'C:\Users\miq\Documents\_dce_memory\recheck';
test_mxl_memory(out,[0 3],'compare');
test_mxl_covariates(fullfile(out,'baseline'),fullfile(out,'covariates'),[0 3]);
test_mxl_memory(fullfile(out,'CH_estimation'),3,'estimate');
```

`tests/measure_mxl_memory.ps1` starts and samples a separate batch session.
Its `-Case` is `CH` or `pooled`; `-FunctionName` selects `LL_mxl`,
`LL_mxl_baseline`, or `LL_mxl_baseline_parfor`. Baseline runs also receive
`-BaselineDir` pointing to the exported HEAD function. A fresh `-OutDir` is
required per series. It rejects output paths inside `replication_package` and
refuses to start if another MATLAB process exists; the operator must also
confirm that other computations are idle. `-SelfTest` checks native counters
without launching MATLAB. The sampler now defaults to the user's R2026b
installation. To reproduce the R2026a performance series in this report, pass
`-Matlab 'C:\Program Files\MATLAB\R2026a\bin\matlab.exe'` explicitly for
every baseline and new run.
