# Extended MXL-family memory allocation report

Date: 2026-10-02. Initial implementation: `6bf56b8`.
Original likelihood reference: `dd6f704`. Branch: `mxl-memory`.

Status: final numerical suites, CH re-estimation, native counter self-check,
full-array benchmark comparisons and local MXL measurements passed. Cluster
outputs were received and verified before this report was finalized. Nine
source/test file hashes matched the validated snapshot. The frozen
replication package and unrelated user changes were not modified.

## Scope

The shared compact choice kernel now handles fixed, normal and lognormal
coefficients (`Dist = -1, 0, 1`), preference/WTP space with one or multiple
mapped costs, diagonal/full coefficient covariance, and missing alternatives
or entire choice tasks. An entirely missing choice task contributes a neutral
probability and no choice score. Function arguments, outputs and existing
`EstimOpt` options remain unchanged.

| Entry | Compact path | Retained paths and limits |
| --- | --- | --- |
| MXL | Value/gradient; respondent, CT and alternative-specific Xm; respondent, CT and alternative-specific Xs | Other distributions, nonlinear/Johnson transformations and analytic Hessian requests keep the existing implementation |
| HMXL | FullCov 0/1; Xm/Xs variants; latent mean interactions; ScaleLV; existing measurement equations | FullCov 2 and other distributions retain the old choice path; global latent normalization and measurement equations are not rewritten |
| LCMXL | All classes use the shared choice kernel, static draw/data cache and mixture-score integration; class-specific distributions/scale | No Xm argument is added to the LCMXL API; unsupported distributions keep the old path |

LCMXL's existing estimation-entry `NumGrad` switches remain unchanged,
including switches for missing alternatives or multiple WTP costs. The direct
likelihood's compact analytic gradient is tested, but the public estimator
may still select numerical gradients for those specifications. This change
does not override existing caller options or constraints.

Invalid/noncost WTP mappings do not enter the compact path. The legacy direct
HMXL two-dimensional mCT Xm representation is reshaped before caching. HMXL
and LCMXL locate their sibling MXL helper directory when only their own folder
has been added to the MATLAB path.

MXL/HMXL preprocessing identifies respondent-constant Xm using available Y
rows, rather than requiring the first alternative to be available. It does
not silently discard invalid Xm values on available alternatives. Respondents
with no available choice rows receive neutral collapsed Xm values.

## Implementation

`mxl_choice` computes one respondent's probabilities and draw-level mean, Xm
and scale scores. With respondent-constant Xm it contracts the design and
probability matrices directly. With CT/alternative-specific Xm it processes
one attribute at a time using alternative-by-draw buffers, rather than
constructing alternative-by-attribute-by-draw derivatives. WTP and lognormal
chain rules are applied before integration, including row-dependent cost
coefficients. Xs scores use each alternative's scaled utility, not only the
first alternative in its task.

`mxl_worker_data` keeps one exact-input/pool-checked `parallel.pool.Constant`
for immutable data and raw draws. A changed dataset or pool replaces that
constant. Without a pool the functions run locally and do not create a pool.
Use `clear mxl_worker_data` to reset the current shared cache. Each worker
has its own resident copy; this is not shared physical memory.

MXL and LCMXL generate coefficients on the worker and contract full covariance
scores into a small coefficient-by-coefficient matrix before selecting
Cholesky entries. LCMXL computes all classes for a respondent and integrates
the mixture score without constructing class-wide coefficient matrices on
the client. Underflow-sensitive reductions retain the legacy multiplication
order where needed.

HMXL preserves global latent-variable normalization and its derivatives.
It combines choice scores with the existing measurement likelihood before
integrating covariance scores. Latent-scale scores include both the utility
scale derivative and its structural-equation contribution. This eliminates
the full covariance-parameter-by-draw tensor and the global coefficient
matrix, including the previous mCT four-dimensional coefficient expansion.
The existing NP-by-NRep-by-NVarA mean-score tensor, latent derivative tensors
and measurement tensors still live on the client. A second respondent loop
integrates covariance scores using the joint choice/measurement posterior;
compact dynamic mean-score slices are transferred to workers again. No
claim is made that all HMXL evaluation transfers have been eliminated.

No GPU backend or respondent blocking is introduced. The
[GPU assessment](GPU_ASSESSMENT.md) remains a separate feasibility report.

## Validation method

`tests/test_mxl_extended.m` defines 92 small deterministic cases (28 MXL,
35 HMXL and 29 LCMXL). It covers preference/WTP with zero, one and two costs,
FullCov 0/1, fixed/lognormal coefficients, Xm/Xs row dependence, missing CTs
and alternatives, empty Xs, HMXL ScaleLV and missing OLS measurements,
flattened/empty HMXL Xm, empty MXL Xm, and a zero-probability latent-class
case. Additional shape checks use two latent variables, three mixture classes,
and NP=1/NRep=1 for MXL/LCMXL. Each case checks value-only against
value/gradient and serial against a three-process pool.

The extension baseline is `6bf56b8`, which retains the original likelihood
for the newly optimized paths. It is exported outside the frozen replication
package. All candidate analytic gradients must agree with central finite
differences within `5e-6` relative error. Where a baseline gradient itself
passes that check, candidate/baseline agreement must also be within `1e-6`.
Values are compared to the baseline within `1e-8` whenever baseline value
evaluation succeeds. Formerly incorrect gradients and reshape exceptions
are recorded as baseline defects, not treated as required behavior. Their
replacement is validated against finite differences.

The original forest/demo HEAD comparisons and full CH estimation are repeated
for the extended MXL implementation. CH uses the published estimates and
must reproduce LL `-4793.8164`; pooled uses NP=8940 and `Weq`. Existing
repository tests and the earlier covariate matrix are also rerun. Numerical
tests use the current default MATLAB R2026b: the small suite ran on D4 and
the core/covariate/CH suite on D3. Results and baselines are saved
outside `replication_package`; the frozen toolbox remains read-only.

| Check | Final result |
| --- | --- |
| 92 central finite-difference cases | Passed; maximum relative gradient error 7.95e-11 |
| 184 serial/three-worker extended comparisons | Passed; serial/parallel value and gradient differences 0; value-only/gradient-branch difference 0 |
| 12 forest/demo HEAD comparisons | Passed; weighted LL difference 0; maximum relative weighted-gradient error 1.64e-12 |
| 102 covariate comparisons | Passed; no candidate exceptions; maximum baseline-comparable gradient error 3.26e-16 |
| Full CH estimation, 3 workers on D3 | Passed; LL -4793.81636367476; relative LL error 1.80e-11; relative parameter error 3.96e-7 |
| Existing repository test, native counter self-check, MATLAB source parsing | test_mdcev_mmdcev passed; counter self-check passed; checkcode found no syntax errors |

All 92 new gradients pass finite differences. The extension baseline has
successful value evaluations for 84 cases and value exceptions for eight;
the maximum comparable relative value error is `2.07e-16`. It has successful
analytic gradients for 33 cases and gradient exceptions for 59; the maximum
comparable relative gradient error is `4.75e-16`. Each of the 33 successful
baseline gradients also passes finite differences. The exceptions include
reshape/indexing failures and an undefined intermediate derivative, not
intentional restrictions preserved by the new code. Counts are per unique
case, not doubled for serial/parallel runs.

| Model | Unique cases | Maximum relative FD error | Baseline value successes/exceptions | Baseline gradient successes/exceptions |
| --- | ---: | ---: | ---: | ---: |
| MXL | 28 | 3.94e-11 | 27 / 1 | 15 / 13 |
| HMXL | 35 | 7.95e-11 | 35 / 0 | 11 / 24 |
| LCMXL | 29 | 3.36e-11 | 22 / 7 | 7 / 22 |

The 102-case covariate matrix now succeeds for all candidate gradients,
including the 32 old reshape exceptions. Those 32 cases pass central finite
differences with maximum relative error `3.03e-11`. All value comparisons
pass with maximum relative error `1.92e-16`.

| Core case | NP | HEAD and extended LL | Maximum relative weighted-gradient error |
| --- | ---: | ---: | ---: |
| Demo preference, diagonal covariance | 789 | -6928.16606386632 | 2.22e-14 |
| Demo preference, full covariance | 789 | -6842.47599520024 | 4.00e-14 |
| Demo WTP, diagonal covariance | 789 | -7051.30135955028 | 3.01e-15 |
| Demo WTP, full covariance | 789 | -7024.18210528086 | 7.63e-16 |
| Forest CH, WTP/full covariance | 644 | -4793.81636376126 | 1.64e-12 |
| Forest pooled, WTP/full covariance, Weq | 8940 | -80435.0119931802 | 2.57e-15 |

The core comparison's maximum absolute per-person gradient difference is
`1.07e-12`; its relative maximum is `7.72e-14`. Recorded final local benchmark
arrays are also compared in full against both actual HEAD and HEAD+parfor.
Their maximum absolute per-person value difference is `3.55e-15`, and their
maximum absolute gradient difference is `1.07e-12`.

The full CH estimation starts from published estimates and follows the
replication trust-region/user-Hessian and quasi-Newton stages. On D3 it took
76.029 seconds with three workers; optimizer time on that machine is not a
local D7 performance benchmark. The result differs from published LL by
`8.65e-8` and rounds to the required `-4793.8164`. This reproduces the
published optimum, not a new multi-start or global-optimum claim.

## Performance

Final MXL CH and pooled value-plus-gradient series use the local Windows 10
computer, i9-13900KS, 192 GB RAM and three process workers, with explicit
MATLAB R2026a. They use the same inputs and native per-PID measurement method
as the initial report: two untimed warm-ups, at least five calls and at least
20 seconds, median evaluation time, sampled simultaneous memory peaks, and
cumulative page-fault deltas sampled at 1 Hz. Benchmarks run only while other
local computations are idle.

The original `dd6f704` gradient runs on the client; its pool workers are idle.
The separate `HEAD+parfor` control changes only two original gradient loops
to parallel loops. Active-worker reduction is evaluated against that clearly
labelled diagnostic, not against idle HEAD workers or the frozen NEE toolbox.
All page faults are counted, not demand-zero faults in isolation. Native
memory peaks include process startup and warm-up; sampled peaks do not.

The unchanged HEAD and HEAD+parfor reference series below are reused from
the initial measurements on the same computer, release, three-worker pool
configuration and saved fixture inputs. Only the extended candidate is
remeasured in this section. MiB means `2^20` bytes.

| Case | Implementation | Calls | Median evaluation (s) | Sampled total peak private (MiB) | Sampled total peak working set (MiB) |
| --- | --- | ---: | ---: | ---: | ---: |
| CH, NP=644 | HEAD | 8 | 2.5066 | 4072.3 | 3952.9 |
| CH, NP=644 | HEAD+parfor diagnostic | 11 | 1.8615 | 6735.4 | 6587.7 |
| CH, NP=644 | Extended, final series | 123 | 0.1618 | 4267.3 | 4132.0 |
| Pooled, NP=8940 | HEAD | 5 | 36.8380 | 9738.7 | 9589.3 |
| Pooled, NP=8940 | HEAD+parfor diagnostic | 5 | 24.2897 | 44618.7 | 42759.4 |
| Pooled, NP=8940 | Extended, final series | 5 | 4.3353 | 8836.3 | 8692.9 |
| Pooled, NP=8940 | Extended, control repetition | 5 | 4.3417 | 8845.4 | 8701.0 |

For CH, the extended median is 15.49 times shorter than HEAD and 11.50 times
shorter than HEAD+parfor. Sampled total private memory is 36.6% lower than
the parallel diagnostic, but 4.8% higher than serial-gradient HEAD because
the extended implementation retains static inputs on active workers.
The final benchmark returned LL `-4793.81636376126`.

For pooled, the final-series median is 8.50 times shorter than HEAD and 5.60
times shorter than HEAD+parfor. Sampled total private memory is 9.3% lower
than HEAD and 80.2% lower than the parallel diagnostic. The final benchmark
returned LL `-80435.01199318023`, identical to both reference series.

The previous candidate series, before the final empty-Xm API guard, measured
0.1627 s for CH and 2.1882 s for pooled. The final pooled series recorded
4.303-4.366 s per call. Both results are retained rather than selecting the
faster series. A fresh control repetition of the final code measured a median
of 4.3417 s, confirming a final-series median range of 4.335-4.342 s. The
reason for the difference from the earlier provisional 2.1882 s is unknown;
it is not attributed to an unverified code or background-process cause.

GoodSync was observed using approximately one core after the first final
series, creating background-activity uncertainty for that run. Before the
control repetition a five-second idle check passed and document writes were
held. GoodSync then consumed 1.28125 CPU seconds across the 57.986-second
controller interval, including preparation/warm-up (less than 0.1% of whole
machine CPU capacity). No other MATLAB calculation was running locally.
This does not establish a causal explanation for the slower timing. The
report uses the reproducible final range, not the faster provisional result.
The memory/faults-per-call comparison below does not depend on selecting the
faster wall-clock result.

The following are per-process means over the steady-state window. Peak fault
rates are the largest sampled rates, not instantaneous peaks. Memory columns
are private commit in MiB. Simultaneous total peaks above are not sums of
independent per-process peaks.

| Case | Implementation | Process | Mean faults/s | Peak sampled faults/s | Sampled peak private (MiB) | Native lifetime peak private (MiB) |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| CH | HEAD | Client | 837440 | 869667 | 1437.1 | 1501.6 |
| CH | HEAD | Worker 1 | 8.4 | 163.4 | 891.8 | 902.8 |
| CH | HEAD | Worker 2 | 7.5 | 136.7 | 905.0 | 913.4 |
| CH | HEAD | Worker 3 | 8.1 | 157.5 | 838.4 | 847.3 |
| CH | HEAD+parfor | Client | 326228 | 624329 | 2522.4 | 2730.4 |
| CH | HEAD+parfor | Worker 1 | 489941 | 625854 | 1416.8 | 1557.6 |
| CH | HEAD+parfor | Worker 2 | 485020 | 611609 | 1410.4 | 1557.9 |
| CH | HEAD+parfor | Worker 3 | 485818 | 618171 | 1414.8 | 1559.7 |
| CH | Extended | Client | 147.3 | 932.7 | 1184.9 | 1328.5 |
| CH | Extended | Worker 1 | 340.4 | 6863.8 | 1021.7 | 1167.2 |
| CH | Extended | Worker 2 | 339.3 | 6847.3 | 1031.9 | 1115.9 |
| CH | Extended | Worker 3 | 356.2 | 7211.3 | 1033.8 | 1117.2 |
| Pooled | HEAD | Client | 789874 | 1265213 | 7049.3 | 7049.4 |
| Pooled | HEAD | Worker 1 | 12.8 | 2167.9 | 898.9 | 905.0 |
| Pooled | HEAD | Worker 2 | 12.3 | 2094.4 | 897.1 | 902.8 |
| Pooled | HEAD | Worker 3 | 15.2 | 2179.8 | 893.5 | 893.5 |
| Pooled | HEAD+parfor | Client | 354563 | 1281173 | 21012.8 | 21203.2 |
| Pooled | HEAD+parfor | Worker 1 | 527420 | 1424904 | 9749.2 | 9749.3 |
| Pooled | HEAD+parfor | Worker 2 | 527401 | 1067428 | 8074.7 | 8245.9 |
| Pooled | HEAD+parfor | Worker 3 | 527289 | 1211891 | 9750.2 | 9750.4 |
| Pooled | Extended | Client | 5586.7 | 18710.2 | 2404.9 | 4768.0 |
| Pooled | Extended | Worker 1 | 1941.1 | 6904.8 | 2141.6 | 4522.3 |
| Pooled | Extended | Worker 2 | 1941.7 | 6908.6 | 2139.1 | 4465.5 |
| Pooled | Extended | Worker 3 | 1942.9 | 6909.6 | 2152.6 | 4537.1 |
| Pooled | Control repetition | Client | 5608.0 | 18706.6 | 2402.4 | 4766.2 |
| Pooled | Control repetition | Worker 1 | 1944.0 | 6911.3 | 2154.1 | 4533.5 |
| Pooled | Control repetition | Worker 2 | 1945.2 | 6919.3 | 2146.6 | 4533.0 |
| Pooled | Control repetition | Worker 3 | 1946.6 | 6912.4 | 2142.6 | 4524.8 |

The three CH active-worker mean fault rates fall by 1439, 1429 and 1364 times
against HEAD+parfor, respectively. Thus every active CH worker exceeds the
requested fivefold reduction without an evaluation-time increase. Worker
faults per completed evaluation are 58.5, 58.3 and 61.2; these divide total
sampled faults by calls and are not direct allocation counts. The pooled
active-worker rates fall by 271.7, 271.6 and 271.4 times; faults per completed
evaluation are 8622.0, 8622.8 and 8628.2. Every active worker in both final
MXL series exceeds the fivefold reduction target. Higher native lifetime
peaks include pool startup and initial constant transfer and are not
steady-state memory requirements.

The pooled control repetition has approximately 1944-1947 faults/s per
worker and 8608-8618 faults/evaluation, consistent with the final series.
The initial HEAD+parfor control has approximately 12.86 million worker
faults/evaluation, so the pooled fault-count reduction is roughly 1490 times
per call as well as more than 270 times per second. These counts include all
page-fault types and a short handshake interval; they do not isolate
demand-zero allocations.

There are no new large-case HMXL or LCMXL performance measurements. Their
allocation/transfer changes are described above, but no measured speedup or
fivefold worker-fault reduction is claimed for those entries. Initial MXL
timings belong to `6bf56b8` and are retained unchanged in the
[historical report](MEMORY_REPORT.md).

## Results and reproduction

Received cluster validation outputs are under
`C:\Users\miq\Documents\_dce_memory\extended_cluster_20261002\output`:

- `extended_final/extended.csv`, `.mat` and `fixtures.mat`: 184 comparisons,
  full baseline/candidate outputs and finite-difference gradients.
- `core_final/correctness.csv` and `.mat`, plus both supplemental CSV/MAT
  files: complete requested demo/forest HEAD comparisons.
- `covariates_final/covariates.csv` and `.mat`: 102 comparisons, including
  baseline exceptions and new finite-difference checks.
- `CH_estimation_final/CH_estimation.csv` and `.mat`: published and final
  estimates, options, and optimizer results.

Final local performance outputs are under
`C:\Users\miq\Documents\_dce_memory\performance\CH_new_extended_final_A` and
`C:\Users\miq\Documents\_dce_memory\performance\pooled_new_extended_final_A`.
`pooled_new_extended_confirm_A` contains the final-code control repetition.
Each series contains `summary.json`, evaluation results, per-PID samples,
process summaries, MATLAB log and result MAT files. The full-array reference
comparison is `performance/extended_comparison.csv`.
The earlier `CH_new_extended_A` and `pooled_new_extended_A` series retain the
pre-guard candidate results for timing-variation context.
All reported completion states require received, nonempty result files and
successful numerical assertions, not just process termination.

```matlab
addpath('C:\Users\miq\Documents\MATLAB\GitHub\DCE\tests');
out = 'C:\Users\miq\Documents\_dce_memory\extended_recheck';
test_mxl_extended(out,[0 3],'6bf56b8');
```

Local performance reproduction uses `tests/measure_mxl_memory.ps1` in a
dedicated otherwise-idle session, with an explicit
`-Matlab 'C:\Program Files\MATLAB\R2026a\bin\matlab.exe'` and a distinct
output directory for each series. Do not write test logs or results into
`replication_package`.
