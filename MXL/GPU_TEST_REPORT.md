# MXL GPU experiment

Historical prototype report. The subsequent production integration, automatic
selection and new measurements are in the [auto-GPU report](GPU_AUTO_REPORT.md).
The helpers described below now live in `MXL`, not `tests`.

Date: 2026-10-02. CPU reference: master `36626de`. Experiment branch:
`mxl-gpu-test`. MATLAB `26.2.0.3386108 (R2026b)`, Windows 10, D7,
i9-13900KS, 192 GB RAM, NVIDIA GeForce RTX 5060 (8151 MiB reported by
`nvidia-smi`), graphics driver 616.92, compute capability 12.0, WDDM.

## Result

The blocked, double-precision MXL prototype is useful on D7: resident-data
value-plus-gradient evaluations including output gathering are 5.57 times
faster for CH and 5.51 times faster for pooled than the current three-worker
CPU implementation in the same MATLAB release. A synthetic CH specification
with varying Xm/Xs is 3.31 times faster. These are likelihood evaluation
speedups, not complete production-estimator or cluster-throughput speedups.

The production `MXL`, `HMXL`, `LCMXL`, their likelihoods and `EstimOpt` API are
unchanged. The prototype lives under `tests`, not on the estimator's call path.
No files in `replication_package`, unrelated user changes or cluster
configuration were modified. All measurements used otherwise-free D7;
the local sessions ran sequentially. NEE queues were inspected read-only.

## Prototype

`mxl_gpu_prepare(C)` uploads immutable choices, design, raw draws and covariates
once. `LL_mxl_gpu(data,b,blockSize)` uploads the small parameter vector and
computes respondent blocks using native `gpuArray` and `pagemtimes`.
It returns device arrays with the same per-person f/g layout as CPU MXL.
There is no GPU `parfor`, custom CUDA code, dependency or public backend option.

Matrix contractions avoid alternative-by-attribute-by-draw expansions.
For varying Xm, the prototype processes attributes separately and preserves
row-dependent WTP/lognormal chain rules. Xs derivatives use each row's scaled
utility. Missing alternatives are masked, and missing tasks have neutral
panel factors. Full covariance gradients contract the draw scores before
selecting Cholesky entries; rare subnormal reductions preserve CPU order.

The tested scope is valid choice-coded data with finite available utilities:
fixed/normal/lognormal coefficients, preference/WTP space (including two
cost coefficients), FullCov 0/1, constant/CT/alternative-specific Xm and Xs, and
value-only/value-plus-gradient calls. Missing alternatives, tasks and a
whole respondent, singleton dimensions, partial blocks, two Xm/two Xs,
RealMin 0/1 and subnormal covariance reductions are tested. Other
distributions, Johnson/NLT/ExpB, analytic Hessians and HMXL/LCMXL are outside
this prototype. Unsupported transformations are rejected, not silently ignored.

## Runtime and method

`validateGPU` passed every reported platform, driver, library, device,
allocation and kernel-launch test. A double `pagemtimes` smoke test matched
CPU exactly. The device reports `SingleDoubleRatio=64`; despite that hardware
constraint, this workload benefits from batching and native contractions.

All CPU/GPU controls use identical saved inputs, b and Sobol draws. CH has
NP=644; pooled has NP=8940 and Weq weights. Both have 17 coefficient rows,
1000 draws, full covariance and WTP. CPU controls use the current optimized
implementation, not the old allocation-heavy HEAD or historical R2026a times.

Each steady-state series has two untimed f/g warmups, at least five calls and
at least 20 seconds. Times below are medians. GPU work is synchronized before
tic and before toc. `resident device` excludes gathering; `resident + gather`
includes both per-person outputs. `transfer-inclusive` includes preparing and
uploading static data every call, plus evaluation and gathering; it uses three
calls in an already-initialized GPU context, not cold driver initialization.
Pool startup and initial constant transfer are outside steady-state CPU times.

See MathWorks on [GPU timing and transfer costs](https://www.mathworks.com/help/parallel-computing/measure-and-improve-gpu-performance.html)
and [GPU synchronization](https://www.mathworks.com/help/parallel-computing/parallel.gpu.gpudevice.wait.html).

## Timings

Seconds per value-plus-gradient evaluation:

| Specification | CPU serial | CPU, 3 workers | GPU resident + gather, block 128 | GPU with static upload each call | GPU speedup vs CPU 3 |
|---|---:|---:|---:|---:|---:|
| CH, published specification | 0.614986 | 0.321210 | 0.057663 | 0.070673 | 5.57x |
| Pooled, Weq | 8.146001 | 4.301399 | 0.780255 | 0.913494 | 5.51x |
| CH, synthetic varying Xm/Xs | Not measured | 0.990353 | 0.299110 | 0.324719 | 3.31x |

The last row retains CH choices/design/draws but adds two deterministic mean
and two scale covariates, varying by respondent, task and alternative, and
small additional coefficients. It is not the published CH specification and
does not establish pooled-Xm/Xs performance. Its LL is -4802.83158728908.

All measured block choices are retained below; no best-case-only selection:

| Sample | Block | Resident device (s) | Resident + gather (s) | Transfer-inclusive (s) | Steady-state calls, device/gather |
|---|---:|---:|---:|---:|---:|
| CH | 32 | 0.064025 | 0.064665 | 0.085588 | 313 / 310 |
| CH | 128 | 0.057320 | 0.057663 | 0.070673 | 348 / 346 |
| CH | 512 | 0.054963 | 0.055160 | 0.067596 | 363 / 362 |
| Pooled | 32 | 0.850226 | 0.852366 | 1.075574 | 24 / 24 |
| Pooled | 128 | 0.777975 | 0.780255 | 0.913494 | 26 / 26 |
| Pooled | 512 | 0.763826 | 0.767690 | 0.930512 | 27 / 26 |

CPU serial/three-worker call counts are 34/68 for CH and 5/5 for pooled.
Xm/Xs counts are 68/67 for GPU device/gather and 22 for CPU three-worker.
Block 128 is the conservative choice: block 512 improves resident/gather
times by only 4.3% for CH and 1.6% for pooled, while increasing memory use.
There is only one invocation per block, not a confidence interval across
independent repeated runs. CH/pooled static payloads are 86.67/1203.13 MiB;
initial preparation/upload took 0.037-0.039/0.176-0.284 seconds across blocks.

## Memory

The external sampler sleeps 500 ms between reads. The table uses samples
whose UTC timestamps fall inside the timed steady-state windows. Host private
commit is the simultaneous sum across this invocation's MATLAB processes.
Device-wide used VRAM is read by `nvidia-smi`; it includes the desktop, other
display contexts and MATLAB allocator caching. It is not MATLAB-only memory
or an exact transient allocation peak. Host and device memory are separate.

| Sample/mode | Sampled peak host private (MiB) | Sampled peak device used (MiB) |
|---|---:|---:|
| CH CPU serial | 882.6 | Not a GPU calculation |
| CH CPU, 3 workers | 3698.0 | Not a GPU calculation |
| CH GPU block 128, resident + gather | 2525.4 | 1988 |
| CH GPU block 512, resident + gather | 3691.2 | 3123 |
| Pooled CPU serial | 2108.2 | Not a GPU calculation |
| Pooled CPU, 3 workers | 8282.0 | Not a GPU calculation |
| Pooled GPU block 128, resident + gather | 5148.9 | 3355 |
| Pooled GPU block 512, resident + gather | 6745.4 | 4806 |
| CH Xm/Xs GPU block 128, resident + gather | 2794.1 | 2691 |

`gpuDevice.AvailableMemory` snapshots are saved separately as allocator
availability, not peak memory. No GPU page-fault reduction is claimed here;
the earlier CPU per-process measurements remain in the memory reports.

## Numerical checks

Errors are max absolute differences divided by max(1, max absolute reference),
over the complete relevant array/vector. They are not per-component relative
guarantees. Weighted LL and weighted full-gradient sums are also checked.
Thresholds are 1e-8 for values/LL and 1e-6 for gradients.

| Check | Result |
|---|---|
| 34 deterministic fixtures, 102 block comparisons | Passed; max LL relative error 1.95e-16, per-person g 6.38e-16, summed g 1.04e-15 |
| Demo preference diagonal/full, WTP diagonal/full | All passed; max per-person value error 1.97e-16 and g error 1.16e-14; value-only also agrees |
| CH, all three block sizes | LL -4793.81636376126; weighted LL difference 0; max per-person g error 4.56e-14, weighted g 9.30e-13 |
| Pooled, all three block sizes | LL -80435.0119931802; weighted LL difference 0; max per-person g error 2.20e-14, weighted g 1.99e-15 |
| CH synthetic Xm/Xs | Weighted LL difference 0; per-person g error 1.15e-14; weighted g 1.02e-14 |
| Existing test_mdcev_mmdcev; GPU source parsing | Passed on R2026b |

The CH optimizer check uses the GPU objective/analytic gradient, starts at the
published parameters and preserves the saved tolerances with quasi-Newton.
It reproduces LL -4793.81636372669 (relative published-LL error 7.21e-12)
and parameters within 2.30e-7 relative error. Exitflag is 2 (step tolerance),
not a new global-optimum result. The final CPU evaluation at these same
parameters matches LL exactly and weighted gradient within 1.72e-12.
Two successful optimizer checks took 0.808/0.530 seconds; this warm-started
likelihood-only check is not a complete MXL run with Hessians/standard errors.

## Reproduction and next boundary

Saved results, raw per-person arrays, sampler records, logs, demo checks and
the synthetic specification are under
`C:\Users\miq\Documents\_dce_memory\gpu_20261002`.
`correctness_v3` is the final 102-comparison suite;
`CH_gpu_estimation_confirm` is the successful optimizer/demo invocation.
The earlier failed runs retain two harness-only issues, fixed before final
validation: an assert identifier and empty-struct report initialization.
No numerical likelihood defect was hidden by those retries.

```matlab
addpath('C:\Users\miq\Documents\MATLAB\GitHub\DCE\tests');
test_mxl_gpu('C:\Users\miq\Documents\_dce_memory\gpu_recheck');
bench_mxl_gpu('CH','C:\Users\miq\Documents\_dce_memory\gpu_CH_recheck',128,'gpu');
test_mxl_gpu_estimate('C:\Users\miq\Documents\_dce_memory\gpu_estimate_recheck',128,120);
```

The tests expect the saved CPU fixtures from the memory-analysis work; the
correctness helper documents how to regenerate the small extended fixtures.
Run each physical benchmark in a fresh otherwise-idle session and use fresh
output directories. The private `run_gpu_benchmark.ps1` and
`summarize_gpu.ps1` record and align the device/process samples.

This evidence justifies considering an optional GPU path for the tested MXL
subset, retaining the CPU path for unsupported options and without sharing
one GPU across CPU workers. It does not justify silently enabling GPU across
the entire package. HMXL's joint choice/measurement posterior and LCMXL's
class mixture need separate integration and tests. No new timings or runtime
validation were performed on D5/D6 GPUs; the historical hardware assessment
is not a benchmark of those nodes. Single precision was not used or proposed
as a substitute for the required accuracy.
