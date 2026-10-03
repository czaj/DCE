# Automatic GPU selection in MXL

Date: 2026-10-03. CPU reference: `36626de`; GPU prototype: `a7ddc13`.
This integrates the tested native `gpuArray`/`pagemtimes` kernel into `LL_mxl`.
HMXL and LCMXL retain their previously optimized CPU paths.

## Selection And Fallback

No option is required: `EstimOpt.GPU = 'auto'` is the default. Optional `cpu`
disables GPU; `gpu` skips the speed comparison but still falls back safely for
unsupported models, unavailable hardware or GPU resource failures. Logical
false/true are aliases for cpu/gpu. Function inputs and CPU-double outputs do
not change.

The client checks the Parallel Computing Toolbox license, native GPU support,
the selected/default device's double-precision capability and available memory.
There is no GPU model-name whitelist or extra CUDA dependency. Immutable choices,
design, draws and covariates are uploaded once. A conservative memory budget
selects blocks of at most 128 respondents; insufficient capacity uses CPU.

Auto warms both paths, times synchronized GPU computation including gathering,
and retains GPU only when the per-person values/gradients agree within 1e-8/1e-6
and GPU takes less than 90% of the CPU time. This is a one-time calibration per
dataset, existing pool and output mode, not a promise about a single cold call.
The benchmark first-call costs below include preparation and calibration.

Parameter updates reuse resident inputs. Exact changes to any input/options,
including same-shaped data changes, or pool replacement invalidate the cache.
Value-only and value-plus-gradient decisions are separate. GPU work runs on the
client, never individually on CPU parfor workers. The device is not reset.

GPU driver/device/allocation exceptions disable that dataset's GPU path and
recompute on CPU. A device switch invalidates resident inputs. An external reset
of the same device conservatively causes CPU fallback rather than stale results;
`clear mxl_gpu_auto` permits a fresh attempt. An unavailable device detected at
initialization is likewise retried after clearing. Non-GPU programming errors
are not swallowed.

## Supported Scope

GPU supports nargout 1/2, FullCov 0/1, normal/lognormal/fixed coefficients,
preference/WTP space, respondent-specific and task/alternative-specific Xm/Xs,
and missing alternatives, tasks or entire respondents. Missing rows are masked
exactly as in the CPU kernel. Legacy zero/two-chosen coding and nonfinite
available utilities preserve CPU output semantics rather than being repaired.

Three-output analytical-Hessian requests, other distributions, nonlinear/Johnson
transforms, ExpB, nonstandard Cholesky derivative order and unsupported input
types retain CPU. Standard FullCov value-only requests do not require derivative
indices. Outer `LL_mxl_MATlike` weight, penalty and RealMin repairs remain unchanged.

## New Local Measurements

Otherwise-free D7, R2026b `26.2.0.3386108`, RTX 5060, driver 616.92.
Saved identical CH/pooled inputs: 1000 Sobol draws, 17 coefficients, full
covariance, WTP, pooled Weq weights. Each timed series has at least five calls
and 20 seconds; medians exclude startup/calibration. An existing three-worker
pool remains open in both compared paths. Auto returns gathered CPU arrays.

| Model | NP | CPU, 3 workers (s) | Default auto (s) | Selected | Speedup | First auto call (s) |
|---|---:|---:|---:|---|---:|---:|
| CH | 644 | 0.164738 | 0.059138 | GPU | 2.79x | 1.052 |
| Pooled, Weq | 8940 | 2.146910 | 0.800524 | GPU | 2.68x | 6.714 |

These new paired production timings supersede prototype timing comparisons;
CPU medians differ from the historical measurements. The historical report
retains its sampled RAM/VRAM table. No new production peak-memory or page-fault
claim is inferred from timings or allocator availability.

## Validation

| Check | Result |
|---|---|
| Direct GPU kernel: 34 fixtures, blocks 1/3/5 | All 102 comparisons passed |
| Public routing/cache/fallback: serial, 3-worker pool, worker call | All 76 comparisons passed; 55 GPU results |
| dd6f704 vs public default: four demo variants, CH, pooled; serial/3 workers | All 12 comparisons passed |
| Xm/Xs baseline and numerical-gradient regressions | All 102 comparisons passed |
| Existing MXL/HMXL/LCMXL extended regressions | All 184 comparisons passed |
| Existing test_mdcev_mmdcev | Passed |
| Simulated unavailable-GPU discovery, serial/3 workers | All 76 comparisons passed; every result used CPU |
| Missing choices/tasks/respondents, multiple Xm/Xs, WTP/covariance variants, subnormal probabilities | Passed |
| Changed b and same-shaped Y/X/draws/Xm/Xs/options | Passed |
| ExpB/spike/custom derivative order and Hessian CPU bypass | Passed |
| Legacy choice coding; Inf/NaN masks | Passed |
| Dedicated-session GPU-context reset, then clear dispatcher | GPU -> CPU -> GPU; fallback exactly matches CPU |
| CH published-point weighted LL | -4793.81636376126 |
| Pooled weighted LL | -80435.0119931802 |
| Full production CH MXL, trust-region/BHHH then quasi-newton | LL -4793.81636367476; 16.807 s |

Direct-kernel maximum scaled per-person value/gradient differences are
1.56e-16/6.38e-16; summed gradient difference is 1.04e-15. Routing tests reach
1.31e-16/3.18e-16. The larger CH/pooled checks reach 4.79e-16/4.56e-14;
the largest weighted-gradient difference is 9.30e-13.
Against the original dd6f704 implementation, the maximum relative LL,
per-person gradient and weighted-gradient differences are
2.66e-16, 7.72e-14 and 2.01e-12, respectively.
Full CH relative LL/parameter differences from publication are
1.80e-11/3.96e-7, within the requested tolerances. This is the complete public
MXL workflow, not only the prototype's bounded optimizer check.

Maximum finite-difference errors in the covariate/extended suites are
3.61e-11/7.98e-11. Both private cluster validation queues were stopped while all
four tasks were still pending because every enabled slot was occupied by peers.
No remote test started. The remaining short regressions ran sequentially on
idle D7 instead. GPU absence was simulated by shadowing only `gpuDeviceCount`
in a dedicated session, not by changing drivers, hardware or production code.
Actual unavailable-GPU cluster execution, deliberate GPU OOM and system-driver
failure injection were not performed; the real context-reset fallback was tested.

## Reproduction

Use a dedicated, otherwise-idle MATLAB session; all outputs must remain outside
`replication_package`. From the repository, add `tests` to the MATLAB path:

```matlab
test_mxl_gpu(outDir,fixtureFile,[1 3 5]);
test_mxl_gpu_auto(outDir,fixtureFile,3);
bench_mxl_gpu('CH',outDir,128,'auto');
bench_mxl_gpu('pooled',outDir,128,'auto');
test_mxl_memory(outDir,0,'estimate',privateReferenceCopy);
```

Each call needs a distinct output directory. Fixtures come from
`test_mxl_extended`; CH/pooled benchmarks use the saved identical inputs from
the memory-analysis harness. Raw CSV/MAT/JSON/log evidence and the tested-source
hash manifest are outside the repository under
`Documents/_dce_memory/gpu_auto_20261003`. No frozen paper files, unrelated user
changes or shared cluster preferences were modified.
The supplemental context-reset check and its CSV/MAT/log are in the sibling
`gpu_auto_local_20261003` directory; no system-driver reset was used.
