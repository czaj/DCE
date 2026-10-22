# Hybrid memory and GPU validation — 2026-10-22

Baseline: DCE `30805de902f70dff8e4c2d0ea6d9530c66691769`, Tools
`77879894d366c3cef6cbf399a52a3e51ad51fc11`. Tests use MATLAB R2026b.
Research inputs and existing mining results were copied for validation; the
replication package and running research jobs were not modified.

## CPU performance

Mining reference HMXL: 1,105 respondents, five choice coefficients, five
latent variables, 12 ordered indicators, 90 parameters. D4: Intel i9-7900X,
20 logical processors, 128 GB RAM. Each fresh process uses three process
workers, two computational threads per worker and two on the client.
The old and new implementations run sequentially on the same idle node.

Each series has two untimed warm-ups, at least five measured calls and at
least 20 seconds of measured work. The native Windows sampler reads all
client/worker PIDs every 200 ms. Memory figures are simultaneous sampled
private-commit peaks after warm-up, not sums of unrelated lifetime maxima.
Faults include every page-fault type, not only demand-zero allocations.

| Draws/respondent | Version | Median value + gradient (s) | Total private peak (MiB) | Total working-set peak (MiB) | Total faults/evaluation |
|---:|---|---:|---:|---:|---:|
| 2,000 | Baseline | 6.67024 | 6,768.1 | 6,886.0 | 3,691,066 |
| 2,000 | Blocked CPU | 2.23534 | 4,126.1 | 4,312.2 | 997,668 |
| 10,000 | Baseline | 33.01863 | 19,914.9 | 20,006.9 | 20,909,727 |
| 10,000 | Blocked CPU | 10.31507 | 6,881.9 | 7,038.4 | 5,018,767 |

At 2,000 draws this is 2.98x faster, 39.0% less peak private memory and 73.0%
fewer total faults per evaluation. At 10,000 draws it is 3.20x faster, 65.4%
less peak private memory and 76.0% fewer total faults per evaluation.
More of the work now runs on workers: worker fault rates alone increase
relative to the old partly client-side measurement calculation. The reduction
claimed here is for the whole job, not each worker separately.

Full respondent outputs were compared, not just the sum of log likelihoods:

| Draws | Maximum absolute value difference | Maximum absolute gradient difference | Relative gradient difference |
|---:|---:|---:|---:|
| 2,000 | 5.33e-15 | 1.07e-12 | 1.26e-14 |
| 10,000 | 2.66e-14 | 7.85e-12 | 1.61e-13 |

These are likelihood-evaluation measurements, not full-estimation speedups.

## Numerical validation

- The existing 92-case MXL/HMXL/LCMXL suite passes 184 serial/three-worker
  comparisons, including central finite differences.
- The new 72-case hybrid suite covers seven measurement types, missing items,
  measurement covariates, multiple latent variables, WTP and latent scales,
  FullCov 2, binary categories, one-respondent blocks, NP=1 and NRep=1. It checks all analytical
  derivatives against central differences and block/serial/pool agreement.
- HLC checks preserve global normalization and category layouts in singleton
  blocks, for both OLS and ordered measurements, on CPU and forced GPU.
- Both mining specifications pass all-parameter finite differences using 64
  draws on all 1,105 respondents. Relative errors are 9.74e-7 (90 parameters)
  and 5.99e-7 (130 parameters). Weighted wrapper scores, BHHH, Hessian symmetry
  and inactive-parameter rows/columns are also checked.
- In the 130-parameter mining check, two respondents hit `ProbMin=1e-300`.
  Their old maximum absolute score was 1.89718 despite a locally constant
  objective. The new score is zero and passes finite differences. Other
  gradients agree to 5.25e-15 relatively; values agree to 9.87e-16.
- Public HMXL, HMNL and MIMIC estimators pass preprocessing, optimization and
  `HessEstFix=5` checks. MIMIC is tested with pseudo-random, Latin hypercube
  and scrambled Sobol draws. These tiny fixtures stop after one optimizer
  iteration and test integration, not convergence or identification. The
  existing MDCEV/MM-DCEV regression test also passes.
- On D7's RTX 5060 (driver 617.42, compute capability 12.0), 144 forced-GPU
  hybrid comparisons pass, including value-only/gradient agreement. The
  refactored GPU MXL passes 68 comparisons from its existing 34 fixtures.
  CPU fallback was exercised on D3, whose driver 582.78 cannot run GPU
  computation under this MATLAB release. Drivers were not changed.
- Routing checks cover repeated calls, changed coefficients/data, automatic
  selection, forced GPU, CPU fallback and invalid options. Native memory
  counters pass their existing allocation/wraparound self-check.

## GPU measurement

The same mining reference at 2,000 draws takes a median **0.80825 s** for
value plus gradient on the RTX 5060, over 25 measured calls after warm-up.
The benchmark asserts that the GPU path was used. Against the D4 CPU
result, maximum absolute differences are 3.55e-15 for values and 6.61e-13
for gradients (relative gradient error 7.76e-15).

This is a standalone GPU measurement, not a controlled CPU/GPU speed ratio:
D7 had a separate CPU research job, and its hardware differs from D4.
The Windows process-memory sampler does not measure GPU device allocations.

## Reproduction

Add both repositories and their test directories to the path. Run in a
dedicated MATLAB session:

```matlab
test_hybrid_memory('output/hybrid',3);
test_hlc_memory;
test_hybrid_auto;
test_hybrid_estimators('output/estimators');
test_hybrid_gpu('output/gpu','output/hybrid/hybrid.mat');
test_hlc_memory(true); % GPU-equipped machine
```

`test_mxl_extended` exports the requested Git baseline or accepts an existing
directory of `baseline_LL_*` aliases. `test_mxl_gpu` uses its saved fixtures.
`test_hybrid_mining(fitFile,outDir)` uses a saved mining fit containing `fit`.
`bench_hybrid_memory(fitFile,outDir,mode,workers,nrep,repeats)` compares
`legacy`, `cpu` and `gpu` in separate fresh processes; its native sampler
attaches to the requesting MATLAB PID. A forced-GPU benchmark asserts that
GPU was actually used. Each benchmark needs its own new output directory.

Local artifacts, inputs and archived baselines are retained outside the
research package in `C:\Users\miq\Documents\_dce_memory\hybrid_20261022`.
Cluster results are under `cluster/output`, including `extended02`,
`tests_release`, `estimators05`, `mining_reference`, `mining_adjusted_final`
and the four `bench_*` directories. GPU validation is under `gpu_release`
and `mxl_gpu_final`, with the mining measurement under `gpu_bench2k_02`.
Research data and raw benchmark artifacts are not committed.
