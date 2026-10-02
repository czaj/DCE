# GPU assessment for DCE/MXL

Date: 2026-10-02, before the D7 driver update. This historical feasibility
assessment used read-only hardware inventory through the shared mcluster
module and one short local `gpuDevice`
preflight. Full `validateGPU` diagnostics, GPU benchmarks, GPU implementation,
driver updates and cluster configuration changes were not performed.
CPU optimization was the recommended first step. The subsequent driver update,
validated D7 runtime and measured GPU prototype are documented in the
[GPU test report](GPU_TEST_REPORT.md); the old-driver failure below is historical.

## Available hardware

VRAM and driver versions were read with `nvidia-smi`; device names were also
checked with `Win32_VideoController`. Eligibility below concerns compute
capability (CC), not a validated driver/runtime combination.

| Node | GPU | VRAM (MiB) | Driver | CC | MATLAB target | CC eligibility |
|---|---|---:|---|---:|---|---|
| D7 | GeForce RTX 5060 | 8151 | 577.00 | 12.0 | R2026b local default | Yes by CC; runtime unavailable (old driver) |
| D4 | GeForce GTX 760 | 4096 | 456.71 | 3.0 | R2026b | No in R2026a or R2026b |
| D3 | GeForce GTX 1060 6GB | 6144 | 560.94 | 6.1 | R2026b | No; eligible in R2026a |
| D2 | GeForce GTX 1070 | 8192 | 566.03 | 6.1 | R2026b | No; eligible in R2026a |
| D6 | GeForce RTX 2070 SUPER | 8192 | 591.86 | 7.5 | R2026b | Yes |
| D5 | GeForce RTX 2070 SUPER | 8192 | 591.86 | 7.5 | R2026b | Yes |

D4, D3, D2, D6 and D5 are enabled in `nodes.json`. D7 remains disabled for
cluster scheduling and has R2026b configured there. R2026b is also the current
local default. The CPU before/after benchmark series started under explicitly
selected R2026a; completed entries remain R2026a results, not R2026b performance
evidence. D7 had 5571 MiB of free VRAM at the inventory snapshot. Free VRAM
changes with display/application use.

MathWorks supports CC 5.0-12.x in
[R2026a](https://www.mathworks.com/help/releases/R2026a/parallel-computing/gpu-computing-requirements.html)
and CC 7.5-12.x in
[R2026b](https://www.mathworks.com/help/parallel-computing/gpu-computing-requirements.html).
The device generations are documented by NVIDIA for
[current](https://developer.nvidia.com/cuda/gpus) and
[legacy](https://developer.nvidia.com/cuda/gpus/legacy) GPUs.

MathWorks lists CUDA 13.1 for
[R2026b](https://www.mathworks.com/help/parallel-computing/run-mex-functions-containing-cuda-code.html),
while NVIDIA lists a minimum driver branch of 580 for
[CUDA 13.x minor-version compatibility](https://docs.nvidia.com/cuda/archive/13.1.0/cuda-toolkit-release-notes/index.html#cuda-driver).
MathWorks recommends a current NVIDIA driver and
[`validateGPU` for diagnosis](https://www.mathworks.com/help/parallel-computing/gpu-computing-requirements.html).
Architecture eligibility alone must not be reported as runtime validation.

## Local R2026b preflight before the driver update

After CPU validation, `gpuDevice` was attempted on D7 using MATLAB
`26.2.0.3386108 (R2026b)`. It failed with
`parallel:gpu:device:OldDriver`:

> GPU computing in MATLAB requires a newer graphics driver. Download and
> install the latest graphics driver for your GPU from NVIDIA.

GPU computing on D7 was unavailable under R2026b with driver 577.00.
A newer NVIDIA driver was needed before a GPU prototype could be tested. This
short initialization check is not a full `validateGPU` run or a GPU benchmark.
Its saved record is `C:\Users\miq\Documents\_dce_memory\gpu_probe_R2026b.mat`;
the accompanying log is `C:\Users\miq\Documents\_dce_memory\validate_R2026b.log`.
No driver was updated and no GPU speedup was measured in that assessment.
The later test used driver 616.92 and passed `validateGPU`.

## Limits relevant to this workload

Consumer GPU arithmetic is a constraint for the required double precision:
NVIDIA lists native FP64 instruction throughput at 1/64 of FP32 for CC 12.0
and 1/32 for CC 7.5 and 6.1. These ratios do not predict MXL elapsed time;
matrix sizes, memory traffic, exponentials and kernel launch overhead also
matter. See the [NVIDIA arithmetic throughput table](https://docs.nvidia.com/cuda/archive/12.9.1/cuda-c-programming-guide/index.html#arithmetic-instructions).
Single precision is not a substitute for the requested LL/gradient tolerances
without separate numerical evidence.

For the pooled case, `17*1000*8940*8 = 1215840000` bytes of static draws is
1.216 GB (1.132 GiB). It can fit once on an 8 GB GPU. However, one full old
expanded buffer, `36*16*1000*8940*8 = 41195520000` bytes, needs 41.2 GB
(38.37 GiB). Processing all respondents together would exceed these GPUs.
Small respondent blocks and contractions that avoid this expansion remain
necessary. Copying draws or coefficients over PCIe on every evaluation would
also undermine the expected benefit.

## Boundary of a useful future experiment

The initial optimized CPU pooled case (`6bf56b8`, `NP=8940`, three workers)
evaluated value plus gradient in a median 2.04 seconds under R2026a. The
extended implementation measures about 4.34 seconds in two final series on
the same inputs; see the
[extended report](EXTENDED_MEMORY_REPORT.md). That is the current measured
CPU reference for deciding whether further GPU work pays off. After a driver
update, only a small measured prototype is justified before any broader port.

1. Keep the CPU implementation and API as the reference. Test a small internal
   double-precision GPU helper for the common MXL case before considering a
   broader DCE port or any backend/API addition.
2. Keep constant `err`, `XXa` and `YY` resident on one GPU between evaluations.
   Update only model parameters, process respondent blocks, and use matrix
   products or `pagemtimes` for gradient contractions. Gather final `f,g`
   together; preserve the existing per-respondent outputs and weights.
3. Use one computation per GPU. Three CPU workers sharing one GPU duplicate
   contexts/data and compete for VRAM. Current mcluster CPU slot counts are
   not GPU slot counts. MathWorks recommends
   [a separate GPU per worker](https://www.mathworks.com/help/parallel-computing/parallel.gpu.gpudevice.html).
4. Preserve model transformations explicitly. Choice-task/alternative-varying
   `Xm` and lognormal costs require row-dependent coefficients and chain-rule
   gradients; `Xs` adds scale calculations and their derivatives. Blindly
   applying `gpuArray` to the current `LL_mxl` does not establish correctness
   or efficient execution for these branches.

[`pagemtimes` supports double and GPU arrays](https://www.mathworks.com/help/matlab/ref/pagemtimes.html).
A future comparison must first validate the GPU runtime, then check identical
inputs/parameters against CPU LL, the whole gradient and estimates. Compare
CPU and GPU timings in the same MATLAB release, synchronize GPU work, include
the intended transfer boundary, and run on an otherwise free computer.
MathWorks documents
[GPU timing and transfer costs](https://www.mathworks.com/help/parallel-computing/measure-and-improve-gpu-performance.html).
No GPU speedup or suitability for the complete DCE package is established here.
