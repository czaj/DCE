# Simulation draws

Update **DCE and Tools together**. Sobol estimation defaults to `Draws=6`
(scrambled Sobol), `NRep=1024`, `HaltonSkip=0`, `HaltonLeap=0` and
`SobolLayout='respondent'`. Explicit counts such as 1000 or 5000 are honored
without rounding. Halton retains skip 1; non-Sobol rules default to 1000 draws.

Powers of two preserve Sobol net balance when points are not skipped or thinned.
They are unrelated to respondent count. Balance does not establish adequate
integration accuracy: check parameters, gradients, uncertainty and scrambles.

## Layout and reproducibility

Tools' `dceQuasiDraws` supplies QMC estimators and blocked demand simulations.
For scrambled Sobol, respondent position i uses substream i of a seeded
`mrg32k3a` stream to construct its own MATLAB `MatousekAffineOwen` scrambling.
Its first R points stay unchanged when R increases, NP increases at the end,
or evaluation is split into blocks. The caller's random stream is restored.
Keep seed, dimensions, skip/leap and respondent ordering fixed for paired
comparisons. These are positional indices, not respondent identifiers;
reordering/subsetting data can change assignments. Explicitly unscrambled
Sobol (`Draws=5`) uses a common deterministic net across respondents.

Output retains R consecutive rows per respondent and the existing dimension
order. Exact uniform endpoints become `eps/2` and `1-eps/2` before inverse-CDF
transforms; interior points are unchanged. This keeps the first Sobol point
instead of dropping it to avoid infinity. The finite endpoint convention
particularly matters for unscrambled Sobol; scrambling is the default.

The new grid changes the simulated objective even at the same seed and R.
To reproduce an old estimation grid, supply its settings explicitly:

```matlab
E.Draws = 6;
E.NRep = 1000;             % actual old count
E.HaltonSkip = 1;          % actual old skip
E.HaltonLeap = 0;
E.SobolLayout = 'global';  % one scramble, consecutive respondent blocks
E.Seed1 = 179424673;       % actual old seed
```

Global blocks are not nested per respondent as R changes. Saved MMDCEV fits
without a layout field retain this historical taste-draw layout if their draws
must be reconstructed for demand. Keep original draws and package versions
for full historical replication.

MDCEV demand defaults to 512 simulations. MMDCEV defaults to the largest power
of two not exceeding `min(512,NRep)`; explicit NSim is unchanged. Demand errors
use observation/task positions and an independent seed: `Seed2` modulo 2^32,
or `Seed1+104729` modulo 2^32 when Seed2 is absent. This avoids reusing taste
scrambles. Simulations remain blocked in memory; the new demand randomization
changes predictions relative to the earlier global-RNG implementation.

Hybrid latent normalization still depends on sample/draw moments despite
nested draws. Use `hybrid_rebase_start` where applicable when moving a fitted
point to a new grid. Nesting removes neither this mapping nor sensitivity checks.

Tests: `Tools/tests/test_dceQuasiDraws` checks defaults, balance, prefix/block
invariance, seeds, legacy equality and endpoints. `tests/test_sobol_draws`
checks likelihood gradients, public hybrid defaults and demand budgets/blocks.

Validated on MATLAB R2026b (2026-10-27): the generator and new draw tests above,
plus `test_hybrid_estimators`, `test_hybrid_public_options`,
`test_hlc_cutoff_output` and `test_vce_estimators` passed. These small synthetic
fixtures check numerical paths and option handling, not statistical identification
or sufficient simulation accuracy for an empirical application.

References: [MathWorks sobolset](https://www.mathworks.com/help/stats/sobolset.html),
[Owen, On dropping the first Sobol' point](https://arxiv.org/abs/2008.08051).
