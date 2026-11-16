# Simulation draws

Update **DCE and Tools together**. Sobol estimation defaults to `Draws=6`
(scrambled Sobol), `NRep=1024`, `HaltonSkip=0`, `HaltonLeap=0` and
`SobolLayout='respondent'`. Explicit counts such as 1000 or 5000 are honored
without rounding. Halton retains skip 1; non-Sobol rules default to 1000 draws.

Powers of two preserve Sobol net balance when points are not skipped or thinned.
They are unrelated to respondent count. Balance does not establish adequate
integration accuracy: check parameters, gradients, uncertainty and scrambles.

## How many draws are enough

Two things are reported after every simulated estimation, in
`Results.diagnostics.Simulator`, under the Hessian rows of the printed output and
in the output file. They answer different questions and both are cheap.

**Effective draws per respondent.** Respondent i's simulated likelihood is the
average of P_i over draws, so the draws enter it with weights proportional to
P_i(beta_ir). Kish's effective sample size `1/sum(w^2)` says how many of them
that average rests on. In MXL it comes from `mxl_drawESS` (normal and lognormal
parameters, diagonal or full covariance, preference or WTP space; it says so
instead of guessing elsewhere). For quasi-random draws it is conservative: it is
the independent-draw formula and ignores the balance of the net. Measured on the
DEFRA farmers data (263 respondents, 12 choice tasks, 13 attributes, MXL_d, 1000
scrambled Sobol draws, at the estimates of the 14.11.2026 standard-error
factorial):

| model | median | 5th pct. | min | below 10 draws | median SE(BHHH)/SE(H) |
|---|---|---|---|---|---|
| preference space, normal cost | 130 | 39 | 3.4 | 1.1% | 1.03 |
| preference space, lognormal cost | 152 | 48 | 9.6 | 0.4% | 1.09 |
| WTP space, normal cost | 129 | 12 | 1.6 | 3.8% | 1.23 |
| WTP space, lognormal cost | 182 | 5.0 | 1.3 | 10.3% | 4.22 |

A panel mixed logit uses 13-18% of its draws whatever the specification: twelve
choices per respondent make the individual posterior much tighter than the
mixing distribution the draws come from. What separates the specification whose
standard errors failed from the ones that passed is therefore not the median but
the tail - a tenth of the respondents have fewer than ten effective draws there,
and the worst has 1.3, one draw carrying 87% of that respondent's likelihood.
The default warning is at more than 5% of respondents below 10 effective draws.

**The draw-set check.** The respondents' draw blocks are rotated by half the
sample and the gradient is re-evaluated at the reported estimate; the implied
Newton step says how far the estimate would move, in its own standard errors.
The rotation is a permutation of blocks between respondents, so the pooled set
of draws, their number and anything normalised on draw moments are untouched:
only which part of the net integrates which respondent changes, and that should
not move the estimate. It costs one gradient evaluation. The linearised
prediction was checked against the thing it predicts - re-estimating the same
model with four draw sequences:

| model | draws | predicted shift | largest actual | median actual |
|---|---|---|---|---|
| demo MXL_d, preference space | 500 | 0.8 - 1.8 SE | 2.9 SE | 0.26 SE |
| demo MXL_d, preference space | 2000 | 0.2 - 0.6 SE | 0.5 SE | 0.15 SE |
| DEFRA, WTP space, lognormal cost | 1000 | 1.3 - 7.9 SE | 24.0 SE | 1.8 SE |
| DEFRA, WTP space, lognormal cost | 5000 | 0.8 - 1.6 SE | 14.7 SE | 1.1 SE |

The default warning is at 1 standard error for the worst parameter, which
separates these cases; the median over parameters is reported next to it,
because a single badly identified standard deviation can carry the maximum. A
shift below the threshold still matters for an individual parameter. The
prediction tracks the actual movement where the model is sound and understates
it by an order of magnitude where the Hessian it linearises around is itself
corrupted, as in the last row: read it as a lower bound on the problem, never
as a reassuring upper one. The demo model needs 2000 draws for its estimates to
be a property of the data rather than of the draw sequence. The WTP-space model
with a lognormal cost coefficient is not one at any draw count tried - there it
is the estimates, not only their standard errors, that belong to the draw set.

**Importance sampling.** Re-centring each respondent's draws on the mode of that
respondent's posterior and re-weighting them (a Laplace proposal, as in the
multivariate Poisson-lognormal of the Natural England work, where it took the
effective draws from 7 of 1000 to 937 of 2048) does not transfer to mixed logit
as it stands, and is not offered here. Measured at the DEFRA estimates with the
same 1000 draws, against a 200,000-draw benchmark, as the root mean squared
error of log L_i over respondents:

| simulator | WTP space, lognormal | preference space, lognormal |
|---|---|---|
| draws as they come | 0.254 | 0.060 |
| Laplace proposal, inflation 1.15 | 1.575 | 0.058 |
| Laplace proposal, inflation 2 | 1.336 | 0.389 |
| half prior, half Laplace, inflation 1.5 | 0.151 | 0.052 |
| three quarters Laplace, inflation 1.15 | 0.158 | 0.027 |

A plain Laplace proposal is six times worse than doing nothing in WTP space: the
posterior of z is not approximately normal there, because the index is bilinear
in the willingness-to-pay parameters and the scale, and the scale is lognormal
with a log standard deviation near 3, so the proposal misses mass and the
weights concentrate on the few draws that do not. A defensive mixture - part of
the draws kept from the prior, the proposal density evaluated as the mixture, so
the weights are bounded by the inverse of the prior's share - is uniformly
better than doing nothing, by a factor of 1.7 to 2.2 in root mean squared error,
which is worth 3 to 5 times the draws. That is an efficiency gain, not a repair:
the error of the simulated likelihood scales with the inverse square root of the
effective draws, so moving a 24-standard-error draw-set shift below one standard
error needs about 600 times the effective draws, not 5. Where the individual
posterior is far tighter than the prior - counts over many sites, or latent
variables measured by many indicators - the gain is orders of magnitude instead,
and importance sampling is worth the per-respondent optimisation it costs. The
measurements and the prototype are in
`MATLAB/.ai/_work/2026-11-16_importance-sampling`.

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
