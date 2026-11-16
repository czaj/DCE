# Simulation draws

Update **DCE and Tools together**. Sobol estimation defaults to `Draws=6`
(scrambled Sobol), `NRep=1024`, `HaltonSkip=0`, `HaltonLeap=0` and
`SobolLayout='respondent'`. Explicit counts such as 1000 or 5000 are honored
without rounding. Halton retains skip 1; non-Sobol rules default to 1000 draws.

Powers of two preserve Sobol net balance when points are not skipped or thinned.
They are unrelated to respondent count. Balance does not establish adequate
integration accuracy: check parameters, gradients, uncertainty and scrambles.

## How many draws are enough

This continues Czajkowski and Budzinski (2019), "Simulation error in maximum
likelihood estimation of discrete choice models", Journal of Choice Modelling 31
(https://doi.org/10.1016/j.jocm.2019.04.003), which measured the same error from
the outside: an MXL data generating process over 27 experimental conditions, 100
to 1,000,000 draws of four kinds, each model re-estimated up to 1000 times. It
called a number of draws sufficient when a different draw set of the same kind
would not move the log-likelihood by half the critical value of a
one-degree-of-freedom likelihood ratio test (1.9207) and would keep the
estimates within 5% of their true values, and its answer was more than 2000
scrambled Sobol draws with five attributes and more than 20,000 with ten - far
above what applications were using, and above this toolbox's own default of
1024 for anything with many attributes. Sobol won against Halton, MLHS and
pseudo-random draws throughout, which is why it is the default here.

The two diagnostics below do not improve on that experiment; they ask a cheaper
question. The paper needs a known data generating process and a thousand
re-estimations to say how many draws a *kind* of model needs. These say whether
*this* model, on this data, at this number of draws, is reporting a maximum of
the likelihood or a maximum of its draw set - for the price of one gradient, with
no true parameters to compare against. That also lets them cover the
specifications the paper closed by listing as untested: correlated parameters,
WTP space, non-normal distributions, latent class and hybrid choice models, and
distributions wide relative to their means. All of those are measured below, and
all of them need more draws than the paper's settings did. The thresholds used
here were calibrated against the paper's own method, re-estimation with several
draw sequences.

Two things are reported after every simulated estimation, in
`Results.diagnostics.Simulator`, under the Hessian rows of the printed output and
in the output file. They answer different questions and both are cheap. A third
row, `What to do`, says which remedy they point to: the number of draws the
square-root law implies for a shift below the threshold (`DrawsSuggested`,
rounded up to a power of two, and a lower bound, since the shift it comes from
is one too); or, when that number is past what machine time buys, that the
specification rather than the simulator is the thing to change; or, only when
`ISCeiling` = draws over effective draws reaches the order of a hundred, that a
better simulator is worth considering.

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

One more check, against the other method entirely: on the toolbox demo MXL_d
with six attributes at 500 draws, the row described below asks for at least 2048
draws, while Czajkowski and Budzinski (2019) ask for more than 2000 scrambled
Sobol draws at five attributes. A per-model gradient and a thousand
re-estimations over artificial data land in the same place.

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

**Hybrid models.** Here the draws are weighted by the choices and by what the
respondent answered to the attitudinal items, so the posterior that carries them
is tighter than in a plain mixed logit - the structure in which importance
sampling transformed the Natural England multivariate Poisson-lognormal, where
the effective draws went from 7 of 1000 with prior draws to 937 of 2048 with a
Laplace proposal. Measured on the toolbox demo survey (789 respondents, 12 choice
tasks, 3 alternatives, 6 attributes, 500 scrambled Sobol draws, one latent
variable, and five 5-point items of one battery as ordered-probit indicators):

| model | draw dimensions | effective draws | below 10 | shift, worst | shift, median | SE(BHHH)/SE(H) |
|---|---|---|---|---|---|---|
| MXL_d, no indicators | 6 | 43 of 500 | 3.9% | 1.80 SE | 0.41 SE | 0.94 |
| HMNL, one indicator | 1 | 125 of 500 | 0% | 0.03 SE | 0.006 SE | 0.93 |
| HMNL, five indicators | 1 | 125 of 500 | 0% | 0.03 SE | 0.001 SE | 1.12 |
| HMNL, five indicators, WTP space | 1 | 75 of 500 | 0.4% | 0.12 SE | 0.04 SE | 1.07 |
| HMXL, five indicators | 7 | 23 of 500 | 16.7% | 1.29 SE | 0.20 SE | 0.99 |
| HMXL, five indicators, WTP space | 7 | 23 of 500 | 18.5% | 3.02 SE | 0.38 SE | 0.98 |
| HMXL, WTP space, lognormal cost | 7 | 23 of 500 | 19.1% | 3.87 SE | 0.33 SE | 0.96 |

The dimension of the integral decides this, not the hybrid structure as such. The
latent variable is one dimension, and 500 quasi-random draws integrate one
dimension so well that rotating the draw blocks moves the worst parameter by 0.03
standard errors; adding four indicators to one changes nothing. There is nothing
in HMNL for importance sampling to repair, whatever its effective draw count
says - which is also the clearest case of that count being conservative for a
balanced net. The random taste parameters are the opposite: the same data as a
plain MXL_d uses 43 of its 500 draws and moves by 1.8 standard errors, and HMXL,
which integrates six tastes and the latent variable together, is down to 23 draws
with a sixth of its respondents below ten, and to a shift of 3 to 4 standard
errors in WTP space. That is the largest ceiling for a better proposal
anywhere in this toolbox - 500 over 23, about 22 times - and in preference space
it is also where a Laplace proposal is well founded: with the scale not driven by
the latent variable, the integrand is log-concave in the draws (a logit
likelihood, ordered-probit indicators and a normal prior), so the mode and
curvature it needs exist and mean what they should. In WTP space that is lost
again, for the same reason as above - the taste coefficients are multiplied by
the cost coefficient. What a proposal actually collects from that ceiling is
measured below, and it is little: 2.4 times the draws in preference space and
nothing at all in WTP space.

The last row also narrows down what breaks the standard errors. A lognormal cost
coefficient in WTP space is not enough by itself: here its log standard deviation
is 0.70 and `SE(BHHH)/SE(H)` stays at 0.96, while the estimates still move by
almost four standard errors between draw sets. It took a log standard deviation
of 2.84 - a mean 56 times the median - for the covariance matrix itself to fail
in the DEFRA model. The two diagnostics therefore catch different things, and the
draw-set check is the one that fires first.

With one latent variable the integral is one-dimensional, so the question can be
settled exactly rather than simulated: the integrand was evaluated on a 321-node
grid at the estimate, which gives the likelihood by quadrature, each respondent's
exact posterior, and the asymptotic efficiency of any proposal (`hyb_is.m` in the
work folder). Three things came out of it. First, the indicators pin the latent
variable down to a posterior standard deviation of 0.29 against the prior's 1
(0.34 in WTP space), and 500 draws already give each respondent's log-likelihood
to 0.005; the asymptotic effective-draw fraction of the prior draws, 0.250,
agrees to three digits with the count the module reports (125 of 500), which
checks the two computations against each other. Second, a plain Laplace proposal
makes the effective-draw count look splendid and the answer worse - 446 of 500
draws, and a root mean squared error in log L_i of 0.023 against 0.005 with the
draws as they come. Its weights have infinite variance, which is arithmetic
rather than a feature of this data set: a Gaussian proposal narrower than the
prior by more than the square root of two makes the integral of phi^2/q diverge,
and 0.29 is a factor of 3.4. The grid puts its asymptotic effective-draw fraction
at 0. Third, the defensive mixture is the only variant with a finite variance
(asymptotic fraction 0.61 with half the draws kept from the prior, 0.77 with a
quarter) and the only one that beats doing nothing: 0.0027 against 0.0048 in
preference space, and 0.0085 against 0.0498 in WTP space, where it is worth some
thirty times the draws. None of that changes a conclusion here, because the error
it improves on is already negligible - but it is the rule to carry to the cases
where it is not: never a plain Laplace proposal, always a mixture with the
distribution the draws came from.

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
weights concentrate on the few draws that do not. A defensive mixture - part of the draws kept from the prior, the proposal
density evaluated as the mixture, so the weights are bounded by the inverse of
the prior's share - is uniformly better than doing nothing, by a factor of 1.7
to 2.2 in root mean squared error, which is worth 3 to 5 times the draws. That
is an efficiency gain, not a repair: the error of the simulated likelihood
scales with the inverse square root of the effective draws, so moving a
24-standard-error draw-set shift below one standard error needs about 600 times
the effective draws, not 5.

**The ceiling, and what a proposal collects from it.** The effective-draw count
fixes how much any importance sampler could possibly gain. If the draws as they
come are worth ESS of R, a proposal proportional to the integrand would be worth
R, so nothing can gain more than R/ESS: 4 times in the one-dimensional HMNL, 5.5
to 6.6 in the 13-dimensional DEFRA MXL_d, 11.6 in the demo MXL_d, 21.6 in the
demo HMXL, and 143 in the multivariate Poisson-lognormal of the Natural England
work at 7 of 1000. Measured against that ceiling, with the posterior mean and
the diagonal posterior standard deviations of the draw vector read off the
weights of a first pass and a share of the draws kept from the prior - two
likelihood passes in all - as the root mean squared error of log L_i against a
large-draw benchmark, averaged over three independent scrambles:

| model | as they come | best proposal | draw-equivalent | passes | net |
|---|---|---|---|---|---|
| DEFRA MXL_d, preference space | 152 of 1000 | 433 of 1000 | 4.4x | 2 | 2.2x |
| DEFRA MXL_d, WTP space | 182 of 1000 | 379 of 1000 | 1.8x | 2 | 0.9x |
| demo HMXL, preference space | 23 of 500 | 132 of 500 | 2.4x | 2 | 1.2x |
| demo HMXL, WTP space | 23 of 500 | 147 of 500 | 0.9x | 2 | 0.5x |
| demo HMXL, WTP space, lognormal cost | 23 of 500 | 170 of 500 | 1.0x | 2 | 0.5x |

Two things follow, and the second is a warning about the count itself. Raising
the number of draws is the better bargain: four times the draws buys 5.2 to 5.4
times the accuracy here - slightly better than the inverse square root, because
a scrambled Sobol net converges faster than pseudo-random draws - while two
passes of importance sampling buy 4.4 times at best and usually less. And the
effective-draw count overstates what a proposal achieves, by a wide margin: in
HMXL in WTP space the mixture raises it from 23 to 147, six times over, and
improves the root mean squared error of log L_i not at all (0.326 against
0.308). For draws from the prior the two move together - four times the draws
raises the effective draws 3.7 times and the accuracy 4.3 times - for draws from
a proposal they do not. The weights become more even at the points that were
drawn, which is not the same as the proposal covering where the integrand has
its mass. The count is therefore a measure of how hard the simulator is working,
not a target to optimise, and "the model effectively uses only a few draws" is
not by itself a case for importance sampling.

What would make it a case is a ceiling of two orders of magnitude. That is where
one error per site in the multivariate Poisson-lognormal sits, and where a model
with many random parameters and several latent variables would sit, because the
effective-draw fraction falls geometrically in the dimension of the integral:
0.250 at one dimension, 0.086 at six and 0.046 at seven on the demo data, about
0.65 per dimension added. Writing one would also need exact derivatives. The
Fisher-identity score treats the proposal as fixed although it moves with the
parameters, so it is not the derivative of the criterion being maximised, and
`fminunc` stops after three iterations with it; freezing the proposal between
stages - estimate, refit the proposal at the new point, freeze it, estimate
again - restores ordinary derivatives and is the construction to start from. The
measurements and the prototype are in
`MATLAB/.ai/_work/2026-11-16_importance-sampling`, and the whole argument,
including the decision rule, in
`MATLAB/outputs/DCE_losowania-i-importance-sampling_v1_2026-11-16.pdf`.

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
