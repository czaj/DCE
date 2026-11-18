# Simulation draws

Update **DCE and Tools together**. Sobol estimation defaults to `Draws=6`
(scrambled Sobol), `NRep=1024`, `HaltonSkip=0`, `HaltonLeap=0` and
`SobolLayout='respondent'`. Explicit counts such as 1000 or 5000 are honored
without rounding. Halton retains skip 1; non-Sobol rules default to 1000 draws.

Powers of two preserve Sobol net balance when points are not skipped or thinned.
They are unrelated to respondent count. Balance does not establish adequate
integration accuracy: check parameters, gradients, uncertainty and scrambles.

## Choosing a draw budget

Czajkowski and Budziński compared repeated fits across simulators, dimensions and draw budgets. Their thresholds concerned LL and coefficient accuracy under the designs studied. Estimating model-based standard errors also requires accurate curvature. Choose the budget by comparing the quantities you plan to report across independently randomized full fits. Record the simulator, seed, start, draw count and elapsed time. [Czajkowski and Budziński (2019)](https://doi.org/10.1016/j.jocm.2019.04.003)

## Reading the diagnostics

1. Check the original score, the Newton step in units of SE and the free-parameter Hessian of the negative LL. A positive-definite observed Hessian and a small step support local convergence on the current objective. Try multiple starts on the same draws and compare their LLs to find better local solutions.
2. Kish ESS, `1/sum(w.^2)` for normalized respondent weights, measures concentration. The default warning flags more than 5% of all respondents with ESS below 10. The output separately counts respondents without a finite ESS. The legacy field `ISCeiling = NRep/MedianESS` bounds the possible increase in this particular weight statistic. Integration error, gradients, SEs and runtime depend on the integrand and the algorithm. [Elvira, Martino and Robert](https://arxiv.org/abs/1809.04129)
3. Rotating respondent draw blocks changes which draws integrate each panel. `RotatedStepSE` uses the rotated score and original curvature. `PairedStepSE` uses the difference between rotated and original scores when both are available. The flag requires a paired score change. Without the original score only the unpaired step is reported. Both are local approximations. An unchanged block assignment supplies no comparison. A large change calls for independent refits. A small change in one rotation can coexist with severely understated SEs. `DrawsSuggested` remains present as NaN for compatibility. The convergence rate, the draw budget and any claim of model misspecification require separate evidence.
4. In DCE estimators using the shared Hessian diagnostics, indefinite or numerically singular free curvature leaves the point estimates and Hessian available and suppresses model-based inference. A free coefficient whose Hessian column is exactly zero carries no information at the estimate. It keeps its value, is reported without an SE and still counts among the estimated parameters; the SEs of the other coefficients come from the remaining block. Coefficients fixed by the user (`BActive = 0`) stay outside the free block. The median SE(BHHH)/SE(H) compares the outer product of scores with curvature. The interval [0.8, 1.25] is descriptive. Disagreement can reflect sampling variation, misspecification, lack of convergence or integration error. Agreement in one fit can coexist with unstable parameters across seeds. For reporting model-based SEs in ordinary MXL, request `RobustStd=0` and `HessEstFix=5` (fresh central differences of the analytical gradient).
5. Compare full fits on independent draw sequences, including coefficients, economically relevant contrasts and SEs. For different simulators, score candidate solutions on a common, sufficiently accurate integration rule before comparing LL. Independent seeds complement multiple starts on the same objective.

## Optional importance sampling after ordinary MXL

```matlab
E.FullCov = 1;                 % correlated normal/lognormal coefficients
E.ImportanceSampling = 'fixed'; % one refit; reuse the pilot points by default
% E.ISNRep = 4096;             % optional separate refit budget
% E.ISSeed = 71629;            % optional independent refit scramble
E.RobustStd = 0;
E.HessEstFix = 5;
E.EffectiveMoments = 0;        % compare both fits on the underlying normal scale
R = MXL(INPUT,Results_old,E,OptimOpt);
if R.IS.success
    disp([R.DetailsA(:,3), R.IS.DetailsA(:,3)])       % SEs of means
    disp([R.DetailsV(:,3), R.IS.DetailsV(:,3)])       % SEs of marginal SDs
    disp([R.LatentCorrelations(:,3:4), R.IS.LatentCorrelations(:,3:4)]) % rho and SE
end
```

By default MXL fits ordinary SML only. Enabling this option retains the ordinary fit in the usual fields and adds a separate fit in `R.IS`. With display/export enabled it writes a second native workbook. An explicit `OutputFilePrefix` gets an `IS_` prefix for this file. `Drawskeep=1` is unnecessary. The IS table gives underlying normal means and positive SDs, including the log-scale parameters of lognormal tastes. Any pilot effective moments or robust SE settings apply to the pilot. The IS SEs use the inverse of its own final observed Hessian. The refit checks its final score against finite differences of the frozen criterion and records the discrepancy as `diagnostics.ScoreCheckStepSE`.

In `'fixed'` mode the proposal is the fitted latent population normal. At the pilot's draw count and seed the refit reuses its actual points, including signs of Cholesky columns. With a different count or IS seed it constructs a scrambled Sobol net using the pilot layout, skip and leap. The refit keeps coefficient values and conditional panel log probabilities fixed while estimating the population mean and covariance using density ratios. The refit caches latent points and panel log probabilities, with memory proportional to respondents × draws × random coefficients. It fits positive diagonal scales in log coordinates and returns the observed Hessian in the DCE parameter coordinates. Full covariance uses column-wise lower-triangular order. `DetailsVcov` and `DetailsVcor` are latent covariance and correlation. For full MXL, `LatentCorrelations` has columns [row, column, latent correlation, delta SE, two-sided p for zero correlation]. The results table and workbook include those pairs with their SEs, p-values and significance stars when inference is available. `DetailsV` reports marginal latent SDs with delta-method SEs, including the covariance among factor entries, and their Wald p-values and significance stars. Compare full runs with independent `Seed1` values and check the SE of the intended contrast.

Both modes support normal/lognormal random coefficients, respondent weights, whole omitted tasks, and preference or single-cost WTP space. Every respondent must have an observed task. Fixed coefficients, constraints, covariates of means or scale, partial choice sets and extra utility transformations raise `MXL:ISSpecification`. A numerical failure of the refit retains the pilot and records `R.IS.success=false` with a message.

`ImportanceSampling=1` and `ImportanceSampling='adaptive'` select the published diagonal option. Its default IS seed is 1234. It follows Wiktor Budziński's five-round adaptive proposal. Within each round it fixes coefficient draws and conditional panel log probabilities and differentiates the population density in the importance weights. Each respondent's proposal mixes a population normal component (about 20%) and a Student t component with five degrees of freedom centered on the weighted draws. The actual mixture proportions account for rounding. The four adaptations after the first optimization update the proposals using the previous round's weights. `R.IS.Rounds` records each LL, optimizer exit flag, iterations and median ESS. Round LLs use different integration points. `R.IS.ESS` gives final respondent ESS. `RefreshStepSE` gives the SE-scaled change in score on new draws from the last proposal, using the final curvature. Repeat the full run with a different `ISSeed` to check both adaptation and fit stability. Changing `Seed1` additionally varies the SML pilot.

### Correlated WTP example with a fixed population proposal

The example has 180 respondents, 12 binary choices per respondent, a normal WTP coefficient and a lognormal cost scale. The true correlation between WTP and **log** scale is 0.65, and the log-scale SD is 2.4. Each method uses the same simulated choices (data seed 918367) and starting coefficients. Each table cell containing two SEs gives results for independent respondent-specific Sobol scrambles, seeds 17431 and 71629, in that order.

| Estimation strategy | SE of mean WTP | SE of WTP SD | SE of latent correlation | Total seconds |
|---|---:|---:|---:|---:|
| Ordinary MXL, 512 draws | 0.0653 / 0.0722 | 0.0469 / 0.0408 | 0.1093 / 0.1124 | 1.81 |
| Pilot 512 + fixed IS 512 | 0.1581 / 0.1514 | 0.1518 / 0.1309 | 0.1454 / 0.1406 | 2.10 |
| Ordinary MXL, 4096 draws | 0.1828 / 0.1737 | 0.1289 / 0.1489 | 0.1723 / 0.1459 | 5.84 |
| Pilot 512 + fixed IS 4096 | 0.1690 / 0.1721 | 0.1414 / 0.1388 | 0.1507 / 0.1575 | 3.08 |
| Independent quadrature refit | 0.1727 | 0.1405 | 0.1575 | — |

All SEs use the inverse of the final observed Hessian. The marginal SD and correlation use the delta method with the full covariance of their Cholesky entries. The refit uses positive Cholesky diagonals. A pilot can have negative diagonals and represent the same covariance, so compare marginal SDs and correlations directly. The reference maximizes the panel likelihood by tensor Gauss–Hermite quadrature. Raising its order from 160 to 224 per axis and rotating the axes changed each of the three reference SEs by less than 0.12%.

Times are medians of six complete, warmed public `MXL` calls per strategy (three runs per seed), measured on one Intel i9-7900X with a two-worker thread pool, MATLAB R2026b, CPU, workbook export off. Both IS times include their 512-draw pilot, final Hessian and fresh-draw check. The IS stage alone took 0.32 s at 512 draws and 1.28 s at 4096. On this sample the 4096-draw fixed refit gave all three SEs within 5% of the reference on both seeds, at 53% of the total time of ordinary MXL 4096. Repeating the benchmark on a second computer reproduced the SEs and gave 3.11 s versus 6.19 s for these two strategies. Two scrambles on one dataset measure numerical sensitivity. Coverage requires repeated datasets.

With DCE, its tests and Tools on the path, run `comparison = test_mxl_frozen_draws(tempname)` to reproduce the timing and the SEs of all five distribution summaries. Calling the test without requesting an output runs only the regression checks.

### Wide-scale WTP example with adaptive IS

The example has 300 respondents, 12 tasks, four independent normal WTP parameters and a lognormal scale with log-SD 2.5. The historical global Sobol pilot uses 1000 draws. The following same-machine timings include the common pilot, its verification and export. Both second stages start at the same pilot solution.

| Fit | SE of mean opt-out constant | Range of WTP SEs relative to reference | Elapsed |
|---|---:|---:|---:|
| Ordinary SML, 1000 | 0.0404 | 20–38% | 20.1 s |
| Ordinary SML, 32768 | 0.1020 | 25–57% | 219.5 s |
| Adaptive IS, 5 × 1000 | 0.2016 | 98–105% | 29.3 s |

The independent reference uses 32768 fixed coefficient draws from one population proposal and gives 0.2064. Freezing just the original 1000 coefficient draws and reweighting their density already gives 0.2018. Current respondent-specific Sobol draws give ordinary-SML SEs of 0.0454 at 1024 and 0.0689 at 32768 for this constant. The small rotated step (0.23 SE) in the historical 32768-draw fit missed the curvature error.

Across 50 independent datasets from this wide-scale design, nominal 95% model-based intervals covered the eight WTP parameters in 146/400 cases with ordinary SML 1000 and 376/400 with adaptive IS. Intervals within each dataset share data. At log-SD 0.5, changing to IS in the control dataset changed the ten SEs by −7% to +3%. In a 13-dimensional empirical DEFRA fit, two 4000-draw adaptive runs gave median BHHH/H ratios of 0.90 and 1.01, yet coefficients differed by a median 2.52 and a maximum 37.38 times the larger SE. Check independent refits before reporting inference. Run `tests/test_mxl_importance_sampling` to reproduce the wide-scale example. `tests/test_mxl_frozen_draws` checks the public refit in both spaces with diagonal and full covariance against independently computed panel likelihoods, scores and finite-difference Hessians.

## Earlier simulator comparisons

The tables below report separate empirical checks. Their comparison metrics differ from coverage of model-based SEs.

Respondent weight concentration in DEFRA (263 respondents, 13 attributes, 1000 scrambled Sobol draws):

| model | median | 5th pct. | min | below 10 draws | median SE(BHHH)/SE(H) |
|---|---|---|---|---|---|
| preference space, normal cost | 130 | 39 | 3.4 | 1.1% | 1.03 |
| preference space, lognormal cost | 152 | 48 | 9.6 | 0.4% | 1.09 |
| WTP space, normal cost | 129 | 12 | 1.6 | 3.8% | 1.23 |
| WTP space, lognormal cost | 182 | 5.0 | 1.3 | 10.3% | 4.22 |

Local rotated-block steps and actual differences between full fits on different sequences, expressed in SEs:

| model | draws | predicted shift | largest actual | median actual |
|---|---|---|---|---|
| demo MXL_d, preference space | 500 | 0.8 - 1.8 SE | 2.9 SE | 0.26 SE |
| demo MXL_d, preference space | 2000 | 0.2 - 0.6 SE | 0.5 SE | 0.15 SE |
| DEFRA, WTP space, lognormal cost | 1000 | 1.3 - 7.9 SE | 24.0 SE | 1.8 SE |
| DEFRA, WTP space, lognormal cost | 5000 | 0.8 - 1.6 SE | 14.7 SE | 1.1 SE |

Hybrid-model checks on the toolbox demo data (500 draws):

| model | draw dimensions | effective draws | below 10 | shift, worst | shift, median | SE(BHHH)/SE(H) |
|---|---|---|---|---|---|---|
| MXL_d, no indicators | 6 | 43 of 500 | 3.9% | 1.80 SE | 0.41 SE | 0.94 |
| HMNL, one indicator | 1 | 125 of 500 | 0% | 0.03 SE | 0.006 SE | 0.93 |
| HMNL, five indicators | 1 | 125 of 500 | 0% | 0.03 SE | 0.001 SE | 1.12 |
| HMNL, five indicators, WTP space | 1 | 75 of 500 | 0.4% | 0.12 SE | 0.04 SE | 1.07 |
| HMXL, five indicators | 7 | 23 of 500 | 16.7% | 1.29 SE | 0.20 SE | 0.99 |
| HMXL, five indicators, WTP space | 7 | 23 of 500 | 18.5% | 3.02 SE | 0.38 SE | 0.98 |
| HMXL, WTP space, lognormal cost | 7 | 23 of 500 | 19.1% | 3.87 SE | 0.33 SE | 0.96 |

Earlier fixed-parameter DEFRA comparisons used 1000 points. Entries give the RMSE of individual log likelihoods against a 200,000-draw benchmark.

| simulator | WTP space, lognormal | preference space, lognormal |
|---|---|---|
| draws as they come | 0.254 | 0.060 |
| Laplace proposal, inflation 1.15 | 1.575 | 0.058 |
| Laplace proposal, inflation 2 | 1.336 | 0.389 |
| half prior, half Laplace, inflation 1.5 | 0.151 | 0.052 |
| three quarters Laplace, inflation 1.15 | 0.158 | 0.027 |

A separate moment-proposal comparison averaged individual log-likelihood RMSE across three scrambles. The final column divides the draw-equivalent gain by an assumed number of likelihood passes. This can change the ranking. In the early DEFRA preference-space lognormal example, the measured-cost gain was 3.44 for ordinary Sobol 4000 versus 1.78 for moment IS 1000 (individual log-likelihood RMSE 0.02549 and 0.02825). Compare elapsed time including the pilot, re-estimation and observed Hessian when choosing an estimator.

| model | as they come | best proposal | draw-equivalent | passes | net |
|---|---|---|---|---|---|
| DEFRA MXL_d, preference space | 152 of 1000 | 433 of 1000 | 4.4x | 2 | 2.2x |
| DEFRA MXL_d, WTP space | 182 of 1000 | 379 of 1000 | 1.8x | 2 | 0.9x |
| demo HMXL, preference space | 23 of 500 | 132 of 500 | 2.4x | 2 | 1.2x |
| demo HMXL, WTP space | 23 of 500 | 147 of 500 | 0.9x | 2 | 0.5x |
| demo HMXL, WTP space, lognormal cost | 23 of 500 | 170 of 500 | 1.0x | 2 | 0.5x |

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
