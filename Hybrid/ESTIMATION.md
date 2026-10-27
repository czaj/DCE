# Estimating hybrid models reliably

These are practical recommendations for HMNL/HMXL and their MIMIC starting
models. They are a workflow, not a guarantee of a global optimum or a fixed
number of draws that works for every data set. Benchmark evidence is recorded
separately in the [validation report](OPTIMIZATION_VALIDATION.md).

## What to improve first

The current hybrid engine already has analytical scores, bounded respondent
blocks, cached inputs and checked automatic CPU/GPU selection. Improving an
individual likelihood evaluation does not remove weak identification, bad
scaling, local maxima or simulation sensitivity. Diagnose those separately.

1. Establish an identified measurement model with sensible loadings, missing
   item handling and anchors. Check direct item effects/DIF where justified.
2. Construct informed starts, then search several genuinely different basins
   at modest integration precision.
3. Continue promising solutions using curvature when available, rather than
   repeatedly restarting BFGS with identity curvature.
4. Increase integration precision, repeat across scrambles, and recompute
   inference only for the retained solutions.

Sequential estimation is useful for development and starts; the final hybrid
estimation must jointly estimate the chosen measurement and choice model.
Treating fitted factor scores as perfectly observed regressors is not full
information estimation. See [Bierlaire's hybrid estimation tutorial](https://transp-or.epfl.ch/pythonbiogeme/documentation/latent/latent.html).

DCE centers each LV globally. An unrestricted structural intercept is therefore
cancelled by that normalization and cannot be identified from the normalized
LV. Do not add it as another free parameter. A one-indicator factor also needs
explicit identifying assumptions about measurement error/reliability; naming an
observed item an LV does not create information. Deterministically related ranks
should not be treated as conditionally independent ordinary ordinal items.
Multiple factors with unrestricted cross-loadings can also retain rotational
ambiguity; an identified loading pattern is needed in addition to sign and
variance conventions. More draws or a different optimizer cannot repair these
identification issues.

## Starting values

| Source | Use | Main limitation |
|---|---|---|
| MIMIC + interaction MXL | Preferred first joint HMXL seed: measurement/structural parameters from MIMIC; tastes, covariance and observed interactions from MXL | Must align sample, coding, normalization, parameter order and units |
| Previous smaller HCM | Add one construct or a justified set of links while retaining existing parameters | Preserve the measurement meaning and identification; do not blindly copy by position |
| HMNL | A cheaper bridge when the joint measurement model is difficult | HMNL may allocate heterogeneity differently; re-estimate everything in HMXL |
| Two-stage factor-score choice model | Approximate sign and magnitude of new links for initialization | Measurement error remains; do not report its naive inference as final HCM inference |
| Small, reproducible perturbations | Explore alternatives around informed seeds | All-zero loadings/variance starts can give symmetry, flat directions or underflow |
| Large unrestricted random vectors | Usually poor value for a large HCM | Frequent extreme probabilities, tiny variances and meaningless scales |

For an initial 2–3 LV investigation, a reasonable **starting budget**, not a
theorem, is 6–12 short runs: two or three informed constructions, modest
perturbations of links/loadings/structural coefficients, all within the same
specification. Alternative link patterns are separate model specifications, not
additional starts of that model. Keep the draw seed fixed across comparable
runs. Changing the integration scramble tests a different source of uncertainty.
Do not count continuing one checkpoint as a new independent start. Flipping
an entire factor and all its linked coefficients may simply relabel the same
model; align signs before comparing results.

Preserve anchors, fixed parameters and ordered cutoffs. Perturb the *raw log
increments*, not cumulative thresholds. Preserve the positive cost convention
and preference/WTP units. Avoid initializing random taste variation at exactly
zero when its local derivative would vanish.

Observed attitude interactions in an MXL are not numerically interchangeable
with LV interactions. Their units, variance and measurement error differ.
Use them to guide seeds, with a scale mapping or a temporary factor-score
regression, then release the joint parameters. Direct demographic effects and
LV structural effects should not inadvertently duplicate unidentified paths.

## Cheap screening and exact fixed parameters

```matlab
E.ConstVarActive = 0;       % DataCleanDCE creates fminunc options
E.HybridOptimize = 1;      % optimize only BActive==1 coordinates
E.HybridFitOnly = 1;       % screening: no final Hessian/posteriors/native tables
E.NumGrad = 0;
E.Draws = 6;               % scrambled Sobol
E.NRep = 256;
E.Drawskeep = 1;
[I,old,E,O] = DataCleanDCE(I,E);
O.Algorithm = 'quasi-newton';
O.SpecifyObjectiveGradient = true;
O.MaxIterations = 300;     % example screening budget; a capped fit is partial
% bStart and active are reviewed full-length vectors for this specification.
startFields = {'HMXL_d','HMXL','HMXL_e'};
old.(startFields{E.FullCov+1}).b0 = bStart(:);
E.BActive = active(:)';
fit = HMXL(I,old,E,O);
```

`HybridOptimize` and `HybridFitOnly` default to **0**, preserving the existing
workflow unless explicitly requested. Screening returns the parameter vector,
LL, gradient, exit flag, optimizer diagnostics and prepared inputs/options.
It deliberately returns no standard errors or native inferential workbook.
For HMNL and MIMIC the corresponding explicit fields are `old.HMNL.b0` and
`old.MIMIC.b0`. Start each public call from the original input/preprocessing
workflow; `fit.INPUT` contains already prepared arrays for the likelihood, not
a replacement for the original long-format input.
Save each screening result under its own name. For retained solutions use
`HybridFitOnly=0`, `HessEstFix=5` and `RobustStd=0`, with the final draw count.
This uses inverse observed-information standard errors. Robust/sandwich
covariance remains an explicit option, not the default recommendation; it
does not correct coefficient bias from misspecification.

The active-coordinate adapter keeps `BActive==0` parameters exactly at their
starting values, including nonzero anchors. Callbacks receive the **full original
parameter vector**, so existing checkpoint formats and parameter guards remain
usable. Callback gradients/search directions are mapped back too; scalar step
sizes and solver optimality measures still refer to solver coordinates.
`output.firstorderopt` is recomputed in original units, and
`output.hybridCoordinateFirstOrderOpt` retains the solver value. Free coordinates
are centered at the starting vector even without supplied curvature. This can
change line-search steps and relative stopping behavior compared with the
legacy full-vector call; identical likelihoods need not follow identical paths.

`HybridOptimize=1` requires fminunc options and `ConstVarActive=0`. If genuine
constraints require fmincon, retain the legacy path. The built-in equality-only
path now fixes coefficients to their supplied initial values, not always zero.
`LVMatrix` now addresses the actual LV-link block after observed interactions
and the selected covariance parameterization. A zero entry fixes the matching
link at its starting value: supply a zero start if the intended path is absent.

## Units and computing resources

Use sensible units for cost and continuous covariates before constructing starts.
Center/standardize structural covariates when appropriate; preserve the coding
and any transformations needed to report effects or WTP in original units.
Inspect conditioning after diagonal scaling as well as in raw units. Merely
setting `TypicalX` does not provide the curvature transformation below.

Parallelize independent starts across machines with explicit resource budgets.
Avoid multiplying several outer starts by unrestricted inner process pools.
At small draw counts, worker or GPU overhead can exceed the saved computation;
retain checked automatic routing or benchmark CPU serial versus an existing
small pool. Reuse a pool within a batch instead of rebuilding it for every fit.

## Continue with available curvature

```matlab
E.HybridOptimize = 1;
E.HybridStartHessian = previous.hess;  % observed negative-LL Hessian, NOT ihess
E.HybridFitOnly = 0;
E.HessEstFix = 5;
E.RobustStd = 0;
% Set b0 = previous.bhat; same parameter layout and active mask.
```

The adapter uses the active positive-definite Hessian as a linear coordinate
transformation, `b_free = b0_free + chol(H) \ z`. This gives BFGS a useful local
metric without changing the likelihood or inventing a new optimizer. It rejects
indefinite, asymmetric, incompatible or severely ill-conditioned matrices; it
does not repair a model by adding arbitrary ridge terms. It accepts either a
full matrix or an active-only matrix in the same parameter order.

A Hessian from a nearby unfinished point can be useful for optimization even
when that point has invalid inference. This does not validate that point.
Curvature is **only a preconditioner**: calculate a fresh observed Hessian and
inverse observed-information covariance at the final point. Do not pass a covariance as
`HybridStartHessian`. Count the time needed to obtain curvature when assessing
speed. If parameters/sample/specification change, do not blindly reuse it.

A cheaper alternative is a positive-definite BHHH starting metric from the
respondent scores: `Jw = J .* I.W; E.HybridStartHessian = Jw' * Jw;`. Here `J`
is the unweighted respondent Jacobian at exactly the proposed start, data and
draws, with the same full parameter order. The active submatrix must pass the
same definiteness/conditioning checks. This uses BHHH only to initialize the
coordinate metric for BFGS; it is distinct from supplying BHHH at every
trust-region iteration. It is not an observed Hessian or final covariance.
The benchmark shows a useful cheap continuation but also a worse solution
from a more distant start, so retain other start/method combinations.

The adapter supports ordinary BFGS, L-BFGS, and trust-region through fminunc.
User Hessian-vector callbacks are explicitly unsupported because their mapping
would need a separate contract. Numeric vectors `TypicalX` and `FiniteDifferenceStepSize` are indexed to free coordinates;
with whitening it still describes solver coordinates. Prefer analytical scores.

## Draw counts and normalization

Suggested adaptive ladder for these models:

| Stage | Initial draw budget | Decision |
|---|---:|---|
| Coding, gradients and nesting | 32–64 | Verify equations and parameter masks, not scientific estimates |
| First exploration | 128–256 | Search several informed starts cheaply |
| Candidate comparison | 512–1024 | Re-optimize distinct plausible basins; inspect loadings and substantive effects |
| Final precision | 2048–4096 initially | Increase further, including 5k/10k if needed, until conclusions and uncertainty are stable |

These numbers are recommendations from the development workflow, not universal
accuracy thresholds. Strong measurements can concentrate each respondent's
posterior and make unconditional simulation inefficient. More latent dimensions,
correlated tastes, extreme responses and sparse indicators can require more
draws. A converged optimizer at 128 draws can be a precise solution to an
insufficiently accurate simulated objective.

At each precision level, hold the scramble fixed while comparing starts. Keep
more than one basin when their substantive implications differ. At the final
level check at least two independent scrambles, with re-optimization. Compare
key coefficients, scenario probabilities and their uncertainty, not just stars
or the sign of one ASC. Record all outcomes. LL at different draw counts is not
a nesting test. Simulation-assisted estimation and variance-reduction principles
are discussed in [Train, chapters 9–10](https://eml.berkeley.edu/books/choice2.html).

The default Sobol layout now keeps each respondent's prefix as `NRep` grows;
see [draw defaults and legacy replication](../DRAWS.md). The explicit legacy
`SobolLayout='global'` layout is not nested by respondent. The normalized LV
still uses global sample/draw moments. Directly
copying raw coefficients can therefore change the implied starting indices.
For an unchanged sample/specification with OLS/ordered measurements and no LV
scale links, use:

```matlab
args = {I.YY,I.XXa,I.Xm,I.Xs,I.Xstr,I.Xmea,I.Xmea_exp,oldDraws,E,b};
[newB,mapping] = hybrid_rebase_start('HMXL',args,newDraws);
% Set E.NRep and b0 consistently with newDraws before estimating.
```

The helper adjusts means, LV links, measurement loadings and first cutoffs
(or OLS intercepts), keeping the structural coefficients/covariance unchanged.
It preserves indices at the **same underlying latent values**; it cannot make
two different simulated integrals equal. It rejects unsupported measurement
families, LV scale effects, or a mapping that would change a fixed coefficient.
Sample changes and adding factors require a reviewed model-specific mapping.
For independent latent errors, one can first rebase a MIMIC seed to the planned
HCM latent-draw rows, then assemble the joint starting vector. FullCov=2 changes
the latent-error covariance as well; copying only the last draw rows is not a
general mapping for that case.

## Which optimizer?

| Method | Practical role |
|---|---|
| Analytical BFGS, full or free coordinates | First baseline for moderate-size HCMs; compare both paths on difficult starts rather than assuming that dropping fixed dimensions is faster |
| BFGS with a previous observed Hessian as metric | Particularly useful for continuing a nearly converged fit |
| BFGS with a BHHH starting metric | Cheap alternative when observed curvature is unavailable; check rank and keep other basins |
| L-BFGS | Worth testing for very many parameters; not automatically better for 100–200 parameters with strong correlations |
| SQP/interior-point | Use when constraints require them; equality-only fixed coefficients do not inherently require them |
| Trust-region with BHHH | A possible search metric; BHHH is not the observed Hessian and is not final robust inference |
| Trust-region with numerical observed Hessian each iteration | Potentially reliable but can cost roughly two gradient evaluations per free parameter per iteration |
| Derivative-free/global population methods | Usually expensive first choices for a smooth HCM with validated gradients |

MATLAB's BFGS/L-BFGS and trust-region options have different Hessian behavior.
In particular, L-BFGS can return an empty Hessian; select `HessEstFix=5` for final
inference. `UseParallel` concerns numerical differentiation, not an automatic
speed-up of analytical-gradient optimization. See the official
[fminunc documentation](https://www.mathworks.com/help/optim/ug/fminunc.html).
For the interpretation of the solver's returned matrix, consult
[Hessian Output](https://www.mathworks.com/help/optim/ug/hessian.html).

Do not change algorithm and draw count simultaneously when diagnosing a fit.
Compare methods from identical starts and draws, on the same hardware, with
common final-objective and original-coordinate gradient checks. Coordinate
scaling changes the meaning of solver stopping tolerances. Time caps and maximum
iterations define partial runs, not valid estimates. Check likelihood and gradient
on the returned vector independently of the optimizer's exit message.
A useful local check is the Newton correction `H \ g` relative to parameter
standard errors, alongside raw and scaled gradients; it does not prove global
optimality or adequate integration.

## Final inference and identification

Inspect the active observed Hessian and its inverse covariance. Numerical
positive definiteness alone is insufficient: also inspect
scaled conditioning, gradient, parameter guards, loadings, boundary parameters,
missingness and agreement across starts/draws. Do not silently accept a pseudoinverse
or use optimizer curvature as a substitute for a validated observed Hessian.
Tools' `hessianFromGradient` computes a fresh central-difference Hessian from
analytical scores. Weighted likelihoods require weights in the observed Hessian. These standard
errors do not by themselves cover error from insufficient simulation draws.
If sandwich inference is needed, retain one score row per independent respondent,
not per choice card or per simulation draw.

Ordered-probit exported thresholds are cumulative transformations: their delta
gradient begins with **1**, not the first threshold's value. The covariance slice
must include the complete choice/structural parameter offset. The public hybrid
tables now apply that correction. HLC output also separates its structural LV
blocks from choice-class blocks when their counts differ. Raw `bhat` remains in
its original parameterization.

Use *** p<0.01, ** p<0.05, * p<0.10 only with valid inference. Retain nonsignificant
and invalid fits; do not select a lower maximum for a preferable sign or stars.
Boundary variance tests need special justification. Measurement fit, choice fit
and substantive usefulness are distinct; a new construct may be useful without
improving choice prediction. Neither HCM estimation nor disappearance of a direct
effect establishes causality or mediation.

## Further improvements worth considering, but not included here

For one or two latent dimensions in HMNL, quadrature may be competitive; a tensor
grid becomes expensive with many random tastes and factors. Importance sampling
around individual posteriors might substantially improve integration efficiency,
but requires correct weights and derivative/normalization checks. It is not a
drop-in replacement for unconditional draws. Adaptive draws must not change
unannounced during optimization.

Alternating measurement/choice blocks can construct starts but may stall along
strongly correlated directions; release all identified parameters for the final
joint fit. EM or Bayesian estimation would require a larger, separately validated
implementation, not just switching a MATLAB optimizer. A fully log-domain joint
kernel could address extreme underflow; investigate floored observations before
assuming an optimizer change can solve it. Do not loosen parameter guards to
force apparent convergence.
