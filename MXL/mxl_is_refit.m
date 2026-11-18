function [fit,tail] = mxl_is_refit(pilot,E,err)
% Frozen population proposal or Wiktor Budzinski's adaptive diagonal MXL refit.
% Each round fixes both the coefficient draws and their proposal density.
started = tic;
previousRng = rng;
restoreRng = onCleanup(@() rng(previousRng));
rng(E.ISSeed,'twister');
d = E.NVarA;
mu = pilot.bhat(1:d);
if E.FullCov
    mask = tril(true(d));
else
    mask = logical(eye(d));
end
raw = zeros(d);
raw(mask) = pilot.bhat(d+1:end);
if any(~isfinite(pilot.bhat)) || any(diag(raw)==0)
    error('MXL:ISFit','The pilot needs finite coefficients and a nonsingular mixing factor.');
end
% Keep the pilot's signed factor for its actual points; canonicalize columns for inference.
L0 = raw.*sign(diag(raw))';
scaleIndex = d+find(ismember(find(mask),1:d+1:d^2));
x = [mu;L0(mask)];
x(scaleIndex) = log(x(scaleIndex));
sd = diag(L0);
fit.Method = E.ISMethod;
fit.Seed = E.ISSeed;
fit.RefreshSeed = mod(double(E.ISSeed)+97531,2^32);
fit.NRep = E.ISNRep;
fixed = strcmp(fit.Method,'fixed');
if fixed && ~(E.ISNRep==E.NRep && E.ISSeed==E.Seed1)
    err = sobolNormals(E,fit.Seed);
end
[Z,logQ] = populationDraws(mu,raw,err,E.NP);
logL = conditionalLogL(Z,pilot.INPUT,E);
[~,weights] = normalize(logL);
fit.PilotESS = 1./sum(weights.^2,1)';
options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true, ...
    'Display','off','MaxIterations',1000,'OptimalityTolerance',1e-8,'StepTolerance',1e-10);
W = pilot.INPUT.W(:);
rounds = 1;
if ~fixed, rounds = 5; end
for k = 1:rounds
    if ~fixed
        proposal = posteriorProposal(Z,weights,mu,sd);
        [Z,logQ] = drawProposal(proposal,E.ISNRep);
        logL = conditionalLogL(Z,pilot.INPUT,E);
    end
    base = logL-logQ;
    objective = @(t) densityObjective(t,Z,base,W,mask,scaleIndex,E.FullCov);
    [startValue,startScore] = objective(x);
    if k == 1, fit.StartLL = -startValue; end
    if ~isfinite(startValue) || any(~isfinite(startScore))
        error('MXL:ISFit','Nonfinite initial criterion or score in IS round %d.',k);
    end
    [x,~,flag,output] = fminunc(objective,x,options);
    if flag <= 0 || any(~isfinite(x))
        error('MXL:ISFit','IS round %d stopped with exit flag %d.',k,flag);
    end
    [value,score,~,weights] = objective(x);
    if ~isfinite(value) || any(~isfinite(score))
        error('MXL:ISFit','Nonfinite final criterion or score in IS round %d.',k);
    end
    mu = x(1:d);
    sd = exp(x(scaleIndex));
    fit.Rounds(k) = struct('LL',-value,'ExitFlag',flag,'Iterations',output.iterations, ...
        'MedianESS',median(1./sum(weights.^2,1)));
end

[~,g,J,~,nll] = objective(x);
H = hessianFromGradient(@(t) gradientAt(objective,t),x);
% Transform derivatives from log-SDs to the public (means, positive SDs) coordinates.
npar = numel(x);
scale = ones(npar,1);
scale(scaleIndex) = sd;
fit.bhat = x;
fit.bhat(scaleIndex) = sd;
fit.g = g./scale;
fit.jacobian = J./scale';
correction = zeros(npar,1);
correction(scaleIndex) = g(scaleIndex);
fit.hess = (H-diag(correction))./(scale*scale');
if any(~isfinite(fit.hess(:))) || any(~isfinite(fit.g))
    error('MXL:ISFit','The final IS score or observed Hessian is nonfinite.');
end
[C,flagH] = chol(fit.hess);
if flagH ~= 0 || rcond(fit.hess) <= eps
    error('MXL:ISFit','The final IS observed Hessian is indefinite or numerically singular.');
end
fit.ihess = C\(C'\eye(npar));
fit.std = sqrt(diag(fit.ihess));
[fit.diagnostics,tail] = hessianDiagnostics(fit.hess,fit.ihess,fit.g,fit.jacobian,ones(npar,1),'MXL IS');
if ~isfinite(fit.diagnostics.MaxNewtonStepSE) || fit.diagnostics.MaxNewtonStepSE > 0.01
    error('MXL:ISFit','The own-draw Newton step is %.3g SE. Check convergence before inference.',fit.diagnostics.MaxNewtonStepSE);
end
% A tiny factor can give a plausible positive Hessian at an inaccurate numerical score.
delta = eps^(1/3)*max(1,abs(x));
finiteScore = zeros(npar,1);
for j = 1:npar
    up = x; down = x;
    up(j) = up(j)+delta(j); down(j) = down(j)-delta(j);
    finiteScore(j) = (objective(up)-objective(down))/(2*delta(j));
end
fit.diagnostics.ScoreCheckStepSE = max(abs(fit.ihess*((finiteScore-g)./scale))./fit.std);
if any(~isfinite(finiteScore)) || ~isfinite(fit.diagnostics.ScoreCheckStepSE) || fit.diagnostics.ScoreCheckStepSE>0.01
    error('MXL:ISFit','Finite differences of the final criterion disagree with its score (%.3g SE).',fit.diagnostics.ScoreCheckStepSE);
end
fit.ESS = 1./sum(weights.^2,1)';
clear objective Z base logL weights err
rng(fit.RefreshSeed,'twister');
if fixed
    [fresh,logQ] = populationDraws(pilot.bhat(1:d),raw,sobolNormals(E,fit.RefreshSeed),E.NP);
else
    [fresh,logQ] = drawProposal(proposal,E.ISNRep);
end
[freshValue,freshG] = densityObjective(x,fresh,conditionalLogL(fresh,pilot.INPUT,E)-logQ,W,mask,scaleIndex,E.FullCov);
fit.RefreshStepSE = abs(fit.ihess*((freshG-g)./scale))./fit.std;
if ~isfinite(freshValue) || any(~isfinite(fit.RefreshStepSE))
    error('MXL:ISFit','The independent draw check returned a nonfinite score.');
end
if max(fit.RefreshStepSE) > 1
    warning('MXL:ISDrawSensitivity','Fresh IS draws change the score by a Newton step of %.3g SE. Compare full refits with independent ISSeed values.',max(fit.RefreshStepSE));
end
tail = [tail; {'Effective IS draws',sprintf('median %.1f of %d, minimum %.1f',median(fit.ESS),E.ISNRep,min(fit.ESS)); ...
    'Fresh IS draws',sprintf('paired score step: max %.3g SE, median %.3g SE (seed %u)',max(fit.RefreshStepSE),median(fit.RefreshStepSE),fit.RefreshSeed)}];
fit.exitf = flag;
fit.output = output;
fit.LLdetailed = W.*nll;
fit.LL = -sum(fit.LLdetailed);
fit.Dist = E.Dist(:);
details = [fit.bhat,zeros(npar,1),fit.std,pv(fit.bhat,fit.std)];
fit.DetailsA = details(1:d,:);
if E.FullCov
    factor = zeros(d);
    factor(mask) = fit.bhat(d+1:end);
    fit.DetailsVcov = factor*factor';
    fit.DetailsVcor = corrcov(fit.DetailsVcov);
    [block,fit.LatentCorrelations] = sdtri(factor(mask),fit.ihess(d+1:end,d+1:end),E);
    fit.DetailsV = [block(:,1),zeros(d,1),block(:,2:3)];
    pairs = fit.LatentCorrelations;
    fit.DetailsCor = [pairs(:,3),zeros(size(pairs,1),1),pairs(:,4:5)];
    fit.chol = details(d+1:end,[1 3 4]);
    fit.chol(scaleIndex-d,3) = NaN;
else
    fit.DetailsV = details(d+1:end,:);
    fit.DetailsV(:,4) = NaN; % zero mixing variance requires a boundary test
end
fit.R = [fit.DetailsA,fit.DetailsV];
fit.stats = pilot.stats;
fit.stats(1) = fit.LL;
fit.stats(3:6) = [1-fit.LL/fit.stats(2); mean(exp(-fit.LLdetailed./E.NCTMiss(:))); ...
    (2*npar-2*fit.LL)/E.NObs; (log(E.NObs)*npar-2*fit.LL)/E.NObs];
fit.EstimationSeconds = toc(started);
fit.success = true;
end

function err = sobolNormals(E,seed)
E.Draws = 6;
E.Seed1 = seed;
E.NRep = E.ISNRep;
err = norminv(dceQuasiDraws(E,E.NVarA))';
end

function [Z,logQ] = populationDraws(mu,factor,err,N)
d = numel(mu);
Z = reshape(mu+factor*err,d,[],N);
logQ = -reshape(sum(err.^2,1),[],N)/2-sum(log(abs(diag(factor))))-d/2*log(2*pi);
end

function logL = conditionalLogL(Z,I,E)
logL = zeros(size(Z,2),E.NP);
for n = 1:E.NP
    B = Z(:,:,n);
    B(E.Dist==1,:) = exp(B(E.Dist==1,:));
    if E.WTP_space
        B(1:end-1,:) = B(1:end-1,:).*B(E.WTP_matrix,:);
    end
    keep = repelem(~E.MissingCT(:,n),E.NAlt);
    logL(:,n) = mxl_log_panel(I.XXa(keep,:,n),I.YY(keep,n)==1,B,E.NAlt)';
end
if any(isnan(logL(:))) || any(~isfinite(max(logL,[],1)))
    error('MXL:ISFit','Nonfinite utilities or no finite panel draw for a respondent.');
end
end

function p = posteriorProposal(Z,w,mu,sd)
[d,~,N] = size(Z);
p.mu = mu;
p.sd = sd;
p.center = zeros(d,N);
p.chol = zeros(d,d,N);
for n = 1:N
    wn = w(:,n);
    ess = 1/sum(wn.^2);
    p.center(:,n) = Z(:,:,n)*wn;
    delta = Z(:,:,n)-p.center(:,n);
    cov = (ess*((delta.*wn')*delta') + d*diag(sd.^2))/(ess+d);
    [p.chol(:,:,n),flag] = chol((cov+cov')/2 + 1e-8*eye(d),'lower');
    if flag, error('MXL:ISFit','Singular proposal for respondent %d.',n); end
end
end

function [Z,logQ] = drawProposal(p,R)
[d,N] = size(p.center);
nPrior = round(0.2*R);
alpha = nPrior/R;  % density of the actual stratified mixture, including rounding
nu = 5;
normal = -d/2*log(2*pi)-sum(log(p.sd));
student = gammaln((nu+d)/2)-gammaln(nu/2)-d/2*log(nu*pi);
Z = zeros(d,R,N);
logQ = zeros(R,N);
for n = 1:N
    L = p.chol(:,:,n);
    a = p.mu+p.sd.*randn(d,nPrior);
    radial = sqrt(nu./(2*randg(nu/2,1,R-nPrior)));
    b = p.center(:,n)+(L*randn(d,R-nPrior)).*radial;
    z = [a,b];
    lp = log(alpha)+normal-sum(((z-p.mu)./p.sd).^2,1)/2;
    lt = log1p(-alpha)+student-sum(log(diag(L))) ...
        -(nu+d)/2*log1p(sum((L\(z-p.center(:,n))).^2,1)/nu);
    m = max(lp,lt);
    logQ(:,n) = (m+log(exp(lp-m)+exp(lt-m)))';
    Z(:,:,n) = z;
end
end

function [f,g,J,w,nll] = densityObjective(x,Z,base,W,mask,scaleIndex,fullCov)
[d,R,N] = size(Z);
sd = exp(x(scaleIndex));
if fullCov
    b = x;
    b(scaleIndex) = sd;
    L = zeros(d);
    L(mask) = b(d+1:end);
    if any(~isfinite(L),'all') || any(sd==0) || rcond(L)<=eps
        f = Inf; g = NaN(size(x)); J = []; w = []; nll = [];
        return
    end
    e = reshape(L\reshape(Z-x(1:d),d,[]),d,R,N);
else
    e = (Z-x(1:d))./sd;
end
logW = base - reshape(sum(e.^2,1)/2,size(base)) - sum(log(sd)) - d/2*log(2*pi);
if any(~isfinite(max(logW,[],1)))
    f = Inf; g = NaN(size(x)); J = []; w = []; nll = [];
    return
end
[logMean,w] = normalize(logW);
nll = -logMean(:);
f = W'*nll;
if nargout > 1
    weights = reshape(w,1,size(Z,2),[]);
    first = reshape(sum(e.*weights,2),d,N);
    if fullCov
        gm = L'\first;
        second = pagemtimes(e.*weights,'none',e,'transpose')-eye(d);
        factorScore = reshape(L'\reshape(second,d,[]),d*d,N);
        gs = factorScore(mask(:),:);
        gs(scaleIndex-d,:) = gs(scaleIndex-d,:).*sd;
    else
        gm = first./sd;
        gs = reshape(sum((e.^2-1).*weights,2),d,N);
    end
    J = -[gm;gs]'.*W;
    g = sum(J,1)';
end
end

function [logMean,w] = normalize(logWeights)
m = max(logWeights,[],1);
w = exp(logWeights-m);
s = sum(w,1);
logMean = m+log(s/size(w,1));
w = w./s;
end

function g = gradientAt(fun,x)
[~,g] = fun(x);
end
