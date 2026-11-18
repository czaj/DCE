function [fit,tail] = mxl_is_refit(pilot,E,err)
% Adaptive importance-sampled refit following Wiktor Budzinski's MXL experiment.
% Each round fixes both the coefficient draws and their proposal density.
started = tic;
previousRng = rng;
restoreRng = onCleanup(@() rng(previousRng));
rng(E.ISSeed,'twister');
d = E.NVarA;
mu = pilot.bhat(1:d);
pilotSD = pilot.bhat(d+1:end);
if any(~isfinite(pilot.bhat)) || any(pilotSD == 0)
    error('MXL:ISFit','The pilot needs finite coefficients and nonzero mixing SDs.');
end
sd = abs(pilotSD);
x = [mu; log(sd)];
Z = mu + pilotSD.*reshape(err,d,E.NRep,E.NP);
logL = conditionalLogL(Z,pilot.INPUT,E);
[~,weights] = normalize(logL);
fit.PilotESS = 1./sum(weights.^2,1)';
fit.Seed = E.ISSeed;
fit.RefreshSeed = mod(double(E.ISSeed)+97531,2^32);
fit.NRep = E.ISNRep;
options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true, ...
    'Display','off','MaxIterations',1000,'OptimalityTolerance',1e-8,'StepTolerance',1e-10);
W = pilot.INPUT.W(:);
for k = 1:5
    % Retain this proposal for an independent draw check after the final fit.
    proposal = posteriorProposal(Z,weights,mu,sd);
    [Z,logQ] = drawProposal(proposal,E.ISNRep);
    base = conditionalLogL(Z,pilot.INPUT,E) - logQ;
    objective = @(t) densityObjective(t,Z,base,W);
    [startValue,startScore] = objective(x);
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
    sd = exp(x(d+1:end));
    fit.Rounds(k) = struct('LL',-value,'ExitFlag',flag,'Iterations',output.iterations, ...
        'MedianESS',median(1./sum(weights.^2,1)));
end

[~,g,J,~,nll] = objective(x);
H = hessianFromGradient(@(t) gradientAt(objective,t),x);
% Transform derivatives from log-SDs to the public (means, positive SDs) coordinates.
scale = [ones(d,1);sd];
fit.bhat = [mu;sd];
fit.g = g./scale;
fit.jacobian = J./scale';
fit.hess = (H-diag([zeros(d,1);g(d+1:end)]))./(scale*scale');
if any(~isfinite(fit.hess(:))) || any(~isfinite(fit.g))
    error('MXL:ISFit','The final IS score or observed Hessian is nonfinite.');
end
[C,flagH] = chol(fit.hess);
if flagH ~= 0 || rcond(fit.hess) <= eps
    error('MXL:ISFit','The final IS observed Hessian is indefinite or numerically singular.');
end
fit.ihess = C\(C'\eye(2*d));
fit.std = sqrt(diag(fit.ihess));
[fit.diagnostics,tail] = hessianDiagnostics(fit.hess,fit.ihess,fit.g,fit.jacobian,ones(2*d,1),'MXL IS');
if ~isfinite(fit.diagnostics.MaxNewtonStepSE) || fit.diagnostics.MaxNewtonStepSE > 0.01
    error('MXL:ISFit','The own-draw Newton step is %.3g SE. Check convergence before inference.',fit.diagnostics.MaxNewtonStepSE);
end
fit.ESS = 1./sum(weights.^2,1)';
rng(fit.RefreshSeed,'twister');
[fresh,logQ] = drawProposal(proposal,E.ISNRep);
[freshValue,freshG] = densityObjective(x,fresh,conditionalLogL(fresh,pilot.INPUT,E)-logQ,W);
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
details = [fit.bhat,zeros(2*d,1),fit.std,pv(fit.bhat,fit.std)];
fit.DetailsA = details(1:d,:);
fit.DetailsV = details(d+1:end,:);
fit.R = [fit.DetailsA,fit.DetailsV];
fit.stats = pilot.stats;
fit.stats(1) = fit.LL;
fit.stats(3:6) = [1-fit.LL/fit.stats(2); mean(exp(-fit.LLdetailed/E.NCT)); ...
    (2*(2*d)-2*fit.LL)/E.NObs; (log(E.NObs)*(2*d)-2*fit.LL)/E.NObs];
fit.EstimationSeconds = toc(started);
fit.success = true;
end

function logL = conditionalLogL(Z,I,E)
logL = zeros(size(Z,2),E.NP);
for n = 1:E.NP
    B = Z(:,:,n);
    B(E.Dist==1,:) = exp(B(E.Dist==1,:));
    if E.WTP_space
        B(1:end-1,:) = B(1:end-1,:).*B(E.WTP_matrix,:);
    end
    logL(:,n) = mxl_log_panel(I.XXa(:,:,n),I.YY(:,n)==1,B,E.NAlt)';
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

function [f,g,J,w,nll] = densityObjective(x,Z,base,W)
d = size(Z,1);
sd = exp(x(d+1:end));
e = (Z-x(1:d))./sd;
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
    gm = reshape(sum(e.*weights,2),d,[])./sd;
    gs = reshape(sum((e.^2-1).*weights,2),d,[]);
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
