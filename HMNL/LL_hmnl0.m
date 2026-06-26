function f = LL_hmnl0(Xmea,EstimOpt,b)
% LL_hmnl0_missing null measurement-equation likelihood for HMNL/MIMIC starts.
%
% Robust version:
%   - missing measurement items are skipped;
%   - the objective is accumulated in log space, so products do not underflow;
%   - ordered-probit upper-tail probabilities use normcdf(-x), avoiding
%     numerical cancellation from 1 - normcdf(x);
%   - each probability contribution is floored in log space by ProbMin;
%   - scale/rate exponentials are guarded against numerical overflow.

if isfield(EstimOpt,'ProbMin') && ~isempty(EstimOpt.ProbMin)
    ProbMin = EstimOpt.ProbMin;
else
    ProbMin = 1e-300;
end
ProbMin = max(realmin,min(ProbMin,1e-6));
LogProbMin = log(ProbMin);

if isfield(EstimOpt,'SigmaMin') && ~isempty(EstimOpt.SigmaMin)
    SigmaMin = EstimOpt.SigmaMin;
else
    SigmaMin = 1e-8;
end

if isfield(EstimOpt,'MaxExp') && ~isempty(EstimOpt.MaxExp)
    MaxExp = EstimOpt.MaxExp;
else
    MaxExp = 50;
end

if isfield(EstimOpt,'MissingIndMea') && ~isempty(EstimOpt.MissingIndMea)
    MissingIndMea = EstimOpt.MissingIndMea;
else
    MissingIndMea = zeros(size(Xmea));
end
if any(size(MissingIndMea) ~= size(Xmea))
    if isequal(size(MissingIndMea'),size(Xmea))
        MissingIndMea = MissingIndMea';
    else
        error('Incorrect size of EstimOpt.MissingIndMea in LL_hmnl0_missing')
    end
end

if isfield(EstimOpt,'MissingInd_tmp') && ~isempty(EstimOpt.MissingInd_tmp)
    MissingIndTmp = EstimOpt.MissingInd_tmp(:) ~= 0;
else
    MissingIndTmp = false(size(Xmea,1),1);
end

if isfield(EstimOpt,'CountCut') && ~isempty(EstimOpt.CountCut)
    CountCut = EstimOpt.CountCut;
else
    CountCut = 80;
end

l = 0;
ll = zeros(size(Xmea,1),1);

for i = 1:size(Xmea,2)
    idxMea = (MissingIndMea(:,i) == 0) & ~MissingIndTmp & isfinite(Xmea(:,i));
    y = Xmea(idxMea,i);

    if EstimOpt.MeaSpecMatrix(i) == 0 % OLS: constant + log(sigma)
        bx = b(l+1:l+2);
        if any(idxMea)
            sigma = max(exp(min(max(bx(2),-MaxExp),MaxExp)),SigmaMin);
            logL = -0.5*log(2*pi) - log(sigma) - 0.5*((y-bx(1))./sigma).^2;
            ll(idxMea,:) = ll(idxMea,:) + max(LogProbMin,logL);
        end
        l = l + 2;

    elseif EstimOpt.MeaSpecMatrix(i) == 1 % MNL: constants for levels 2..J
        UniqueMea = unique(y);
        k = length(UniqueMea) - 1;
        bx = [0;b(l+1:l+k)];
        if any(idxMea)
            bx = bx - max(bx);
            logP = bx - log(sum(exp(bx)));
            logL = LogProbMin*ones(size(y));
            for j = 1:length(UniqueMea)
                logL(y == UniqueMea(j)) = max(LogProbMin,logP(j));
            end
            ll(idxMea,:) = ll(idxMea,:) + logL;
        end
        l = l + k;

    elseif EstimOpt.MeaSpecMatrix(i) == 2 % Ordered probit: thresholds only
        UniqueMea = unique(y);
        k = length(UniqueMea) - 1;
        bx = b(l+1:l+k);
        if any(idxMea)
            if k == 0
                logL = zeros(size(y));
            else
                alpha = cumsum([bx(1);exp(min(max(bx(2:end),-MaxExp),MaxExp))]);
                P = zeros(size(y));
                P(y == UniqueMea(1)) = normcdf(alpha(1));
                P(y == UniqueMea(end)) = normcdf(-alpha(end));
                for j = 2:k
                    P(y == UniqueMea(j)) = normcdf(alpha(j)) - normcdf(alpha(j-1));
                end
                logL = log(max(ProbMin,P));
            end
            ll(idxMea,:) = ll(idxMea,:) + logL;
        end
        l = l + k;

    elseif EstimOpt.MeaSpecMatrix(i) == 3 % Poisson: constant only
        bx = b(l+1);
        if any(idxMea)
            y(y > CountCut) = CountCut;
            loglam = min(max(bx,-MaxExp),MaxExp);
            lam = exp(loglam);
            logL = loglam.*y - lam - gammaln(y+1);
            ll(idxMea,:) = ll(idxMea,:) + max(LogProbMin,logL);
        end
        l = l + 1;

    elseif EstimOpt.MeaSpecMatrix(i) == 4 % Negative binomial: constant + log(theta)
        bx = b(l+1);
        theta = exp(min(max(b(l+2),-MaxExp),MaxExp));
        if any(idxMea)
            y(y > CountCut) = CountCut;
            loglam = min(max(bx,-MaxExp),MaxExp);
            lam = exp(loglam);
            u = theta./(theta+lam);
            logL = gammaln(theta+y) - gammaln(theta) - gammaln(y+1) + theta.*log(u) + y.*log(1-u);
            ll(idxMea,:) = ll(idxMea,:) + max(LogProbMin,logL);
        end
        l = l + 2;

    elseif EstimOpt.MeaSpecMatrix(i) == 5 % ZIP: zero-inflation logit + Poisson constants
        bx = b(l+1:l+2);
        if any(idxMea)
            y(y > CountCut) = CountCut;
            p = 1./(1+exp(-min(max(bx(1),-MaxExp),MaxExp)));
            loglam = min(max(bx(2),-MaxExp),MaxExp);
            lam = exp(loglam);
            logL = zeros(size(y));
            IndxZIP = y == 0;
            P0 = p + (1-p).*exp(-lam);
            logL(IndxZIP) = log(max(ProbMin,P0));
            logL(~IndxZIP) = log(max(ProbMin,1-p)) + loglam.*y(~IndxZIP) - lam - gammaln(y(~IndxZIP)+1);
            ll(idxMea,:) = ll(idxMea,:) + max(LogProbMin,logL);
        end
        l = l + 2;

    elseif EstimOpt.MeaSpecMatrix(i) == 6 % ZINB: zero-inflation logit + NB constant + log(theta)
        bx = b(l+1:l+2);
        theta = exp(min(max(b(l+3),-MaxExp),MaxExp));
        if any(idxMea)
            y(y > CountCut) = CountCut;
            p = 1./(1+exp(-min(max(bx(1),-MaxExp),MaxExp)));
            loglam = min(max(bx(2),-MaxExp),MaxExp);
            lam = exp(loglam);
            u = theta./(theta+lam);
            logL = zeros(size(y));
            IndxZIP = y == 0;
            P0 = p + (1-p).*(u.^theta);
            logL(IndxZIP) = log(max(ProbMin,P0));
            logL(~IndxZIP) = log(max(ProbMin,1-p)) + gammaln(theta+y(~IndxZIP)) - gammaln(theta) - gammaln(y(~IndxZIP)+1) + theta.*log(u) + y(~IndxZIP).*log(1-u);
            ll(idxMea,:) = ll(idxMea,:) + max(LogProbMin,logL);
        end
        l = l + 3;
    end
end

f = -sum(ll);
if ~isfinite(f)
    f = 1e100;
end
end
