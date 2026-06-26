function b0 = MIMIC0_starting_values(Xmea,EstimOpt)
% MIMIC0_starting_values creates robust data-based starting values for the
% null measurement model used before HMNL/MIMIC estimation.
%
% This avoids the fragile 0.001*ones(...) initial vector. In particular, for
% ordered-probit indicators with many levels, equal 0.001 starts imply very
% high thresholds and can make upper-tail probabilities numerically zero.

if isfield(EstimOpt,'MissingIndMea') && ~isempty(EstimOpt.MissingIndMea)
    MissingIndMea = EstimOpt.MissingIndMea;
else
    MissingIndMea = zeros(size(Xmea));
end
if any(size(MissingIndMea) ~= size(Xmea))
    if isequal(size(MissingIndMea'),size(Xmea))
        MissingIndMea = MissingIndMea';
    else
        MissingIndMea = zeros(size(Xmea));
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

if isfield(EstimOpt,'SigmaMin') && ~isempty(EstimOpt.SigmaMin)
    SigmaMin = EstimOpt.SigmaMin;
else
    SigmaMin = 1e-8;
end

MinGap = 1e-4;
b0 = [];

for i = 1:size(Xmea,2)
    idxMea = (MissingIndMea(:,i) == 0) & ~MissingIndTmp & isfinite(Xmea(:,i));
    y = Xmea(idxMea,i);

    if isempty(y)
        b0 = [b0; local_default_start(EstimOpt.MeaSpecMatrix(i))]; %#ok<AGROW>
        continue
    end

    if EstimOpt.MeaSpecMatrix(i) == 0 % OLS: constant + log(sigma)
        mu = mean(y);
        sigma = std(y);
        if ~isfinite(mu)
            mu = 0;
        end
        if ~isfinite(sigma) || sigma < SigmaMin
            sigma = 1;
        end
        b0 = [b0;mu;log(sigma)]; %#ok<AGROW>

    elseif EstimOpt.MeaSpecMatrix(i) == 1 % MNL: constants for levels 2..J
        UniqueMea = unique(y);
        counts = zeros(length(UniqueMea),1);
        for j = 1:length(UniqueMea)
            counts(j) = sum(y == UniqueMea(j));
        end
        counts = counts + 0.5; % small smoothing for rare levels
        p = counts./sum(counts);
        b0 = [b0;log(p(2:end)./p(1))]; %#ok<AGROW>

    elseif EstimOpt.MeaSpecMatrix(i) == 2 % Ordered probit: thresholds only
        UniqueMea = unique(y);
        counts = zeros(length(UniqueMea),1);
        for j = 1:length(UniqueMea)
            counts(j) = sum(y == UniqueMea(j));
        end
        if length(UniqueMea) > 1
            pmin = min(1e-4,0.5/sum(counts));
            cumP = cumsum(counts);
            cumP = cumP(1:end-1)./sum(counts);
            cumP = min(max(cumP,pmin),1-pmin);
            alpha = norminv(cumP);
            for j = 2:length(alpha)
                if alpha(j) <= alpha(j-1) + MinGap
                    alpha(j) = alpha(j-1) + MinGap;
                end
            end
            b0 = [b0;alpha(1);log(max(diff(alpha),MinGap))]; %#ok<AGROW>
        end

    elseif EstimOpt.MeaSpecMatrix(i) == 3 % Poisson: constant only
        y = min(y,CountCut);
        lam = max(mean(y),1e-4);
        b0 = [b0;log(lam)]; %#ok<AGROW>

    elseif EstimOpt.MeaSpecMatrix(i) == 4 % Negative binomial: constant + log(theta)
        y = min(y,CountCut);
        lam = max(mean(y),1e-4);
        b0 = [b0;log(lam);0]; %#ok<AGROW>

    elseif EstimOpt.MeaSpecMatrix(i) == 5 % ZIP: zero-inflation logit + Poisson constants
        y = min(y,CountCut);
        positive_y = y(y > 0);
        if isempty(positive_y)
            lam = max(mean(y),1e-4);
        else
            lam = max(mean(positive_y),1e-4);
        end
        p0 = mean(y == 0);
        pzip = min(max(p0 - exp(-lam),0.01),0.99);
        b0 = [b0;local_logit(pzip);log(lam)]; %#ok<AGROW>

    elseif EstimOpt.MeaSpecMatrix(i) == 6 % ZINB: zero-inflation logit + NB constant + log(theta)
        y = min(y,CountCut);
        positive_y = y(y > 0);
        if isempty(positive_y)
            lam = max(mean(y),1e-4);
        else
            lam = max(mean(positive_y),1e-4);
        end
        p0 = mean(y == 0);
        pzip = min(max(p0 - 1/(1+lam),0.01),0.99);
        b0 = [b0;local_logit(pzip);log(lam);0]; %#ok<AGROW>
    end
end

b0(~isfinite(b0)) = 0;
if isfield(EstimOpt,'NVarcut0') && length(b0) ~= EstimOpt.NVarcut0
    % Last-resort guard. This should not normally be used; it avoids a hard
    % crash if old code counted categories differently from the missing-safe
    % code above.
    btmp = 0.001*ones(EstimOpt.NVarcut0,1);
    btmp(1:min(length(b0),EstimOpt.NVarcut0)) = b0(1:min(length(b0),EstimOpt.NVarcut0));
    b0 = btmp;
end
end

function x = local_logit(p)
p = min(max(p,1e-6),1-1e-6);
x = log(p./(1-p));
end

function b = local_default_start(spec)
if spec == 0
    b = [0;0];
elseif spec == 1
    b = [];
elseif spec == 2
    b = [];
elseif spec == 3
    b = 0;
elseif spec == 4
    b = [0;0];
elseif spec == 5
    b = [0;0];
elseif spec == 6
    b = [0;0;0];
else
    b = [];
end
end
