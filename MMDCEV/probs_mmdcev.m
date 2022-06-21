function [f] = probs_mmdcev(data, EstimOpt, variables)
% Function calculates probabilities of the MDCEV model.
% It returns log(probabilities)
% See:
% Bhat, Chandra R. 2008. 
% ‘The Multiple Discrete-Continuous Extreme Value (MDCEV) Model: Role of 
% Utility Function Parameters, Identification Considerations, and Model 
% Extensions’. Transportation Research Part B: Methodological 42 (3): 
% 274–303. https://doi.org/10.1016/j.trb.2007.06.002.
%
% Estimating probabilities as defined in eq. (19) or (20) in the article.
%
%
% EstimOpt.Profile -- sets the utility profile
% 1 -- "alpha" profile
% 2 -- "gamma" profile
%
% variables -- includes variables to be estimated:
%   1. betas    -- alternative attributes parameters (NVarA)
%   2. alphas   -- satiation parameters (NAlt)
%      or
%   2. gammas   -- translation parameters (enable corner solutions) (NAlt)
%   3. sigma    -- scale parameter (1)
%
% data contains:
%   Y           -- dependent variables (quantities demanded); (NAltxN)
%   Xa          -- covariates (alternatives attributes) (NxNVar)
%   priceMat    -- matrix of prices of each alternative (NAltxN)


% save LL_mnl
% return

y = data.Y; % dependent variables (quantities demands); size: NAltxN
err = data.err;
FullCov = EstimOpt.FullCov;
NVarA = EstimOpt.NVarA; % Number of attributes
NVarP = EstimOpt.NVarP; % Number of parameters for alpha/gamma profile
NVarM = EstimOpt.NVarM;
NRep = EstimOpt.NRep;
Profile = EstimOpt.Profile; % Utility function version
SpecProfile = EstimOpt.SpecProfile;
Dist = EstimOpt.Dist;
NAlt = EstimOpt.NAlt; % Number of alternatives
NP = EstimOpt.NP;
NCT = EstimOpt.NCT;
RealMin = EstimOpt.RealMin;
N = size(y,2); % Number of decisions

% define variables to be optimized
betas = variables(1:NVarA); % betas
if FullCov == 0
    b0v = variables(NVarA+1:2*NVarA);
    VC = diag(b0v);
    l = 2*NVarA;
else
    b0v = variables(NVarA+1:NVarA+sum(1:NVarA));
    VC = tril(ones(NVarA));
    VC(VC == 1) = b0v;
    l = NVarA+sum(1:NVarA);
end
b_mtx = betas + VC*err; 
if sum(Dist == 1) > 0 % Log-normal
    b_mtx(Dist == 1,:) = exp(b_mtx(Dist == 1,:));
end
b_mtx = reshape(b_mtx,[NVarA,NRep,NP]);

b_profile = variables(l+1:l+NVarP*(1+NVarM));
scale = exp(variables(l+NVarP*(1+NVarM)+1)); % scale parameter (sigma); one for all dataset
% scale = exp(variables(end)); % scale parameter (sigma); one for all dataset
% exp() for ensuring > 0

Xm = data.Xm; % covariates (socio-demographic) (NxNVarM)
% alphas and gammas
if Profile == 1 % alpha profile
    if NVarM == 0
        alphas = b_profile;
        alphas = 1 - exp(-alphas); % between 0 and 1
    else
        alpha = reshape(b_profile, [NVarM+1, NAlt]);
        alphas = 1 - exp(-Xm*alpha)'; % NAlt x N
    end
    gammas = ones(size(alphas));
elseif Profile == 2 % gamma profile
    if NVarM == 0
        gammas = b_profile;
        gammas = exp(gammas); % greater than 0
    else
        gamma = reshape(b_profile, [NVarM+1, NAlt]);
        gammas= exp(Xm*gamma)'; % NAlt x N
    end
    alphas = zeros(size(gammas));
elseif Profile == 3 % SpecProfile
    indx = length(unique(SpecProfile(1,SpecProfile(1,:) ~= 0))); % unique alphas
    if NVarM == 0
        a = b_profile(1:indx);
        a = 1 - exp(-a);
        a = a(SpecProfile(1,SpecProfile(1,:) ~= 0)); 
        alphas = zeros(NAlt,1);
        alphas(SpecProfile(1,:) ~= 0) = a;
    
        g = b_profile(indx+1:end);
        g = exp(g); % greater than 0
        g = g(SpecProfile(2,SpecProfile(2,:) ~= 0)); 
        gammas = ones(NAlt,1);
        gammas(SpecProfile(2,:) ~= 0) = g;
    else

        a = reshape(b_profile(1:indx*(NVarM+1)), [NVarM+1, indx]);
        alpha = 1 - exp(-Xm*a)';
        alpha = alpha(SpecProfile(1,SpecProfile(1,:) ~= 0),:); 
        alphas = zeros(NAlt,N);
        alphas(SpecProfile(1,:) ~= 0,:) = alpha;
    
        g = reshape(b_profile(indx*(NVarM+1)+1:end), [NVarM+1, NVarP - indx]);
        gamma = exp(Xm*g)'; % greater than 0
        gamma = gamma(SpecProfile(2,SpecProfile(2,:) ~= 0),:); 
        gammas = ones(NAlt,N);
        gammas(SpecProfile(2,:) ~= 0,:) = gamma;
    end
end
 
% define other variables

Xa = data.Xa; % covariates (alternatives attributes) (NxNVar)
priceMat = data.priceMat; % matrix of prices of each alternative (NAltxN)

% budget = data.budget; % available budget for spending on alternatives (income)

isChosen = (y ~= 0); % if the alternative was chosen or not

M = sum(isChosen, 1); % number of consumed goods in each decision

% computing utility levels
betasZ = zeros(NAlt*NCT, NRep, NP);
for i = 1:NP
    betasZ(:,:,i) = Xa(:,:,i)*b_mtx(:,:,i);
end
betasZ = reshape(permute(betasZ, [1, 3, 2]), [NAlt, N, NRep]);
% betasZ = reshape(Xa*betas(1:NVarA), [NAlt, N]);
% %betasZ = reshape(Xa*betas(1:NVarA), [NAlt, N]); % beta*X part of utility

%% logarithms
% logf_i = log(1 - alphas) - log(y + gammas); % (NAltxN)
logc_i = log(1 - alphas) - log(y + gammas) - log(priceMat); % (NAltxN)
V = betasZ + (alphas - 1) .* log(y ./ gammas + 1) - log(priceMat);
logV = V / scale; % NAlt x N x NRep

% size(M)
% size(isChosen)
% size(priceMat)
% size(logf_i)
% size(logV)

% priceRatio = priceMat ./ exp(logf_i);
priceRatio = exp(-logc_i);
sumV = sum(exp(logV), 1);

logprobs = (1 - M) .* log(scale) + ...
    sum(isChosen .* logc_i, 1) + ...
    log(sum(isChosen .* priceRatio, 1)) + ...
    (sum(isChosen .* logV, 1) - M .* log(sumV)) + ...
    gammaln(M); % log(factorial(M-1)) = gammaln(M)

probs = reshape(exp(logprobs), [NCT, NP, NRep]);
probs = mean(prod(probs,1),3);

if RealMin == 1
    f = log(max(probs',realmin));
else
    f = log(probs');
end
end
