function [f] = probs_mdcev(data, EstimOpt, variables)
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
NVarA = EstimOpt.NVarA; % Number of attributes
Profile = EstimOpt.Profile; % Utility function version
NAlt = EstimOpt.NAlt; % Number of alternatives
N = size(y,2); % Number of decisions

% define variables to be optimized
betas = variables(1:NVarA); % betas
scale = exp(variables(NVarA+NAlt+1)); % scale parameter (sigma); one for all dataset
% scale = exp(variables(end)); % scale parameter (sigma); one for all dataset
% exp() for ensuring > 0

% alphas and gammas
if Profile == 1 % alpha profile
    alphas = variables(NVarA+1:NVarA+NAlt);
    gammas = ones(size(alphas));
    
    % another parametrisation accounting for constraints
    alphas = 1 - exp(-alphas); % between 0 and 1
elseif Profile == 2 % gamma profile
    gammas = variables(NVarA+1:NVarA+NAlt);
    alphas = zeros(size(gammas));
    
    % another parametrisation accounting for constraints
    gammas = exp(gammas); % greater than 0
end
 
% define other variables

Xa = data.Xa; % covariates (alternatives attributes) (NxNVar)
priceMat = data.priceMat; % matrix of prices of each alternative (NAltxN)
% budget = data.budget; % available budget for spending on alternatives (income)

isChosen = (y ~= 0); % if the alternative was chosen or not

M = sum(isChosen, 1); % number of consumed goods in each decision

% computing utility levels
betasZ = reshape(Xa*betas(1:NVarA), [NAlt, N]); % beta*X part of utility

%% logarithms
logf_i = log(1 - alphas) - log(y + gammas); % (NAltxN)
V = betasZ + (alphas - 1) .* log(y ./ gammas + 1) - log(priceMat);
logV = V / scale; % log(exp(V_i / scale)

% size(M)
% size(isChosen)
% size(priceMat)
% size(logf_i)
% size(logV)

priceRatio = priceMat ./ exp(logf_i);
sumV = sum(exp(logV), 1);

logprobs = (1 - M) .* log(scale) + ...
    sum(isChosen .* logf_i, 1) + ...
    log(sum(isChosen .* priceRatio, 1)) + ...
    (sum(isChosen .* logV, 1) - M .* log(sumV)) + ...
    gammaln(M); % log(factorial(M-1)) = gammaln(M)

f = logprobs';

end
