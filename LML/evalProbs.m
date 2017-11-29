function [P, M] = evalProbs(Z, GridMat, EstimOpt, bhat, iHess)

% NGrid = EstimOpt.NGrid;

P = mnlquick(Z, bhat);

h = @(b) sum(mnlquick(Z,b).*GridMat,2); % Calculates Mean
H = jacobianest(h, bhat);
M.Mean = [h(bhat), zeros(EstimOpt.NVarA,1), sqrt(diag(H*iHess*H')), pv(h(bhat), sqrt(diag(H*iHess*H')))];  % Using Delta Method

h = @(b) sqrt(sum(mnlquick(Z,b).*(GridMat.^2),2) - sum(mnlquick(Z,b).*GridMat,2).^2); % Calculates Std. Dev
H = jacobianest(h, bhat);
M.Std = [h(bhat), zeros(EstimOpt.NVarA,1), sqrt(diag(H*iHess*H')), pv(h(bhat), sqrt(diag(H*iHess*H')))];

Pcumsum = cumsum(P,2);
P_tmp = find(Pcumsum > 0.1);
M.Quantile = [GridMat(:, P_tmp(1))];
P_tmp = find(Pcumsum > 0.25);
M.Quantile = [M.Quantile, GridMat(:, P_tmp(1))];
P_tmp = find(Pcumsum > 0.5);
M.Quantile = [M.Quantile, GridMat(:, P_tmp(1))];
P_tmp = find(Pcumsum > 0.75);
M.Quantile = [M.Quantile, GridMat(:, P_tmp(1))];
P_tmp = find(Pcumsum > 0.9);
M.Quantile = [M.Quantile, GridMat(:, P_tmp(1))];

end

function PX = mnlquick(Z, bhat)
    Fit = Z'*bhat; % NGrid x 1
    Fit = exp(Fit - max(Fit));
    Fit_sum = sum(Fit);
    PX = (Fit./Fit_sum)';
    
end
