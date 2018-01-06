function [f,j,h] = LL_lml(GridProbs,b_mtx,EstimOpt,b0)

% save tmp_LL_lml
% return

NP = EstimOpt.NP;
NRep = EstimOpt.NRep;
NVar = size(b_mtx,1);

if nargout == 1
    
    Fit_exp = reshape(b_mtx'*b0,[NRep,NP]); % NRep x NP
    Fit_exp = exp(Fit_exp - max(Fit_exp,[],1)); % NRep x NP
    Fit_sum = sum(Fit_exp,1); % 1 x NP
    Fit_exp = Fit_exp./Fit_sum; % NRep x NP
    P_sum = sum(GridProbs.*Fit_exp',2); % NP x 1
    
elseif nargout == 2
    
    Fit = reshape(b_mtx'*b0,[NRep,NP]); % NRep x NP
    Fit_exp = exp(Fit - max(Fit,[],1)); % NRep x NP
    Fit_sum = sum(Fit_exp,1); % 1 x NP
    Fit_exp = Fit_exp./Fit_sum; % NRep x NP
    P = GridProbs.*Fit_exp'; % NP x NRep
    P_sum = sum(P,2); % NP x 1
    
    P_share = P./P_sum; % NP x NRep
    j = (P_share - Fit_exp').*permute(reshape(b_mtx,[NVar,NRep,NP]),[3,2,1]); % NP x NRep x NVar
    j = -reshape(sum(j,2),[NP,NVar]); % NP x NVar
    
elseif nargout == 3
    
    h = [];
    error('Analytical Hessian missing')
    
end

if EstimOpt.RealMin == 1
    f = -log(max(P_sum,realmin));
else
    f = -log(P_sum);
end