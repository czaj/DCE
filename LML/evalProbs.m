function P = evalProbs(Z, EstimOpt, bhat)

NP = EstimOpt.NP;
NRep = EstimOpt.NRep;

Fit = reshape(Z'*bhat,[NRep, NP]); % NRep x NP x 1
Fit = Fit(:,1);
Fit = exp(Fit - max(Fit));
Fit_sum = sum(Fit);
P = Fit./Fit_sum;

