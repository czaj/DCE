function f = LL_lmxl(Q,Z,EstimOpt,b0)


NP = EstimOpt.NP;
NRep = EstimOpt.NRep;

Fit = reshape(Z'*b0,[NRep, NP]); % NRep x NP x 1
Fit = exp(Fit - max(Fit));
Fit_sum = sum(Fit);
Fit = Fit./Fit_sum;

P = sum(Q.*Fit',2); 
f = -log(P);

