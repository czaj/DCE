function f = LL_lml(GridProbs,b_mtx,EstimOpt,b0)

% save tmp1

NP = EstimOpt.NP;
NRep = EstimOpt.NRep;

Fit = reshape(b_mtx'*b0,[NRep,NP]); % NRep x NP x 1
Fit = exp(Fit - max(Fit));
Fit_sum = sum(Fit);
Fit = Fit./Fit_sum;

P = sum(GridProbs.*Fit',2);

if EstimOpt.RealMin == 1
    f = -log(max(P,realmin));
else
    f = -log(P);
end
