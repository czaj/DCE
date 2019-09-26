function CV = LL_mnl_bs(Y, XXa, EstimOpt, b0)

NP = size(Y,2);
%Y = reshape(Y, 
Xa = reshape(XXa,EstimOpt.NAlt*EstimOpt.NCT*NP, EstimOpt.NVarA);
betaX = Xa*b0;

v = reshape(betaX,EstimOpt.NAlt,EstimOpt.NCT, NP);


maxv = max(v,[],1);
evdiff = exp(v- maxv(ones(EstimOpt.NAlt,1),:,:)); %clear v maxv
sum_evdiff = sum(evdiff,1); %1 by N

P = evdiff./sum_evdiff(ones(EstimOpt.NAlt,1),:,:); % NAlt x N

probs= reshape(P(Y == 1), EstimOpt.NCT, NP);

CV = prod(probs,1).^(1/EstimOpt.NCT);

CV = -CV';






