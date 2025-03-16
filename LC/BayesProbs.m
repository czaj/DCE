function Scores = BayesProbs(YY,Xa,Xc,Xs, MissingInd,EstimOpt,B)

% save tmp_LL_lc
% return

if sum(EstimOpt.BActiveClass == 0,1) == 0
    Bclass = reshape(B(1:EstimOpt.NClass*EstimOpt.NVarA),[EstimOpt.NVarA,EstimOpt.NClass]);
    l = EstimOpt.NClass*EstimOpt.NVarA;
else
   Bclass = B(1:EstimOpt.NVarA)*ones(1,EstimOpt.NClass);
   for i = 1:(EstimOpt.NClass- 1)
       Bclass(EstimOpt.BActiveClass == 1,i+1) = B(EstimOpt.NVarA+(i-1)*sum(EstimOpt.BActiveClass,1)+1:EstimOpt.NVarA+i*sum(EstimOpt.BActiveClass,1));
   end
   l = (EstimOpt.NClass-1)*sum(EstimOpt.BActiveClass,1)+EstimOpt.NVarA;
end
if EstimOpt.WTP_space > 0
    Bclass(1:end-EstimOpt.WTP_space,:) = Bclass(1:end-EstimOpt.WTP_space,:).*Bclass(EstimOpt.WTP_matrix,:);
end
if EstimOpt.NVarS > 0
    bs = reshape(B(l+1:l+EstimOpt.NClass*EstimOpt.NVarS),[EstimOpt.NVarS,EstimOpt.NClass]);
    Scale = exp(Xs*bs);
    U = exp(reshape((Xa*Bclass).*Scale,[EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP*EstimOpt.NClass])); % NAlt*NCT*NP x NClass
    l = l+EstimOpt.NClass*EstimOpt.NVarS;
else
    U = exp(reshape(Xa*Bclass,[EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP*EstimOpt.NClass])); % NAlt*NCT*NP x NClass
end
U(MissingInd == 1) = 0;% do not include alternatives which were not available

P = reshape(sum(YY.*U./sum(U,1),1),[EstimOpt.NCT,EstimOpt.NP*EstimOpt.NClass]); % NCT x NP*NClass
P(reshape(MissingInd(1:EstimOpt.NAlt:end),[EstimOpt.NCT,EstimOpt.NP*EstimOpt.NClass]) == 1) = 1; % do not include choice tasks which were not completed
probs = prod(P,1);
probs = reshape(probs,[EstimOpt.NP,EstimOpt.NClass]);

Pclass = exp(Xc*reshape([B(l+1:end);zeros(EstimOpt.NVarC,1)],[EstimOpt.NVarC,EstimOpt.NClass]));
Pclass = Pclass./sum(Pclass,2);
Scores = probs.*Pclass./sum(probs.*Pclass,2);

end