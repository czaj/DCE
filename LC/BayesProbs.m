function Scores = BayesProbs(YY,Xa,Xc,MissingInd,EstimOpt,B)

% save tmp_LL_lc
% return

if sum(EstimOpt.BActiveClass == 0,1) == 0
    Bclass = reshape(B(1:EstimOpt.NClass*EstimOpt.NVarA), EstimOpt.NVarA, EstimOpt.NClass);
else
   Bclass = B(1:EstimOpt.NVarA)*ones(1,EstimOpt.NClass);
   for i = 1:(EstimOpt.NClass - 1)
       Bclass(EstimOpt.BActiveClass == 1,i+1) = B(EstimOpt.NVarA + (i-1)*sum(EstimOpt.BActiveClass,1)+1:EstimOpt.NVarA + i*sum(EstimOpt.BActiveClass,1));
   end
end
if EstimOpt.WTP_space > 0
    Bclass(1:end - EstimOpt.WTP_space,:) = Bclass(1:end - EstimOpt.WTP_space,:).*Bclass(EstimOpt.WTP_matrix,:);
end
U = exp(reshape(Xa * Bclass,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP*EstimOpt.NClass)); ...% NAlt*NCT*NP x NClass
U(MissingInd==1) = 0;... % do not include alternatives which were not available
U_sum = sum(U,1);... 
P = reshape(sum(YY .* U ./ U_sum(ones(EstimOpt.NAlt,1),:),1),EstimOpt.NCT,EstimOpt.NP*EstimOpt.NClass);... % NCT x NP*NClass
P(reshape(MissingInd(1:EstimOpt.NAlt:end),EstimOpt.NCT,EstimOpt.NP*EstimOpt.NClass)==1) = 1;... % do not include choice tasks which were not completed
probs = prod(P,1); ...
probs = reshape(probs,EstimOpt.NP,EstimOpt.NClass); ...
if sum(EstimOpt.BActiveClass == 0,1) == 0
    Pclass = exp(Xc*reshape([B(EstimOpt.NClass*EstimOpt.NVarA+1:end);zeros(EstimOpt.NVarC,1)],EstimOpt.NVarC, EstimOpt.NClass));...
else
    Pclass = exp(Xc*reshape([B((EstimOpt.NClass-1)*sum(EstimOpt.BActiveClass,1)+EstimOpt.NVarA+1:end);zeros(EstimOpt.NVarC,1)],EstimOpt.NVarC, EstimOpt.NClass));...
end
Pclass_sum = sum(Pclass,2);...
Pclass = Pclass./Pclass_sum(:,ones(EstimOpt.NClass,1)); ...

fx = sum(probs.*Pclass,2);

Scores = probs.*Pclass./fx(:,ones(EstimOpt.NClass,1));
% f = -log(max(sum(probs.*Pclass,2),realmin));

end