function f = sdtriHe2(B1,EstimOpt,Corr)

if Corr == -1
    EstimOpt.NLatent = 0;
    BActive = EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2); 
    B1 = B1(BActive);
end    
VC = tril(ones(EstimOpt.NLatent+sum(EstimOpt.Dist ~= -1)));
VCtmp = diag(ones(EstimOpt.NLatent+sum(EstimOpt.Dist ~= -1),1));
VCtmp(1:sum(EstimOpt.Dist ~= -1),1:sum(EstimOpt.Dist ~= -1)) = 0;
VC(VCtmp == 1) = 0;
VCX = VC;
VC(VC == 1) = B1;
VC(VCtmp == 1) = 1;
if Corr > -1
    VCtmp2 = sqrt(sum(VC(sum(EstimOpt.Dist ~= -1)+1:end,:).^2,2));
    VC(sum(EstimOpt.Dist ~= -1)+1:end,:) = VC(sum(EstimOpt.Dist ~= -1)+1:end,:)./VCtmp2(:,ones(1, EstimOpt.NLatent+sum(EstimOpt.Dist ~= -1)));
end
COV = VC*VC';
if Corr == 0
    f = COV(VCX == 1);
elseif Corr == 1
    VCtmp = diag(ones(EstimOpt.NLatent+sum(EstimOpt.Dist ~= -1),1));
    VCtmp(sum(EstimOpt.Dist ~= -1)+1:end,sum(EstimOpt.Dist ~= -1)+1:end) = 0;
    VCX(VCtmp == 1) = 0;
    COV = corrcov(COV);
    f = COV(VCX == 1);
elseif Corr == 2 % only standard deviations
    f = sqrt(diag(COV(1:sum(EstimOpt.Dist ~= -1),1:sum(EstimOpt.Dist ~= -1))));
elseif Corr == -1
    f = sqrt(diag(COV));
end





