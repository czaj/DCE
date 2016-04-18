function f = sdtriHe2(B1,EstimOpt, Corr)

VC = tril(ones(EstimOpt.NLatent+sum(EstimOpt.Dist ~= -1)));
VCtmp = diag(ones(EstimOpt.NLatent+sum(EstimOpt.Dist ~= -1),1));
VCtmp(1:sum(EstimOpt.Dist ~= -1), 1:sum(EstimOpt.Dist ~= -1)) = 0;
VC(VCtmp == 1) = 0;
VCX = VC;
VC(VC == 1)= B1;
VC(VCtmp == 1) = 1;
VCtmp2 = sqrt(sum(VC(sum(EstimOpt.Dist ~= -1)+1:end,:).^2,2));
VC(sum(EstimOpt.Dist ~= -1)+1:end,:) = VC(sum(EstimOpt.Dist ~= -1)+1:end,:)./VCtmp2(:, ones(1, EstimOpt.NLatent+sum(EstimOpt.Dist ~= -1)));

COV = VC*VC';
if Corr == 0
    f = COV(VCX == 1);
else
    VCtmp = diag(ones(EstimOpt.NLatent+sum(EstimOpt.Dist ~= -1),1));
    VCtmp(sum(EstimOpt.Dist ~= -1)+1:end, sum(EstimOpt.Dist ~= -1)+1:end) = 0;
    VCX(VCtmp == 1) = 0;
    COV = corrcov(COV);
    f = COV(VCX == 1);
    
end





