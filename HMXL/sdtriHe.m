function [COV, CORR] = sdtriHe(B1, VARB1, EstimOpt)

EstimOpt.NVarA = EstimOpt.NVarA - sum(EstimOpt.Dist == -1);
m1 = length(B1);

[m2,m3] = size(VARB1);

if m1 ~= m2 || m1 ~= m3 || m2 ~= m3
    error('Dimentions of B, and VARB must agree.');
end;

%m0 = (sqrt(8*m1 + 1) - 1)/2;

%D = zeros(EstimOpt.NSdSim,m0);

PRM = mvnrnd(B1,VARB1,EstimOpt.NSdSim);
% if isfield(EstimOpt,'BActive') == 1 && sum(EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA + sum(1:EstimOpt.NVarA+EstimOpt.NLatent))) > 0
%     PRM(:,EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent)'==0) = B1(EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent)==0,ones(EstimOpt.NSdSim,1))';
% end

PRMX = zeros(EstimOpt.NSdSim, sum(1:EstimOpt.NVarA+EstimOpt.NLatent));
VCtmp = diag(ones(EstimOpt.NVarA+EstimOpt.NLatent,1));
VCtmp(1:EstimOpt.NVarA, 1:EstimOpt.NVarA) = 0;

Indx = VCtmp(tril(ones(EstimOpt.NVarA+EstimOpt.NLatent))==1);
PRMX(:, Indx == 0) = PRM;
PRMX(:, Indx == 1) = 1;

Vt = tril(ones(EstimOpt.NVarA+EstimOpt.NLatent));
Vt_t = Vt;

COVX = zeros(EstimOpt.NVarA+EstimOpt.NLatent,EstimOpt.NVarA+EstimOpt.NLatent,EstimOpt.NSdSim);
CORRX = zeros(EstimOpt.NVarA+EstimOpt.NLatent,EstimOpt.NVarA+EstimOpt.NLatent,EstimOpt.NSdSim);

for i = 1:EstimOpt.NSdSim
    Vt_t(Vt==1) = PRMX(i,:); ...
    tmp = sqrt(sum(Vt_t(EstimOpt.NVarA+1:end,:).^2,2));
    Vt_t(EstimOpt.NVarA+1:end,:) = Vt_t(EstimOpt.NVarA+1:end,:)./tmp(:, ones(1,EstimOpt.NVarA+EstimOpt.NLatent));
    StdX = sqrt(diag(Vt_t*Vt_t'))';
    A = (StdX')*StdX;
    
    COVX(:,:,i) = Vt_t*Vt_t';
    CORRX(:,:,i) = COVX(:,:,i)./A;
end

COV = zeros(EstimOpt.NVarA+EstimOpt.NLatent,EstimOpt.NVarA+EstimOpt.NLatent,3);
COV(:,:,1) = mean(COVX,3);
COV(:,:,2) = prctile(COVX,2.5,3);
COV(:,:,3) = prctile(COVX,97.5,3);

CORR = zeros(EstimOpt.NVarA+EstimOpt.NLatent,EstimOpt.NVarA+EstimOpt.NLatent,3);
CORR(:,:,1) = mean(CORRX,3);
CORR(:,:,2) = prctile(CORRX,2.5,3);
CORR(:,:,3) = prctile(CORRX,97.5,3);


