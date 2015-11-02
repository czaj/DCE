function C = sdtri(B1, VARB1, EstimOpt)

m1 = length(B1);
[m2,m3] = size(VARB1);

if m1 ~= m2 || m1 ~= m3 || m2 ~= m3
    error('Dimentions of B, and VARB must agree.');
end;

m0 = (sqrt(8*m1 + 1) - 1)/2;

D = zeros(EstimOpt.NSdSim,m0);

PRM = mvnrnd(B1,VARB1,EstimOpt.NSdSim);
if isfield(EstimOpt,'BActive') == 1 && sum(EstimOpt.BActive(EstimOpt.NVarA+1:end)) > 0
    PRM(:,EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA))'==0) = B1(EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA))==0,ones(EstimOpt.NSdSim,1))';
end
if sum(EstimOpt.Dist(2:end) >= 3 & EstimOpt.Dist(2:end) <= 5) > 0
    PRM(:, EstimOpt.DiagIndex(EstimOpt.Dist(2:end) >= 3 & EstimOpt.Dist(2:end) <= 5)) = 1;
end

Vt = tril(ones(m0));
Vt_t = Vt;
for i = 1:EstimOpt.NSdSim
    Vt_t(Vt==1) = PRM(i,:); ...
    if sum(EstimOpt.Dist(2:end) >= 3 & EstimOpt.Dist(2:end) <= 5) > 0
        tmp = sqrt(sum(Vt_t(EstimOpt.Dist(2:end) >= 3 & EstimOpt.Dist(2:end) <= 5,:).^2,2));
        Vt_t(EstimOpt.Dist(2:end) >= 3 & EstimOpt.Dist(2:end) <= 5,:) = Vt_t(EstimOpt.Dist(2:end) >= 3 & EstimOpt.Dist(2:end) <= 5,:)./tmp(:, ones(1,m0));
    end
    D(i,:) = sqrt(diag(Vt_t*Vt_t'))';
end

C = [mean(D); std(D)];
C = [C; pv(mean(D),std(D))]';


