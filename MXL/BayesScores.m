function Score = BayesScores(YY,XXa,XXm,Xs,err,EstimOpt,b0)

% save LL_mxl
% return

NAlt = EstimOpt.NAlt;
NCT = EstimOpt.NCT;
NP = EstimOpt.NP;
NRep = EstimOpt.NRep;
NVarA = EstimOpt.NVarA;
NVarM = EstimOpt.NVarM;
NVarS = EstimOpt.NVarS;
Dist = EstimOpt.Dist;
WTP_space = EstimOpt.WTP_space;
WTP_matrix = EstimOpt.WTP_matrix;
DiagIndex = EstimOpt.DiagIndex;
Triang = EstimOpt.Triang;
Johnson = EstimOpt.Johnson;

b0a = b0(1:NVarA);
if any(Dist == 3)
    b0triag_c = exp(b0a(Dist == 3)) + Triang';
    b0a(Dist == 3) = 0;
end
if any(Dist == 4)
    b0weibA = exp(b0a(Dist == 4));
    b0a(Dist == 4) = 0;
end
if any(Dist == 5)
    b0sinhA = b0a(Dist == 5);
    b0a(Dist == 5) = 0;
end

if any(Dist == 9)   % Uni-Log distribution
    b0UniLogA = exp(b0a(Dist == 9));
    b0a(Dist == 9) = 0;
end

if any(Dist == 10)   % Pareto distribution
    b0ParetoA = exp(b0a(Dist == 10));
    b0a(Dist == 10) = 0;
end

if any(Dist == 11)   % Lomax distribution
    b0LomaxA = exp(b0a(Dist == 11));
    b0a(Dist == 11) = 0;
end

if any(Dist == 12)   % Logistic distribution
    b0LogA = exp(b0a(Dist == 12));
    b0a(Dist == 12) = 0;
end

if any(Dist == 13)   % Log-Logistic distribution
    b0LogLA = exp(b0a(Dist == 13));
    b0a(Dist == 13) = 0;
end

if any(Dist == 14)   % Gumbel distribution
    b0GumA = exp(b0a(Dist == 14));
    b0a(Dist == 14) = 0;
end

if any(Dist == 15)   % Cauchy distribution
    b0CaucA = exp(b0a(Dist == 15));
    b0a(Dist == 15) = 0;
end

if any(Dist == 16)   % Rayleigh distribution
    b0RayA = exp(b0a(Dist == 16));
    b0a(Dist == 16) = 0;
end

if any(Dist == 17)   % Exponential distribution
    b0ExpA = exp(b0a(Dist == 17));
    b0a(Dist == 17) = 0;
end


if EstimOpt.FullCov == 0
    b0v = (b0(NVarA+1:NVarA*2));
    if any(Dist == 3)
        b0triag_b = exp(b0v(Dist == 3)) + b0triag_c;
        b0v(Dist == 3) = 1;
    end
    if any(Dist == 4)
        b0weibB = exp(-b0v(Dist == 4));
        b0v(Dist == 4) = 1;
    end
    if any(Dist == 5)
        b0sinhB = b0v(Dist == 5).^2;
        b0v(Dist == 5) = 1;
    end
    
    if any(Dist == 9)   % Uni-Log
        b0UniLogB = exp(b0v(Dist == 9));
        b0v(Dist == 9) = 1;
    end
    if any(Dist == 10)   % Pareto
        b0ParetoB = exp(b0v(Dist == 10));
        b0v(Dist == 10) = 1;
    end
    if any(Dist == 11)   % Lomax
        b0LomaxB = exp(b0v(Dist == 11));
        b0v(Dist == 11) = 1;
    end
    if any(Dist == 12)   % Logistic
        b0LogB = exp(b0v(Dist == 12));
        b0v(Dist == 12) = 1;
    end    
    if any(Dist == 13)   % Log-Logistic
        b0LogLB = exp(b0v(Dist == 13));
        b0v(Dist == 13) = 1;
    end     
    if any(Dist == 14)   % Gumbel
        b0GumB = exp(b0v(Dist == 14));
        b0v(Dist == 14) = 1;
    end     
    if any(Dist == 15)   % Cauchy
        b0CaucB = exp(b0v(Dist == 15));
        b0v(Dist == 15) = 1;
    end 
    if any(Dist == 16)   % Rayleigh
        % b0RayB = exp(b0v(Dist == 16));
        b0v(Dist == 16) = 1;
    end
    if any(Dist == 17)   % Exponential
        % b0ExpB = exp(b0v(Dist == 17));
        b0v(Dist == 17) = 1;
    end    
    
    b0v = b0v.^2;
    VC = diag(b0v);
    b0m = b0(NVarA*2+1:NVarA*(NVarM+2));
    b0m = reshape(b0m,[NVarA,NVarM]);
    b0s = b0(NVarA*(NVarM+2)+1:NVarA*(NVarM+2)+NVarS);
    b0j = b0(NVarA*(NVarM+2)+NVarS+1:NVarA*(NVarM+2)+NVarS+2*Johnson);
else
    b0v = b0(NVarA+1:NVarA+sum(1:NVarA));
    tmp = b0v(DiagIndex);
    b0v(DiagIndex(Dist >=3)) = 1;
    if any(Dist == 3)
        b0triag_b = exp(tmp(Dist == 3))+b0triag_c;
    end
    if any(Dist == 4)
        b0weibB = exp(-tmp(Dist == 4));
    end
    if any(Dist == 5)
        b0sinhB = tmp(Dist == 5).^2;
    end
    
    if any(Dist == 9)
        b0UniLogB = exp(tmp(Dist == 9));
    end
    if any(Dist == 10)
        b0ParetoB = exp(tmp(Dist == 10));
    end
    if any(Dist == 11)
        b0LomaxB = exp(tmp(Dist == 11));
    end
    if any(Dist == 12)
        b0LogB = exp(tmp(Dist == 12));
    end
    if any(Dist == 13)
        b0LogLB = exp(tmp(Dist == 13));
    end
    if any(Dist == 14)
        b0GumB = exp(tmp(Dist == 14));
    end     
    if any(Dist == 15)
        b0CaucB = exp(tmp(Dist == 15));
    end
    if any(Dist == 16)
        % b0RayB = exp(tmp(Dist == 16));
    end
    if any(Dist == 17)
        % b0ExpB = exp(tmp(Dist == 17));
    end
    
    VC = tril(ones(NVarA));
    VC(VC==1) = b0v;
    
    if any(Dist >= 3 & Dist <= 5 | Dist == 9 | Dist == 10 | Dist == 11 | Dist == 12 | Dist == 13 | Dist == 14 | Dist == 15 | Dist == 16 | Dist == 17)
        tmp = sqrt(sum(VC(Dist >= 3 & Dist <= 5 | Dist == 9 | Dist == 10 | Dist == 11 | Dist == 12 | Dist == 13 | Dist == 14 | Dist == 15 | Dist == 16 | Dist == 17,:).^2,2));
        VC(Dist >= 3 & Dist <= 5 | Dist == 9 | Dist == 10 | Dist == 11 | Dist == 12 | Dist == 13 | Dist == 14 | Dist == 15 | Dist == 16 | Dist == 17,:) = VC(Dist >= 3 & Dist <= 5 | Dist == 9 | Dist == 10 | Dist == 11 | Dist == 12 | Dist == 13 | Dist == 14 | Dist == 15 | Dist == 16 | Dist == 17,:)./tmp(:,ones(1,NVarA));
    end
    
    b0m = b0(NVarA*(NVarA/2+1.5)+1:NVarA*(NVarA/2+1.5+NVarM));
    b0m = reshape(b0m,[NVarA,NVarM]);
    b0s = b0(NVarA*(NVarA/2+1.5+NVarM)+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS);
    b0j = b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS+2*Johnson);
end

b0n = b0a + b0m*XXm;
b0n = reshape(b0n((1:size(b0n,1))'*ones(1,NRep),(1:size(b0n,2))'),[NVarA,NRep*NP]); % NVarA x NRep*NP

b_mtx_n = b0n + VC*err; % NVarA x NRep*NP

if sum(Dist == 1) > 0 % Log - normal
    b_mtx_n(Dist == 1,:) = exp(b_mtx_n(Dist == 1,:));
end
if sum(Dist == 2) > 0 % Spike
    b_mtx_n(Dist == 2,:) = max(b_mtx_n(Dist == 2,:),0);
end
if sum(Dist == 3) > 0 % Triangular
    tmp = normcdf(b_mtx_n(Dist == 3,:));
    Triang = Triang(ones(NRep*NP,1),:)';
    b0triag_c = b0triag_c(:,ones(NRep*NP,1));
    b0triag_b = b0triag_b(:,ones(NRep*NP,1));
    Ftriang =  (b0triag_c - Triang)./(b0triag_b - Triang);
    bmtx_triang = zeros(size(tmp));
    tmp2 = (b0triag_b - Triang).*(b0triag_c - Triang);
    bmtx_triang(tmp < Ftriang) = Triang(tmp < Ftriang) + sqrt(tmp(tmp < Ftriang).*tmp2(tmp < Ftriang));
    tmp2 = (b0triag_b - Triang).*(b0triag_b - b0triag_c);
    bmtx_triang(tmp >= Ftriang) = b0triag_b(tmp >= Ftriang)- sqrt((1-tmp(tmp >= Ftriang)).*tmp2(tmp >= Ftriang));
    b_mtx_n(Dist == 3,:) = bmtx_triang;
end
if sum(Dist == 4) > 0 % Weibull
    tmp = -log(1 - normcdf(b_mtx_n(Dist == 4,:)));
    b_mtx_n(Dist == 4,:) = b0weibA.*(tmp.^b0weibB);
end
if sum(Dist >= 5) > 0 % Johnson
    if sum(Dist == 5) > 0 % Sinh-Arcsinh
        b_mtx_n(Dist == 5,:) = b0sinhA + b0sinhB.*asinh(b_mtx_n(Dist == 5,:));
        b_mtx_n(Dist == 5,:) = b0j(1:Johnson,ones(NRep*NP,1)) + exp(b0j(Johnson+1:end,ones(NRep*NP,1))).*sinh(b_mtx_n(Dist == 5,:));
    end
    if sum(Dist == 6) > 0 % Johnson Sb
        tmp = exp(b_mtx_n(Dist == 6,:));
        b_mtx_n(Dist == 6,:) = tmp./(1+tmp);
        b_mtx_n(Dist == 6,:) = b0j(1:Johnson,ones(NRep*NP,1)) + exp(b0j(Johnson+1:end,ones(NRep*NP,1))).*b_mtx_n(Dist == 6,:);
    end
    if sum(Dist == 7) > 0 % Johnson Su
        b_mtx_n(Dist == 7,:) = sinh(b_mtx_n(Dist==7,:));
        b_mtx_n(Dist == 7,:) = b0j(1:Johnson,ones(NRep*NP,1)) + exp(b0j(Johnson+1:end,ones(NRep*NP,1))).*b_mtx_n(Dist == 7,:);
    end
end

if sum(Dist == 9) > 0 % Uni-Log
    tmp = normcdf(b_mtx(Dist == 9,:));
    b_mtx_n(Dist == 9,:) = exp(log(b0UniLogA)+tmp.*(log(b0UniLogB)-log(b0UniLogA)));   % inverse CDF function
end

if sum(Dist == 10) > 0 % Pareto
    tmp = normcdf(b_mtx(Dist == 10,:));
    b_mtx_n(Dist == 10,:) = b0ParetoA./((1-tmp).^(1./b0ParetoB));   % inverse CDF function
end

if sum(Dist == 11) > 0 % Lomax
    tmp = normcdf(b_mtx(Dist == 11,:));
    b_mtx_n(Dist == 11,:) = b0LomaxB.*(((1./(1-tmp)).^(1./b0LomaxA))-1);   % inverse CDF function
end

if sum(Dist == 12) > 0 % Logistic
    tmp = normcdf(b_mtx(Dist == 12,:));
    b_mtx_n(Dist == 12,:) = b0LogA+b0LogB.*log(tmp./(1-tmp));   % inverse CDF function
end

if sum(Dist == 13) > 0 % Log-Logistic
    tmp = normcdf(b_mtx(Dist == 13,:));
    b_mtx_n(Dist == 13,:) = b0LogLA.*(1./((1./tmp)-1)).^(1./b0LogLB);   % inverse CDF function
end

if sum(Dist == 14) > 0 % Gumbel
    tmp = normcdf(b_mtx(Dist == 14,:));
    b_mtx_n(Dist == 14,:) = b0GumA+b0GumB.*log(1./log(1./tmp));   % inverse CDF function
end

if sum(Dist == 15) > 0 % Cauchy
    tmp = normcdf(b_mtx(Dist == 15,:));
    b_mtx_n(Dist == 15,:) = b0CaucA+b0CaucB.*tan(pi.*(tmp-1/2));   % inverse CDF function
end

if sum(Dist == 16) > 0 % Rayleigh
    tmp = normcdf(b_mtx(Dist == 16,:));
    b_mtx_n(Dist == 16,:) = (2.*(b0RayA.^2).*log(1./(1-tmp))).^(1/2);   % inverse CDF function
end

if sum(Dist == 17) > 0 % Exponential
    tmp = normcdf(b_mtx(Dist == 17,:));
    b_mtx_n(Dist == 17,:) = (log(1./(1-tmp))).*(1./b0ExpA);   % inverse CDF function
end

b_score = b_mtx_n;
if WTP_space > 0
    b_mtx_n(1:end-WTP_space,:) = b_mtx_n(1:end-WTP_space,:).*b_mtx_n(WTP_matrix,:);
end

cs = reshape(exp(Xs*b0s),[NAlt*NCT,1,NP]);
XXa_n = XXa.*cs;

b_mtx_n = reshape(b_mtx_n,[NVarA,NRep,NP]);
p0 = zeros(NP,NRep);

if any(isnan(XXa(:))) == 0 % faster version for complete dataset
    for n = 1:NP
        U = reshape(XXa_n(:,:,n)*b_mtx_n(:,:,n),[NAlt,NCT,NRep]);
        U = exp(U - max(U,[],1)); % rescale utility to avoid exploding
        U_sum = reshape(sum(U,1),[NCT,NRep]);
        U_selected = reshape(U(YY(:,n*ones(NRep,1)) == 1),[NCT,NRep]);
        p0(n,:) = prod(U_selected./U_sum,1);
    end
else
    for n = 1:NP
        U = reshape(XXa_n(~isnan(YY(:,n)),:,n)*b_mtx_n(:,:,n),[NAlt,NCT-sum(isnan(YY(1:NAlt:end,n))),NRep]);
        U = exp(U - max(U,[],1));
        U_sum = reshape(nansum(U,1),[NCT-sum(isnan(YY(1:NAlt:end,n))),NRep]);
        U_selected = reshape(U(YY(~isnan(YY(:,n)),n*ones(NRep,1)) == 1),[NCT-sum(isnan(YY(1:NAlt:end,n))),NRep]);
        p0(n,:) = prod(U_selected./U_sum,1);
    end
end

fx = mean(p0,2); % NP x 1
Score = zeros(NP,NVarA);
if WTP_space > 0
    b_score = reshape(b_score,[NVarA,NRep,NP]);
    for j = 1:NVarA
        bx = squeeze(b_score(j,:,:))'; % NP x NRep
        Score(:,j) = mean(p0.*bx,2)./fx;
    end
else
    b_score = reshape(b_score,[NVarA,NRep,NP]);
    fee = squeeze(b_score(end,:,:))';
    for j = 1:NVarA-1
        bx = squeeze(b_score(j,:,:))'; % NP x NRep
        Score(:,j) = mean(p0.*bx./fee,2)./fx;
    end
    Score(:,end) = mean(p0.*fee,2)./fx;
end