function Score = BayesScores(YY,XXa,XXm,Xs,err,EstimOpt,b0)

% save LL_mxl
% return

% NAlt = EstimOpt.NAlt; ...
% NCT = EstimOpt.NCT; ...
% NP = EstimOpt.NP; ...
% NRep = EstimOpt.NRep; ...
% NVarA = EstimOpt.NVarA; ...
% NVarM = EstimOpt.NVarM; ...
% NVarS = EstimOpt.NVarS; ...
% Dist = EstimOpt.Dist; ...
% WTP_space = EstimOpt.WTP_space; ...
% WTP_matrix = EstimOpt.WTP_matrix; ...
% DiagIndex = EstimOpt.DiagIndex; ...
% Triang = EstimOpt.Triang; ...
% Johnson = EstimOpt.Johnson; ...
% 
% b0a = b0(1:NVarA); ...
% if any(Dist(2:end) == 3)
%     b0triag_c = exp(b0a(Dist(2:end) == 3)) + Triang';
%     b0a(Dist(2:end) == 3) = 0;
% end
% if any(Dist(2:end) == 4)
%     b0weibA = exp(b0a(Dist(2:end) == 4));
%     b0a(Dist(2:end) == 4) = 0;
% end
% if any(Dist(2:end) == 5)
%     b0sinhA = b0a(Dist(2:end) == 5);
%     b0a(Dist(2:end) == 5) = 0;
% end
% 
% if EstimOpt.FullCov == 0 ...
%     b0v = (b0(NVarA+1:NVarA*2)); ...
%     if any(Dist(2:end) == 3)
%         b0triag_b = exp(b0v(Dist(2:end) == 3)) + b0triag_c;
%         b0v(Dist(2:end) == 3) = 1;
%     end
%     if any(Dist(2:end) == 4)
%         b0weibB = exp(-b0v(Dist(2:end) == 4));
%         b0v(Dist(2:end) == 4) = 1;
%     end
%     if any(Dist(2:end) == 5)
%         b0sinhB = b0v(Dist(2:end) == 5).^2;
%         b0v(Dist(2:end) == 5) = 1;
%     end
%     b0v = b0v.^2;
%     VC = diag(b0v); ...
%     b0m = b0(NVarA*2+1:NVarA*(NVarM+2)); ...    
%     b0m = reshape(b0m,NVarA, NVarM); ...
%     b0s = b0(NVarA*(NVarM+2)+1:NVarA*(NVarM+2)+NVarS); ...
%     b0j = b0(NVarA*(NVarM+2)+NVarS+1:NVarA*(NVarM+2)+NVarS+2*Johnson); ...
% 
% else ...
% 	b0v = b0(NVarA+1:NVarA+sum(1:NVarA)); ...
%     tmp = b0v(DiagIndex);
%     b0v(DiagIndex(Dist(2:end) >=3)) = 1;
%     if any(Dist(2:end) == 3)
%         b0triag_b = exp(tmp(Dist(2:end) == 3))+b0triag_c;
%     end    
%     if any(Dist(2:end) == 4)
%         b0weibB = exp(-tmp(Dist(2:end) == 4));
%     end    
%     if any(Dist(2:end) == 5)
%         b0sinhB = tmp(Dist(2:end) == 5).^2;
%     end  
%     VC = tril(ones(NVarA)); ...
%     VC(VC==1) = b0v; ...
%     if any(Dist(2:end) >= 3 & Dist(2:end) <= 5)
%         tmp = sqrt(sum(VC(Dist(2:end) >= 3 & Dist(2:end) <= 5,:).^2,2));
%         VC(Dist(2:end) >= 3 & Dist(2:end) <= 5,:) = VC(Dist(2:end) >= 3 & Dist(2:end) <= 5,:)./tmp(:, ones(1,NVarA));
%     end
%     b0m = b0(NVarA*(NVarA/2+1.5)+1:NVarA*(NVarA/2+1.5+NVarM)); ...
%     b0m = reshape(b0m,NVarA, NVarM); ...
%     b0s = b0(NVarA*(NVarA/2+1.5+NVarM)+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS); ...
%     b0j = b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS+2*Johnson); ...
% end ...
% 
% b0n = b0a(:,ones(NP,1)) + b0m*XXm; ...
% b0n = reshape(b0n([1:size(b0n,1)]'*ones(1,NRep),[1:size(b0n,2)]'),NVarA,NRep*NP); ... % NVarA x NRep*NP
% 
% b_mtx_n = b0n + VC*err; ... % NVarA x NRep*NP
% 
% if sum(Dist(2:end)==1) > 0 % Log - normal
%     b_mtx_n(Dist(2:end)==1,:) = exp(b_mtx_n(Dist(2:end)==1,:)); ...
% end
% if sum(Dist(2:end)==2) > 0 % Spike
%     b_mtx_n(Dist(2:end)==2,:) = max(b_mtx_n(Dist(2:end)==2,:),0); ...
% end
% if sum(Dist(2:end) ==3) > 0 % Triangular
%     tmp = normcdf(b_mtx_n(Dist(2:end)==3,:)); ...
%     Triang = Triang(ones(NRep*NP,1),:)';
%     b0triag_c = b0triag_c(:, ones(NRep*NP,1));
%     b0triag_b = b0triag_b(:, ones(NRep*NP,1));
%     Ftriang =  (b0triag_c - Triang)./(b0triag_b- Triang);
%     bmtx_triang = zeros(size(tmp));
%     tmp2 = (b0triag_b- Triang).*(b0triag_c - Triang);
%     bmtx_triang(tmp < Ftriang) = Triang(tmp < Ftriang)+ sqrt(tmp(tmp < Ftriang).*tmp2(tmp < Ftriang));
%     tmp2 = (b0triag_b- Triang).*(b0triag_b-b0triag_c);
%     bmtx_triang(tmp >= Ftriang) = b0triag_b(tmp >= Ftriang)- sqrt((1-tmp(tmp >= Ftriang)).*tmp2(tmp >= Ftriang));  
%     b_mtx_n(Dist(2:end)==3,:) = bmtx_triang;
% end
% if sum(Dist(2:end) ==4) > 0 % Weibull
%     tmp = -log(1-normcdf(b_mtx_n(Dist(2:end)==4,:))); ...
%     b_mtx_n(Dist(2:end)==4,:) = b0weibA(:, ones(1,NP*NRep,1)).*(tmp.^b0weibB(:, ones(1,NP*NRep,1)));
% end
% if sum(Dist(2:end)>=5) > 0 % Johnson
%     if sum(Dist(2:end)==5) > 0 % Sinh-Arcsinh
%         b_mtx_n(Dist(2:end)==5,:) = b0sinhA(:,ones(NRep*NP,1))+ b0sinhB(:,ones(NRep*NP,1)).*asinh(b_mtx_n(Dist(2:end)==5,:));
%         b_mtx_n(Dist(2:end)==5,:) = b0j(1:Johnson, ones(NRep*NP,1)) + exp(b0j(Johnson+1:end, ones(NRep*NP,1))).*sinh(b_mtx_n(Dist(2:end)==5,:));  
%     end
%     if sum(Dist(2:end)==6) > 0 % Johnson Sb
%         tmp = exp(b_mtx_n(Dist(2:end) ==6,:));
%         b_mtx_n(Dist(2:end)==6,:) = tmp./(1+tmp); ...
%         b_mtx_n(Dist(2:end) == 6,:) =b0j(1:Johnson, ones(NRep*NP,1)) + exp(b0j(Johnson+1:end, ones(NRep*NP,1))).*b_mtx_n(Dist(2:end)==6,:);
%     end
%     if sum(Dist(2:end)==7) > 0 % Johnson Su
%         b_mtx_n(Dist(2:end)==7,:) = sinh(b_mtx_n(Dist(2:end)==7,:)); ...
%         b_mtx_n(Dist(2:end) ==7,:) =b0j(1:Johnson, ones(NRep*NP,1)) + exp(b0j(Johnson+1:end, ones(NRep*NP,1))).*b_mtx_n(Dist(2:end)==7,:);
%     end
% 
% end
% 
% 
% b_score = b_mtx_n;
% if WTP_space > 0
%     b_mtx_n(1:end-WTP_space,:) = b_mtx_n(1:end-WTP_space,:).*b_mtx_n(WTP_matrix,:); ...
% end
% 
% 
% cs = reshape(exp(Xs*b0s),NAlt*NCT,1,NP); ...
% XXa_n = XXa .* cs(:,ones(1,NVarA,1),:); ...
% 
% b_mtx_n = reshape(b_mtx_n,NVarA,NRep,NP); ...
% 
% p0 = zeros(NP,NRep); ...

% if any(isnan(XXa(:))) == 0 ... % faster version for complete dataset
% 
%     for n = 1:NP ...
%         U = reshape(XXa_n(:,:,n)*b_mtx_n(:,:,n),NAlt,NCT,NRep); ...
%         U_max = max(U); ...
%         U = exp(U - U_max(ones(NAlt,1),:,:)); ... % rescale utility to avoid exploding
%         U_sum = reshape(sum(U,1),NCT,NRep); ...
%         U_selected = reshape(U(YY(:,n*ones(NRep,1))==1),NCT,NRep); ...   
%         p0(n,:) = prod(U_selected ./ U_sum,1);
%     end; ...
% else ...
%     for n = 1:NP ...
%         save tmp1
%         U = reshape(XXa_n(~isnan(YY(:,n)),:,n)*b_mtx_n(:,:,n),NAlt,NCT-sum(isnan(YY(1:NAlt:end,n))),NRep); ...
%         U_max = max(U); ...
%         U = exp(U - U_max(ones(NAlt,1),:,:)); ...
%         U_sum = reshape(nansum(U,1),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep); ...
%         U_selected = reshape(U(YY(~isnan(YY(:,n)),n*ones(NRep,1))==1),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep); ...   
%         p0(n,:) = prod(U_selected ./ U_sum,1);
%     end; ...
% end ...


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
FullCov = EstimOpt.FullCov;
DiagIndex = EstimOpt.DiagIndex;
Triang = EstimOpt.Triang;
NVarNLT = EstimOpt.NVarNLT;
NLTVariables = EstimOpt.NLTVariables;
NLTType = EstimOpt.NLTType;
Johnson = EstimOpt.Johnson;
NCTMiss = EstimOpt.NCTMiss;
NAltMiss = EstimOpt.NAltMiss;
NAltMissInd = EstimOpt.NAltMissInd;
NAltMissIndExp = EstimOpt.NAltMissIndExp;
%NCTMissIndExp = EstimOpt.NCTMissIndExp;
MissingCT = EstimOpt.MissingCT;
RealMin = EstimOpt.RealMin;
ExpB = EstimOpt.ExpB;
mCT = EstimOpt.mCT; % = 1 if Xm is CT or Alt specific

% if nargout == 3
%    XXX = permute(mmx('square',permute(XXa,[2,4,1,3]),[]),[3,1,2,4]);
% %    VCx = EstimOpt.VCx;
% end
if FullCov == 1 && nargout > 1
    indx1 = EstimOpt.indx1;
    indx2 = EstimOpt.indx2;
else
    indx1 = [];
    indx2 = [];
end

if ~isempty(ExpB)
    b0(ExpB) = exp(b0(ExpB));
end

%% First parameters of the distributions
b0a = b0(1:NVarA);  % first parameters of the distributions

% Additional variables for distributions other than normal
% Introduction of new variables requires zeroing the corresponding elements of the parameters vector
if any(Dist == 3)
    b0triang_c = exp(b0a(Dist == 3)) + Triang';
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


%% Second parameters of the distributions (std devs as default)
if FullCov == 0
    b0v = (b0(NVarA+1:NVarA*2));
    if any(Dist == 3)
        b0triang_b = exp(b0v(Dist == 3)) + b0triang_c;
        b0v(Dist == 3) = 1;
    end
    if any(Dist == 4)
        b0weibB = exp(-b0v(Dist == 4));
        b0v(Dist == 4) = 1;
    else
        b0weibA = [];
        b0weibB = [];
    end
    if any(Dist == 5)
        b0sinhB = b0v(Dist == 5).^2;
        b0v(Dist == 5) = 1;
    end
    
    
    if any(Dist == 9)   % Uni-Log
        b0UniLogB = exp(b0v(Dist == 9));
        b0v(Dist == 9) = 1;
    else
        b0UniLogA = [];  % initialize variable for parfor loop
        b0UniLogB = []; % initialize variable for parfor loop
    end
        
    if any(Dist == 10)   % Pareto
        b0ParetoB = exp(b0v(Dist == 10));
        b0v(Dist == 10) = 1;
    else
        b0ParetoA = [];  % initialize variable for parfor loop
        b0ParetoB = []; % initialize variable for parfor loop
    end
    
    if any(Dist == 11)   % Lomax
        b0LomaxB = exp(b0v(Dist == 11));
        b0v(Dist == 11) = 1;
    else
        b0LomaxA = [];  % initialize variable for parfor loop
        b0LomaxB = []; % initialize variable for parfor loop
    end
    
    if any(Dist == 12)   % Logistic
        b0LogB = exp(b0v(Dist == 12));
        b0v(Dist == 12) = 1;
    else
        b0LogA = [];  % initialize variable for parfor loop
        b0LogB = []; % initialize variable for parfor loop
    end
        
    if any(Dist == 13)   % Log-Logistic
        b0LogLB = exp(b0v(Dist == 13));
        b0v(Dist == 13) = 1;
    else
        b0LogLA = [];  % initialize variable for parfor loop
        b0LogLB = []; % initialize variable for parfor loop
    end

    if any(Dist == 14)   % Gumbel
        b0GumB = exp(b0v(Dist == 14));
        b0v(Dist == 14) = 1;
    else
        b0GumA = [];  % initialize variable for parfor loop
        b0GumB = []; % initialize variable for parfor loop
    end
    
    if any(Dist == 15)   % Cauchy
        b0CaucB = exp(b0v(Dist == 15));
        b0v(Dist == 15) = 1;
    else
        b0CaucA = [];  % initialize variable for parfor loop
        b0CaucB = []; % initialize variable for parfor loop
    end
        
    if any(Dist == 16)   % Rayleigh
        % b0RayB = exp(b0v(Dist == 16));
        b0v(Dist == 16) = 1;
    else
        b0RayA = [];  % initialize variable for parfor loop
        % b0RayB = []; % initialize variable for parfor loop
    end
    
    if any(Dist == 17)   % Exponential
        % b0ExpB = exp(b0v(Dist == 17));
        b0v(Dist == 17) = 1;
    else
        b0ExpA = [];  % initialize variable for parfor loop
        % b0ExpB = []; % initialize variable for parfor loop
    end
    
    
    %     b0v = b0v.^2;
    VC = diag(b0v);
    b0m = b0(NVarA*2+1:NVarA*(NVarM+2));
    b0m = reshape(b0m,[NVarA,NVarM]);
    b0s = b0(NVarA*(NVarM+2)+1:NVarA*(NVarM+2)+NVarS);
    b0t = b0(NVarA*(NVarM+2)+NVarS+1:NVarA*(NVarM+2)+NVarS+NVarNLT);
    b0j = b0(NVarA*(NVarM+2)+NVarS+NVarNLT+1:NVarA*(NVarM+2)+NVarS+NVarNLT+2*Johnson);
else
    b0v = b0(NVarA+1:NVarA+sum(1:NVarA));
    tmp = b0v(DiagIndex);
    
    b0v(DiagIndex(Dist >=3 & Dist <=5 | Dist == 9 | Dist == 10 | Dist == 11 | Dist == 12 | Dist == 13 | Dist == 14 | Dist == 15 | Dist == 16 | Dist == 17)) = 1;    % triangular, weibull, sinh-arcsinh, uni-log(reciprocal), pareto, lomax, logistic, log-logistic, gumbel, cauchy, rayleigh, exponential distributions

    if any(Dist == 3)
        b0triang_b = exp(tmp(Dist == 3)) + b0triang_c;
    end
    if any(Dist == 4)
        b0weibB = exp(-tmp(Dist == 4));
    else
        b0weibA = [];
        b0weibB = [];        
    end
    if any(Dist == 5)
        b0sinhB = tmp(Dist == 5).^2;
    end
    
    if any(Dist == 9)
        b0UniLogB = exp(tmp(Dist == 9));
    else
        b0UniLogA = [];
        b0UniLogB = [];
    end
    if any(Dist == 10)
        b0ParetoB = exp(tmp(Dist == 10));
    else
        b0ParetoA = [];
        b0ParetoB = [];
    end
    if any(Dist == 11)
        b0LomaxB = exp(tmp(Dist == 11));
    else
        b0LomaxA = [];
        b0LomaxB = [];
    end
    if any(Dist == 12)
        b0LogB = exp(tmp(Dist == 12));
    else
        b0LogA = [];
        b0LogB = [];
    end  
    if any(Dist == 13)
        b0LogLB = exp(tmp(Dist == 13));
    else
        b0LogLA = [];
        b0LogLB = [];
    end
    if any(Dist == 14)
        b0GumB = exp(tmp(Dist == 14));
    else
        b0GumA = [];
        b0GumB = [];
    end    
    if any(Dist == 15)
        b0CaucB = exp(tmp(Dist == 15));
    else
        b0CaucA = [];
        b0CaucB = [];
    end
    if any(Dist == 16)
        % b0RayB = exp(tmp(Dist == 16));
    else
        b0RayA = [];
        % b0RayB = [];
    end
    if any(Dist == 17)
        % b0ExpB = exp(tmp(Dist == 17));
    else
        b0ExpA = [];
        % b0ExpB = [];
    end
    
    
    VC = tril(ones(NVarA));
    VC(VC == 1) = b0v;
    
    if any(Dist >= 3 & Dist <= 5 | Dist == 9 | Dist == 10 | Dist == 11 | Dist == 12 | Dist == 13 | Dist == 14 | Dist == 15 | Dist == 16 | Dist == 17)
        tmp = sqrt(sum(VC(Dist >= 3 & Dist <= 5 | Dist == 9 | Dist == 10 | Dist == 11 | Dist == 12 | Dist == 13 | Dist == 14 | Dist == 15 | Dist == 16 | Dist == 17,:).^2,2));
        VC(Dist >= 3 & Dist <= 5 | Dist == 9 | Dist == 10 | Dist == 11 | Dist == 12 | Dist == 13 | Dist == 14 | Dist == 15 | Dist == 16 | Dist == 17,:) = VC(Dist >= 3 & Dist <= 5 | Dist == 9 | Dist == 10 | Dist == 11 | Dist == 12 | Dist == 13 | Dist == 14 | Dist == 15 | Dist == 16 | Dist == 17,:)./tmp;
    end
    
    b0m = b0(NVarA*(NVarA/2+1.5)+1:NVarA*(NVarA/2+1.5+NVarM));
    b0m = reshape(b0m,[NVarA,NVarM]);
    b0s = b0(NVarA*(NVarA/2+1.5+NVarM)+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS);
    b0t = b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarNLT);
    b0j = b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarNLT+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarNLT+2*Johnson);
end

%% Nonlinear transformations
if NVarNLT > 0
    % IndTransNon0 = (abs(bt) > 0.00001)';
    IndTransNon0 = (abs(b0t) > eps)';
    Xt = XXa(:,NLTVariables,:);
    %     bt_tmp = permute(b0t(:, ones(NAlt*NCT,1), ones(NP,1)), [2 1 3]);
    bt_tmp = b0t(:,ones(NAlt*NCT,1))';
    bt_tmp = bt_tmp(:,:,ones(1,1,NP));
    
    % Transform variables with the chosen transformation type (Box-Cox, Yeo-Johnson)    
    if NLTType == 1 % BC
        Xt(:,IndTransNon0,:) = -(Xt(:,IndTransNon0,:).^bt_tmp(:,IndTransNon0,:) - 1)./bt_tmp(:,IndTransNon0,:);
        Xt(:,~IndTransNon0,:) = -log(Xt(:,~IndTransNon0,:));
    elseif NLTType == 2 % YJ
        IndXtNon0 = (Xt >= 0);
        IndXtCase1 = IndXtNon0 & IndTransNon0; % X >= 0, lam ~= 0
        IndXtCase2 = IndXtNon0 & ~IndTransNon0; % X >= 0, lam = 0
        %     IndTransNon2 = (abs(bt - 2) < 0.00001)';
        IndTransNon2 = (abs(b0t - 2) > eps)';
        IndXtCase3 = ~IndXtNon0 & IndTransNon2;  % X < 0, lam ~= 2
        IndXtCase4 = ~IndXtNon0 & ~IndTransNon2; % X < 0, lam = 2
        %         bt_tmp = b0t(:,ones(size(XXa,1),1))';
        Xt(IndXtCase1) = ((Xt(IndXtCase1) + 1).^bt_tmp(IndXtCase1) - 1)./bt_tmp(IndXtCase1);
        Xt(IndXtCase2) = log(Xt(IndXtCase2) + 1);
        Xt(IndXtCase3) = -((-Xt(IndXtCase3) + 1).^(2 - bt_tmp(IndXtCase3)) - 1)./(2 - bt_tmp(IndXtCase3));
        Xt(IndXtCase4) = -log(-Xt(IndXtCase4) + 1);
    end
    
    if EstimOpt.NumGrad == 0 %
        if NLTType == 1 % BC
            XXt = XXa(:,NLTVariables,:);
            XXt(:,IndTransNon0,:) = -(XXt(:,IndTransNon0,:).^bt_tmp(:,IndTransNon0,:).*(bt_tmp(:,IndTransNon0,:).*log(XXt(:,IndTransNon0,:))-1)+1)./(bt_tmp(:,IndTransNon0,:).^2);
            XXt(:,IndTransNon0 == 0,:) = -0.5*log(XXt(:,IndTransNon0 == 0)).^2;
        elseif NLTType == 2 % YJ
            XXt = XXa(:,NLTVariables,:);
            XXt(IndXtCase1) = ((XXt(IndXtCase1)+1).^bt_tmp(IndXtCase1).*(bt_tmp(IndXtCase1).*log(XXt(IndXtCase1)+1)-1)+1)./(bt_tmp(IndXtCase1).^2);% X >= 0, lam ~= 0
            XXt(IndXtCase2) = 0.5*log(XXt(IndXtCase2)+1).^2;% X >= 0, lam == 0
            XXt(IndXtCase3) = -((-XXt(IndXtCase3)+1).^(2-bt_tmp(IndXtCase3)).*(1-(2-bt_tmp(IndXtCase3)).*log(-XXt(IndXtCase3)+1))-1)./((2-bt_tmp(IndXtCase3)).^2);% X < 0, lam ~= 2
            XXt(IndXtCase4) = -0.5*log(-XXt(IndXtCase4)+1).^2;% X < 0, lam == 2
        end
    end
    XXa(:,NLTVariables,:) = Xt;
else
    %     if EstimOpt.NumGrad == 0
    XXt = zeros(0,0,NP);
    %     end
end

%% Transformation from standard normal to normal(mu, sig2) draws
if mCT == 0
    b0n = b0a + b0m*XXm;
    b0n = reshape(b0n((1:size(b0n,1))'*ones(1,NRep),(1:size(b0n,2))'),[NVarA,NRep*NP]);  % NVarA x NRep*NP
    b_mtx = b0n + VC*err;  % NVarA x NRep*NP
    % draws for distributions 3, 4 and 5 remain standard normal in b_mtx
    % and are adjusted below
else
    Xmfit = reshape(b0m*XXm, [NVarA, 1, NAlt*NCT,NP]); 
    b_mtx = reshape(b0a + VC*err,[NVarA, NRep, 1, NP]);  % NVarA x NRep*NP
    b_mtx = reshape(b_mtx + Xmfit, [NVarA, NRep*NAlt*NCT*NP]);
    
end
%% Transformations for other distributions
if sum(Dist == 1) > 0 % Log-normal
    b_mtx(Dist == 1,:) = exp(b_mtx(Dist == 1,:));
end
if sum(Dist == 2) > 0 % Spike
    b_mtx(Dist == 2,:) = max(b_mtx(Dist == 2,:),0);
end
if sum(Dist == 3) > 0 % Triangular
    tmp = normcdf(b_mtx(Dist == 3,:));
    Triang = Triang(ones(NRep*NP,1),:)';
    Ftriang = (b0triang_c - Triang)./(b0triang_b - Triang);
    bmtx_triang = zeros(size(tmp));
    tmp2 = (b0triang_b - Triang).*(b0triang_c - Triang);
    bmtx_triang(tmp < Ftriang) = Triang(tmp < Ftriang)+ sqrt(tmp(tmp < Ftriang).*tmp2(tmp < Ftriang));
    tmp2 = (b0triang_b - Triang).*(b0triang_b-b0triang_c);
    %bmtx_triang(tmp >= Ftriang) = b0triang_b(tmp >= Ftriang)- sqrt((1-tmp(tmp >= Ftriang)).*tmp2(tmp >= Ftriang));
    bmtx_triang(tmp >= Ftriang) = b0triang_b- sqrt((1-tmp(tmp >= Ftriang)).*tmp2(tmp >= Ftriang));
    b_mtx(Dist == 3,:) = bmtx_triang;
end
if sum(Dist == 4) > 0 % Weibull
    tmpWeib = -log(1-normcdf(b_mtx(Dist == 4,:)));
    b_mtx(Dist == 4,:) = b0weibA.*(tmpWeib.^b0weibB);   % inverse CDF function
    
    if EstimOpt.NumGrad == 0
        tmpWeib = reshape(tmpWeib, [], NRep, NP);
    end
else
    tmpWeib = double.empty(0, 0, NP); % initialize variable for parfor loop (even thoough it is not used)
end
if sum(Dist >= 5) > 0 % Johnson
    if sum(Dist == 5) > 0 % Sinh-Arcsinh
        b_mtx(Dist == 5,:) = b0sinhA + b0sinhB.*asinh(b_mtx(Dist == 5,:));
        b_mtx(Dist == 5,:) = b0j(1:Johnson,ones(NRep*NP,1)) + exp(b0j(Johnson+1:end,ones(NRep*NP,1))).*sinh(b_mtx(Dist == 5,:));
        b_mtx(Dist == 5,:) = b0j(1:Johnson,:) + exp(b0j(Johnson+1:end,:)).*sinh(b_mtx(Dist == 5,:));
    end
    if sum(Dist == 6) > 0 % Johnson Sb
        tmp = exp(b_mtx(Dist == 6,:));
        b_mtx(Dist == 6,:) = tmp./(1+tmp);
        b_mtx(Dist == 6,:) = b0j(1:Johnson,:) + exp(b0j(Johnson+1:end,:)).*b_mtx(Dist == 6,:);
    end
    if sum(Dist == 7) > 0 % Johnson Su
        b_mtx(Dist == 7,:) = sinh(b_mtx(Dist == 7,:));
        b_mtx(Dist == 7,:) = b0j(1:Johnson,:) + exp(b0j(Johnson+1:end,:)).*b_mtx(Dist == 7,:);
    end
end

if sum(Dist == 9) > 0 % Uni-Log
    tmpUniLog = normcdf(b_mtx(Dist == 9,:));
    b_mtx(Dist == 9,:) = exp(log(b0UniLogA)+tmpUniLog.*(log(b0UniLogB)-log(b0UniLogA)));   % inverse CDF function (dziala)
    if EstimOpt.NumGrad == 0
        tmpUniLog = reshape(tmpUniLog, [], NRep, NP);
    end
else
    tmpUniLog = double.empty(0, 0, NP); % initialize variable for parfor loop (even thoough it is not used)
end

if sum(Dist == 10) > 0 % Pareto
    tmpPareto = normcdf(b_mtx(Dist == 10,:));
    b_mtx(Dist == 10,:) = b0ParetoA./((1-tmpPareto).^(1./b0ParetoB));   % inverse CDF function
    if EstimOpt.NumGrad == 0
        tmpPareto = reshape(tmpPareto, [], NRep, NP);
    end
else
    tmpPareto = double.empty(0, 0, NP); % initialize variable for parfor loop (even thoough it is not used)
end

if sum(Dist == 11) > 0 % Lomax
    tmpLomax = normcdf(b_mtx(Dist == 11,:));
    b_mtx(Dist == 11,:) = b0LomaxB.*(((1./(1-tmpLomax)).^(1./b0LomaxA))-1);   % inverse CDF function
    if EstimOpt.NumGrad == 0
        tmpLomax = reshape(tmpLomax, [], NRep, NP);
    end
else
    tmpLomax = double.empty(0, 0, NP); % initialize variable for parfor loop (even thoough it is not used)
end

if sum(Dist == 12) > 0 % Logistic
    tmpLog = normcdf(b_mtx(Dist == 12,:));
    b_mtx(Dist == 12,:) = b0LogA+b0LogB.*log(tmpLog./(1-tmpLog));   % inverse CDF function
    if EstimOpt.NumGrad == 0
        tmpLog = reshape(tmpLog, [], NRep, NP);
    end
else
    tmpLog = double.empty(0, 0, NP); % initialize variable for parfor loop (even thoough it is not used)
end

if sum(Dist == 13) > 0 % Log-Logistic
    tmpLogL = normcdf(b_mtx(Dist == 13,:));
    b_mtx(Dist == 13,:) = b0LogLA.*(1./((1./tmpLogL)-1)).^(1./b0LogLB);   % inverse CDF function
    if EstimOpt.NumGrad == 0
        tmpLogL = reshape(tmpLogL, [], NRep, NP);
    end
else
    tmpLogL = double.empty(0, 0, NP); % initialize variable for parfor loop (even thoough it is not used)
end

if sum(Dist == 14) > 0 % Gumbel
    tmpGum = normcdf(b_mtx(Dist == 14,:));
    b_mtx(Dist == 14,:) = b0GumA+b0GumB.*log(1./log(1./tmpGum));   % inverse CDF function
    if EstimOpt.NumGrad == 0
        tmpGum = reshape(tmpGum, [], NRep, NP);
    end
else
    tmpGum = double.empty(0, 0, NP); % initialize variable for parfor loop (even thoough it is not used)
end

if sum(Dist == 15) > 0 % Cauchy
    tmpCauc = normcdf(b_mtx(Dist == 15,:));
    b_mtx(Dist == 15,:) = b0CaucA+b0CaucB.*tan(pi.*(tmpCauc-1/2));   % inverse CDF function
    if EstimOpt.NumGrad == 0
        tmpCauc = reshape(tmpCauc, [], NRep, NP);
    end
else
    tmpCauc = double.empty(0, 0, NP); % initialize variable for parfor loop (even thoough it is not used)
end

if sum(Dist == 16) > 0 % Rayleigh
    tmpRay = normcdf(b_mtx(Dist == 16,:));
    b_mtx(Dist == 16,:) = (2.*(b0RayA.^2).*log(1./(1-tmpRay))).^(1/2);   % inverse CDF function
    if EstimOpt.NumGrad == 0
        tmpRay = reshape(tmpRay, [], NRep, NP);
    end
else
    tmpRay = double.empty(0, 0, NP); % initialize variable for parfor loop (even thoough it is not used)
end

if sum(Dist == 17) > 0 % Exponential
    tmpExp = normcdf(b_mtx(Dist == 17,:));
    b_mtx(Dist == 17,:) = (log(1./(1-tmpExp))).*(1./b0ExpA);   % inverse CDF function
    if EstimOpt.NumGrad == 0
        tmpExp = reshape(tmpExp, [], NRep, NP);
    end
else
    tmpExp = double.empty(0, 0, NP); % initialize variable for parfor loop (even thoough it is not used)
end

% save tmp1
b_score = b_mtx;

if WTP_space > 0
    if mCT == 0
        b_mtx_grad = reshape(b_mtx,[NVarA,NRep,NP]); % needed for gradient calculation in WTP_space
    else
        b_mtx_grad = reshape(b_mtx,[NVarA,NRep,NAlt*NCT, NP]); % needed for gradient calculation in WTP_space
        b_mtx_grad = permute(b_mtx_grad, [3 1 2 4]);
    end
    b_mtx(1:end-WTP_space,:) = b_mtx(1:end-WTP_space,:).*b_mtx(WTP_matrix,:);
else
    b_mtx_grad = zeros(0,0,NP);
end

if NVarS > 0
    cs = reshape(exp(Xs*b0s),[NAlt*NCT,1,NP]);
    XXa = XXa .* cs;
end
if mCT == 0
    b_mtx = reshape(b_mtx,[NVarA,NRep,NP]);
else
    b_mtx = permute(reshape(b_mtx,[NVarA,NRep,NAlt*NCT, NP]), [3, 1, 2, 4]); % NAlt*NCT x NvarA x NRep x NP
end
% p0 = zeros(NP,1);
p0 = zeros(NP,NRep);


    if ~any(isnan(XXa(:))) % faster version for complete dataset
        % chosen alternative indicator        
        YYy = YY == 1;

        if mCT == 0 
            % calculation of the chosen alternative probability estimator for
            % each individual        
            parfor n = 1:NP   % for each person
                % utility function    
                %if mCT == 0
                U = reshape(XXa(:,:,n)*b_mtx(:,:,n),[NAlt,NCT,NRep]);

                U = exp(U - max(U,[],1)); % rescale utility to avoid exploding
                % denominator term of the conditional probability fraction
                % (denominator of the logit probability (sum over alternatives))            
                U_sum = reshape(sum(U,1),[NCT,NRep]);
                YYy_n = YYy(:,n);
                % numerator term of the conditional probability fraction            
                U_selected = reshape(U(YYy_n(:,ones(NRep,1))),[NCT,NRep]);
                % calculate probability estimator for the chosen alternative
                % prod for panel data, mean for simulated probability            
                p0(n,:) = prod(U_selected./U_sum,1);
            end            
        else
            % calculation of the chosen alternative probability estimator for
            % each individual        
            parfor n = 1:NP   % for each person
                % utility function    
                %if mCT ~= 0
                U = reshape(sum(XXa(:,:,n).*b_mtx(:,:,:,n),2), [NAlt, NCT, NRep]);
                
                U = exp(U - max(U,[],1)); % rescale utility to avoid exploding
                % denominator term of the conditional probability fraction
                % (denominator of the logit probability (sum over alternatives))            
                U_sum = reshape(sum(U,1),[NCT,NRep]);
                YYy_n = YYy(:,n);
                % numerator term of the conditional probability fraction            
                U_selected = reshape(U(YYy_n(:,ones(NRep,1))),[NCT,NRep]);
                % calculate probability estimator for the chosen alternative
                % prod for panel data, mean for simulated probability            
                % p0(n) = mean(prod(U_selected./U_sum,1),2);
                p0(n,:) = prod(U_selected./U_sum,1);
            end 
        end

    else
        if mCT == 0
            parfor n = 1:NP
                YnanInd = ~isnan(YY(:,n));
                XXa_n = XXa(:,:,n);
                NAltMissIndExp_n = NAltMissIndExp(:,n);
                NAltMissIndExp_n = NAltMissIndExp_n(YnanInd);
                if var(NAltMissIndExp_n(NAltMissIndExp_n > 0)) == 0 % if NAlt is constant per individual (but can vary between individuals)
                    %if mCT == 0
                    U = reshape(XXa_n(YnanInd,:)*b_mtx(:,:,n),[NAltMiss(n),NCTMiss(n),NRep]);

                    U = exp(U - max(U,[],1));
                    U_sum = reshape(sum(U,1),[NCTMiss(n),NRep]);
                else
                    NAltMissInd_n = NAltMissInd(:,n);
                    NAltMissInd_n = NAltMissInd_n(MissingCT(:,n) == 0);
                    %if mCT == 0
                    U = XXa_n(YnanInd,:)*b_mtx(:,:,n);

                    Uniq = unique(NAltMissIndExp_n);
                    U_sum = zeros(NCTMiss(n),NRep);
                    if length(Uniq) == 2
                        U_tmp = U(NAltMissIndExp_n == Uniq(1),:);
                        U_tmp = reshape(U_tmp,[Uniq(1),size(U_tmp,1)/Uniq(1),NRep]);
                        U_tmp = exp(U_tmp - max(U_tmp));
                        U_sum(NAltMissInd_n == Uniq(1),:) = reshape(sum(U_tmp,1),[size(U_tmp,2),NRep]);
                        U(NAltMissIndExp_n == Uniq(1),:) = reshape(U_tmp,[size(U_tmp,2)*Uniq(1),NRep]);
                        U_tmp = U(NAltMissIndExp_n == Uniq(2),:);
                        U_tmp = reshape(U_tmp,[Uniq(2),size(U_tmp,1)/Uniq(2),NRep]);
                        U_tmp = exp(U_tmp - max(U_tmp));
                        U_sum(NAltMissInd_n == Uniq(2),:) = reshape(sum(U_tmp,1),[size(U_tmp,2),NRep]);
                        U(NAltMissIndExp_n == Uniq(2),:) = reshape(U_tmp,[size(U_tmp,2)*Uniq(2),NRep]);
                    else
                        for i = 1:length(Uniq)
                            U_tmp = U(NAltMissIndExp_n == Uniq(i),:);
                            U_tmp = reshape(U_tmp,[Uniq(i),size(U_tmp,1)/Uniq(i),NRep]);
                            U_tmp = exp(U_tmp - max(U_tmp));
                            U_sum(NAltMissInd_n == Uniq(i),:) = reshape(sum(U_tmp,1),[size(U_tmp,2),NRep]);
                            U(NAltMissIndExp_n == Uniq(i),:) = reshape(U_tmp,[size(U_tmp,2)*Uniq(i),NRep]);
                        end
                    end
                end
                YYy_n = YY(:,n)==1;
                U_selected = reshape(U(YYy_n(YnanInd,ones(NRep,1))),[NCTMiss(n),NRep]);
                % p0(n) = mean(prod(U_selected./U_sum,1));
                p0(n,:) = prod(U_selected./U_sum,1);
            end            
        else
            parfor n = 1:NP
                YnanInd = ~isnan(YY(:,n));
                XXa_n = XXa(:,:,n);
                NAltMissIndExp_n = NAltMissIndExp(:,n);
                NAltMissIndExp_n = NAltMissIndExp_n(YnanInd);
                if var(NAltMissIndExp_n(NAltMissIndExp_n > 0)) == 0 % if NAlt is constant per individual (but can vary between individuals)
                    %if mCT ~= 0
                    U = reshape(sum(XXa_n(YnanInd,:).*b_mtx(YnanInd,:,:,n),2), [NAltMiss(n),NCTMiss(n), NRep]);

                    U = exp(U - max(U,[],1));
                    U_sum = reshape(sum(U,1),[NCTMiss(n),NRep]);
                else
                    NAltMissInd_n = NAltMissInd(:,n);
                    NAltMissInd_n = NAltMissInd_n(MissingCT(:,n) == 0);
                    %if mCT ~= 0
                    U = reshape(sum(XXa_n(YnanInd,:).*b_mtx(YnanInd,:,:,n),2), [NAltMiss(n)*NCTMiss(n), NRep]);

                    Uniq = unique(NAltMissIndExp_n);
                    U_sum = zeros(NCTMiss(n),NRep);
                    if length(Uniq) == 2
                        U_tmp = U(NAltMissIndExp_n == Uniq(1),:);
                        U_tmp = reshape(U_tmp,[Uniq(1),size(U_tmp,1)/Uniq(1),NRep]);
                        U_tmp = exp(U_tmp - max(U_tmp));
                        U_sum(NAltMissInd_n == Uniq(1),:) = reshape(sum(U_tmp,1),[size(U_tmp,2),NRep]);
                        U(NAltMissIndExp_n == Uniq(1),:) = reshape(U_tmp,[size(U_tmp,2)*Uniq(1),NRep]);
                        U_tmp = U(NAltMissIndExp_n == Uniq(2),:);
                        U_tmp = reshape(U_tmp,[Uniq(2),size(U_tmp,1)/Uniq(2),NRep]);
                        U_tmp = exp(U_tmp - max(U_tmp));
                        U_sum(NAltMissInd_n == Uniq(2),:) = reshape(sum(U_tmp,1),[size(U_tmp,2),NRep]);
                        U(NAltMissIndExp_n == Uniq(2),:) = reshape(U_tmp,[size(U_tmp,2)*Uniq(2),NRep]);
                    else
                        for i = 1:length(Uniq)
                            U_tmp = U(NAltMissIndExp_n == Uniq(i),:);
                            U_tmp = reshape(U_tmp,[Uniq(i),size(U_tmp,1)/Uniq(i),NRep]);
                            U_tmp = exp(U_tmp - max(U_tmp));
                            U_sum(NAltMissInd_n == Uniq(i),:) = reshape(sum(U_tmp,1),[size(U_tmp,2),NRep]);
                            U(NAltMissIndExp_n == Uniq(i),:) = reshape(U_tmp,[size(U_tmp,2)*Uniq(i),NRep]);
                        end
                    end
                end
                YYy_n = YY(:,n)==1;
                U_selected = reshape(U(YYy_n(YnanInd,ones(NRep,1))),[NCTMiss(n),NRep]);
                % p0(n) = mean(prod(U_selected./U_sum,1));
                p0(n,:) = prod(U_selected./U_sum,1);
            end              
        end
    end

    % save tmp1


fx = mean(p0,2); % NP x 1
Score = zeros(NP, NVarA);
if WTP_space > 0
    b_score = reshape(b_score, NVarA,NRep,NP); ...
    for j = 1:NVarA
        bx = squeeze(b_score(j,:,:))'; % NP x NRep
        Score(:,j) = mean(p0.*bx,2)./fx;
    end
else
    b_score = reshape(b_score, NVarA,NRep,NP); ...
% This calculates individual scores in terms of WTP, even for PS models. 
%     fee = squeeze(b_score(end,:,:))';
%     for j = 1:NVarA-1
%         bx = squeeze(b_score(j,:,:))'; % NP x NRep
%         Score(:,j) = mean(p0.*bx./fee,2)./fx;
%     end

% This keeps individual scores in PS: 
    for j = 1:NVarA
        bx = squeeze(b_score(j,:,:))'; % NP x NRep
        Score(:,j) = mean(p0.*bx,2)./fx;
    end

end