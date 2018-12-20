function Score = BayesScores(YY,XXa,XXm,Xs,err,EstimOpt,b0,varargin)

if nargin < 8
    model_idx = 0;
else
    model_idx = varargin{1};
    Xstr = varargin{2};
end

save BayesScores
% return

if model_idx == 1 % HMXL
    NAlt = EstimOpt.NAlt;
    NCT = EstimOpt.NCT;
    NP = EstimOpt.NP;
    NRep = EstimOpt.NRep;
    NVarA = EstimOpt.NVarA;
    WTP_space = EstimOpt.WTP_space;
    WTP_matrix = EstimOpt.WTP_matrix;
    NCTMiss = EstimOpt.NCTMiss;
    NAltMiss = EstimOpt.NAltMiss;
    NLatent = EstimOpt.NLatent;
    NVarStr = EstimOpt.NVarStr;
    MeaMatrix = EstimOpt.MeaMatrix;
    MeaSpecMatrix = EstimOpt.MeaSpecMatrix;
    NVarMeaExp = EstimOpt.NVarMeaExp;
    MeaExpMatrix = EstimOpt.MeaExpMatrix;
    Dist = EstimOpt.Dist;
    FullCov = EstimOpt.FullCov;
    RealMin = EstimOpt.RealMin;
    ScaleLV = EstimOpt.ScaleLV;
    NVarS = EstimOpt.NVarS;
    NVarM = EstimOpt.NVarM;
    b_mtx_grad = [];
    MissingIndMea = EstimOpt.MissingIndMea;
    NAltMissIndExp = EstimOpt.NAltMissIndExp;
    NAltMissInd = EstimOpt.NAltMissInd;
    MissingCT = EstimOpt.MissingCT;
    %     VC = tril(ones(NVarA));
    %     VC(VC == 1) = (1:(NVarA*(NVarA-1)/2+NVarA))';
    %     EstimOpt.DiagIndex = diag(VC);
    if isfield(EstimOpt,'Triang') == 0 || length(EstimOpt.Triang) ~= sum(EstimOpt.Dist == 3,2) % Needed only if any parameter has triangular distribution
        EstimOpt.Triang = zeros(1,sum(EstimOpt.Dist == 3,2));
    elseif length(EstimOpt.Triang) == 1
        EstimOpt.Triang = EstimOpt.Triang*ones(1,sum(EstimOpt.Dist == 3,2));
    else
        EstimOpt.Triang = EstimOpt.Triang(:)';
    end
    EstimOpt.Johnson = sum(EstimOpt.Dist >= 5);
else
    DiagIndex = EstimOpt.DiagIndex;
end


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

if model_idx == 0 % MXL
    
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
        VC = tril(ones(NVarA));
        VC(VC==1) = b0v;
        if any(Dist >= 3 & Dist <= 5)
            tmp = sqrt(sum(VC(Dist >= 3 & Dist <= 5,:).^2,2));
            VC(Dist >= 3 & Dist <= 5,:) = VC(Dist >= 3 & Dist <= 5,:)./tmp(:,ones(1,NVarA));
        end
        b0m = b0(NVarA*(NVarA/2+1.5)+1:NVarA*(NVarA/2+1.5+NVarM));
        b0m = reshape(b0m,[NVarA,NVarM]);
        b0s = b0(NVarA*(NVarA/2+1.5+NVarM)+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS);
        b0j = b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS+2*Johnson);
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
    
elseif model_idx == 1 % HMXL
    
    if FullCov == 0
        ba = b0(1:NVarA); % b atrybutów
        bv = b0(NVarA+1:2*NVarA);
        %VC = diag(bv.^2);
        VC = diag(bv);
        bm = reshape(b0(2*NVarA+1:NVarA*(NVarM+2)),[NVarA,NVarM]); % b mean covariates
        bl = reshape(b0((2+NVarM)*NVarA+1:NVarA*(NLatent+2+NVarM)),[NVarA,NLatent]); % b interakcji z LV
        bs = b0(NVarA*(NLatent+NVarM+2)+1:NVarA*(NLatent+NVarM+2)+NVarS); % b scale
        bstr = reshape(b0(NVarA*(NLatent+2+NVarM)+NVarS+1:(NVarA+NVarStr)*NLatent+(2+NVarM)*NVarA+NVarS),[NVarStr,NLatent]); % b równania struktury
        bmea = b0((NVarA+NVarStr)*NLatent+(2+NVarM)*NVarA+NVarS+1:end); % b measurement
    elseif FullCov == 1
        ba = b0(1:NVarA); % b atrybutów
        bv = b0(NVarA+1:NVarA+sum(1:NVarA,2));
        VC = tril(ones(NVarA));
        VC(VC == 1) = bv;
        bm = reshape(b0(NVarA+sum(1:NVarA,2)+1:NVarA*(NVarM+1)+sum(1:NVarA,2)),[NVarA,NVarM]); % b mean covariates
        bl = reshape(b0(NVarA*(1+NVarM)+sum(1:NVarA,2)+1:NVarA*(NLatent+1+NVarM)+sum(1:NVarA,2)),[NVarA,NLatent]); % b interakcji z LV
        bs = b0(NVarA*(NLatent+1+NVarM)+sum(1:NVarA,2)+1:NVarA*(NLatent+1+NVarM)+sum(1:NVarA,2)+NVarS); % b scale
        bstr = reshape(b0(NVarA*(NLatent+NVarM+1)+sum(1:NVarA,2)+NVarS+1:(NVarA+NVarStr)*NLatent+NVarA*(1+NVarM)+sum(1:NVarA,2)+NVarS),[NVarStr,NLatent]); % b równania struktury
        bmea = b0((NVarA+NVarStr)*NLatent+NVarA*(1+NVarM)+sum(1:NVarA,2)+NVarS+1:end); % b measurement
    elseif FullCov == 2
        ba = b0(1:NVarA); % b atrybutów
        bv = b0(NVarA+1:NVarA+sum(1:NVarA+NLatent,2)-NLatent);
        VC = tril(ones(NVarA+NLatent));
        VCtmp2 = diag(ones(NLatent+NVarA,1));
        VCtmp2(1:NVarA,1:NVarA) = 0;
        VC(VCtmp2 == 1) = 0;
        VC(VC == 1) = bv;
        VC(VCtmp2 == 1) = 1;
        VCGrad = VC; % before normalization, for gradient
        NormVar = sqrt(sum(VC(NVarA+1:end,:).^2,2));
        VC(NVarA+1:end,:) = VC(NVarA+1:end,:)./NormVar;
        bm = reshape(b0(NVarA+sum(1:NVarA+NLatent,2)-NLatent+1:NVarA*(NVarM+1)+sum(1:NVarA+NLatent,2)-NLatent),[NVarA,NVarM]); % b mean covariates
        bl = reshape(b0(NVarA*(1+NVarM)+sum(1:NVarA+NLatent,2)-NLatent+1:NVarA*(NLatent+1+NVarM)+sum(1:NVarA+NLatent,2)-NLatent),[NVarA,NLatent]); % b interakcji z LV
        bs = b0(NVarA*(NLatent+1+NVarM)+sum(1:NVarA+NLatent,2)-NLatent+1:NVarA*(NLatent+1+NVarM)+sum(1:NVarA+NLatent,2)-NLatent+NVarS); % b scale
        bstr = reshape(b0(NVarA*(NLatent+NVarM+1)+sum(1:NVarA+NLatent,2)-NLatent+NVarS+1:(NVarA+NVarStr)*NLatent+NVarA*(1+NVarM)+sum(1:NVarA+NLatent,2)-NLatent+NVarS),[NVarStr,NLatent]); % b równania struktury
        bmea = b0((NVarA+NVarStr)*NLatent+NVarA*(1+NVarM)+sum(1:NVarA+NLatent,2)-NLatent+NVarS+1:end); % b measurement
        err2 = VC*err;
    end
    if ScaleLV == 1
        bsLV = bs(NVarS-NLatent+1:end)';
        bs = bs(1:NVarS-NLatent);
        NVarS = NVarS - NLatent;
    else
        bsLV = zeros(1,0);
    end
    
    LV_tmp = Xstr*bstr; % NP x Latent
    LV_tmp = reshape(permute(LV_tmp(:,:,ones(NRep,1)),[2 3 1]),[NLatent,NRep*NP]); % Latent*NRep*NP
    if FullCov == 2
        LV_tmp = LV_tmp + err2(NVarA+1:end,:); % Latent x NRep*NP
    else
        LV_tmp = LV_tmp + err(NVarA+1:end,:); % Latent x NRep*NP
    end
    sLV = std(LV_tmp,0,2);
    LV = (LV_tmp - mean(LV_tmp,2))./sLV; % normalilzing for 0 mean and std
    
    if EstimOpt.NVarM > 0
        ba = ba(:,ones(NP,1)) + bm*Xm; % NVarA x NP
        ba = reshape(permute(ba(:,:,ones(NRep,1)),[1 3 2]),[NVarA,NRep*NP]);
    else
        ba = ba(:,ones(NP*NRep,1));
    end
    
    if FullCov == 2
        b_mtx = ba + bl*LV + err2(1:NVarA,:); % NVarA x NRep*NP
    else
        b_mtx = ba + bl*LV + VC*err(1:NVarA,:); % NVarA x NRep*NP
    end
    
    if sum(Dist == 1) > 0 % Log - normal
        b_mtx(Dist == 1,:) = exp(b_mtx(Dist == 1,:));
    elseif sum(Dist == 2) > 0 % Spike
        b_mtx(Dist == 2,:) = max(b_mtx(Dist == 2,:),0);
    end
    
    if WTP_space > 0
        if EstimOpt.NumGrad == 0
            b_mtx_grad = reshape(b_mtx,[NVarA,NRep,NP]);% needed for gradient calculation in WTP_space
        end
        b_mtx(1:end-WTP_space,:) = b_mtx(1:end-WTP_space,:).*b_mtx(WTP_matrix,:);
    elseif nargout == 2
        if ScaleLV == 1 && any(Dist == 1)
            b_mtx_grad = reshape(b_mtx,[NVarA,NRep,NP]);% needed for gradient calculation in WTP_space
        else
            b_mtx_grad = zeros(0,0,NP);
        end
    end
    if ScaleLV == 1
        ScaleLVX = exp(bsLV*LV);
        b_mtx = ScaleLVX.*b_mtx;
        ScaleLVX = permute(reshape(ScaleLVX,[1,NRep,NP]),[1 3 2]);
    else
        ScaleLVX = zeros(0,NP,0);
    end
    
    b_mtx_sliced = reshape(b_mtx,[NVarA,NRep,NP]); % NVarA x NRep x NP
    if NVarS > 0
        Scale = reshape(exp(Xs*bs),[NAlt*NCT,1,NP]);
        Xa = Xa.*Scale;
    end
    
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

p0 = zeros(NP,NRep);

if nargout == 1 % function value only
    if any(isnan(Xa(:))) == 0 % faster version for complete dataset
        parfor n = 1:NP
            Yy_n = Y(:,n) == 1;
            Xa_n = Xa(:,:,n);
            U = reshape(Xa_n*b_mtx_sliced(:,:,n),[NAlt,NCT,NRep]);
            U = exp(U - max(U)); % rescale utility to avoid exploding
            U_sum = reshape(sum(U,1),[NCT,NRep]);
            U_selected = reshape(U(Yy_n(:,ones(NRep,1))),[NCT,NRep]);
            p0(n,:) = prod(U_selected./U_sum,1);
        end
    else
        parfor n = 1:NP
            YnanInd = ~isnan(Y(:,n));
            Yy_n = Y(:,n) == 1;
            Xa_n = Xa(:,:,n);
            U = reshape(Xa_n(YnanInd,:)*b_mtx_sliced(:,:,n),NAltMiss(n),NCTMiss(n),NRep);
            U = exp(U - max(U)); % rescale utility to avoid exploding
            U_sum = reshape(sum(U,1),[NCTMiss(n),NRep]);
            U_selected = reshape(U(Yy_n(YnanInd,ones(NRep,1))),[NCTMiss(n),NRep]);
            p0(n,:) = prod(U_selected./U_sum,1);
        end
    end
    
    %     L_mea = ones(NP,NRep);
    %     l = 0;
    %
    %     if NVarMeaExp > 0
    %         Xmea_exp = reshape(Xmea_exp,[1,NP,NVarMeaExp]);
    %         Xmea_exp = reshape(Xmea_exp(ones(NRep,1),:,:),[NP*NRep,NVarMeaExp]);
    %     end
    %     for i = 1:size(Xmea,2)
    %         if MeaSpecMatrix(i) == 0 % OLS
    %             if MeaExpMatrix(i) == 0
    %                 X = [ones(NRep*NP,1),LV(MeaMatrix(:,i)' == 1,:)'];
    %             else
    %                 X = [ones(NRep*NP,1),LV(MeaMatrix(:,i)' == 1,:)',Xmea_exp];
    %             end
    %             b = bmea(l+1:l+size(X,2)+1);
    %             X = reshape(X*b(1:end-1),[NRep,NP])';
    %             L_mea(MissingIndMea(:,i) == 0,:) = L_mea(MissingIndMea(:,i) == 0,:).*normpdf(Xmea(MissingIndMea(:,i) == 0,i*ones(NRep,1)),X(MissingIndMea(:,i) == 0,:),exp(b(end)));
    %             l = l + size(X,2) + 1;
    %         elseif MeaSpecMatrix(i) == 1 % MNL
    %             UniqueMea = unique(Xmea(:,i));
    %             k = length(UniqueMea) - 1;
    %             if MeaExpMatrix(i) == 0
    %                 X = [ones(NRep*NP,1),LV(MeaMatrix(:,i)' == 1,:)'];
    %             else
    %                 X = [ones(NRep*NP,1),LV(MeaMatrix(:,i)' == 1,:)',Xmea_exp];
    %             end
    %             V = exp(X*reshape([zeros(size(X,2),1);bmea(l+1:l+size(X,2)*k)],[size(X,2),k+1])); % NRep*NP x unique values of attitude
    %             Vsum = sum(V,2);
    %             V = reshape(V./Vsum,[NRep,NP,k+1]); %NRep x NP x unique
    %             V = permute(V,[2 1 3]); % NP x NRep x unique
    %             L = zeros(NP,NRep);
    %             for j = 1:length(UniqueMea)
    %                 L(Xmea(:,i) == UniqueMea(j),:) = V(Xmea(:,i) == UniqueMea(j),:,j);
    %             end
    %             L_mea = L_mea.*L;
    %             l = l + size(X,2)*k;
    %
    %         elseif MeaSpecMatrix(i) == 2 % Ordered Probit
    %             UniqueMea = unique(Xmea(EstimOpt.MissingIndMea(:,i) == 0,i));
    %             k = length(UniqueMea) - 1;
    %             if MeaExpMatrix(i) == 0
    %                 X = LV(MeaMatrix(:,i)' == 1,:)';
    %             else
    %                 X = [LV(MeaMatrix(:,i)' == 1,:)',Xmea_exp];
    %             end
    %             tmp = (MeaExpMatrix(i) ~= 0)*NVarMeaExp;
    %             b = bmea(l+1:l+k+size(X,2));
    %             Xb = reshape(X*b(1:sum(MeaMatrix(:,i),1)+tmp),[NRep,NP])'; % NP x NRep
    %             Xb = Xb(MissingIndMea(:,i) == 0,:);
    %             alpha = cumsum([b(sum(MeaMatrix(:,i))+tmp+1);exp(b(sum(MeaMatrix(:,i))+tmp+2:end))]);
    %             L = zeros(sum(MissingIndMea(:,i) == 0),NRep);
    %             L(Xmea(MissingIndMea(:,i) == 0,i) == min(UniqueMea),:) = normcdf(alpha(1)-Xb(Xmea(MissingIndMea(:,i) == 0,i) == min(UniqueMea),:));
    %             L(Xmea(MissingIndMea(:,i) == 0,i) == max(UniqueMea),:) = 1 - normcdf(alpha(end)-Xb(Xmea(MissingIndMea(:,i) == 0,i) == max(UniqueMea),:));
    %             for j = 2:k
    %                 L(Xmea(MissingIndMea(:,i) == 0,i) == UniqueMea(j),:) = normcdf(alpha(j)-Xb(Xmea(MissingIndMea(:,i) == 0,i) == UniqueMea(j),:)) - normcdf(alpha(j-1)-Xb(Xmea(MissingIndMea(:,i) == 0,i) == UniqueMea(j),:));
    %             end
    %             L_mea(MissingIndMea(:,i) == 0,:) = L_mea(MissingIndMea(:,i) == 0,:).*L;
    %             l = l + k + size(X,2);
    %         elseif MeaSpecMatrix(i) == 3 % Poisson
    %             if MeaExpMatrix(i) == 0
    %                 X = [ones(NRep*NP,1),LV(MeaMatrix(:,i)' == 1,:)'];
    %             else
    %                 X = [ones(NRep*NP,1),LV(MeaMatrix(:,i)' == 1,:)',Xmea_exp];
    %             end
    %             b = bmea(l+1:l+size(X,2));
    %             fit = reshape(X*b,[NRep,NP])';
    %             lam = exp(fit);
    %             if RealMin == 1
    %                 L = exp(fit.*Xmea(:,i*ones(NRep,1))-lam)./min(gamma(Xmea(:,i*ones(NRep,1))+1),realmax);
    %             else
    %                 %L = exp(fit.*Xmea(:,i*ones(NRep,1))-lam)./ gamma(Xmea(:,i*ones(NRep,1))+1);
    %                 L = exp(fit.*Xmea(:,i*ones(NRep,1))-lam-gammaln(Xmea(:,i*ones(NRep,1))+1));
    %             end
    %             L_mea = L_mea.*L;
    %             l = l + size(X,2);
    %         elseif MeaSpecMatrix(i) == 4 % NB
    %             if MeaExpMatrix(i) == 0
    %                 X = [ones(NRep*NP,1),LV(MeaMatrix(:,i)' == 1,:)'];
    %             else
    %                 X = [ones(NRep*NP,1),LV(MeaMatrix(:,i)' == 1,:)',Xmea_exp];
    %             end
    %             b = bmea(l+1:l+size(X,2));
    %             fit = reshape(X*b,[NRep,NP])';
    %             lam = exp(fit);
    %             theta = exp(bmea(l+size(X,2)+1));
    %             u = theta./(theta+lam);
    %             if RealMin == 1
    %                 L = min(gamma(theta+Xmea(:,i*ones(NRep,1))),realmax)./(gamma(theta).*min(gamma(Xmea(:,i*ones(NRep,1))+1),realmax));
    %             else
    %                 L = exp(gammaln(theta+Xmea(:,i*ones(NRep,1))) - gammaln(theta) - gammaln(Xmea(:,i*ones(NRep,1))+1));
    %             end
    %             L = L.*(u.^theta).*((1-u).^Xmea(:,i*ones(NRep,1)));
    %             L_mea = L_mea.*L;
    %             l = l + size(X,2) + 1;
    %         elseif MeaSpecMatrix(i) == 5 % ZIP
    %             if MeaExpMatrix(i) == 0
    %                 X = [ones(NRep*NP,1),LV(MeaMatrix(:,i)' == 1,:)'];
    %             else
    %                 X = [ones(NRep*NP,1),LV(MeaMatrix(:,i)' == 1,:)',Xmea_exp];
    %             end
    %             bzip = bmea(l+1:l+size(X,2));
    %             bpoiss = bmea(l+size(X,2)+1:l+2*size(X,2));
    %             fit = reshape(X*bpoiss,[NRep,NP])';
    %             pzip = reshape(exp(X*bzip),[NRep,NP])';
    %             pzip = pzip./(1+pzip);
    %             L = zeros(NP,NRep);
    %             lam = exp(fit);
    %             IndxZIP = Xmea(:,i) == 0;
    %             L(IndxZIP,:) = pzip(IndxZIP,:) + (1 - pzip(IndxZIP,:)).*exp(-lam(IndxZIP,:));
    %             if RealMin == 1
    %                 L(~IndxZIP,:) = (1 - pzip(~IndxZIP,:)).*exp(fit(~IndxZIP,:).*Xmea(~IndxZIP,i*ones(NRep,1))-lam(~IndxZIP,:))./min(gamma(Xmea(~IndxZIP,i*ones(NRep,1))+1),realmax);
    %             else
    %                 L(~IndxZIP,:) = (1 - pzip(~IndxZIP,:)).*exp(fit(~IndxZIP,:).*Xmea(~IndxZIP,i*ones(NRep,1))-lam(~IndxZIP,:)-gammaln(Xmea(~IndxZIP,i*ones(NRep,1))+1));
    %             end
    %             L_mea = L_mea.*L;
    %             l = l + 2*size(X,2);
    %         elseif MeaSpecMatrix(i) == 6 % ZINB
    %             if MeaExpMatrix(i) == 0
    %                 X = [ones(NRep*NP,1),LV(MeaMatrix(:,i)' == 1,:)'];
    %             else
    %                 X = [ones(NRep*NP,1),LV(MeaMatrix(:,i)' == 1,:)',Xmea_exp];
    %             end
    %             bzip = bmea(l+1:l+size(X,2));
    %             bpoiss = bmea(l+size(X,2)+1:l+2*size(X,2));
    %             fit = reshape(X*bpoiss,[NRep,NP])';
    %             theta = exp(bmea(l+2*size(X,2)+1));
    %             pzip = reshape(exp(X*bzip),[NRep,NP])';
    %             pzip = pzip./(1+pzip);
    %             L = zeros(NP,NRep);
    %             lam = exp(fit);
    %             u = theta./(theta+lam);
    %             IndxZIP = Xmea(:,i) == 0;
    %             L(IndxZIP,:) = pzip(IndxZIP,:) + (1 - pzip(IndxZIP,:)).*(u(IndxZIP,:).^theta);
    %             if RealMin == 1
    %                 L(~IndxZIP,:) = (1 - pzip(~IndxZIP,:)).*min(gamma(theta+Xmea(~IndxZIP,i*ones(NRep,1))),realmax)./(gamma(theta).*min(gamma(Xmea(~IndxZIP,i*ones(NRep,1))+1),realmax));
    %             else
    %                 L(~IndxZIP,:) = (1 - pzip(~IndxZIP,:)).*exp(gammaln(theta+Xmea(~IndxZIP,i*ones(NRep,1))) - gammaln(Xmea(~IndxZIP,i*ones(NRep,1))+1)-gammaln(theta));
    %             end
    %             L(~IndxZIP,:) = L(~IndxZIP,:).*(u(~IndxZIP,:).^theta).*((1-u(~IndxZIP,:)).^Xmea(~IndxZIP,i*ones(NRep,1)));
    %             L_mea = L_mea.*L;
    %             l = l + 2*size(X,2) + 1;
    %         end
    %     end
    %     if RealMin == 0
    %         f = -log(mean(p0.*L_mea,2));
    %     else
    %         f = -log(max(realmin,mean(p0.*L_mea,2)));
    %     end
    
    
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
    
end

