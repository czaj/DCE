function Score = BayesScoresHMXL(Y,Xa,Xm,Xs,Xstr,Xmea,Xmea_exp,err_sliced,EstimOpt,B)

save BayesScoresHMXL
% return

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

if FullCov == 0
    ba = B(1:NVarA); % b atrybutów
    bv = B(NVarA+1:2*NVarA);
    %VC = diag(bv.^2);
    VC = diag(bv);
    bm = reshape(B(2*NVarA+1:NVarA*(NVarM+2)),[NVarA,NVarM]); % b mean covariates
    bl = reshape(B((2+NVarM)*NVarA+1:NVarA*(NLatent+2+NVarM)),[NVarA,NLatent]); % b interakcji z LV
    bs = B(NVarA*(NLatent+NVarM+2)+1:NVarA*(NLatent+NVarM+2)+NVarS); % b scale
    bstr = reshape(B(NVarA*(NLatent+2+NVarM)+NVarS+1:(NVarA+NVarStr)*NLatent+(2+NVarM)*NVarA+NVarS),[NVarStr,NLatent]); % b równania struktury
    bmea = B((NVarA+NVarStr)*NLatent+(2+NVarM)*NVarA+NVarS+1:end); % b measurement    
elseif FullCov == 1
    ba = B(1:NVarA); % b atrybutów
    bv = B(NVarA+1:NVarA+sum(1:NVarA,2));
    VC = tril(ones(NVarA));
    VC(VC == 1) = bv;
    bm = reshape(B(NVarA+sum(1:NVarA,2)+1:NVarA*(NVarM+1)+sum(1:NVarA,2)),[NVarA,NVarM]); % b mean covariates
    bl = reshape(B(NVarA*(1+NVarM)+sum(1:NVarA,2)+1:NVarA*(NLatent+1+NVarM)+sum(1:NVarA,2)),[NVarA,NLatent]); % b interakcji z LV
    bs = B(NVarA*(NLatent+1+NVarM)+sum(1:NVarA,2)+1:NVarA*(NLatent+1+NVarM)+sum(1:NVarA,2)+NVarS); % b scale
    bstr = reshape(B(NVarA*(NLatent+NVarM+1)+sum(1:NVarA,2)+NVarS+1:(NVarA+NVarStr)*NLatent+NVarA*(1+NVarM)+sum(1:NVarA,2)+NVarS),[NVarStr,NLatent]); % b równania struktury
    bmea = B((NVarA+NVarStr)*NLatent+NVarA*(1+NVarM)+sum(1:NVarA,2)+NVarS+1:end); % b measurement
elseif FullCov == 2
    ba = B(1:NVarA); % b atrybutów
    bv = B(NVarA+1:NVarA+sum(1:NVarA+NLatent,2)-NLatent);
    VC = tril(ones(NVarA+NLatent));
    VCtmp2 = diag(ones(NLatent+NVarA,1));
    VCtmp2(1:NVarA,1:NVarA) = 0;
    VC(VCtmp2 == 1) = 0;
    VC(VC == 1) = bv;
    VC(VCtmp2 == 1) = 1;
    VCGrad = VC; % before normalization, for gradient
    NormVar = sqrt(sum(VC(NVarA+1:end,:).^2,2)); 
    VC(NVarA+1:end,:) = VC(NVarA+1:end,:)./NormVar;
    bm = reshape(B(NVarA+sum(1:NVarA+NLatent,2)-NLatent+1:NVarA*(NVarM+1)+sum(1:NVarA+NLatent,2)-NLatent),[NVarA,NVarM]); % b mean covariates
    bl = reshape(B(NVarA*(1+NVarM)+sum(1:NVarA+NLatent,2)-NLatent+1:NVarA*(NLatent+1+NVarM)+sum(1:NVarA+NLatent,2)-NLatent),[NVarA,NLatent]); % b interakcji z LV
    bs = B(NVarA*(NLatent+1+NVarM)+sum(1:NVarA+NLatent,2)-NLatent+1:NVarA*(NLatent+1+NVarM)+sum(1:NVarA+NLatent,2)-NLatent+NVarS); % b scale
    bstr = reshape(B(NVarA*(NLatent+NVarM+1)+sum(1:NVarA+NLatent,2)-NLatent+NVarS+1:(NVarA+NVarStr)*NLatent+NVarA*(1+NVarM)+sum(1:NVarA+NLatent,2)-NLatent+NVarS),[NVarStr,NLatent]); % b równania struktury
    bmea = B((NVarA+NVarStr)*NLatent+NVarA*(1+NVarM)+sum(1:NVarA+NLatent,2)-NLatent+NVarS+1:end); % b measurement
    err = VC*err_sliced;
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
    LV_tmp = LV_tmp + err(NVarA+1:end,:); % Latent x NRep*NP
else
    LV_tmp = LV_tmp + err_sliced(NVarA+1:end,:); % Latent x NRep*NP
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
    b_mtx = ba + bl*LV + err(1:NVarA,:); % NVarA x NRep*NP  
else
    b_mtx = ba + bl*LV + VC*err_sliced(1:NVarA,:); % NVarA x NRep*NP
end  

if sum(Dist == 1) > 0 % Log - normal
    b_mtx(Dist == 1,:) = exp(b_mtx(Dist == 1,:));
elseif sum(Dist == 2) > 0 % Spike       
    b_mtx(Dist == 2,:) = max(b_mtx(Dist == 2,:),0);
end

b_score = b_mtx;
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

probs = zeros(NP,NRep);

    if any(isnan(Xa(:))) == 0 % faster version for complete dataset
        parfor n = 1:NP
            Yy_n = Y(:,n) == 1;
            Xa_n = Xa(:,:,n);
            U = reshape(Xa_n*b_mtx_sliced(:,:,n),[NAlt,NCT,NRep]);
            U = exp(U - max(U)); % rescale utility to avoid exploding
            U_sum = reshape(sum(U,1),[NCT,NRep]);
            U_selected = reshape(U(Yy_n(:,ones(NRep,1))),[NCT,NRep]);   
            probs(n,:) = prod(U_selected./U_sum,1);
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
            probs(n,:) = prod(U_selected./U_sum,1);
        end
    end 
    
fx = mean(probs,2); % NP x 1
Score = zeros(NP,NVarA);
if WTP_space > 0
    b_score = reshape(b_score,[NVarA,NRep,NP]);
    for j = 1:NVarA
        bx = squeeze(b_score(j,:,:))'; % NP x NRep
        Score(:,j) = mean(probs.*bx,2)./fx;
    end
else
    b_score = reshape(b_score,[NVarA,NRep,NP]);
    fee = squeeze(b_score(end,:,:))';
    for j = 1:NVarA-1
        bx = squeeze(b_score(j,:,:))'; % NP x NRep
        Score(:,j) = mean(probs.*bx./fee,2)./fx;
    end
    Score(:,end) = mean(probs.*fee,2)./fx;
end    

