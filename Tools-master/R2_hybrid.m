function R2 = R2_hybrid(Y,Xa,X_str,Xc,Xm,Xs,MissingInd,err_sliced,EstimOpt,B,type)

NP = EstimOpt.NP;
NRep = EstimOpt.NRep;
NAlt = EstimOpt.NAlt;
NCT = EstimOpt.NCT;
NVarA = EstimOpt.NVarA;
NLatent = EstimOpt.NLatent;
NVarStr = EstimOpt.NVarStr;
NCTMiss = EstimOpt.NCTMiss;

if type == 0 || type == 2
    NVarM = EstimOpt.NVarM;
    NVarS = EstimOpt.NVarS;
    NAltMiss = EstimOpt.NAltMiss;
    NAltMissIndExp = EstimOpt.NAltMissIndExp;
    ScaleLV = EstimOpt.ScaleLV;
end

if type == 0 % HMNL
    ba = B(1:NVarA); % b atrybutów
    bm = reshape(B(NVarA+1:NVarA*(1+NVarM)),[NVarA,NVarM]);
    bl = reshape(B(NVarA*(1+NVarM)+1:NVarA*(NLatent+NVarM+1)),[NVarA,NLatent]); % b interakcji z LV
    bs = B(NVarA*(NLatent+NVarM+1)+1:NVarA*(NLatent+NVarM+1)+NVarS); % b równania struktury
    bstr = reshape(B(NVarA*(NLatent+NVarM+1)+NVarS+1:(NVarA+NVarStr)*NLatent+NVarA*(1+NVarM)+NVarS),[NVarStr,NLatent]); % b równania struktury

    if ScaleLV == 1
        bsLV = bs(NVarS-NLatent+1:end)';
        bs = bs(1:NVarS-NLatent); 
        NVarS = NVarS - NLatent;
    end
    
    LV_tmp = X_str*bstr; % NP x NLatent
    LV_tmp = reshape(permute(LV_tmp(:,:, ones(NRep,1)),[2 3 1]),[NLatent,NRep*NP]); % NLatent*NRep*NP

    LV = LV_tmp + err_sliced; % NLatent x NRep*NP
    LV = (LV - mean(LV,2))./std(LV,0,2); % normalilzing for 0 mean and std
    if NVarM > 0
        ba = ba + bm*Xm'; % NVarA x NP
        ba = reshape(permute(ba(:,:,ones(EstimOpt.NRep,1)),[1 3 2]),[NVarA,NRep*NP]);
    else
        ba = ba(:,ones(EstimOpt.NP*EstimOpt.NRep,1)); 
    end
    b_mtx = ba + bl*LV;  % NVarA x NRep*NP  

    if EstimOpt.WTP_space > 0
        b_mtx(1:end-EstimOpt.WTP_space,:) = b_mtx(1:end-EstimOpt.WTP_space,:).*b_mtx(EstimOpt.WTP_matrix,:);
    end
    if any(EstimOpt.Dist ~= 0)
        b_mtx(EstimOpt.Dist == 1,:) = exp(b_mtx(EstimOpt.Dist == 1,:)); 
        b_mtx(EstimOpt.Dist == 2,:) = max(b_mtx(EstimOpt.Dist == 2,:),0);
    end
    if ScaleLV == 1
        ScaleLVX = exp(bsLV*LV);
        b_mtx = ScaleLVX.*b_mtx; 
    end

    if NVarS > 0
       Scale = reshape(exp(Xs*bs),[EstimOpt.NAlt*EstimOpt.NCT,1,EstimOpt.NP]);
       Xa = Xa.*Scale;
    end
    probs = zeros(EstimOpt.NP,EstimOpt.NRep);
    b_mtx = reshape(b_mtx,[NVarA,EstimOpt.NRep,NP]);
    
    if any(isnan(Xa(:))) == 0 % faster version for complete dataset
        parfor n = 1:EstimOpt.NP
            b_mtx_n = b_mtx(:,:,n);
            Xa_n = Xa(:,:,n);
            Yy_n = Y(:,n) == 1;
            U = reshape(Xa_n*b_mtx_n,[NAlt,NCT,NRep]);
            U = exp(U - max(U,[],1)); % rescale utility to avoid exploding
            U_sum = reshape(sum(U,1),[1,NCT,NRep]);
            U_prob = U./U_sum; % NAlt x NCT x NRep
            probs(n,:) = prod(reshape(U_prob(Yy_n(:,ones(NRep,1))),[NCT,NRep]),1); % 1 x NRep
        end
    else
        parfor n = 1:NP
            Xa_n = Xa(:,:,n);
            Yy_n = Y(:,n) == 1;
            YnanInd = ~isnan(Y(:,n));
            b_mtx_n = b_mtx(:,:,n);
            NAltMissIndExp_n = NAltMissIndExp(:,n);
            NAltMissIndExp_n = NAltMissIndExp_n(YnanInd);
            if var(NAltMissIndExp_n(NAltMissIndExp_n > 0)) == 0 
                U = reshape(Xa_n(YnanInd,:)*b_mtx_n,[NAltMiss(n),NCTMiss(n),NRep]);
                U = exp(U - max(U,[],1)); % NAlt x NCT - NaNs x NRep
                U_sum = reshape(sum(U,1),[1,NCTMiss(n),NRep]);
                U_prob = U./U_sum; % NAlt x NCT x NRep
                probs(n,:) = prod(reshape(U_prob(Yy_n(YnanInd,ones(NRep,1))),[NCTMiss(n),NRep]),1); % 1 x NRep
            else
                U = Xa_n(YnanInd,:)*b_mtx_n;
                Uniq = unique(NAltMissIndExp_n);
                U_prob = zeros(size(U,1),1,NRep);

                for i = 1:length(Uniq)
                    U_tmp = U(NAltMissIndExp_n == Uniq(i),:);
                    U_tmp = reshape(U_tmp,[Uniq(i),size(U_tmp,1)/Uniq(i),NRep]);
                    U_tmp = exp(U_tmp - max(U_tmp));
                    U_sum = reshape(sum(U_tmp,1),[1,size(U_tmp,2),NRep]);
                    U_tmp = U_tmp./U_sum;
                    U_prob(NAltMissIndExp_n == Uniq(i),:,:) = reshape(U_tmp,[size(U_tmp,2)*Uniq(i),1,NRep]);
                end
                probs(n,:) = prod(reshape(U_prob(Yy_n(YnanInd,ones(NRep,1))),[NCTMiss(n),NRep]),1);  % 1 x NRep
            end
        end
    end
    
    if any(MissingInd == 1) % In case of some missing data
        R2 = mean(mean(probs,2).^(1./NCTMiss),1);
    else
        R2 = mean(mean(probs,2).^(1/NCT),1);
    end
    
elseif type == 1 % HLC
    
    NClass = EstimOpt.NClass;
    NVarC = EstimOpt.NVarC;
    
    beta = reshape(B(1:NClass*NVarA),[NVarA,NClass]);
    if EstimOpt.WTP_space > 0
        beta(1:end - EstimOpt.WTP_space,:) = beta(1:end - EstimOpt.WTP_space,:).*beta(EstimOpt.WTP_matrix,:);
    end
    U = reshape(Xa*beta,[NAlt,NCT*NP*NClass]);
    U = exp(U - max(U));% NAlt*NCT*NP x NClass U(isnan(U)) = 0;... % do not include alternatives which were not available
    P = reshape(sum(Y .* U ./ sum(U,1),1),[NCT,NP*NClass]); % NCT x NP*NClass
    P(isnan(reshape(Y(1,:),[NCT,NP*NClass]))) = 1;
    probs = prod(P,1); % 1 x NP*NClass
    probs = reshape(probs,[NP,NClass]);
    probs = probs(:,:,ones(NRep,1));
    probs = permute(probs, [3 2 1]); %NRep x NClass x NP

    bNClass = [B(NClass*NVarA+1:NClass*NVarA + (NVarC+NLatent)*(NClass-1));zeros(NVarC + NLatent,1)];
    bNClass = reshape(bNClass,[NVarC+NLatent,NClass]);
    bstr = reshape(B(NClass*NVarA + (NVarC+NLatent)*(NClass-1)+1:NClass*NVarA + (NVarC+NLatent)*(NClass-1)+NLatent*NVarStr),[NVarStr,NLatent]);

    LV_tmp = X_str*bstr; % NP x NLatent
    LV_tmp = reshape(permute(LV_tmp(:,:, ones(EstimOpt.NRep,1)),[2 3 1]),[NLatent,NRep*NP]);
    LV_tmp = LV_tmp + err_sliced; % NLatent x NRep*NP

    LV = (LV_tmp - mean(LV_tmp,2))./std(LV_tmp,0,2); % normalilzing for 0 mean and std

    p = zeros(NP, NRep);
    parfor i = 1:NP
        Xc_i = Xc(i,:);
        XXc = [Xc_i(ones(NRep,1),:), LV(:,(i-1)*NRep+1:i*NRep)'];
        PNClass = exp(XXc*bNClass); % NRep x NClass
        PNClass_sum = sum(PNClass,2);
        PNClass = PNClass./PNClass_sum; % NRep x NClass
        p(i,:) = sum(PNClass.*probs(:,:,i),2)'; 

    end
    if any(MissingInd == 1) % In case of some missing data
        R2 = mean(mean(p,2).^(1./NCTMiss),1);
    else
        R2 = mean(mean(p,2).^(1/NCT),1);
    end
    
elseif type == 2 % HMXL
    
    if EstimOpt.FullCov == 0
        ba = B(1:NVarA); % b atrybutów
        bv = B(NVarA+1:2*NVarA);
        VC = diag(bv);
        bm = reshape(B(2*NVarA+1:NVarA*(NVarM+2)),[NVarA,NVarM]); % b mean covariates
        bl = reshape(B((2+NVarM)*NVarA+1:NVarA*(NLatent+2+NVarM)),[NVarA,NLatent]); % b interakcji z LV
        bs = B(NVarA*(NLatent+NVarM+2)+1:NVarA*(NLatent+NVarM+2)+NVarS); % b scale
        bstr = reshape(B(NVarA*(NLatent+2+NVarM)+NVarS+1:(NVarA+NVarStr)*NLatent+(2+NVarM)*NVarA+NVarS),[NVarStr,NLatent]); % b równania struktury
        %bstr = reshape(B(EstimOpt.NVarA*(EstimOpt.NLatent+2+EstimOpt.NVarM)+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA), EstimOpt.NVarStr, EstimOpt.NLatent); % b równania struktury
    else
        ba = B(1:NVarA); % b atrybutów
        bv = B(NVarA+1:NVarA+sum(1:NVarA,2));
        VC = tril(ones(NVarA));
        VC(VC == 1) = bv;
        bm = reshape(B(NVarA+sum(1:NVarA,2)+1:NVarA*(NVarM+1)+sum(1:NVarA,2)),[NVarA,NVarM]); % b mean covariates
        bl = reshape(B(NVarA*(1+NVarM)+sum(1:NVarA,2)+1:NVarA*(NLatent+1+NVarM)+sum(1:NVarA,2)),[NVarA,NLatent]); % b interakcji z LV
        bs = B(NVarA*(NLatent+1+NVarM)+sum(1:NVarA,2)+1:NVarA*(NLatent+1+NVarM)+sum(1:NVarA,2)+NVarS); % b scale
        bstr = reshape(B(NVarA*(NLatent+NVarM+1)+sum(1:NVarA,2)+NVarS+1:(NVarA+NVarStr)*NLatent+NVarA*(1+NVarM)+sum(1:NVarA,2)+NVarS),NVarStr,NLatent); % b równania struktury
        %bstr = reshape(B(EstimOpt.NVarA*(EstimOpt.NLatent++EstimOpt.NVarM+1)+sum(1EstimOpt.NVarA,2)+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+EstimOpt.NVarA*(1+EstimOpt.NVarM)+sum(1:EstimOpt.NVarA,2)), EstimOpt.NVarStr, EstimOpt.NLatent); % b równania struktury
    end
    if ScaleLV == 1
       bsLV = bs(NVarS-NLatent+1:end)';
       bs = bs(1:NVarS-NLatent); 
       NVarS = NVarS - NLatent;
    end
    LV_tmp = X_str*bstr; % NP x NLatent
    LV_tmp = reshape(permute(LV_tmp(:,:, ones(NRep,1)),[2 3 1]),[NLatent,NRep*NP]); % NLatent*NRep*NP

    LV = LV_tmp + err_sliced(NVarA+1:end,:); % NLatent x NRep*NP
    LV = (LV - mean(LV,2))./std(LV,0,2); % normalilzing for 0 mean and std

    if NVarM > 0
        ba = ba + bm*Xm; % NVarA x NP
        ba = reshape(permute(ba(:,:, ones(NRep,1)), [1 3 2]),[NVarA,NRep*NP]);
    else
        ba = ba(:, ones(NP*NRep,1));
    end

    b_mtx = ba + bl*LV + VC*err_sliced(1:NVarA,:); % NVarA x NRep*NP
    if sum(EstimOpt.Dist==1) > 0 % Log - normal
        b_mtx(EstimOpt.Dist==1,:) = exp(b_mtx(EstimOpt.Dist == 1,:));
    elseif sum(EstimOpt.Dist==2) > 0 % Spike       
        b_mtx(EstimOpt.Dist==2,:) = max(b_mtx(EstimOpt.Dist==2,:),0);
    end

    if EstimOpt.WTP_space > 0
        b_mtx(1:end-EstimOpt.WTP_space,:) = b_mtx(1:end-EstimOpt.WTP_space,:).*b_mtx(EstimOpt.WTP_matrix,:);
    end
    if ScaleLV == 1
        ScaleLVX = exp(bsLV*LV);
        b_mtx = ScaleLVX.*b_mtx; 
    end
    if NVarS > 0
       Scale = reshape(exp(Xs*bs),[NAlt*NCT,1,NP]);
       Xa = Xa.*Scale;
    end
    probs = zeros(NP,NRep);
    b_mtx = reshape(b_mtx,[NVarA,NRep,NP]);
    
    if any(isnan(Xa(:))) == 0 % faster version for complete dataset

        parfor n = 1:NP
            b_mtx_n = b_mtx(:,:,n);
            Xa_n = Xa(:,:,n);
            Yy_n = Y(:,n) == 1;
            U = reshape(Xa_n*b_mtx_n,[NAlt,NCT,NRep]);
            U = exp(U - max(U,[],1)); % rescale utility to avoid exploding
            U_sum = reshape(sum(U,1),[1,NCT,NRep]);
            U_prob = U./U_sum; % NAlt x NCT x NRep
            probs(n,:) = prod(reshape(U_prob(Yy_n(:,ones(NRep,1))),[NCT,NRep]),1); % 1 x NRep
        end
    else
        parfor n = 1:EstimOpt.NP
            Xa_n = Xa(:,:,n);
            Yy_n = Y(:,n) == 1;
            YnanInd = ~isnan(Y(:,n));
            b_mtx_n = b_mtx(:,:,n);
            NAltMissIndExp_n = NAltMissIndExp(:,n);
            NAltMissIndExp_n = NAltMissIndExp_n(YnanInd);
            if var(NAltMissIndExp_n(NAltMissIndExp_n > 0)) == 0 
                U = reshape(Xa_n(YnanInd,:)*b_mtx_n,[NAltMiss(n),NCTMiss(n),NRep]);
                U = exp(U - max(U,[],1)); % NAlt x NCT - NaNs x NRep
                U_sum = reshape(sum(U,1),[1,NCTMiss(n),NRep]);
                U_prob = U./U_sum; % NAlt x NCT x NRep
                probs(n,:) = prod(reshape(U_prob(Yy_n(YnanInd,ones(NRep,1))),[NCTMiss(n),NRep]),1); % 1 x NRep
            else
                U = Xa_n(YnanInd,:)*b_mtx_n;
                Uniq = unique(NAltMissIndExp_n);
                U_prob = zeros(size(U,1),1,NRep);

                for i = 1:length(Uniq)
                    U_tmp = U(NAltMissIndExp_n == Uniq(i),:);
                    U_tmp = reshape(U_tmp,[Uniq(i),size(U_tmp,1)/Uniq(i),NRep]);
                    U_tmp = exp(U_tmp - max(U_tmp));
                    U_sum = reshape(sum(U_tmp,1),[1,size(U_tmp,2),NRep]);
                    U_tmp = U_tmp./U_sum;
                    U_prob(NAltMissIndExp_n == Uniq(i),:,:) = reshape(U_tmp,[size(U_tmp,2)*Uniq(i),1,NRep]);
                end
                probs(n,:) = prod(reshape(U_prob(Yy_n(YnanInd,ones(NRep,1))),[NCTMiss(n),NRep]),1);  % 1 x NRep
            end
        end
    end
    
    if any(MissingInd == 1) % In case of some missing data
        R2 = mean(mean(probs,2).^(1./NCTMiss),1);
    else
        R2 = mean(mean(probs,2).^(1/NCT),1);
    end
end
