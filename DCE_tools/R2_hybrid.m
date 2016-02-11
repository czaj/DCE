function R2 = R2_hybrid(Y,Xa,X_str,Xc,Xm,MissingInd,err_sliced,EstimOpt,B,type)

NAltMiss = EstimOpt.NAltMiss; 

if any(MissingInd == 1) % In case of some missing data
   idx = sum(reshape(MissingInd,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP)) == EstimOpt.NAlt;
   idx = sum(reshape(idx, EstimOpt.NCT, EstimOpt.NP),1)'; % no. of missing NCT for every respondent
   idx = EstimOpt.NCT - idx;
end

if type == 0 % HMNL
    ba = B(1:EstimOpt.NVarA); % b atrybutów
    bl = reshape(B(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NLatent+1)), EstimOpt.NVarA, EstimOpt.NLatent); % b interakcji z LV
    bstr = reshape(B(EstimOpt.NVarA*(EstimOpt.NLatent+1)+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+EstimOpt.NVarA), EstimOpt.NVarStr, EstimOpt.NLatent); % b równania struktury
   % bmea = B((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+EstimOpt.NVarA+1:end); % b measurement

    LV_tmp = X_str*bstr; % NP x NLatent
    LV_tmp = reshape(permute(LV_tmp(:,:, ones(EstimOpt.NRep,1)),[2 3 1]), EstimOpt.NLatent, EstimOpt.NRep*EstimOpt.NP); % NLatent*NRep*NP

    LV = LV_tmp + err_sliced; % NLatent x NRep*NP
    mLV = mean(LV,2);
    sLV = std(LV,0,2);
    % LV = LV - mLV(:,ones(1,size(LV,2))); % normalilzing for 0 mean
    LV = (LV - mLV(:,ones(1,size(LV,2))))./sLV(:,ones(1,size(LV,2))); % normalilzing for 0 mean and std

    b_mtx = ba(:,ones(EstimOpt.NRep*EstimOpt.NP,1)) + bl*LV; % NVarA x NRep*NP
    if EstimOpt.WTP_space > 0
        b_mtx(1:end-EstimOpt.WTP_space,:) = b_mtx(1:end-EstimOpt.WTP_space,:).*b_mtx(EstimOpt.WTP_matrix,:);
    end

    probs = zeros(EstimOpt.NP,EstimOpt.NRep);

    if any(isnan(Xa(:))) == 0 % faster version for complete dataset
        parfor n = 1:EstimOpt.NP
            U = reshape(Xa(:,:,n)*b_mtx(:,((n-1)*EstimOpt.NRep+1):n*EstimOpt.NRep),EstimOpt.NAlt,EstimOpt.NCT, EstimOpt.NRep);
            U_max = max(U);
            U = exp(U - U_max(ones(EstimOpt.NAlt,1),:,:)); % rescale utility to avoid exploding
            U_sum = reshape(sum(U,1),EstimOpt.NCT,EstimOpt.NRep);
            U_seleNCTed = reshape(U(Y(:,n*ones(EstimOpt.NRep,1))==1),EstimOpt.NCT,EstimOpt.NRep);   
            probs(n,:) = prod(U_seleNCTed ./ U_sum,1);
        end;
    else
        parfor n = 1:EstimOpt.NP
%             U = reshape(Xa(~isnan(Y(:,n)),:,n)*b_mtx(:,((n-1)*EstimOpt.NRep+1):n*EstimOpt.NRep),EstimOpt.NAlt,EstimOpt.NCT-sum(isnan(Y(1:EstimOpt.NAlt:end,n))),EstimOpt.NRep);
            U = reshape(Xa(~isnan(Y(:,n)),:,n)*b_mtx(:,((n-1)*EstimOpt.NRep+1):n*EstimOpt.NRep),numel(Y(~isnan(Y(:,n))))./(EstimOpt.NCT-sum(isnan(Y(1:EstimOpt.NAlt:end,n)))),EstimOpt.NCT-sum(isnan(Y(1:EstimOpt.NAlt:end,n))),EstimOpt.NRep);
            U_max = max(U);
%             U = exp(U - U_max(ones(EstimOpt.NAlt,1),:,:)); % NAlt x NCT - NaNs x NRep
            U = exp(U - U_max(ones(numel(Y(~isnan(Y(:,n))))./(EstimOpt.NCT-sum(isnan(Y(1:EstimOpt.NAlt:end,n)))),1),:,:)); % NAlt x NCT - NaNs x NRep
            U_sum = reshape(nansum(U,1),EstimOpt.NCT-sum(isnan(Y(1:EstimOpt.NAlt:end,n))),EstimOpt.NRep);
            U_seleNCTed = reshape(U(Y(~isnan(Y(:,n)),n*ones(EstimOpt.NRep,1))==1),EstimOpt.NCT-sum(isnan(Y(1:EstimOpt.NAlt:end,n))),EstimOpt.NRep);   
            probs(n,:) = prod(U_seleNCTed ./ U_sum,1);
        end;
    end
    
    if any(MissingInd == 1) % In case of some missing data
        R2 = mean(mean(probs,2).^(1./idx),1);
    else
        R2 = mean(mean(probs,2).^(1/EstimOpt.NCT),1);
    end
    
elseif type == 1 % HLC
    
    beta = reshape(B(1:EstimOpt.NClass*EstimOpt.NVarA), EstimOpt.NVarA, EstimOpt.NClass);
    if EstimOpt.WTP_space > 0
        beta(1:end - EstimOpt.WTP_space,:) = beta(1:end - EstimOpt.WTP_space,:).*beta(EstimOpt.WTP_matrix,:);
    end
    U = reshape(Xa*beta,EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP*EstimOpt.NClass);
    maxU = max(U);
    U = exp(U - maxU(ones(EstimOpt.NAlt,1),:));% NAlt*NCT*NP x NClass U(isnan(U)) = 0;... % do not include alternatives which were not available
    U_sum = nansum(U,1);... 
    P = reshape(sum(Y .* U ./ U_sum(ones(EstimOpt.NAlt,1),:),1),EstimOpt.NCT,EstimOpt.NP*EstimOpt.NClass);... % NCT x NP*NClass
    P(isnan(reshape(Y(1,:),EstimOpt.NCT, EstimOpt.NP*EstimOpt.NClass))) = 1;
    probs = prod(P,1); % 1 x NP*NClass
    probs = reshape(probs,EstimOpt.NP,EstimOpt.NClass);
    probs = probs(:,:, ones(EstimOpt.NRep,1));
    probs = permute(probs, [3 2 1]); %NRep x NClass x NP

    bNClass = [B(EstimOpt.NClass*EstimOpt.NVarA+1:EstimOpt.NClass*EstimOpt.NVarA + (EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1));zeros(EstimOpt.NVarC + EstimOpt.NLatent,1)];
    bNClass = reshape(bNClass, EstimOpt.NVarC+EstimOpt.NLatent, EstimOpt.NClass);
    bstr = reshape(B(EstimOpt.NClass*EstimOpt.NVarA + (EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1)+1:EstimOpt.NClass*EstimOpt.NVarA + (EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1)+EstimOpt.NLatent*EstimOpt.NVarStr),EstimOpt.NVarStr, EstimOpt.NLatent);

    LV_tmp = X_str*bstr; % NP x NLatent
    LV_tmp = reshape(permute(LV_tmp(:,:, ones(EstimOpt.NRep,1)),[2 3 1]), EstimOpt.NLatent, EstimOpt.NRep*EstimOpt.NP);
    LV_tmp = LV_tmp + err_sliced; % NLatent x NRep*NP

    mLV = mean(LV_tmp,2);
    sLV = std(LV_tmp,0,2);
    LV = (LV_tmp - mLV(:,ones(1,size(LV_tmp,2))))./sLV(:,ones(1,size(LV_tmp,2))); % normalilzing for 0 mean and std

    p = zeros(EstimOpt.NP, EstimOpt.NRep);
    for i = 1:EstimOpt.NP
        Xc_i = Xc(i,:);
        XXc = [Xc_i(ones(EstimOpt.NRep,1),:), LV(:,(i-1)*EstimOpt.NRep+1:i*EstimOpt.NRep)'];
        PNClass = exp(XXc*bNClass); % NRep x NClass
        PNClass_sum = sum(PNClass,2);...
        PNClass = PNClass./PNClass_sum(:,ones(EstimOpt.NClass,1)); % NRep x NClass
        p(i,:) = sum(PNClass.*probs(:,:,i),2)'; 

    end
    if any(MissingInd == 1) % In case of some missing data
        R2 = mean(mean(p,2).^(1./idx),1);
    else
        R2 = mean(mean(p,2).^(1/EstimOpt.NCT),1);
    end
    
elseif type ==2 % HMXL
    
    if EstimOpt.FullCov == 0
        ba = B(1:EstimOpt.NVarA); % b atrybutów
        bv = B(EstimOpt.NVarA+1:2*EstimOpt.NVarA);
        VC = diag(bv.^2);
        bm = reshape(B(2*EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarM+2)), EstimOpt.NVarA, EstimOpt.NVarM); % b mean covariates
        bl = reshape(B((2+EstimOpt.NVarM)*EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NLatent+2+EstimOpt.NVarM)), EstimOpt.NVarA, EstimOpt.NLatent); % b interakcji z LV
        bstr = reshape(B(EstimOpt.NVarA*(EstimOpt.NLatent+2+EstimOpt.NVarM)+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA), EstimOpt.NVarStr, EstimOpt.NLatent); % b równania struktury
    else
        ba = B(1:EstimOpt.NVarA); % b atrybutów
        bv = B(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA,2));
        VC = tril(ones(EstimOpt.NVarA));
        VC(VC==1) = bv;
        bm = reshape(B(EstimOpt.NVarA+sum(1:EstimOpt.NVarA,2)+1:EstimOpt.NVarA*(EstimOpt.NVarM+1)+sum(1:EstimOpt.NVarA,2)), EstimOpt.NVarA, EstimOpt.NVarM); % b mean covariates
        bl = reshape(B(EstimOpt.NVarA*(1+EstimOpt.NVarM)+sum(1:EstimOpt.NVarA,2)+1:EstimOpt.NVarA*(EstimOpt.NLatent+1+EstimOpt.NVarM)+sum(1:EstimOpt.NVarA,2)), EstimOpt.NVarA, EstimOpt.NLatent); % b interakcji z LV
        bstr = reshape(B(EstimOpt.NVarA*(EstimOpt.NLatent++EstimOpt.NVarM+1)+sum(1:EstimOpt.NVarA,2)+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+EstimOpt.NVarA*(1+EstimOpt.NVarM)+sum(1:EstimOpt.NVarA,2)), EstimOpt.NVarStr, EstimOpt.NLatent); % b równania struktury
    end

    LV_tmp = X_str*bstr; % NP x NLatent
    LV_tmp = reshape(permute(LV_tmp(:,:, ones(EstimOpt.NRep,1)),[2 3 1]), EstimOpt.NLatent, EstimOpt.NRep*EstimOpt.NP); % NLatent*NRep*NP


    LV = LV_tmp + err_sliced(EstimOpt.NVarA+1:end,:); % NLatent x NRep*NP
    mLV = mean(LV,2);
    sLV = std(LV,0,2);
    % LV = LV - mLV(:,ones(1,size(LV,2))); % normalilzing for 0 mean
    LV = (LV - mLV(:,ones(1,size(LV,2))))./sLV(:,ones(1,size(LV,2))); % normalilzing for 0 mean and std

    if EstimOpt.NVarM > 0
        ba = ba(:,ones(EstimOpt.NP,1))+bm*Xm; % NVarA x NP
        ba = reshape(permute(ba(:,:, ones(EstimOpt.NRep,1)), [1 3 2]), EstimOpt.NVarA, EstimOpt.NRep*EstimOpt.NP);
    else
        ba = ba(:, ones(EstimOpt.NP*EstimOpt.NRep,1));
    end

    b_mtx = ba + bl*LV + VC*err_sliced(1:EstimOpt.NVarA,:); % NVarA x NRep*NP
    if sum(EstimOpt.Dist==1) > 0; % Log - normal
        b_mtx(EstimOpt.Dist==1,:) = exp(b_mtx(EstimOpt.Dist == 1,:));
    elseif sum(EstimOpt.Dist==2) > 0; % Spike       
        b_mtx(EstimOpt.Dist==2,:) = max(b_mtx(EstimOpt.Dist==2,:),0);
    end

    if EstimOpt.WTP_space > 0
        b_mtx(1:end-EstimOpt.WTP_space,:) = b_mtx(1:end-EstimOpt.WTP_space,:).*b_mtx(EstimOpt.WTP_matrix,:);
    end

    probs = zeros(EstimOpt.NP,EstimOpt.NRep);

    if any(isnan(Xa(:))) == 0 % faster version for complete dataset

        parfor n = 1:EstimOpt.NP

            U = reshape(Xa(:,:,n)*b_mtx(:,((n-1)*EstimOpt.NRep+1):n*EstimOpt.NRep),EstimOpt.NAlt,EstimOpt.NCT, EstimOpt.NRep);
            U_max = max(U);
            U = exp(U - U_max(ones(EstimOpt.NAlt,1),:,:)); % rescale utility to avoid exploding
            U_sum = reshape(sum(U,1),EstimOpt.NCT,EstimOpt.NRep);
            U_seleNCTed = reshape(U(Y(:,n*ones(EstimOpt.NRep,1))==1),EstimOpt.NCT,EstimOpt.NRep);   
            probs(n,:) = prod(U_seleNCTed ./ U_sum,1);
        end;
    else
        parfor n = 1:EstimOpt.NP
%             U = reshape(Xa(~isnan(Y(:,n)),:,n)*b_mtx(:,((n-1)*EstimOpt.NRep+1):n*EstimOpt.NRep),EstimOpt.NAlt,EstimOpt.NCT-sum(isnan(Y(1:EstimOpt.NAlt:end,n))),EstimOpt.NRep);
            U = reshape(Xa(~isnan(Y(:,n)),:,n)*b_mtx(:,((n-1)*EstimOpt.NRep+1):n*EstimOpt.NRep),NAltMiss(n),EstimOpt.NCT-sum(isnan(Y(1:EstimOpt.NAlt:end,n))),EstimOpt.NRep);
            U_max = max(U);
%             U = exp(U - U_max(ones(EstimOpt.NAlt,1),:,:)); % NAlt x NCT - NaNs x NRep
            U = exp(U - U_max(ones(NAltMiss(n),1),:,:));
            U_sum = reshape(nansum(U,1),EstimOpt.NCT-sum(isnan(Y(1:EstimOpt.NAlt:end,n))),EstimOpt.NRep);
            U_seleNCTed = reshape(U(Y(~isnan(Y(:,n)),n*ones(EstimOpt.NRep,1))==1),EstimOpt.NCT-sum(isnan(Y(1:EstimOpt.NAlt:end,n))),EstimOpt.NRep);   
            probs(n,:) = prod(U_seleNCTed ./ U_sum,1);
        end;
    end
    
    if any(MissingInd == 1) % In case of some missing data
        R2 = mean(mean(probs,2).^(1./idx),1);
    else
        R2 = mean(mean(probs,2).^(1/EstimOpt.NCT),1);
    end
end
