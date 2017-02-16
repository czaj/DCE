function [f,g] = LL_hmnl(Y,Xa,Xm,Xs,Xstr,Xmea,Xmea_exp,err_sliced,EstimOpt,B)

% save tmp_LL_hmnl
% return

NAlt = EstimOpt.NAlt;
NCT = EstimOpt.NCT;
NP = EstimOpt.NP;
NRep = EstimOpt.NRep;
NVarA = EstimOpt.NVarA;
NVarS = EstimOpt.NVarS;
NVarM = EstimOpt.NVarM;
WTP_space = EstimOpt.WTP_space;
WTP_matrix = EstimOpt.WTP_matrix;
NCTMiss = EstimOpt.NCTMiss;
NAltMiss = EstimOpt.NAltMiss;
Dist = EstimOpt.Dist;
NLatent = EstimOpt.NLatent;
NVarStr = EstimOpt.NVarStr;
% NVarMea = EstimOpt.NVarMea;
MeaMatrix = EstimOpt.MeaMatrix;
MeaSpecMatrix = EstimOpt.MeaSpecMatrix;
NVarMeaExp = EstimOpt.NVarMeaExp;
MeaExpMatrix = EstimOpt.MeaExpMatrix;
RealMin = EstimOpt.RealMin;
ScaleLV = EstimOpt.ScaleLV;

ba = B(1:NVarA);  % b atrybutów
bm = reshape(B(NVarA+1:NVarA*(1+NVarM)), NVarA, NVarM);
bl = reshape(B(NVarA*(1+NVarM)+1:NVarA*(NLatent+NVarM+1)), NVarA, NLatent); % b interakcji z LV
bs = B(NVarA*(NLatent+NVarM+1)+1:NVarA*(NLatent+NVarM+1)+NVarS); % b równania struktury
bstr = reshape(B(NVarA*(NLatent+NVarM+1)+NVarS+1:(NVarA+NVarStr)*NLatent+NVarA*(1+NVarM)+NVarS), NVarStr, NLatent); % b równania struktury
bmea = B((NVarA+NVarStr)*NLatent+NVarA*(1+NVarM)+NVarS+1:end); % b measurement

if ScaleLV == 1
    bsLV = bs(NVarS - NLatent+1:end)';
    bs = bs(1:NVarS - NLatent);
    NVarS = NVarS - NLatent;
else
    bsLV = zeros(1,0);
end
LV_tmp = Xstr*bstr; % NP x NLatent
LV_tmp = reshape(permute(LV_tmp(:,:, ones(NRep,1)),[2 3 1]), NLatent, NRep*NP);
LV_tmp = LV_tmp + err_sliced; % NLatent x NRep*NP
mLV = mean(LV_tmp,2);
sLV = std(LV_tmp,0,2);
LV = (LV_tmp - mLV(:,ones(1,size(LV_tmp,2))))./sLV(:,ones(1,size(LV_tmp,2))); % normalilzing for 0 mean and std

if NVarM > 0
    ba = ba(:,ones(NP,1))+bm*Xm'; % NVarA x NP
    ba = reshape(permute(ba(:,:,ones(NRep,1)),[1 3 2]),NVarA,NRep*NP);
else
    ba = ba(:,ones(NP*NRep,1));
end
b_mtx = ba + bl*LV;  % NVarA x NRep*NP

if any(Dist ~= 0)
    %     if sum(Dist == 1) > 0; % Log - normal
    b_mtx(Dist == 1,:) = exp(b_mtx(Dist == 1,:));
    %     elseif sum(Dist == 2) > 0; % Spike
    b_mtx(Dist == 2,:) = max(b_mtx(Dist == 2,:),0);
    %     end
end


if WTP_space > 0
    if nargout == 2
        b_mtx_grad = reshape(b_mtx,NVarA,NRep,NP);% needed for gradient calculation in WTP_space
    end
    b_mtx(1:end-WTP_space,:) = b_mtx(1:end-WTP_space,:).*b_mtx(WTP_matrix,:);
elseif nargout == 2
    if ScaleLV == 1 && any(Dist == 1)
        b_mtx_grad = reshape(b_mtx,NVarA,NRep,NP);% needed for gradient calculation in WTP_space
    else
        b_mtx_grad = zeros(0,0,NP);
    end
end
if ScaleLV == 1
    ScaleLVX =exp(bsLV*LV);
    b_mtx = bsxfun(@times,ScaleLVX, b_mtx);
    ScaleLVX = permute(reshape(ScaleLVX, 1,NRep, NP), [1 3 2]);
else
    ScaleLVX = zeros(0, NP, 0);
end

b_mtx_sliced = reshape(b_mtx,[NVarA,NRep,NP]); % NVarA x NRep x NP
if NVarS > 0
    Scale = reshape(exp(Xs*bs), NAlt*NCT, 1, NP);
    Xa = bsxfun(@times, Xa, Scale);
end


probs = zeros(NP,NRep);

%for n = 1:NP
%   U = reshape(exp(Xa(:,:,n)*b_mtx(:,((n-1)*NRep+1):n*NRep)),NAlt,NCT,NRep);
%U(isnan(U)) = 0; % skip missing ALT

%  U_sum = sum(U,1);
%  P = sum(Y(:,:,:,n).*U ./ U_sum(ones(NAlt,1),:,:),1);
% P(isnan(P)) = 1; %
%  probs(:,n,:) = prod(P,2);
%end;
%probs = squeeze(probs); % taking this out of the loop saves almost 10% of time

if nargout == 1 % function value only
    
    if any(isnan(Xa(:))) == 0  % faster version for complete dataset
        
        parfor n = 1:NP
            b_mtx_n = b_mtx_sliced(:,:,n);
            Xa_n = Xa(:,:,n);
            Yy_n = Y(:,n)==1;
            U = reshape(Xa_n*b_mtx_n,NAlt,NCT, NRep);
            U_max = max(U);
            U = exp(U - U_max(ones(NAlt,1),:,:)); % rescale utility to avoid exploding
            U_sum = reshape(sum(U,1),NCT,NRep);
            %             U_selected = reshape(U(Y(:,n*ones(NRep,1))==1),NCT,NRep);
            U_selected = reshape(U(Yy_n(:,ones(NRep,1))),NCT,NRep);
            probs(n,:) = prod(U_selected ./ U_sum,1);
        end;
    else
        parfor n = 1:NP
            YnanInd = ~isnan(Y(:,n));
            %             U = reshape(Xa(~isnan(Y(:,n)),:,n)*b_mtx(:,((n-1)*NRep+1):n*NRep),NAlt,NCT-sum(isnan(Y(1:NAlt:end,n))),NRep);
            b_mtx_n = b_mtx_sliced(:,:,n);
            Xa_n = Xa(:,:,n);
            Yy_n = Y(:,n)==1;
            U = reshape(Xa_n(YnanInd,:)*b_mtx_n,NAltMiss(n),NCTMiss(n),NRep);
            U_max = max(U);
            %             U = exp(U - U_max(ones(NAlt,1),:,:)); % NAlt x NCT - NaNs x NRep
            U = exp(U - U_max(ones(NAltMiss(n),1),:,:)); % NAlt x NCT - NaNs x NRep
            U_sum = reshape(sum(U,1),NCTMiss(n),NRep);
            U_selected = reshape(U(Yy_n(YnanInd ,ones(NRep,1))),NCTMiss(n),NRep);
            probs(n,:) = prod(U_selected ./ U_sum,1);
        end;
    end
    
    L_mea = ones(NP,NRep);
    l = 0;
    
    if NVarMeaExp > 0
        Xmea_exp = reshape(Xmea_exp,1,NP,NVarMeaExp);
        Xmea_exp = reshape(Xmea_exp(ones(NRep,1),:,:), NP*NRep, NVarMeaExp);
    end
    for i = 1:size(Xmea,2)
        if MeaSpecMatrix(i) == 0 % OLS
            if MeaExpMatrix(i) == 0
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)'];
            else
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)', Xmea_exp];
            end
            b = bmea(l+1:l+size(X,2)+1);
            L_mea = L_mea.*normpdf(Xmea(:,i*ones(NRep,1)),reshape(X*b(1:end-1), NRep, NP)',exp(b(end)));
            l = l+size(X,2)+1;
        elseif MeaSpecMatrix(i) == 1 % MNL
            UniqueMea = unique(Xmea(:,i));
            k = length(UniqueMea)-1;
            if MeaExpMatrix(i) == 0
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)'];
            else
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)', Xmea_exp];
            end
            V = exp(X*reshape([zeros(size(X,2),1); bmea(l+1:l+size(X,2)*k)], size(X,2), k+1)); % NRep*NP x unique values of attitude
            Vsum = sum(V,2);
            V = reshape(V./Vsum(:, ones(k+1,1)), NRep, NP, k+1); % NRep x NP x unique
            V = permute(V, [2 1 3]); % NP x NRep x unique
            L = zeros(NP,NRep);
            for j = 1:length(UniqueMea)
                L(Xmea(:,i) == UniqueMea(j),:) = V(Xmea(:,i) == UniqueMea(j),:,j);
            end
            L_mea = L_mea.*L;
            l = l + size(X,2)*k;
            
        elseif MeaSpecMatrix(i) == 2 % Ordered Probit
            UniqueMea = unique(Xmea(:,i));
            k = length(UniqueMea)-1;
            if MeaExpMatrix(i) == 0
                X = LV(MeaMatrix(:,i)'== 1,:)';
            else
                X = [LV(MeaMatrix(:,i)'== 1,:)', Xmea_exp];
            end
            tmp = (MeaExpMatrix(i)~= 0)*NVarMeaExp;
            b = bmea(l+1:l+ k+ size(X,2));
            Xb = reshape(X*b(1:sum(MeaMatrix(:,i),1)+tmp), NRep, NP)'; % NP x NRep
            alpha = cumsum([b(sum(MeaMatrix(:,i))+tmp+1); exp(b(sum(MeaMatrix(:,i))+tmp+2:end))]);
            L = zeros(NP,NRep);
            L(Xmea(:,i) == min(UniqueMea),:) = normcdf(alpha(1)-Xb(Xmea(:,i) == min(UniqueMea),:));
            L(Xmea(:,i) == max(UniqueMea),:) = 1-normcdf(alpha(end)-Xb(Xmea(:,i) == max(UniqueMea),:));
            for j = 2:k
                L(Xmea(:,i) == UniqueMea(j),:) = normcdf(alpha(j)-Xb(Xmea(:,i) == UniqueMea(j),:)) - normcdf(alpha(j-1)-Xb(Xmea(:,i) == UniqueMea(j),:));
            end
            L_mea = L_mea.*L;
            l = l+k+size(X,2);
        elseif MeaSpecMatrix(i) == 3 % Poisson
            if MeaExpMatrix(i) == 0
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)'];
            else
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)', Xmea_exp];
            end
            b = bmea(l+1:l+size(X,2));
            fit = reshape(X*b, NRep, NP)';
            lam = exp(fit);
            %             L = exp(fit.*Xmea(:,i*ones(NRep,1))-lam)./min(gamma(Xmea(:,i*ones(NRep,1))+1),realmax);
            %             L = exp(fit.*Xmea(:,i*ones(NRep,1))-lam)./gamma(Xmea(:,i*ones(NRep,1))+1);
            if RealMin == 1
                L = exp(fit.*Xmea(:,i*ones(NRep,1))-lam)./min(gamma(Xmea(:,i*ones(NRep,1))+1),realmax);
            else
                L = exp(fit.*Xmea(:,i*ones(NRep,1))-lam-gammaln(Xmea(:,i*ones(NRep,1))+1));
            end
            L_mea = L_mea.*L;
            l = l+size(X,2);
        elseif MeaSpecMatrix(i) == 4 % NB
            if MeaExpMatrix(i) == 0
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)'];
            else
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)', Xmea_exp];
            end
            b = bmea(l+1:l+size(X,2));
            fit = reshape(X*b, NRep, NP)';
            lam = exp(fit);
            theta = exp(bmea(l+size(X,2)+1));
            u = theta./(theta+lam);
            %             L = min(gamma(theta+Xmea(:,i*ones(NRep,1))), realmax)./(gamma(theta).*min(gamma(Xmea(:,i*ones(NRep,1))+1),realmax));
            %             L = gamma(theta+Xmea(:,i*ones(NRep,1)))./(gamma(theta).*gamma(Xmea(:,i*ones(NRep,1))+1));
            if RealMin == 1
                L = min(gamma(theta+Xmea(:,i*ones(NRep,1))), realmax)./(gamma(theta).*min(gamma(Xmea(:,i*ones(NRep,1))+1),realmax));
            else
                L = exp(gammaln(theta+Xmea(:,i*ones(NRep,1))) -gammaln(theta) -gammaln(Xmea(:,i*ones(NRep,1))+1));
            end
            L = L.*(u.^theta).*((1-u).^Xmea(:,i*ones(NRep,1)));
            L_mea = L_mea.*L;
            l = l+size(X,2)+1;
        elseif MeaSpecMatrix(i) == 5 % ZIP
            if MeaExpMatrix(i) == 0
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)'];
            else
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)', Xmea_exp];
            end
            bzip = bmea(l+1:l+size(X,2));
            bpoiss = bmea(l+size(X,2)+1:l+2*size(X,2));
            fit = reshape(X*bpoiss, NRep, NP)';
            p = reshape(exp(X*bzip), NRep, NP)';
            p = p./(1+p);
            L = zeros(NP, NRep);
            lam = exp(fit);
            IndxZIP = Xmea(:,i) == 0;
            L(IndxZIP,:) = p(IndxZIP,:) + (1-p(IndxZIP,:)).*exp(-lam(IndxZIP,:));
            %             L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*exp(fit(~IndxZIP,:).*Xmea(~IndxZIP,i*ones(NRep,1))-lam(~IndxZIP,:))./min(gamma(Xmea(~IndxZIP,i*ones(NRep,1))+1),realmax);
            %             L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*exp(fit(~IndxZIP,:).*Xmea(~IndxZIP,i*ones(NRep,1))-lam(~IndxZIP,:))./ gamma(Xmea(~IndxZIP,i*ones(NRep,1))+1);
            if RealMin == 1
                L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*exp(fit(~IndxZIP,:).*Xmea(~IndxZIP,i*ones(NRep,1))-lam(~IndxZIP,:))./min(gamma(Xmea(~IndxZIP,i*ones(NRep,1))+1),realmax);
            else
                %  L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*exp(fit(~IndxZIP,:).*Xmea(~IndxZIP,i*ones(NRep,1))-lam(~IndxZIP,:))./ gamma(Xmea(~IndxZIP,i*ones(NRep,1))+1);
                L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*exp(fit(~IndxZIP,:).*Xmea(~IndxZIP,i*ones(NRep,1))-lam(~IndxZIP,:)-gammaln(Xmea(~IndxZIP,i*ones(NRep,1))+1));
            end
            L_mea = L_mea.*L;
            l = l+2*size(X,2);
        elseif MeaSpecMatrix(i) == 6 % ZINB
            if MeaExpMatrix(i) == 0
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)']; ...
            else
            X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)', Xmea_exp]; ...
            end
        bzip = bmea(l+1:l+size(X,2)); ...
            bpoiss = bmea(l+size(X,2)+1:l+2*size(X,2));
        fit = reshape(X*bpoiss, NRep, NP)';
        theta = exp(bmea(l+2*size(X,2)+1));
        p = reshape(exp(X*bzip), NRep, NP)';
        p = p./(1+p);
        L = zeros(NP, NRep);
        lam = exp(fit);
        u = theta./(theta+lam);
        IndxZIP = Xmea(:,i) == 0;
        L(IndxZIP,:) = p(IndxZIP,:) + (1-p(IndxZIP,:)).*(u(IndxZIP,:).^theta);
        if RealMin == 1
            L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*min(gamma(theta+Xmea(~IndxZIP,i*ones(NRep,1))), realmax)./(gamma(theta).*min(gamma(Xmea(~IndxZIP,i*ones(NRep,1))+1),realmax));
        else
            L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*exp(gammaln(theta+Xmea(~IndxZIP,i*ones(NRep,1))) - gammaln(Xmea(~IndxZIP,i*ones(NRep,1))+1)-gammaln(theta));
        end
        L(~IndxZIP,:) = L(~IndxZIP,:).*(u(~IndxZIP,:).^theta).*((1-u(~IndxZIP,:)).^Xmea(~IndxZIP,i*ones(NRep,1)));
        L_mea = L_mea.*L; ...
            l = l+2*size(X,2)+1; ...
        end
    end
    
    %     f = -log(max(realmin,mean(probs.*L_mea,2)));
    %     f = -log(mean(probs.*L_mea,2));
    if RealMin == 0
        f = -log(mean(probs.*L_mea,2));
    else
        f = -log(max(realmin,mean(probs.*L_mea,2)));
    end
    
else % function value + gradient
    
    %     if ~exist('b_mtx_grad','var')
    %         b_mtx_grad = [];
    %     end
    %     gmnl = zeros(NP, NRep, NVarA, (1+NLatent)); % gradient for mnl parameters
    gmnl = zeros(NP, NRep, NVarA); % gradient for mnl parameters
    gstr = zeros(NP, NRep, NVarStr,NLatent); % gradient for parameters from structural equations
    gmea = zeros(NP, NRep, size(bmea,1));% gradient for other parameters
    gs = zeros(NP, NRep, NVarS+(ScaleLV == 1)*NLatent);
    % terms from LV normalization
    mXstr = mean(Xstr,1);
    Xstr_expand = reshape(Xstr - mXstr(ones(NP,1),:), 1,NP, NVarStr);
    Xstr_expand = reshape(Xstr_expand(ones(NRep*NLatent,1), :,:), NLatent,NRep*NP, NVarStr);
    LV_tmp = LV_tmp - mLV(:,ones(1,size(LV_tmp,2))); % NLatent x NRep*NP
    LV_std = sum(LV_tmp(:,:,ones(NVarStr,1)).*Xstr_expand,2)/(NRep*NP-1); % NLatent x 1 x NVarstr
    
    Xstr_expand = Xstr_expand./sLV(:, ones(NRep*NP,1),ones(NVarStr,1));
    LV_tmp = LV_tmp./(sLV(:, ones(NRep*NP,1)).^3);
    LV_tmp = LV_tmp(:,:, ones(NVarStr,1));
    
    LV_der = reshape(Xstr_expand - LV_tmp.*LV_std(:,ones(NRep*NP,1),:), NLatent, NRep, NP, NVarStr); % NLatent x NRep x NP x NVarstr
    %LV_der = reshape(permute(LV_der, [3 2 4 1]), NP, NRep, NVarstr*NLatent);
    LV_der = permute(LV_der, [3 2 4 1]); % NP x NRep x NVarstr x NLatent
    LV_expand = permute(reshape(LV',NRep,NP,NLatent), [2 1 3]); % NP x NRep x NLatent
    if NVarS > 0
        Xs = reshape(Xs(1:NAlt:end,:), NCT, NP, NVarS);
        Xs = permute(Xs, [1 3 2]);
    else
        Xs = zeros(0,0, NP);
    end
    if ScaleLV > 0
        bsLV = bsLV';
    end
    if any(isnan(Xa(:))) == 0  % faster version for complete dataset
        
        parfor n = 1:NP
            b_mtx_n = b_mtx_sliced(:,:,n);
            Xa_n = Xa(:,:,n);
            Yy_n = Y(:,n)==1;
            %             U = reshape(Xa(:,:,n)*b_mtx(:,((n-1)*NRep+1):n*NRep),NAlt,NCT,NRep); % NAlt x NCT x NRep
            U = reshape(Xa_n*b_mtx_n,[NAlt,NCT,NRep]); % NAlt x NCT x NRep
            U_max = max(U);
            U = exp(U - U_max(ones(NAlt,1),:,:)); % rescale utility to avoid exploding
            U_sum = reshape(sum(U,1),1,NCT,NRep);
            U_prob = U./U_sum(ones(NAlt,1,1),:,:); % NAlt x NCT x NRep
            probs(n,:) = prod(reshape(U_prob(Yy_n(:,ones(NRep,1))),NCT,NRep),1); % 1 x NRep
            
            % calculations for gradient
            U_prob = reshape(U_prob, NAlt*NCT,1, NRep); % NAlt*NCT x NVarA x NRep
            if WTP_space == 0
                X_hat = sum(reshape(bsxfun(@times,U_prob,Xa_n),[NAlt,NCT,NVarA,NRep]),1);
                F = bsxfun(@minus,Xa_n(Y(:,n)==1,:,:),reshape(X_hat,[NCT,NVarA,NRep])); %NCT x NVarA x NRep
                if NVarS > 0 || ScaleLV == 1
                    bss = reshape(b_mtx_n,1,NVarA, NRep);
                    Fs = sum(bsxfun(@times, F, bss),2); % NCT x 1 x NRep
                    if ScaleLV == 1
                        LV_tmp = permute(LV_expand(n,:,:), [1 3 2]); % 1 x NLatent x NRep
                        FsLV = bsxfun(@times, Fs, LV_tmp);
                        FsX = sum(Fs,1); % 1 x 1 x NRep
                        FsLV = reshape(sum(FsLV,1), NLatent, NRep);
                    end
                    if NVarS > 0
                        Fs = bsxfun(@times, Fs, Xs(:,:,n)); % NCT x NVarS x NRep
                        Fs = reshape(sum(Fs,1), NVarS, NRep);
                    end
                end
                if ScaleLV == 1
                    F = bsxfun(@times, ScaleLVX(:,n,:), F);
                end
                sumFsqueezed = reshape(sum(F,1),[NVarA,NRep]); %NVarA x NRep
                if ScaleLV == 1 && any(Dist == 1)
                    b_mtx_grad_n = b_mtx_grad(:,:,n);
                    sumFsqueezed(Dist==1,:) = sumFsqueezed(Dist==1,:).*b_mtx_grad_n(Dist==1,:);
                else
                    sumFsqueezed(Dist==1,:) = sumFsqueezed(Dist==1,:).*b_mtx_n(Dist==1,:);
                end
                sumFsqueezed_LV = sumFsqueezed'*bl; % NRep x NLatent
            else
                b_mtx_wtp = reshape(b_mtx_n, 1, NVarA, NRep);
                Xalpha = Xa_n(:,1:end-WTP_space,ones(NRep,1)).*b_mtx_wtp(ones(NAlt*NCT,1),WTP_matrix,:);
                % for non-cost variables
                X_hat1 = sum(reshape(U_prob(:,ones(1,NVarA-WTP_space),:).*Xalpha, NAlt, NCT, NVarA-WTP_space, NRep),1);
                F1 = Xalpha(Y(:,n) == 1,:,:) - squeeze(X_hat1); %NCT x NVarA-WTP_space x NRep
                % for cost variables
                b_mtx_grad_n = b_mtx_grad(:,:,n);
                if WTP_space == 1 % without starting the loop
                    Xbmxlfit = Xa_n(:,1:end-WTP_space)*b_mtx_grad_n(1:end-WTP_space,:);
                    pX = squeeze(Xa_n(:,NVarA,ones(NRep,1))) + Xbmxlfit;
                    X_hat2 = sum(reshape(squeeze(U_prob).*pX, NAlt, NCT, WTP_space, NRep),1);
                else
                    pX = zeros(NCT*NAlt, WTP_space, NRep);
                    for i = 1:WTP_space
                        Xbmxlfit = Xa_n(:,WTP_matrix==NVarA-WTP_space+i)*b_mtx_grad_n(WTP_matrix==NVarA-WTP_space+i,:);
                        pX(:,i,:) = squeeze(Xa_n(:,NVarA-WTP_space+i,ones(NRep,1))) + Xbmxlfit;
                    end
                    X_hat2 = sum(reshape(U_prob(:,ones(1,WTP_space),:).* pX, NAlt, NCT, WTP_space, NRep),1);
                end
                F2 = pX(Y(:,n) == 1,:,:) - squeeze(X_hat2); %NCT x WTP_space x NRep
                
                if NVarS > 0 || ScaleLV == 1
                    bss = reshape(b_mtx_n(NVarA-WTP_space+1:end,:),1,WTP_space, NRep);
                    Fs = sum(bsxfun(@times, reshape(F2, NCT, WTP_space, NRep), bss),2); % NCT x 1 x NRep
                    if ScaleLV == 1
                        LV_tmp = permute(LV_expand(n,:,:), [1 3 2]); % 1 x NLatent x NRep
                        FsLV = bsxfun(@times, Fs, LV_tmp);
                        FsX = sum(Fs,1); % 1 x 1 x NRep
                        FsLV = reshape(sum(FsLV,1), NLatent, NRep);
                    end
                    if NVarS > 0
                        Fs = bsxfun(@times, Fs, Xs(:,:,n)); % NCT x NVarS x NRep
                        Fs = reshape(sum(Fs,1), NVarS, NRep);
                    end
                end
                if ScaleLV == 1
                    F2 = bsxfun(@times, ScaleLVX(:,n,:), reshape(F2, NCT, WTP_space, NRep));
                end
                sumFsqueezed = [reshape(sum(F1,1),NVarA - WTP_space, NRep);reshape(sum(F2,1),WTP_space, NRep) ]; %NVarA x NRep
                sumFsqueezed(Dist ==1, :) = sumFsqueezed(Dist ==1, :).*b_mtx_grad_n(Dist==1,:);
                sumFsqueezed_LV = sumFsqueezed'*bl; % NRep x NLatent
            end
            
            %             sumFsqueezed = sumFsqueezed'; % NRep x NVarA
            %             gmnl(n,:,:,:) = sumFsqueezed(:,:,ones(1+NLatent,1));
            gmnl(n,:,:) = sumFsqueezed';
            if NVarS > 0 && ScaleLV == 0
                gs(n,:,:) = Fs';
            elseif NVarS > 0 && ScaleLV == 1
                gs(n,:,:) = [Fs; FsLV]';
            elseif NVarS == 0 && ScaleLV == 1
                gs(n,:,:) = FsLV';
            end
            
            sumFsqueezed_LV = permute(sumFsqueezed_LV(:,:, ones(NVarStr,1)), [1 3 2]);
            gstr(n,:,:,:) = sumFsqueezed_LV;
            if ScaleLV == 1
                bsLVx = bsLV;
                FsLV_tmp = FsX(ones(NLatent,1),:,:).*bsLVx(:,:, ones(1, NRep));
                FsLV_tmp = permute(FsLV_tmp, [4 3 2 1]);
                gstr(n,:,:,:) = bsxfun(@plus, gstr(n,:,:,:), FsLV_tmp);
            end
            
        end;
    else
        if NVarS > 0 || ScaleLV == 1
            YCT = reshape(sum(reshape(~isnan(Y), NAlt, NCT, NP),1) ~= 0, NCT, NP);
        else
            YCT = zeros(0,NP);
        end
        parfor n = 1:NP
            YnanInd = ~isnan(Y(:,n));
            b_mtx_n = b_mtx_sliced(:,:,n);
            Xa_n = Xa(:,:,n);
            Yy_n = Y(:,n)==1;
            U = reshape(Xa_n(YnanInd,:)*b_mtx_n,NAltMiss(n),NCTMiss(n),NRep); % NAlt x NCT x NRep
            U_max = max(U);
            U = exp(U - U_max(ones(NAltMiss(n),1),:,:)); % rescale utility to avoid exploding
            U_sum = reshape(sum(U,1),1,NCTMiss(n),NRep);
            U_prob = U./U_sum(ones(NAltMiss(n),1,1),:,:); % NAlt x NCT x NRep
            probs(n,:) = prod(reshape(U_prob(Yy_n(YnanInd,ones(NRep,1))),NCTMiss(n),NRep),1); % 1 x NRep
            %probs(n,:) = U_prod;
            
            % calculations for gradient
            U_prob = reshape(U_prob, NAltMiss(n)*NCTMiss(n),1, NRep); % NAlt*NCT x NVarA x NRep
            if WTP_space == 0
                X_hat = sum(reshape(U_prob(:,ones(1,NVarA,1),:).* Xa_n(YnanInd,:,ones(NRep,1)), NAltMiss(n), NCTMiss(n), NVarA, NRep),1);
                F = bsxfun(@minus,Xa_n(Yy_n,:,:),reshape(X_hat,[NCTMiss(n),NVarA,NRep])); %NCT x NVarA x NRep
                if NVarS > 0 || ScaleLV == 1
                    bss = reshape(b_mtx_n,1,NVarA, NRep);
                    Fs = sum(bsxfun(@times, F, bss),2); % NCT x 1 x NRep
                    if ScaleLV == 1
                        LV_tmp = permute(LV_expand(n,:,:), [1 3 2]); % 1 x NLatent x NRep
                        FsLV = bsxfun(@times, Fs, LV_tmp);
                        FsX = sum(Fs,1); % 1 x 1 x NRep
                        FsLV = reshape(sum(FsLV,1), NLatent, NRep);
                    end
                    if NVarS > 0
                        Xs_n = Xs(:,:,n);
                        Fs = bsxfun(@times, Fs, Xs_n(YCT(:,n),:)); % NCT x NVarS x NRep
                        Fs = reshape(sum(Fs,1), NVarS, NRep);
                    end
                end
                if ScaleLV == 1
                    F = bsxfun(@times, ScaleLVX(:,n,:), F);
                end
                sumFsqueezed = reshape(sum(F,1),[NVarA,NRep]); %NVarA x NRep
                if ScaleLV == 1 && any(Dist == 1)
                    b_mtx_grad_n = b_mtx_grad(:,:,n);
                    sumFsqueezed(Dist==1,:) = sumFsqueezed(Dist==1,:).*b_mtx_grad_n(Dist==1,:);
                else
                    sumFsqueezed(Dist==1,:) = sumFsqueezed(Dist==1,:).*b_mtx_n(Dist==1,:);
                end
                %sumFsqueezed(Dist ==1, :) = sumFsqueezed(Dist ==1, :).*b_mtx_n(Dist==1,:);
                sumFsqueezed_LV = sumFsqueezed'*bl; % NRep x NLatent
            else
                b_mtx_wtp = reshape(b_mtx_n, 1, NVarA, NRep);
                Xalpha = Xa_n(YnanInd,1:end-WTP_space,ones(NRep,1)).*b_mtx_wtp(ones(NAltMiss(n)*NCTMiss(n),1),WTP_matrix,:);
                % for non-cost variables
                X_hat1 = sum(reshape(U_prob(:,ones(1,NVarA-WTP_space),:).* Xalpha, NAltMiss(n), NCTMiss(n), NVarA-WTP_space, NRep),1);
                F1 = Xalpha(Yy_n(YnanInd),:,:) - reshape(X_hat1, NCTMiss(n), NVarA-WTP_space, NRep); %NCT x NVarA-WTP_space x NRep
                % for cost variables
                b_mtx_grad_n = b_mtx_grad(:,:,n);
                if WTP_space == 1 % without starting the loop
                    Xbmxlfit = Xa_n(YnanInd,1:end-WTP_space)*b_mtx_grad_n(1:end-WTP_space,:);
                    pX = squeeze(Xa_n(YnanInd,NVarA,ones(NRep,1))) + Xbmxlfit;
                    X_hat2 = sum(reshape(squeeze(U_prob).*pX, NAltMiss(n), NCTMiss(n), WTP_space, NRep),1);
                    F2 = pX(Yy_n(YnanInd),:,:) - reshape(X_hat2, NCTMiss(n), NRep); %NCT x NRep
                else
                    pX = zeros(NCTMiss(n)*NAltMiss(n), WTP_space, NRep);
                    for i = 1:WTP_space
                        Xbmxlfit = Xa_n(YnanInd,WTP_matrix==NVarA-WTP_space+i)*b_mtx_grad_n(WTP_matrix==NVarA-WTP_space+i,:);
                        pX(:,i,:) = squeeze(Xa_n(YnanInd,NVarA-WTP_space+i,ones(NRep,1))) + Xbmxlfit;
                    end
                    X_hat2 = sum(reshape(U_prob(:,ones(1,WTP_space),:).*pX, NAltMiss(n), NCTMiss(n), WTP_space, NRep),1);
                    F2 = pX(Yy_n(YnanInd),:,:) - reshape(X_hat2, NCTMiss(n),WTP_space, NRep); %NCT x WTP x NRep
                end
                if NVarS > 0 || ScaleLV == 1
                    bss = reshape(b_mtx_n(NVarA-WTP_space+1:end,:),1,WTP_space, NRep);
                    Fs = sum(bsxfun(@times, reshape(F2, NCTMiss(n), WTP_space, NRep), bss),2); % NCT x 1 x NRep
                    if ScaleLV == 1
                        LV_tmp = permute(LV_expand(n,:,:), [1 3 2]); % 1 x NLatent x NRep
                        FsLV = bsxfun(@times, Fs, LV_tmp);
                        FsX = sum(Fs,1); % 1 x 1 x NRep
                        FsLV = reshape(sum(FsLV,1), NLatent, NRep);
                    end
                    if NVarS > 0
                        Xs_n = Xs(:,:,n);
                        Fs = bsxfun(@times, Fs, Xs_n(YCT(:,n),:)); % NCT x NVarS x NRep
                        Fs = reshape(sum(Fs,1), NVarS, NRep);
                    end
                end
                if ScaleLV == 1
                    F2 = bsxfun(@times, ScaleLVX(:,n,:), reshape(F2, NCTMiss(n), WTP_space, NRep));
                end
                sumFsqueezed = [reshape(sum(F1,1),NVarA - WTP_space, NRep);reshape(sum(F2,1),WTP_space, NRep) ]; %NVarA x NRep
                sumFsqueezed(Dist ==1, :) = sumFsqueezed(Dist ==1, :).*b_mtx_grad_n(Dist==1,:);
                sumFsqueezed_LV = sumFsqueezed'*bl; % NRep x NLatent
            end
            %             sumFsqueezed = sumFsqueezed'; % NRep x NVarA
            %             gmnl(n,:,:,:) = sumFsqueezed(:,:,ones(1+NLatent,1));
            gmnl(n,:,:) = sumFsqueezed';
            if NVarS > 0 && ScaleLV == 0
                gs(n,:,:) = Fs';
            elseif NVarS > 0 && ScaleLV == 1
                gs(n,:,:) = [Fs; FsLV]';
            elseif NVarS == 0 && ScaleLV == 1
                gs(n,:,:) = FsLV';
            end
            sumFsqueezed_LV = permute(sumFsqueezed_LV(:,:, ones(NVarStr,1)), [1 3 2]);
            gstr(n,:,:,:) = sumFsqueezed_LV;
            if ScaleLV == 1
                bsLVx = bsLV;
                FsLV_tmp = FsX(ones(NLatent,1),:,:).*bsLVx(:,:, ones(1, NRep));
                FsLV_tmp = permute(FsLV_tmp, [4 3 2 1]);
                gstr(n,:,:,:) = bsxfun(@plus, gstr(n,:,:,:), FsLV_tmp);
            end
        end;
        
    end
    
    gstr = gstr.*LV_der;
    %     LV_expand = reshape(LV_expand,[NP,NRep,1,NLatent]);
    %     gmnl0(:,:,:,2:end) = gmnl0(:,:,:,2:end).*LV_expand(:,:, ones(NVarA,1),:);
    %     gmnl = cat(4,gmnl,bsxfun(@times,gmnl,LV_expand));
    LV_expand_tmp = cat(3,ones(NP,NRep,1),LV_expand);
    LV_expand_tmp = reshape(LV_expand_tmp,[NP,NRep,1,1+NLatent]);
    gmnl = bsxfun(@times,gmnl,LV_expand_tmp);
    %     LV_expand = squeeze(LV_expand);
    L_mea = ones(NP,NRep);
    l = 0;
    
    if NVarMeaExp > 0
        Xmea_exp = reshape(Xmea_exp,1,NP,NVarMeaExp);
        Xmea_exp = reshape(Xmea_exp(ones(NRep,1),:,:), NP*NRep, NVarMeaExp);
    end
    if ScaleLV == 1
        NVarS = NVarS + NLatent;
    end
    Xmea_exp_expand = permute(reshape(Xmea_exp, NRep, NP,NVarMeaExp), [2 1 3]);
    
    for i = 1:size(Xmea,2)
        if MeaSpecMatrix(i) == 0 % OLS
            if MeaExpMatrix(i) == 0
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)'];
            else
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)', Xmea_exp];
            end
            b = bmea(l+1:l+size(X,2)+1);
            fit = reshape(X*b(1:end-1), NRep, NP)';
            L_mea = L_mea.*normpdf(Xmea(:,i*ones(NRep,1)),fit,exp(b(end)));
            grad_tmp = - (fit - Xmea(:,i*ones(NRep,1)))/exp(2*b(end)); % NP x NRep
            LVindx = find(MeaMatrix(:,i)'== 1);
            if sum(MeaMatrix(:,i)'== 1) > 1
                bx = b(2:1+sum(MeaMatrix(:,i)'== 1));
                bx = permute(bx(:, ones(NP, 1), ones(NRep, 1), ones(NVarStr,1)), [2 3 4 1]);
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp(:,:, ones(NVarStr,1), ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_der(:,:,:,LVindx).*bx;
            else
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp(:,:, ones(NVarStr,1)).*LV_der(:,:,:,LVindx)*b(2);
            end
            
            gmea(:,:,l+1) = grad_tmp; % constant
            %LV_expand = permute(reshape(LV(MeaMatrix(:,i)'== 1,:)',NRep,NP,sum(MeaMatrix(:,i)'== 1)), [2 1 3])  ;
            gmea(:,:,l+2:l+1+sum(MeaMatrix(:,i)'== 1)) = grad_tmp(:,:, ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_expand(:,:,MeaMatrix(:,i)'== 1); % parameters for LV
            if MeaExpMatrix(i) == 0
                gmea(:,:,l+2+sum(MeaMatrix(:,i)'== 1)) = -1+((fit - Xmea(:,i*ones(NRep,1))).^2)/(exp(2*b(end))) ; % variance
            else
                gmea(:,:,l+2+sum(MeaMatrix(:,i)'== 1):l+1+sum(MeaMatrix(:,i)'== 1)+NVarMeaExp) = ...
                    grad_tmp(:,:,ones(NVarMeaExp,1)).*Xmea_exp_expand ;
                gmea(:,:,l+2+sum(MeaMatrix(:,i)'== 1)+NVarMeaExp) = -1+((fit - Xmea(:,i*ones(NRep,1))).^2)/(exp(2*b(end))) ; % variance
                
            end
            l = l+size(X,2)+1;
        elseif MeaSpecMatrix(i) == 1 % MNL
            UniqueMea = unique(Xmea(:,i));
            k = length(UniqueMea)-1;
            if MeaExpMatrix(i) == 0
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)'];
            else
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)', Xmea_exp];
            end
            bx = reshape([zeros(size(X,2),1); bmea(l+1:l+size(X,2)*k)], size(X,2), k+1);
            V = exp(X*bx); % NRep*NP x unique values of attitude
            Vsum = sum(V,2);
            V = reshape(V./Vsum(:, ones(k+1,1)), NRep, NP, k+1); % NRep x NP x unique
            V = permute(V, [2 1 3]); % NP x NRep x unique
            L = zeros(NP,NRep);
            L(Xmea(:,i) == UniqueMea(1),:) = V(Xmea(:,i) == UniqueMea(1),:,1);
            MAT = permute(reshape(X, NRep, NP, size(X,2)), [2 1 3]);
            for j = 2:length(UniqueMea)
                L(Xmea(:,i) == UniqueMea(j),:) = V(Xmea(:,i) == UniqueMea(j),:,j);
                
                gmea(Xmea(:,i) == UniqueMea(j),:,l+1+(j-2)*size(X,2):l+size(X,2)*(j-1)) = (1 - V(Xmea(:,i) == UniqueMea(j),:,j*ones(1, size(X,2)))).*MAT(Xmea(:,i) == UniqueMea(j),:,:);
                gmea(Xmea(:,i) ~= UniqueMea(j),:,l+1+(j-2)*size(X,2):l+size(X,2)*(j-1)) = (  - V(Xmea(:,i) ~= UniqueMea(j),:,j*ones(1, size(X,2)))).*MAT(Xmea(:,i) ~= UniqueMea(j),:,:);
            end
            alpha = bx(2:1+sum(MeaMatrix(:,i),1),:);
            alpha = alpha(:,:, ones(NP,1), ones(NRep,1));
            alphax = reshape(permute(alpha, [2 3 4 1]),(k+1)*NP,NRep, sum(MeaMatrix(:,i)));
            alpha = permute(alpha, [3 4 2 1]);
            LVindx = find(MeaMatrix(:,i)'== 1);
            Xmea_i = reshape(dummyvar(Xmea(:,i))', (k+1)*NP,1);
            if sum(MeaMatrix(:,i)'== 1) > 1
                
                for j = 1:sum(MeaMatrix(:,i)'== 1)
                    grad_tmp = alphax(Xmea_i==1,:,j)  -sum(V.*alpha(:,:,:,j),3); % NP x NRep
                    gstr(:,:,:,LVindx(j)) = gstr(:,:,:,LVindx(j)) +grad_tmp(:,:, ones(NVarStr,1)).*LV_der(:,:,:,LVindx(j));
                end
            else
                grad_tmp = alphax(Xmea_i==1,:)  -sum(V.*alpha,3); % NP x NRep
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx) +grad_tmp(:,:, ones(NVarStr,1)).*LV_der(:,:,:,LVindx);
            end
            L_mea = L_mea.*L;
            l = l + size(X,2)*k;
            
        elseif MeaSpecMatrix(i) == 2 % Ordered Probit
            UniqueMea = unique(Xmea(:,i));
            k = length(UniqueMea)-1;
            if MeaExpMatrix(i) == 0
                X = LV(MeaMatrix(:,i)'== 1,:)';
            else
                X = [LV(MeaMatrix(:,i)'== 1,:)', Xmea_exp];
            end
            tmp = (MeaExpMatrix(i)~= 0)*NVarMeaExp;
            b = bmea(l+1:l+ k+ size(X,2));
            Xb = reshape(X*b(1:sum(MeaMatrix(:,i),1)+tmp), NRep, NP)'; % NP x NRep
            alpha = cumsum([b(sum(MeaMatrix(:,i))+tmp+1); exp(b(sum(MeaMatrix(:,i))+tmp+2:end))]);
            L = zeros(NP,NRep);
            grad_tmp = zeros(NP,NRep);
            L(Xmea(:,i) == min(UniqueMea),:) = normcdf(alpha(1)-Xb(Xmea(:,i) == min(UniqueMea),:));
            L(Xmea(:,i) == max(UniqueMea),:) = 1-normcdf(alpha(end)-Xb(Xmea(:,i) == max(UniqueMea),:));
            %              grad_tmp(Xmea(:,i) == min(UniqueMea),:) = normpdf(alpha(1)-Xb(Xmea(:,i) == min(UniqueMea),:))./max(L(Xmea(:,i) == min(UniqueMea),:),realmin);
            %             grad_tmp(Xmea(:,i) == min(UniqueMea),:) = normpdf(alpha(1)-Xb(Xmea(:,i) == min(UniqueMea),:))./L(Xmea(:,i) == min(UniqueMea),:);
            %              grad_tmp(Xmea(:,i) == max(UniqueMea),:) = -normpdf(alpha(end)-Xb(Xmea(:,i) == max(UniqueMea),:))./max(L(Xmea(:,i) == max(UniqueMea),:), realmin);
            %             grad_tmp(Xmea(:,i) == max(UniqueMea),:) = -normpdf(alpha(end)-Xb(Xmea(:,i) == max(UniqueMea),:))./L(Xmea(:,i) == max(UniqueMea),:);
            if RealMin == 1
                grad_tmp(Xmea(:,i) == min(UniqueMea),:) = normpdf(alpha(1)-Xb(Xmea(:,i) == min(UniqueMea),:))./max(L(Xmea(:,i) == min(UniqueMea),:),realmin);
                grad_tmp(Xmea(:,i) == max(UniqueMea),:) = -normpdf(alpha(end)-Xb(Xmea(:,i) == max(UniqueMea),:))./max(L(Xmea(:,i) == max(UniqueMea),:), realmin);
            else
                grad_tmp(Xmea(:,i) == min(UniqueMea),:) = normpdf(alpha(1)-Xb(Xmea(:,i) == min(UniqueMea),:))./L(Xmea(:,i) == min(UniqueMea),:);
                grad_tmp(Xmea(:,i) == max(UniqueMea),:) = -normpdf(alpha(end)-Xb(Xmea(:,i) == max(UniqueMea),:))./L(Xmea(:,i) == max(UniqueMea),:);
            end
            
            for j = 2:k
                L(Xmea(:,i) == UniqueMea(j),:) = normcdf(alpha(j)-Xb(Xmea(:,i) == UniqueMea(j),:)) - normcdf(alpha(j-1)-Xb(Xmea(:,i) == UniqueMea(j),:));
                grad_tmp(Xmea(:,i) == UniqueMea(j),:) = normpdf(alpha(j)-Xb(Xmea(:,i) == UniqueMea(j),:)) - normpdf(alpha(j-1)-Xb(Xmea(:,i) == UniqueMea(j),:));
                %                 grad_tmp(Xmea(:,i) == UniqueMea(j),:) = grad_tmp(Xmea(:,i) == UniqueMea(j),:)./max(L(Xmea(:,i) == UniqueMea(j),:),realmin);
                %                 grad_tmp(Xmea(:,i) == UniqueMea(j),:) = grad_tmp(Xmea(:,i) == UniqueMea(j),:)./L(Xmea(:,i) == UniqueMea(j),:);
                if RealMin == 1
                    grad_tmp(Xmea(:,i) == UniqueMea(j),:) = grad_tmp(Xmea(:,i) == UniqueMea(j),:)./max(L(Xmea(:,i) == UniqueMea(j),:),realmin);
                else
                    grad_tmp(Xmea(:,i) == UniqueMea(j),:) = grad_tmp(Xmea(:,i) == UniqueMea(j),:)./L(Xmea(:,i) == UniqueMea(j),:);
                end
            end
            
            LVindx = find(MeaMatrix(:,i)'== 1);
            if sum(MeaMatrix(:,i)'== 1) > 1
                bx = b(1:sum(MeaMatrix(:,i)'== 1));
                bx = permute(bx(:, ones(NP, 1), ones(NRep, 1), ones(NVarStr,1)), [2 3 4 1]);
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)- grad_tmp(:,:, ones(NVarStr,1), ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_der(:,:,:,LVindx).*bx;
            else
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)- grad_tmp(:,:, ones(NVarStr,1)).*LV_der(:,:,:,LVindx)*b(1);
            end
            
            gmea(:,:,l+1:l+sum(MeaMatrix(:,i)'== 1)) = -grad_tmp(:,:, ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_expand(:,:,MeaMatrix(:,i)'== 1); % this is for LV parameters
            if MeaExpMatrix(i) ~= 0
                gmea(:,:,l+1+sum(MeaMatrix(:,i)'== 1):l+sum(MeaMatrix(:,i)'== 1)+tmp) = ...
                    -grad_tmp(:,:,ones(NVarMeaExp,1)).*Xmea_exp_expand ; % additional explanatory variables of Xmea
            end
            gmea(:,:,l+1+sum(MeaMatrix(:,i)'== 1)+tmp) = grad_tmp; %first threshold level
            for j = 2:k
                gmea(Xmea(:,i) > UniqueMea(j),:,l+j+sum(MeaMatrix(:,i)'== 1)+tmp) = grad_tmp(Xmea(:,i) > UniqueMea(j),:)*exp(b(sum(MeaMatrix(:,i))+tmp+j)); %other thresholds levels
                %                 gmea(Xmea(:,i) == UniqueMea(j),:,l+j+sum(MeaMatrix(:,i)'== 1)+tmp) = (normpdf(alpha(j)-Xb(Xmea(:,i) == UniqueMea(j),:))./max(L(Xmea(:,i) == UniqueMea(j),:),realmin))*exp(b(sum(MeaMatrix(:,i))+tmp+j)); %other thresholds levels
                %                 gmea(Xmea(:,i) == UniqueMea(j),:,l+j+sum(MeaMatrix(:,i)'== 1)+tmp) = (normpdf(alpha(j)-Xb(Xmea(:,i) == UniqueMea(j),:))./L(Xmea(:,i) == UniqueMea(j),:))*exp(b(sum(MeaMatrix(:,i))+tmp+j)); %other thresholds levels
                if RealMin == 1
                    gmea(Xmea(:,i) == UniqueMea(j),:,l+j+sum(MeaMatrix(:,i)'== 1)+tmp) = (normpdf(alpha(j)-Xb(Xmea(:,i) == UniqueMea(j),:))./max(L(Xmea(:,i) == UniqueMea(j),:),realmin))*exp(b(sum(MeaMatrix(:,i))+tmp+j)); %other thresholds levels
                else
                    gmea(Xmea(:,i) == UniqueMea(j),:,l+j+sum(MeaMatrix(:,i)'== 1)+tmp) = (normpdf(alpha(j)-Xb(Xmea(:,i) == UniqueMea(j),:))./L(Xmea(:,i) == UniqueMea(j),:))*exp(b(sum(MeaMatrix(:,i))+tmp+j)); %other thresholds levels
                end
            end
            L_mea = L_mea.*L;
            l = l+k+size(X,2);
        elseif MeaSpecMatrix(i) == 3 % POISS
            if MeaExpMatrix(i) == 0
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)'];
            else
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)', Xmea_exp];
            end
            b = bmea(l+1:l+size(X,2));
            fit = reshape(X*b, NRep, NP)';
            lam = exp(fit);
            
            %             L = L./min(gamma(Xmea(:,i*ones(NRep,1))+1),realmax);
            %             L = L./gamma(Xmea(:,i*ones(NRep,1))+1);
            if RealMin == 1
                L = exp(fit.*Xmea(:,i*ones(NRep,1))-lam);
                L = L./min(gamma(Xmea(:,i*ones(NRep,1))+1),realmax);
            else
                L =exp(fit.*Xmea(:,i*ones(NRep,1))-lam-gammaln(Xmea(:,i*ones(NRep,1))+1));
            end
            L_mea = L_mea.*L;
            grad_tmp = Xmea(:,i*ones(NRep,1)) - lam;
            LVindx = find(MeaMatrix(:,i)'== 1);
            
            % gradient for structural equation
            if sum(MeaMatrix(:,i)'== 1) > 1
                bx = b(2:1+sum(MeaMatrix(:,i)'== 1)); % parameters for LV
                bx = permute(bx(:, ones(NP, 1), ones(NRep, 1), ones(NVarStr,1)), [2 3 4 1]);
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp(:,:, ones(NVarStr,1), ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_der(:,:,:,LVindx).*bx;
            else
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp(:,:, ones(NVarStr,1)).*LV_der(:,:,:,LVindx)*b(2);
            end
            gmea(:,:,l+1) = grad_tmp; % constant
            gmea(:,:,l+2:l+1+sum(MeaMatrix(:,i)'== 1)) = grad_tmp(:,:, ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_expand(:,:,MeaMatrix(:,i)'== 1); % parameters for LV
            if MeaExpMatrix(i) ~= 0
                gmea(:,:,l+2+sum(MeaMatrix(:,i)'== 1):l+1+sum(MeaMatrix(:,i)'== 1)+NVarMeaExp) = grad_tmp(:,:,ones(NVarMeaExp,1)).*Xmea_exp_expand ;
            end
            l = l+size(X,2);
            
        elseif MeaSpecMatrix(i) == 4 % NB
            if MeaExpMatrix(i) == 0
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)'];
            else
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)', Xmea_exp];
            end
            b = bmea(l+1:l+size(X,2));
            fit = reshape(X*b, NRep, NP)';
            lam = exp(fit);
            theta = exp(bmea(l+size(X,2)+1));
            u = theta./(theta+lam);
            %             L = min(gamma(theta+Xmea(:,i*ones(NRep,1))), realmax)./(gamma(theta).*min(gamma(Xmea(:,i*ones(NRep,1))+1),realmax));
            %             L = gamma(theta+Xmea(:,i*ones(NRep,1))) ./ (gamma(theta) .* gamma(Xmea(:,i*ones(NRep,1))+1));
            if RealMin == 1
                L = min(gamma(theta+Xmea(:,i*ones(NRep,1))), realmax)./(gamma(theta).*min(gamma(Xmea(:,i*ones(NRep,1))+1),realmax));
            else
                L = exp(gammaln(theta+Xmea(:,i*ones(NRep,1)))-gammaln(theta)-gammaln(Xmea(:,i*ones(NRep,1))+1));
            end
            L = L.*(u.^theta).*((1-u).^Xmea(:,i*ones(NRep,1)));
            L_mea = L_mea.*L;
            
            % Calculations for gradient
            grad_tmp = u.*(Xmea(:,i*ones(NRep,1))-  lam);
            LVindx = find(MeaMatrix(:,i)'== 1);
            % gradient for structural equation
            if sum(MeaMatrix(:,i)'== 1) > 1
                bx = b(2:1+sum(MeaMatrix(:,i)'== 1)); % parameters for LV
                bx = permute(bx(:, ones(NP, 1), ones(NRep, 1), ones(NVarStr,1)), [2 3 4 1]);
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp(:,:, ones(NVarStr,1), ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_der(:,:,:,LVindx).*bx;
            else
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp(:,:, ones(NVarStr,1)).*LV_der(:,:,:,LVindx)*b(2);
            end
            
            gmea(:,:,l+1) = grad_tmp; % constant
            gmea(:,:,l+2:l+1+sum(MeaMatrix(:,i)'== 1)) = grad_tmp(:,:, ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_expand(:,:,MeaMatrix(:,i)'== 1); % parameters for LV
            
            if MeaExpMatrix(i) == 0
                gmea(:,:,l+2+sum(MeaMatrix(:,i)'== 1)) = psi(theta+Xmea(:,i*ones(NRep,1)))*theta - psi(theta)*theta+theta*log(u)...
                    + theta*(1-u) - Xmea(:,i*ones(NRep,1)).*u; % theta
            else
                gmea(:,:,l+2+sum(MeaMatrix(:,i)'== 1):l+1+sum(MeaMatrix(:,i)'== 1)+NVarMeaExp) = ...
                    grad_tmp(:,:,ones(NVarMeaExp,1)).*Xmea_exp_expand ;
                gmea(:,:,l+2+sum(MeaMatrix(:,i)'== 1)+NVarMeaExp) = psi(theta+Xmea(:,i*ones(NRep,1)))*theta - psi(theta)*theta+theta*log(u)...
                    + theta*(1-u) - Xmea(:,i*ones(NRep,1)).*u; % theta
            end
            l = l+size(X,2)+1;
        elseif MeaSpecMatrix(i) == 5 % ZIP
            if MeaExpMatrix(i) == 0
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)'];
            else
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)', Xmea_exp];
            end
            bzip = bmea(l+1:l+size(X,2));
            bpoiss = bmea(l+size(X,2)+1:l+2*size(X,2));
            fit = reshape(X*bpoiss, NRep, NP)';
            p = reshape(exp(X*bzip), NRep, NP)';
            p = p./(1+p);
            L = zeros(NP, NRep);
            lam = exp(fit);
            IndxZIP = Xmea(:,i) == 0;
            L(IndxZIP,:) = p(IndxZIP,:) + (1-p(IndxZIP,:)).*exp(-lam(IndxZIP,:));
            %             L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*exp(fit(~IndxZIP,:).*Xmea(~IndxZIP,i*ones(NRep,1))-lam(~IndxZIP,:))./min(gamma(Xmea(~IndxZIP,i*ones(NRep,1))+1),realmax);
            %             L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*exp(fit(~IndxZIP,:).*Xmea(~IndxZIP,i*ones(NRep,1))-lam(~IndxZIP,:))./ gamma(Xmea(~IndxZIP,i*ones(NRep,1))+1);
            if RealMin == 1
                L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*exp(fit(~IndxZIP,:).*Xmea(~IndxZIP,i*ones(NRep,1))-lam(~IndxZIP,:))./min(gamma(Xmea(~IndxZIP,i*ones(NRep,1))+1),realmax);
            else
                L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*exp(fit(~IndxZIP,:).*Xmea(~IndxZIP,i*ones(NRep,1))-lam(~IndxZIP,:)-gammaln(Xmea(~IndxZIP,i*ones(NRep,1))+1));
            end
            L_mea = L_mea.*L;
            
            % Calculations for gradient
            grad_tmp1 = zeros(NP, NRep);
            grad_tmp2 = zeros(NP, NRep);
            %             % For ZIP
            %             grad_tmp1(IndxZIP,:) = (-(1-p(IndxZIP,:)).*p(IndxZIP,:).*exp(-lam(IndxZIP,:))+p(IndxZIP,:)-p(IndxZIP,:).^2)./L(IndxZIP,:);
            %             grad_tmp1(~IndxZIP,:) = -p(~IndxZIP,:);
            %             % For Poiss
            %             grad_tmp2(IndxZIP,:) = (p(IndxZIP,:)-1).*exp(fit(IndxZIP,:)-lam(IndxZIP,:))./L(IndxZIP,:);
            %             grad_tmp2(~IndxZIP,:) = Xmea(~IndxZIP,i*ones(NRep,1)) - lam(~IndxZIP,:);
            if RealMin == 1
                grad_tmp1(IndxZIP,:) = (-(1-p(IndxZIP,:)).*p(IndxZIP,:).*exp(-lam(IndxZIP,:))+p(IndxZIP,:)-p(IndxZIP,:).^2)./max(L(IndxZIP,:),realmin);
                grad_tmp2(IndxZIP,:) = (p(IndxZIP,:)-1).*exp(fit(IndxZIP,:)-lam(IndxZIP,:))./max(L(IndxZIP,:),realmin);
            else
                grad_tmp1(IndxZIP,:) = (-(1-p(IndxZIP,:)).*p(IndxZIP,:).*exp(-lam(IndxZIP,:))+p(IndxZIP,:)-p(IndxZIP,:).^2)./L(IndxZIP,:);
                grad_tmp2(IndxZIP,:) = (p(IndxZIP,:)-1).*exp(fit(IndxZIP,:)-lam(IndxZIP,:))./L(IndxZIP,:);
            end
            grad_tmp1(~IndxZIP,:) = -p(~IndxZIP,:);
            grad_tmp2(~IndxZIP,:) = Xmea(~IndxZIP,i*ones(NRep,1)) - lam(~IndxZIP,:);
            
            LVindx = find(MeaMatrix(:,i)'== 1);
            % gradient for structural equation
            if sum(MeaMatrix(:,i)'== 1) > 1
                bx1 = bzip(2:1+sum(MeaMatrix(:,i)'== 1)); % parameters for LV
                bx1 = permute(bx1(:, ones(NP, 1), ones(NRep, 1), ones(NVarStr,1)), [2 3 4 1]);
                bx2 = bpoiss(2:1+sum(MeaMatrix(:,i)'== 1)); % parameters for LV
                bx2 = permute(bx2(:, ones(NP, 1), ones(NRep, 1), ones(NVarStr,1)), [2 3 4 1]);
                
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp1(:,:, ones(NVarStr,1), ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_der(:,:,:,LVindx).*bx1+...
                    grad_tmp2(:,:, ones(NVarStr,1), ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_der(:,:,:,LVindx).*bx2;
            else
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp1(:,:, ones(NVarStr,1)).*LV_der(:,:,:,LVindx)*bzip(2)+...
                    grad_tmp2(:,:, ones(NVarStr,1)).*LV_der(:,:,:,LVindx)*bpoiss(2);
            end
            gmea(:,:,l+1) = grad_tmp1; % constant
            gmea(:,:,l+2:l+1+sum(MeaMatrix(:,i)'== 1)) = grad_tmp1(:,:, ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_expand(:,:,MeaMatrix(:,i)'== 1); % parameters for LV
            gmea(:,:,l+sum(MeaMatrix(:,i)'== 1)+MeaExpMatrix(i)*NVarMeaExp+2) = grad_tmp2; % constant
            gmea(:,:,l+sum(MeaMatrix(:,i)'== 1)+MeaExpMatrix(i)*NVarMeaExp+3:l+2*sum(MeaMatrix(:,i)'== 1)+MeaExpMatrix(i)*NVarMeaExp+2) = ...
                grad_tmp2(:,:, ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_expand(:,:,MeaMatrix(:,i)'== 1); % parameters for LV
            
            if MeaExpMatrix(i) ~= 0
                gmea(:,:,l+2+sum(MeaMatrix(:,i)'== 1):l+1+sum(MeaMatrix(:,i)'== 1)+NVarMeaExp) = ...
                    grad_tmp1(:,:,ones(NVarMeaExp,1)).*Xmea_exp_expand ;
                gmea(:,:,l+2*sum(MeaMatrix(:,i)'== 1)+MeaExpMatrix(i)*NVarMeaExp+3:l+2*sum(MeaMatrix(:,i)'== 1)+2*MeaExpMatrix(i)*NVarMeaExp+2) = ...
                    grad_tmp2(:,:,ones(NVarMeaExp,1)).*Xmea_exp_expand ;
            end
            l = l+2*size(X,2);
        elseif MeaSpecMatrix(i) == 6 % ZINB
            if MeaExpMatrix(i) == 0
                X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)']; ...
            else
            X = [ones(NRep*NP,1), LV(MeaMatrix(:,i)'== 1,:)', Xmea_exp]; ...
            end
        bzip = bmea(l+1:l+size(X,2)); ...
            bpoiss = bmea(l+size(X,2)+1:l+2*size(X,2));
        fit = reshape(X*bpoiss, NRep, NP)';
        theta = exp(bmea(l+2*size(X,2)+1));
        p = reshape(exp(X*bzip), NRep, NP)';
        p = p./(1+p);
        L = zeros(NP, NRep);
        lam = exp(fit);
        u = theta./(theta+lam);
        IndxZIP = Xmea(:,i) == 0;
        L(IndxZIP,:) = p(IndxZIP,:) + (1-p(IndxZIP,:)).*(u(IndxZIP,:).^theta);
        if RealMin == 1
            L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*min(gamma(theta+Xmea(~IndxZIP,i*ones(NRep,1))), realmax)./(gamma(theta).*min(gamma(Xmea(~IndxZIP,i*ones(NRep,1))+1),realmax));
        else
            L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*exp(gammaln(theta+Xmea(~IndxZIP,i*ones(NRep,1))) - gammaln(Xmea(~IndxZIP,i*ones(NRep,1))+1)-gammaln(theta));
        end
        L(~IndxZIP,:) = L(~IndxZIP,:).*(u(~IndxZIP,:).^theta).*((1-u(~IndxZIP,:)).^Xmea(~IndxZIP,i*ones(NRep,1)));
        L_mea = L_mea.*L; ...
            
    % Calculations for gradient
    grad_tmp1 = zeros(NP, NRep);
    grad_tmp2 = zeros(NP, NRep);
    grad_tmp3 = zeros(NP, NRep); % for theta
    if RealMin == 1
        grad_tmp1(IndxZIP,:) = (-(1-p(IndxZIP,:)).*p(IndxZIP,:).*(u(IndxZIP,:).^theta)+p(IndxZIP,:)-p(IndxZIP,:).^2)./max(L(IndxZIP,:),realmin);
        grad_tmp2(IndxZIP,:) = (p(IndxZIP,:)-1).*lam(IndxZIP,:).*(u(IndxZIP,:).^(theta+1))./max(L(IndxZIP,:),realmin);
        grad_tmp3(IndxZIP,:) = (1-p(IndxZIP,:)).*(u(IndxZIP,:).^theta).*(log(u(IndxZIP,:))+1-u(IndxZIP,:))./max(L(IndxZIP,:),realmin);
    else
        grad_tmp1(IndxZIP,:) = (-(1-p(IndxZIP,:)).*p(IndxZIP,:).*(u(IndxZIP,:).^theta)+p(IndxZIP,:)-p(IndxZIP,:).^2)./L(IndxZIP,:);
        grad_tmp2(IndxZIP,:) = (p(IndxZIP,:)-1).*lam(IndxZIP,:).*(u(IndxZIP,:).^(theta+1))./L(IndxZIP,:);
        grad_tmp3(IndxZIP,:) = theta*(1-p(IndxZIP,:)).*(u(IndxZIP,:).^theta).*(log(u(IndxZIP,:))+1-u(IndxZIP,:))./L(IndxZIP,:);
    end
    % For ZIP
    grad_tmp1(~IndxZIP,:) = -p(~IndxZIP,:);
    % For NB
    grad_tmp2(~IndxZIP,:) = u(~IndxZIP,:).*(Xmea(~IndxZIP,i*ones(NRep,1)) - lam(~IndxZIP,:));
    grad_tmp3(~IndxZIP,:) = psi(theta+Xmea(~IndxZIP,i*ones(NRep,1)))*theta - psi(theta)*theta+theta*log(u(~IndxZIP,:))...
        + theta*(1-u(~IndxZIP,:)) - Xmea(~IndxZIP,i*ones(NRep,1)).*u(~IndxZIP,:); % theta;
    
    LVindx = find(MeaMatrix(:,i)'== 1);
    % gradient for structural equation
    if sum(MeaMatrix(:,i)'== 1) > 1
        bx1 = bzip(2:1+sum(MeaMatrix(:,i)'== 1)); % parameters for LV
        bx1 = permute(bx1(:, ones(NP, 1), ones(NRep, 1), ones(NVarStr,1)), [2 3 4 1]);
        bx2 = bpoiss(2:1+sum(MeaMatrix(:,i)'== 1)); % parameters for LV
        bx2 = permute(bx2(:, ones(NP, 1), ones(NRep, 1), ones(NVarStr,1)), [2 3 4 1]);
        
        gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp1(:,:, ones(NVarStr,1), ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_der(:,:,:,LVindx).*bx1+...
            grad_tmp2(:,:, ones(NVarStr,1), ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_der(:,:,:,LVindx).*bx2;
    else
        gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp1(:,:, ones(NVarStr,1)).*LV_der(:,:,:,LVindx)*bzip(2)+...
            grad_tmp2(:,:, ones(NVarStr,1)).*LV_der(:,:,:,LVindx)*bpoiss(2);
    end
    gmea(:,:,l+1) = grad_tmp1; % constant
    gmea(:,:,l+2:l+1+sum(MeaMatrix(:,i)'== 1)) = grad_tmp1(:,:, ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_expand(:,:,MeaMatrix(:,i)'== 1); % parameters for LV
    gmea(:,:,l+sum(MeaMatrix(:,i)'== 1)+MeaExpMatrix(i)*NVarMeaExp+2) = grad_tmp2; % constant
    gmea(:,:,l+sum(MeaMatrix(:,i)'== 1)+MeaExpMatrix(i)*NVarMeaExp+3:l+2*sum(MeaMatrix(:,i)'== 1)+MeaExpMatrix(i)*NVarMeaExp+2) = ...
        grad_tmp2(:,:, ones(sum(MeaMatrix(:,i)'== 1),1)).*LV_expand(:,:,MeaMatrix(:,i)'== 1); % parameters for LV
    
    if MeaExpMatrix(i) ~= 0
        gmea(:,:,l+2+sum(MeaMatrix(:,i)'== 1):l+1+sum(MeaMatrix(:,i)'== 1)+NVarMeaExp) = ...
            grad_tmp1(:,:,ones(NVarMeaExp,1)).*Xmea_exp_expand ;
        gmea(:,:,l+2*sum(MeaMatrix(:,i)'== 1)+MeaExpMatrix(i)*NVarMeaExp+3:l+2*sum(MeaMatrix(:,i)'== 1)+2*MeaExpMatrix(i)*NVarMeaExp+2) = ...
            grad_tmp2(:,:,ones(NVarMeaExp,1)).*Xmea_exp_expand ;
    end
    gmea(:,:, l+2*sum(MeaMatrix(:,i)'== 1)+2*MeaExpMatrix(i)*NVarMeaExp+3) = grad_tmp3;
    
    l = l+2*size(X,2)+1; ...
        end
    end
    probs = probs.*L_mea;
    %     p = max(realmin,mean(probs,2));
    %     p = mean(probs,2);
    if RealMin == 1
        p = max(realmin,mean(probs,2));
    else
        p = mean(probs,2);
    end
    f = -log(p);
    
    
    gstr = reshape(gstr, NP, NRep, NVarStr*NLatent);
    gmnl = reshape(gmnl, NP, NRep, NVarA*(1+NLatent));
    %     g = squeeze(mean(probs(:,:, ones(NVarA*(1+NLatent),1)).*gmnl,2)); %
    g = reshape(mean(bsxfun(@times,probs,gmnl),2),[NP,NVarA*(1+NLatent)]); %
    %     g3 = squeeze(mean(probs(:,:,ones(length(bmea),1)).*gmea,2)); % NP x NVarmea
    g3 = reshape(mean(bsxfun(@times,probs,gmea),2),[NP,size(bmea,1)]); % NP x NVarmea
    %     g2 = squeeze(mean(probs(:,:, ones(NVarstr*NLatent,1)).*gstr,2)); % NP x NVarstr*NLatent
    g2 = reshape(mean(bsxfun(@times,probs,gstr),2),[NP,NVarStr*NLatent]); % NP x NVarstr*NLatent
    if NVarM > 0
        gm =  g(:,repmat(1:NVarA,1,EstimOpt.NVarM)).*(Xm(:,kron(1:EstimOpt.NVarM,ones(1,NVarA))));
        g = [g(:,1:EstimOpt.NVarA), gm, g(:, EstimOpt.NVarA+1:end),g2,g3];
    else
        g = [g,g2,g3];
    end
    if NVarS > 0
        gss = reshape(mean(bsxfun(@times, gs, probs),2), NP, NVarS);
        g = [g(:,1:NVarA*(1+NVarM+NLatent)), gss, g(:, NVarA*(1+NVarM+NLatent)+1:end)];
    end
    
    g = -g./p(:,ones(1,length(B)));
    
end

% f = -log(mean(probs.*L_mea,2));

