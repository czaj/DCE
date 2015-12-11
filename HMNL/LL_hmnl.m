function [f,g] = LL_hmnl(Y,Xa,X_str,X_mea,Xmea_exp, err_sliced,EstimOpt,B)

% save tmp_LL_hmnl
% return

ba = B(1:EstimOpt.NVarA); ... % b atrybutów
bl = reshape(B(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NLatent+1)), EstimOpt.NVarA, EstimOpt.NLatent); ... % b interakcji z LV
bstr = reshape(B(EstimOpt.NVarA*(EstimOpt.NLatent+1)+1:(EstimOpt.NVarA+EstimOpt.NVarstr)*EstimOpt.NLatent+EstimOpt.NVarA), EstimOpt.NVarstr, EstimOpt.NLatent); ... b równania struktury
bmea = B((EstimOpt.NVarA+EstimOpt.NVarstr)*EstimOpt.NLatent+EstimOpt.NVarA+1:end); ... % b measurement

LV_tmp = X_str*bstr; ... % NP x NLatent
LV_tmp = reshape(permute(LV_tmp(:,:, ones(EstimOpt.NRep,1)),[2 3 1]), EstimOpt.NLatent, EstimOpt.NRep*EstimOpt.NP); ...
LV_tmp = LV_tmp + err_sliced; ... % NLatent x NRep*NP
mLV = mean(LV_tmp,2); ...
sLV = std(LV_tmp,0,2); ...
LV = (LV_tmp - mLV(:,ones(1,size(LV_tmp,2))))./sLV(:,ones(1,size(LV_tmp,2))); ... % normalilzing for 0 mean and std

b_mtx = ba(:,ones(EstimOpt.NRep*EstimOpt.NP,1)) + bl*LV; ... % NVarA x NRep*NP
if EstimOpt.WTP_space > 0
    if EstimOpt.NumGrad == 0
       b_mtx_grad = reshape(b_mtx ,EstimOpt.NVarA,EstimOpt.NRep,EstimOpt.NP); ...% needed for gradient calculation in WTP_space
    end
    b_mtx(1:end-EstimOpt.WTP_space,:) = b_mtx(1:end-EstimOpt.WTP_space,:).*b_mtx(EstimOpt.WTP_matrix,:); ...

end

probs = zeros(EstimOpt.NP,EstimOpt.NRep); ...

%for n = 1:EstimOpt.NP
 %   U = reshape(exp(Xa(:,:,n)*b_mtx(:,((n-1)*EstimOpt.NRep+1):n*EstimOpt.NRep)),EstimOpt.NAlt,EstimOpt.NCT,EstimOpt.NRep); ...
    %U(isnan(U)) = 0; ... % skip missing ALT
    
  %  U_sum = sum(U,1); ...
  %  P = sum(Y(:,:,:,n).*U ./ U_sum(ones(EstimOpt.NAlt,1),:,:),1); ... 
   % P(isnan(P)) = 1; ... % 
  %  probs(:,n,:) = prod(P,2); ...
%end; 
%probs = squeeze(probs); % taking this out of the loop saves almost 10% of time

if EstimOpt.NumGrad == 1 % numerical gradient
    
    if any(isnan(Xa(:))) == 0 ... % faster version for complete dataset

        for n = 1:EstimOpt.NP ...

            U = reshape(Xa(:,:,n)*b_mtx(:,((n-1)*EstimOpt.NRep+1):n*EstimOpt.NRep),EstimOpt.NAlt,EstimOpt.NCT, EstimOpt.NRep); ...
            U_max = max(U); ...
            U = exp(U - U_max(ones(EstimOpt.NAlt,1),:,:)); ... % rescale utility to avoid exploding
            U_sum = reshape(sum(U,1),EstimOpt.NCT,EstimOpt.NRep); ...
            U_selected = reshape(U(Y(:,n*ones(EstimOpt.NRep,1))==1),EstimOpt.NCT,EstimOpt.NRep); ...   
            probs(n,:) = prod(U_selected ./ U_sum,1);
        end; ...
    else ...
        for n = 1:EstimOpt.NP ...

            U = reshape(Xa(~isnan(Y(:,n)),:,n)*b_mtx(:,((n-1)*EstimOpt.NRep+1):n*EstimOpt.NRep),EstimOpt.NAlt,EstimOpt.NCT-sum(isnan(Y(1:EstimOpt.NAlt:end,n))),EstimOpt.NRep); ...
            U_max = max(U); ...
            U = exp(U - U_max(ones(EstimOpt.NAlt,1),:,:)); ... % NAlt x NCT - NaNs x NRep
            U_sum = reshape(nansum(U,1),EstimOpt.NCT-sum(isnan(Y(1:EstimOpt.NAlt:end,n))),EstimOpt.NRep); ...
            U_selected = reshape(U(Y(~isnan(Y(:,n)),n*ones(EstimOpt.NRep,1))==1),EstimOpt.NCT-sum(isnan(Y(1:EstimOpt.NAlt:end,n))),EstimOpt.NRep); ...   
            probs(n,:) = prod(U_selected ./ U_sum,1);
        end; ...
    end ...

    L_mea = ones(EstimOpt.NP,EstimOpt.NRep); ...
    l = 0; ...

    if EstimOpt.NVarmea_exp > 0
        Xmea_exp = reshape(Xmea_exp,1,EstimOpt.NP,EstimOpt.NVarmea_exp);
        Xmea_exp = reshape(Xmea_exp(ones(EstimOpt.NRep,1),:,:), EstimOpt.NP*EstimOpt.NRep, EstimOpt.NVarmea_exp);
    end
    for i = 1:size(X_mea,2)
        if EstimOpt.MeaSpecMatrix(i) == 0 % OLS
            if EstimOpt.MeaExpMatrix(i) == 0
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)']; ...
            else
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)', Xmea_exp]; ...
            end
            b = bmea(l+1:l+size(X,2)+1); ...
            L_mea = L_mea.*normpdf(X_mea(:,i*ones(EstimOpt.NRep,1)),reshape(X*b(1:end-1), EstimOpt.NRep, EstimOpt.NP)',exp(b(end))); ...
            l = l+size(X,2)+1; ...
        elseif EstimOpt.MeaSpecMatrix(i) == 1 % MNL 
            UniqueMea = unique(X_mea(:,i));  ...
            k = length(UniqueMea)-1; ...
            if EstimOpt.MeaExpMatrix(i) == 0
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)']; ...
            else
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)', Xmea_exp]; ...
            end
            V = exp(X*reshape([zeros(size(X,2),1); bmea(l+1:l+size(X,2)*k)], size(X,2), k+1)); ... % NRep*NP x unique values of attitude
            Vsum = sum(V,2); ...
            V = reshape(V./Vsum(:, ones(k+1,1)), EstimOpt.NRep, EstimOpt.NP, k+1); ... NRep x NP x unique
            V = permute(V, [2 1 3]); % NP x NRep x unique
            L = zeros(EstimOpt.NP,EstimOpt.NRep); ...
            for j = 1:length(UniqueMea)
                L(X_mea(:,i) == UniqueMea(j),:) = V(X_mea(:,i) == UniqueMea(j),:,j); ...
            end
            L_mea = L_mea.*L; ...
            l = l + size(X,2)*k; ... 

        elseif EstimOpt.MeaSpecMatrix(i) == 2 % Ordered Probit
            UniqueMea = unique(X_mea(:,i));  ...
            k = length(UniqueMea)-1; ...
            if EstimOpt.MeaExpMatrix(i) == 0
                X = LV(EstimOpt.MeaMatrix(:,i)'== 1,:)'; ...
            else
                X = [LV(EstimOpt.MeaMatrix(:,i)'== 1,:)', Xmea_exp]; ...
            end
            tmp = (EstimOpt.MeaExpMatrix(i)~= 0)*EstimOpt.NVarmea_exp;
            b = bmea(l+1:l+ k+ size(X,2)); ...
            Xb = reshape(X*b(1:sum(EstimOpt.MeaMatrix(:,i),1)+tmp), EstimOpt.NRep, EstimOpt.NP)'; ... % NP x NRep
            alpha = cumsum([b(sum(EstimOpt.MeaMatrix(:,i))+tmp+1); exp(b(sum(EstimOpt.MeaMatrix(:,i))+tmp+2:end))]); ...
            L = zeros(EstimOpt.NP,EstimOpt.NRep); ...
            L(X_mea(:,i) == min(UniqueMea),:) = normcdf(alpha(1)-Xb(X_mea(:,i) == min(UniqueMea),:)); ...
            L(X_mea(:,i) == max(UniqueMea),:) = 1-normcdf(alpha(end)-Xb(X_mea(:,i) == max(UniqueMea),:)); ...
            for j = 2:k
                L(X_mea(:,i) == UniqueMea(j),:) = normcdf(alpha(j)-Xb(X_mea(:,i) == UniqueMea(j),:)) - normcdf(alpha(j-1)-Xb(X_mea(:,i) == UniqueMea(j),:)); ...
            end
            L_mea = L_mea.*L; ...
            l = l+k+size(X,2); ...   
        elseif EstimOpt.MeaSpecMatrix(i) == 3 % Poisson
            if EstimOpt.MeaExpMatrix(i) == 0
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)']; ...
            else
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)', Xmea_exp]; ...
            end
            b = bmea(l+1:l+size(X,2)); ...
            fit = reshape(X*b, EstimOpt.NRep, EstimOpt.NP)';
            lam = exp(fit);
            L = exp(fit.*X_mea(:,i*ones(EstimOpt.NRep,1))-lam)./min(gamma(X_mea(:,i*ones(EstimOpt.NRep,1))+1),realmax);
            L_mea = L_mea.*L; ...
            l = l+size(X,2); ...        
        elseif EstimOpt.MeaSpecMatrix(i) == 4 % NB
            if EstimOpt.MeaExpMatrix(i) == 0
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)']; ...
            else
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)', Xmea_exp]; ...
            end
            b = bmea(l+1:l+size(X,2)); ...
            fit = reshape(X*b, EstimOpt.NRep, EstimOpt.NP)';
            lam = exp(fit);
            theta = exp(bmea(l+size(X,2)+1));
            u = theta./(theta+lam);  
            L = min(gamma(theta+X_mea(:,i*ones(EstimOpt.NRep,1))), realmax)./(gamma(theta).*min(gamma(X_mea(:,i*ones(EstimOpt.NRep,1))+1),realmax));
            L = L.*(u.^theta).*((1-u).^X_mea(:,i*ones(EstimOpt.NRep,1)));
            L_mea = L_mea.*L; ...
            l = l+size(X,2)+1; ...
        elseif EstimOpt.MeaSpecMatrix(i) == 5 % ZIP
            if EstimOpt.MeaExpMatrix(i) == 0
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)']; ...
            else
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)', Xmea_exp]; ...
            end
            bzip = bmea(l+1:l+size(X,2)); ...
            bpoiss = bmea(l+size(X,2)+1:l+2*size(X,2));
            fit = reshape(X*bpoiss, EstimOpt.NRep, EstimOpt.NP)';
            p = reshape(exp(X*bzip), EstimOpt.NRep, EstimOpt.NP)';
            p = p./(1+p);
            L = zeros(EstimOpt.NP, EstimOpt.NRep);
            lam = exp(fit);
            IndxZIP = X_mea(:,i) == 0;
            L(IndxZIP,:) = p(IndxZIP,:) + (1-p(IndxZIP,:)).*exp(-lam(IndxZIP,:));
            L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*exp(fit(~IndxZIP,:).*X_mea(~IndxZIP,i*ones(EstimOpt.NRep,1))-lam(~IndxZIP,:))./min(gamma(X_mea(~IndxZIP,i*ones(EstimOpt.NRep,1))+1),realmax);
            L_mea = L_mea.*L; ...
            l = l+2*size(X,2); ...
        end
    end
    
    f = -log(max(realmin,mean(probs.*L_mea,2)));

else % analitical
    if ~exist('b_mtx_grad','var')
        b_mtx_grad = [];
    end
    gmnl = zeros(EstimOpt.NP, EstimOpt.NRep, EstimOpt.NVarA, (1+EstimOpt.NLatent)); % gradient for mnl parameters
    gstr = zeros(EstimOpt.NP, EstimOpt.NRep, EstimOpt.NVarstr,EstimOpt.NLatent); % gradient for parameters from structural equations
    gmea = zeros(EstimOpt.NP, EstimOpt.NRep, length(bmea));% gradient for other parameters
    
    
    % terms from LV normalization
    mXstr = mean(X_str,1);
    X_str_expand = reshape(X_str- mXstr(ones(EstimOpt.NP,1),:), 1,EstimOpt.NP, EstimOpt.NVarstr);
    X_str_expand = reshape(X_str_expand(ones(EstimOpt.NRep*EstimOpt.NLatent,1), :,:), EstimOpt.NLatent,EstimOpt.NRep*EstimOpt.NP, EstimOpt.NVarstr);
    LV_tmp = LV_tmp - mLV(:,ones(1,size(LV_tmp,2))); %NLatent x NRep*NP
    LV_std = sum(LV_tmp(:,:,ones(EstimOpt.NVarstr,1)).*X_str_expand,2)/(EstimOpt.NRep*EstimOpt.NP-1);% NLatent x 1 x NVarstr
    
    X_str_expand = X_str_expand./sLV(:, ones(EstimOpt.NRep*EstimOpt.NP,1),ones(EstimOpt.NVarstr,1));
    LV_tmp = LV_tmp./(sLV(:, ones(EstimOpt.NRep*EstimOpt.NP,1)).^3);
    LV_tmp = LV_tmp(:,:, ones(EstimOpt.NVarstr,1));

    LV_der = reshape(X_str_expand - LV_tmp.*LV_std(:,ones(EstimOpt.NRep*EstimOpt.NP,1),:), EstimOpt.NLatent, EstimOpt.NRep, EstimOpt.NP, EstimOpt.NVarstr); % NLatent x NRep x NP x NVarstr
    %LV_der = reshape(permute(LV_der, [3 2 4 1]), EstimOpt.NP, EstimOpt.NRep, EstimOpt.NVarstr*EstimOpt.NLatent);
    LV_der = permute(LV_der, [3 2 4 1]); % NP x NRep x NVarstr x NLatent
    LV_expand = permute(reshape(LV',EstimOpt.NRep,EstimOpt.NP,EstimOpt.NLatent), [2 1 3])  ; 
    
    
    if any(isnan(Xa(:))) == 0 ... % faster version for complete dataset
%         b_mtx_grad = [];
         parfor n = 1:EstimOpt.NP
            U = reshape(Xa(:,:,n)*b_mtx(:,((n-1)*EstimOpt.NRep+1):n*EstimOpt.NRep),EstimOpt.NAlt,EstimOpt.NCT,EstimOpt.NRep); ... % NAlt x NCT x NRep
            U_max = max(U); ...
            U = exp(U - U_max(ones(EstimOpt.NAlt,1),:,:)); ... % rescale utility to avoid exploding 
            U_sum = reshape(sum(U,1),1,EstimOpt.NCT,EstimOpt.NRep); ... 
            U_prob = U./U_sum(ones(EstimOpt.NAlt,1,1),:,:); ... % NAlt x NCT x NRep
            probs(n,:) = prod(reshape(U_prob(Y(:,n*ones(EstimOpt.NRep,1))==1),EstimOpt.NCT,EstimOpt.NRep),1); ... % 1 x NRep
            %probs(n,:) = U_prod;    
        
            % calculations for gradient
            U_prob = reshape(U_prob, EstimOpt.NAlt*EstimOpt.NCT,1, EstimOpt.NRep); ... % NAlt*NCT x NVarA x NRep
            if EstimOpt.WTP_space == 0   
                X_hat = sum(reshape(U_prob(:,ones(1,EstimOpt.NVarA,1),:).* Xa(:,:, n*ones(EstimOpt.NRep,1)), EstimOpt.NAlt, EstimOpt.NCT, EstimOpt.NVarA, EstimOpt.NRep),1); ...

                if EstimOpt.NCT ~= 1
                    F = Xa(Y(:,n) == 1,:,n*ones(EstimOpt.NRep,1)) - squeeze(X_hat); ... %NCT x NVarA x NRep 
                    sumFsqueezed = squeeze(sum(F,1)); ... %NVarA x NRep     
                else
                    sumFsqueezed = squeeze(Xa(Y(:,n) == 1,:,n*ones(EstimOpt.NRep,1))) - squeeze(X_hat); %NVarA x NRep 
                end   

                
                sumFsqueezed_LV = sumFsqueezed'*bl; % NRep x NLatent
            else
                b_mtx_wtp = reshape(b_mtx(:,((n-1)*EstimOpt.NRep+1):n*EstimOpt.NRep), 1, EstimOpt.NVarA, EstimOpt.NRep);
                Xalpha = Xa(:,1:end-EstimOpt.WTP_space, n*ones(EstimOpt.NRep,1)).*b_mtx_wtp(ones(EstimOpt.NAlt*EstimOpt.NCT,1),EstimOpt.WTP_matrix,:);
                % for non-cost variables
                X_hat1 = sum(reshape(U_prob(:,ones(1,EstimOpt.NVarA-EstimOpt.WTP_space),:).* Xalpha, EstimOpt.NAlt, EstimOpt.NCT, EstimOpt.NVarA-EstimOpt.WTP_space, EstimOpt.NRep),1); ...
                F1 = Xalpha(Y(:,n) == 1,:,:) - squeeze(X_hat1); ... %NCT x NVarA-WTP_space x NRep   
                % for cost variables
                if EstimOpt.WTP_space == 1 % without starting the loop
                    Xbmxlfit = Xa(:,1:end-EstimOpt.WTP_space,n)*b_mtx_grad(1:end-EstimOpt.WTP_space,:,n);
                    pX = squeeze(Xa(:,EstimOpt.NVarA, n*ones(EstimOpt.NRep,1))) + Xbmxlfit;
                    X_hat2 = sum(reshape(squeeze(U_prob).*pX, EstimOpt.NAlt, EstimOpt.NCT, EstimOpt.WTP_space, EstimOpt.NRep),1); ...

                else
                    pX = zeros(EstimOpt.NCT*EstimOpt.NAlt, EstimOpt.WTP_space, EstimOpt.NRep);
                    for i = 1:EstimOpt.WTP_space
                        Xbmxlfit = Xa(:,EstimOpt.WTP_matrix == EstimOpt.NVarA-EstimOpt.WTP_space+i,n)*b_mtx_grad(EstimOpt.WTP_matrix == EstimOpt.NVarA-EstimOpt.WTP_space+i,:,n);
                        pX(:,i,:) = squeeze(Xa(:,EstimOpt.NVarA-EstimOpt.WTP_space+i, n*ones(EstimOpt.NRep,1))) + Xbmxlfit;

                    end
                    X_hat2 = sum(reshape(U_prob(:,ones(1,EstimOpt.WTP_space),:).* pX, EstimOpt.NAlt, EstimOpt.NCT, EstimOpt.WTP_space, EstimOpt.NRep),1); ...
                end
                F2 = pX(Y(:,n) == 1,:,:) - squeeze(X_hat2); ... %NCT x WTP_space x NRep  
                

                sumFsqueezed = [squeeze(sum(F1,1));squeeze(sum(F2,1)) ]; ... %NVarA x NRep
                sumFsqueezed_LV = sumFsqueezed'*bl; % NRep x NLatent
            end
        
            sumFsqueezed = sumFsqueezed'; % NRep x NVarA
            gmnl(n,:,:,:) = sumFsqueezed(:,:,ones(1+EstimOpt.NLatent,1)); 
            sumFsqueezed_LV = permute(sumFsqueezed_LV(:,:, ones(EstimOpt.NVarstr,1)), [1 3 2]);
            gstr(n,:,:,:) = sumFsqueezed_LV;
            

        end; .....
    else ...
    
        parfor n = 1:EstimOpt.NP
            YnanInd = ~isnan(Y(:,n));
            NCT_n =EstimOpt.NCT-sum(isnan(Y(1:EstimOpt.NAlt:end,n)));
            U = reshape(Xa(YnanInd,:,n)*b_mtx(:,((n-1)*EstimOpt.NRep+1):n*EstimOpt.NRep),EstimOpt.NAlt,NCT_n,EstimOpt.NRep); ... % NAlt x NCT x NRep
            U_max = max(U); ...
            U = exp(U - U_max(ones(EstimOpt.NAlt,1),:,:)); ... % rescale utility to avoid exploding 
            U_sum = reshape(sum(U,1),1,NCT_n,EstimOpt.NRep); ... 
            U_prob = U./U_sum(ones(EstimOpt.NAlt,1,1),:,:); ... % NAlt x NCT x NRep
            probs(n,:) = prod(reshape(U_prob(Y(YnanInd,n*ones(EstimOpt.NRep,1))==1),NCT_n,EstimOpt.NRep),1); ... % 1 x NRep
            %probs(n,:) = U_prod;    
        
            % calculations for gradient
            U_prob = reshape(U_prob, EstimOpt.NAlt*NCT_n,1, EstimOpt.NRep); ... % NAlt*NCT x NVarA x NRep
            if EstimOpt.WTP_space == 0   
                X_hat = sum(reshape(U_prob(:,ones(1,EstimOpt.NVarA,1),:).* Xa(YnanInd,:, n*ones(EstimOpt.NRep,1)), EstimOpt.NAlt, NCT_n, EstimOpt.NVarA, EstimOpt.NRep),1); ...

                if NCT_n ~= 1
                    F = Xa(Y(:,n) == 1,:,n*ones(EstimOpt.NRep,1)) - squeeze(X_hat); ... %NCT x NVarA x NRep 
                    sumFsqueezed = squeeze(sum(F,1)); ... %NVarA x NRep
                    
                else
                    sumFsqueezed = squeeze(Xa(Y(:,n) == 1,:,n*ones(EstimOpt.NRep,1))) - squeeze(X_hat); %NVarA x NRep 
                end   

                sumFsqueezed_LV = sumFsqueezed'*bl; % NRep x NLatent
            else
                b_mtx_wtp = reshape(b_mtx(:,((n-1)*EstimOpt.NRep+1):n*EstimOpt.NRep), 1, EstimOpt.NVarA, EstimOpt.NRep);
                Xalpha = Xa(YnanInd,1:end-EstimOpt.WTP_space, n*ones(EstimOpt.NRep,1)).*b_mtx_wtp(ones(EstimOpt.NAlt*NCT_n,1),EstimOpt.WTP_matrix,:);
                % for non-cost variables
                X_hat1 = sum(reshape(U_prob(:,ones(1,EstimOpt.NVarA-EstimOpt.WTP_space),:).* Xalpha, EstimOpt.NAlt, NCT_n, EstimOpt.NVarA-EstimOpt.WTP_space, EstimOpt.NRep),1); ...
                F1 = Xalpha(Y(YnanInd,n) == 1,:,:) - squeeze(X_hat1); ... %NCT x NVarA-WTP_space x NRep   
                % for cost variables
                if EstimOpt.WTP_space == 1 % without starting the loop
                    Xbmxlfit = Xa(YnanInd,1:end-EstimOpt.WTP_space,n)*b_mtx_grad(1:end-EstimOpt.WTP_space,:,n);
                    pX = squeeze(Xa(YnanInd,EstimOpt.NVarA, n*ones(EstimOpt.NRep,1))) + Xbmxlfit;
                    X_hat2 = sum(reshape(squeeze(U_prob).*pX, EstimOpt.NAlt, NCT_n, EstimOpt.WTP_space, EstimOpt.NRep),1); ...

                else
                    pX = zeros(NCT_n*EstimOpt.NAlt, EstimOpt.WTP_space, EstimOpt.NRep);
                    for i = 1:EstimOpt.WTP_space
                        Xbmxlfit = Xa(YnanInd,EstimOpt.WTP_matrix == EstimOpt.NVarA-EstimOpt.WTP_space+i,n)*b_mtx_grad(EstimOpt.WTP_matrix == EstimOpt.NVarA-EstimOpt.WTP_space+i,:,n);
                        pX(:,i,:) = squeeze(Xa(YnanInd,EstimOpt.NVarA-EstimOpt.WTP_space+i, n*ones(EstimOpt.NRep,1))) + Xbmxlfit;
                        
                    end
                    X_hat2 = sum(reshape(U_prob(:,ones(1,EstimOpt.WTP_space),:).* pX, EstimOpt.NAlt, NCT_n, EstimOpt.WTP_space, EstimOpt.NRep),1); ...
                end
                F2 = pX(Y(YnanInd,n) == 1,:,:) - squeeze(X_hat2); ... %NCT x WTP_space x NRep  
                
                
                sumFsqueezed = [squeeze(sum(F1,1));squeeze(sum(F2,1)) ]; ... %NVarA x NRep
                sumFsqueezed_LV = sumFsqueezed'*bl; % NRep x NLatent
            end
        
            sumFsqueezed = sumFsqueezed'; % NRep x NVarA
            gmnl(n,:,:,:) = sumFsqueezed(:,:,ones(1+EstimOpt.NLatent,1)); 
            sumFsqueezed_LV = permute(sumFsqueezed_LV(:,:, ones(EstimOpt.NVarstr,1)), [1 3 2]);
            gstr(n,:,:,:) = sumFsqueezed_LV;
            

        end; .....
    
    end ...
        
    gstr = gstr.*LV_der;
    LV_expand = reshape(LV_expand, EstimOpt.NP, EstimOpt.NRep,1, EstimOpt.NLatent );
    gmnl(:,:,:,2:end) = gmnl(:,:,:,2:end).*LV_expand(:,:, ones(EstimOpt.NVarA,1),:);
    LV_expand = squeeze(LV_expand);
    L_mea = ones(EstimOpt.NP,EstimOpt.NRep); ...
    l = 0; ...

    if EstimOpt.NVarmea_exp > 0
        Xmea_exp = reshape(Xmea_exp,1,EstimOpt.NP,EstimOpt.NVarmea_exp);
        Xmea_exp = reshape(Xmea_exp(ones(EstimOpt.NRep,1),:,:), EstimOpt.NP*EstimOpt.NRep, EstimOpt.NVarmea_exp);
    end
    Xmea_exp_expand = permute(reshape(Xmea_exp, EstimOpt.NRep, EstimOpt.NP,EstimOpt.NVarmea_exp), [2 1 3]);

    for i = 1:size(X_mea,2)
        if EstimOpt.MeaSpecMatrix(i) == 0 % OLS
            if EstimOpt.MeaExpMatrix(i) == 0
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)']; ...
            else
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)', Xmea_exp]; ...
            end
            b = bmea(l+1:l+size(X,2)+1); ...
            fit = reshape(X*b(1:end-1), EstimOpt.NRep, EstimOpt.NP)';
            L_mea = L_mea.*normpdf(X_mea(:,i*ones(EstimOpt.NRep,1)),fit,exp(b(end))); ...
            grad_tmp = - (fit - X_mea(:,i*ones(EstimOpt.NRep,1)))/exp(2*b(end)); % NP x NRep 
            LVindx = find(EstimOpt.MeaMatrix(:,i)'== 1);
            if sum(EstimOpt.MeaMatrix(:,i)'== 1) > 1
                bx = b(2:1+sum(EstimOpt.MeaMatrix(:,i)'== 1));
                bx = permute(bx(:, ones(EstimOpt.NP, 1), ones(EstimOpt.NRep, 1), ones(EstimOpt.NVarstr,1)), [2 3 4 1]);
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp(:,:, ones(EstimOpt.NVarstr,1), ones(sum(EstimOpt.MeaMatrix(:,i)'== 1),1)).*LV_der(:,:,:,LVindx).*bx;
            else
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp(:,:, ones(EstimOpt.NVarstr,1)).*LV_der(:,:,:,LVindx)*b(2);
            end
            
            gmea(:,:,l+1) = grad_tmp; % constant
            %LV_expand = permute(reshape(LV(EstimOpt.MeaMatrix(:,i)'== 1,:)',EstimOpt.NRep,EstimOpt.NP,sum(EstimOpt.MeaMatrix(:,i)'== 1)), [2 1 3])  ; 
            gmea(:,:,l+2:l+1+sum(EstimOpt.MeaMatrix(:,i)'== 1)) = grad_tmp(:,:, ones(sum(EstimOpt.MeaMatrix(:,i)'== 1),1)).*LV_expand(:,:,EstimOpt.MeaMatrix(:,i)'== 1); % parameters for LV
            if EstimOpt.MeaExpMatrix(i) == 0
                gmea(:,:,l+2+sum(EstimOpt.MeaMatrix(:,i)'== 1)) = -1+((fit - X_mea(:,i*ones(EstimOpt.NRep,1))).^2)/(exp(2*b(end))) ; % variance 
            else
                gmea(:,:,l+2+sum(EstimOpt.MeaMatrix(:,i)'== 1):l+1+sum(EstimOpt.MeaMatrix(:,i)'== 1)+EstimOpt.NVarmea_exp) = ...
                    grad_tmp(:,:,ones(EstimOpt.NVarmea_exp,1)).*Xmea_exp_expand ;
                gmea(:,:,l+2+sum(EstimOpt.MeaMatrix(:,i)'== 1)+EstimOpt.NVarmea_exp) = -1+((fit - X_mea(:,i*ones(EstimOpt.NRep,1))).^2)/(exp(2*b(end))) ; % variance 
                
            end
            l = l+size(X,2)+1; ...
        elseif EstimOpt.MeaSpecMatrix(i) == 1 % MNL 
            UniqueMea = unique(X_mea(:,i));  ...
            k = length(UniqueMea)-1; ...
            if EstimOpt.MeaExpMatrix(i) == 0
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)']; ...
            else
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)', Xmea_exp]; ...
            end
            bx = reshape([zeros(size(X,2),1); bmea(l+1:l+size(X,2)*k)], size(X,2), k+1);
            V = exp(X*bx); ... % NRep*NP x unique values of attitude
            Vsum = sum(V,2); ...
            V = reshape(V./Vsum(:, ones(k+1,1)), EstimOpt.NRep, EstimOpt.NP, k+1); ... NRep x NP x unique
            V = permute(V, [2 1 3]); % NP x NRep x unique
            L = zeros(EstimOpt.NP,EstimOpt.NRep); ...
            L(X_mea(:,i) == UniqueMea(1),:) = V(X_mea(:,i) == UniqueMea(1),:,1); ...
            MAT = permute(reshape(X, EstimOpt.NRep, EstimOpt.NP, size(X,2)), [2 1 3]);
            for j = 2:length(UniqueMea)
                L(X_mea(:,i) == UniqueMea(j),:) = V(X_mea(:,i) == UniqueMea(j),:,j); ...
                
                gmea(X_mea(:,i) == UniqueMea(j),:,l+1+(j-2)*size(X,2):l+size(X,2)*(j-1)) = (1 - V(X_mea(:,i) == UniqueMea(j),:,j*ones(1, size(X,2)))).*MAT(X_mea(:,i) == UniqueMea(j),:,:);
                gmea(X_mea(:,i) ~= UniqueMea(j),:,l+1+(j-2)*size(X,2):l+size(X,2)*(j-1)) = (  - V(X_mea(:,i) ~= UniqueMea(j),:,j*ones(1, size(X,2)))).*MAT(X_mea(:,i) ~= UniqueMea(j),:,:);
            end
            alpha = bx(2:1+sum(EstimOpt.MeaMatrix(:,i),1),:); 
            alpha = alpha(:,:, ones(EstimOpt.NP,1), ones(EstimOpt.NRep,1));
            alphax = reshape(permute(alpha, [2 3 4 1]),(k+1)*EstimOpt.NP,EstimOpt.NRep, sum(EstimOpt.MeaMatrix(:,i)));
            alpha = permute(alpha, [3 4 2 1]);
            
            LVindx = find(EstimOpt.MeaMatrix(:,i)'== 1);
            X_mea_i = reshape(dummyvar(X_mea(:,i))', (k+1)*EstimOpt.NP,1);
            if sum(EstimOpt.MeaMatrix(:,i)'== 1) > 1
                
                for j = 1:sum(EstimOpt.MeaMatrix(:,i)'== 1)
                    grad_tmp = alphax(X_mea_i==1,:,j)  -sum(V.*alpha(:,:,:,j),3); % NP x NRep
                    gstr(:,:,:,LVindx(j)) = gstr(:,:,:,LVindx(j)) +grad_tmp(:,:, ones(EstimOpt.NVarstr,1)).*LV_der(:,:,:,LVindx(j));
                end
            else
                grad_tmp = alphax(X_mea_i==1,:)  -sum(V.*alpha,3); % NP x NRep
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx) +grad_tmp(:,:, ones(EstimOpt.NVarstr,1)).*LV_der(:,:,:,LVindx);
            end
            L_mea = L_mea.*L; ...
            l = l + size(X,2)*k; ... 

        elseif EstimOpt.MeaSpecMatrix(i) == 2 % Ordered Probit
            UniqueMea = unique(X_mea(:,i));  ...
            k = length(UniqueMea)-1; ...
            if EstimOpt.MeaExpMatrix(i) == 0
                X = LV(EstimOpt.MeaMatrix(:,i)'== 1,:)'; ...
            else
                X = [LV(EstimOpt.MeaMatrix(:,i)'== 1,:)', Xmea_exp]; ...
            end
            tmp = (EstimOpt.MeaExpMatrix(i)~= 0)*EstimOpt.NVarmea_exp;
            b = bmea(l+1:l+ k+ size(X,2)); ...
            Xb = reshape(X*b(1:sum(EstimOpt.MeaMatrix(:,i),1)+tmp), EstimOpt.NRep, EstimOpt.NP)'; ... % NP x NRep
            alpha = cumsum([b(sum(EstimOpt.MeaMatrix(:,i))+tmp+1); exp(b(sum(EstimOpt.MeaMatrix(:,i))+tmp+2:end))]); ...
            L = zeros(EstimOpt.NP,EstimOpt.NRep); ...
            grad_tmp = zeros(EstimOpt.NP,EstimOpt.NRep); ...
            L(X_mea(:,i) == min(UniqueMea),:) = normcdf(alpha(1)-Xb(X_mea(:,i) == min(UniqueMea),:)); ...
            L(X_mea(:,i) == max(UniqueMea),:) = 1-normcdf(alpha(end)-Xb(X_mea(:,i) == max(UniqueMea),:)); ...
            grad_tmp(X_mea(:,i) == min(UniqueMea),:) = normpdf(alpha(1)-Xb(X_mea(:,i) == min(UniqueMea),:))./max(L(X_mea(:,i) == min(UniqueMea),:),realmin); ...
            grad_tmp(X_mea(:,i) == max(UniqueMea),:) = -normpdf(alpha(end)-Xb(X_mea(:,i) == max(UniqueMea),:))./max(L(X_mea(:,i) == max(UniqueMea),:), realmin); ...
            
            for j = 2:k
                L(X_mea(:,i) == UniqueMea(j),:) = normcdf(alpha(j)-Xb(X_mea(:,i) == UniqueMea(j),:)) - normcdf(alpha(j-1)-Xb(X_mea(:,i) == UniqueMea(j),:)); ...
                grad_tmp(X_mea(:,i) == UniqueMea(j),:) = normpdf(alpha(j)-Xb(X_mea(:,i) == UniqueMea(j),:)) - normpdf(alpha(j-1)-Xb(X_mea(:,i) == UniqueMea(j),:)); ...   
                grad_tmp(X_mea(:,i) == UniqueMea(j),:) = grad_tmp(X_mea(:,i) == UniqueMea(j),:)./max(L(X_mea(:,i) == UniqueMea(j),:),realmin);
               
            end
            
            LVindx = find(EstimOpt.MeaMatrix(:,i)'== 1);
            if sum(EstimOpt.MeaMatrix(:,i)'== 1) > 1
                bx = b(1:sum(EstimOpt.MeaMatrix(:,i)'== 1));
                bx = permute(bx(:, ones(EstimOpt.NP, 1), ones(EstimOpt.NRep, 1), ones(EstimOpt.NVarstr,1)), [2 3 4 1]);
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)- grad_tmp(:,:, ones(EstimOpt.NVarstr,1), ones(sum(EstimOpt.MeaMatrix(:,i)'== 1),1)).*LV_der(:,:,:,LVindx).*bx;
            else
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)- grad_tmp(:,:, ones(EstimOpt.NVarstr,1)).*LV_der(:,:,:,LVindx)*b(1);
            end
            
            gmea(:,:,l+1:l+sum(EstimOpt.MeaMatrix(:,i)'== 1)) = -grad_tmp(:,:, ones(sum(EstimOpt.MeaMatrix(:,i)'== 1),1)).*LV_expand(:,:,EstimOpt.MeaMatrix(:,i)'== 1); % this is for LV parameters
            if EstimOpt.MeaExpMatrix(i) ~= 0
                gmea(:,:,l+1+sum(EstimOpt.MeaMatrix(:,i)'== 1):l+sum(EstimOpt.MeaMatrix(:,i)'== 1)+tmp) = ...
                    -grad_tmp(:,:,ones(EstimOpt.NVarmea_exp,1)).*Xmea_exp_expand ; % additional explanatory variables of Xmea
            end
            gmea(:,:,l+1+sum(EstimOpt.MeaMatrix(:,i)'== 1)+tmp) = grad_tmp; %first threshold level
            for j = 2:k
                gmea(X_mea(:,i) > UniqueMea(j),:,l+j+sum(EstimOpt.MeaMatrix(:,i)'== 1)+tmp) = grad_tmp(X_mea(:,i) > UniqueMea(j),:)*exp(b(sum(EstimOpt.MeaMatrix(:,i))+tmp+j)); %other thresholds levels
                gmea(X_mea(:,i) == UniqueMea(j),:,l+j+sum(EstimOpt.MeaMatrix(:,i)'== 1)+tmp) = (normpdf(alpha(j)-Xb(X_mea(:,i) == UniqueMea(j),:))./max(L(X_mea(:,i) == UniqueMea(j),:),realmin))*exp(b(sum(EstimOpt.MeaMatrix(:,i))+tmp+j)); %other thresholds levels
            end           
            L_mea = L_mea.*L; ...
            l = l+k+size(X,2); ...   
        elseif EstimOpt.MeaSpecMatrix(i) == 3 % POISS
            if EstimOpt.MeaExpMatrix(i) == 0
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)']; ...
            else
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)', Xmea_exp]; ...
            end
            b = bmea(l+1:l+size(X,2)); ...
            fit = reshape(X*b, EstimOpt.NRep, EstimOpt.NP)';
            lam = exp(fit);
            L = exp(fit.*X_mea(:,i*ones(EstimOpt.NRep,1))-lam);
            L = L./min(gamma(X_mea(:,i*ones(EstimOpt.NRep,1))+1),realmax);
            
            L_mea = L_mea.*L; ...
            grad_tmp = X_mea(:,i*ones(EstimOpt.NRep,1)) - lam;
            LVindx = find(EstimOpt.MeaMatrix(:,i)'== 1);
            % gradient for structural equation
            if sum(EstimOpt.MeaMatrix(:,i)'== 1) > 1
                bx = b(2:1+sum(EstimOpt.MeaMatrix(:,i)'== 1)); % parameters for LV
                bx = permute(bx(:, ones(EstimOpt.NP, 1), ones(EstimOpt.NRep, 1), ones(EstimOpt.NVarstr,1)), [2 3 4 1]);
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp(:,:, ones(EstimOpt.NVarstr,1), ones(sum(EstimOpt.MeaMatrix(:,i)'== 1),1)).*LV_der(:,:,:,LVindx).*bx;
            else
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp(:,:, ones(EstimOpt.NVarstr,1)).*LV_der(:,:,:,LVindx)*b(2);
            end
            
            gmea(:,:,l+1) = grad_tmp; % constant
            gmea(:,:,l+2:l+1+sum(EstimOpt.MeaMatrix(:,i)'== 1)) = grad_tmp(:,:, ones(sum(EstimOpt.MeaMatrix(:,i)'== 1),1)).*LV_expand(:,:,EstimOpt.MeaMatrix(:,i)'== 1); % parameters for LV
            if EstimOpt.MeaExpMatrix(i) ~= 0
                gmea(:,:,l+2+sum(EstimOpt.MeaMatrix(:,i)'== 1):l+1+sum(EstimOpt.MeaMatrix(:,i)'== 1)+EstimOpt.NVarmea_exp) = ...
                    grad_tmp(:,:,ones(EstimOpt.NVarmea_exp,1)).*Xmea_exp_expand ;
            end
            l = l+size(X,2); ...
                
        elseif EstimOpt.MeaSpecMatrix(i) == 4 % NB
            if EstimOpt.MeaExpMatrix(i) == 0
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)']; ...
            else
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)', Xmea_exp]; ...
            end
            b = bmea(l+1:l+size(X,2)); ...
            fit = reshape(X*b, EstimOpt.NRep, EstimOpt.NP)';
            lam = exp(fit);
            theta = exp(bmea(l+size(X,2)+1));
            u = theta./(theta+lam);  
            L = min(gamma(theta+X_mea(:,i*ones(EstimOpt.NRep,1))), realmax)./(gamma(theta).*min(gamma(X_mea(:,i*ones(EstimOpt.NRep,1))+1),realmax));
            L = L.*(u.^theta).*((1-u).^X_mea(:,i*ones(EstimOpt.NRep,1)));
            L_mea = L_mea.*L; ...
            % Calculations for gradient 
                    
            grad_tmp = u.*(X_mea(:,i*ones(EstimOpt.NRep,1))-  lam);
            LVindx = find(EstimOpt.MeaMatrix(:,i)'== 1);
            % gradient for structural equation
            if sum(EstimOpt.MeaMatrix(:,i)'== 1) > 1
                bx = b(2:1+sum(EstimOpt.MeaMatrix(:,i)'== 1)); % parameters for LV
                bx = permute(bx(:, ones(EstimOpt.NP, 1), ones(EstimOpt.NRep, 1), ones(EstimOpt.NVarstr,1)), [2 3 4 1]);
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp(:,:, ones(EstimOpt.NVarstr,1), ones(sum(EstimOpt.MeaMatrix(:,i)'== 1),1)).*LV_der(:,:,:,LVindx).*bx;
            else
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp(:,:, ones(EstimOpt.NVarstr,1)).*LV_der(:,:,:,LVindx)*b(2);
            end
            
            gmea(:,:,l+1) = grad_tmp; % constant
            gmea(:,:,l+2:l+1+sum(EstimOpt.MeaMatrix(:,i)'== 1)) = grad_tmp(:,:, ones(sum(EstimOpt.MeaMatrix(:,i)'== 1),1)).*LV_expand(:,:,EstimOpt.MeaMatrix(:,i)'== 1); % parameters for LV
            
            if EstimOpt.MeaExpMatrix(i) == 0
                gmea(:,:,l+2+sum(EstimOpt.MeaMatrix(:,i)'== 1)) = psi(theta+X_mea(:,i*ones(EstimOpt.NRep,1)))*theta - psi(theta)*theta+theta*log(u)...
                    + theta*(1-u) - X_mea(:,i*ones(EstimOpt.NRep,1)).*u; % theta
            else
                gmea(:,:,l+2+sum(EstimOpt.MeaMatrix(:,i)'== 1):l+1+sum(EstimOpt.MeaMatrix(:,i)'== 1)+EstimOpt.NVarmea_exp) = ...
                    grad_tmp(:,:,ones(EstimOpt.NVarmea_exp,1)).*Xmea_exp_expand ;
                gmea(:,:,l+2+sum(EstimOpt.MeaMatrix(:,i)'== 1)+EstimOpt.NVarmea_exp) = psi(theta+X_mea(:,i*ones(EstimOpt.NRep,1)))*theta - psi(theta)*theta+theta*log(u)...
                    + theta*(1-u) - X_mea(:,i*ones(EstimOpt.NRep,1)).*u; % theta
            end
            l = l+size(X,2)+1; ...
        elseif EstimOpt.MeaSpecMatrix(i) == 5 % ZIP
            if EstimOpt.MeaExpMatrix(i) == 0
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)']; ...
            else
                X = [ones(EstimOpt.NRep*EstimOpt.NP,1), LV(EstimOpt.MeaMatrix(:,i)'== 1,:)', Xmea_exp]; ...
            end
            bzip = bmea(l+1:l+size(X,2)); ...
            bpoiss = bmea(l+size(X,2)+1:l+2*size(X,2));
            fit = reshape(X*bpoiss, EstimOpt.NRep, EstimOpt.NP)';
            p = reshape(exp(X*bzip), EstimOpt.NRep, EstimOpt.NP)';
            p = p./(1+p);
            L = zeros(EstimOpt.NP, EstimOpt.NRep);
            lam = exp(fit);
            IndxZIP = X_mea(:,i) == 0;
            L(IndxZIP,:) = p(IndxZIP,:) + (1-p(IndxZIP,:)).*exp(-lam(IndxZIP,:));
            L(~IndxZIP,:) = (1-p(~IndxZIP,:)).*exp(fit(~IndxZIP,:).*X_mea(~IndxZIP,i*ones(EstimOpt.NRep,1))-lam(~IndxZIP,:))./min(gamma(X_mea(~IndxZIP,i*ones(EstimOpt.NRep,1))+1),realmax);
            L_mea = L_mea.*L; ...

            % Calculations for gradient 
            grad_tmp1 = zeros(EstimOpt.NP, EstimOpt.NRep);
            grad_tmp2 = zeros(EstimOpt.NP, EstimOpt.NRep);
            % For ZIP 
            grad_tmp1(IndxZIP,:) = (-(1-p(IndxZIP,:)).*p(IndxZIP,:).*exp(-lam(IndxZIP,:))+p(IndxZIP,:)-p(IndxZIP,:).^2)./L(IndxZIP,:);
            grad_tmp1(~IndxZIP,:) = -p(~IndxZIP,:);
            % For Poiss
            grad_tmp2(IndxZIP,:) = (p(IndxZIP,:)-1).*exp(fit(IndxZIP,:)-lam(IndxZIP,:))./L(IndxZIP,:);
            grad_tmp2(~IndxZIP,:) = X_mea(~IndxZIP,i*ones(EstimOpt.NRep,1)) - lam(~IndxZIP,:);

            LVindx = find(EstimOpt.MeaMatrix(:,i)'== 1);
            % gradient for structural equation
            if sum(EstimOpt.MeaMatrix(:,i)'== 1) > 1
                bx1 = bzip(2:1+sum(EstimOpt.MeaMatrix(:,i)'== 1)); % parameters for LV
                bx1 = permute(bx1(:, ones(EstimOpt.NP, 1), ones(EstimOpt.NRep, 1), ones(EstimOpt.NVarstr,1)), [2 3 4 1]);
                bx2 = bpoiss(2:1+sum(EstimOpt.MeaMatrix(:,i)'== 1)); % parameters for LV
                bx2 = permute(bx2(:, ones(EstimOpt.NP, 1), ones(EstimOpt.NRep, 1), ones(EstimOpt.NVarstr,1)), [2 3 4 1]);
                
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp1(:,:, ones(EstimOpt.NVarstr,1), ones(sum(EstimOpt.MeaMatrix(:,i)'== 1),1)).*LV_der(:,:,:,LVindx).*bx1+...
                    grad_tmp2(:,:, ones(EstimOpt.NVarstr,1), ones(sum(EstimOpt.MeaMatrix(:,i)'== 1),1)).*LV_der(:,:,:,LVindx).*bx2;
            else
                gstr(:,:,:,LVindx) = gstr(:,:,:,LVindx)+ grad_tmp1(:,:, ones(EstimOpt.NVarstr,1)).*LV_der(:,:,:,LVindx)*bzip(2)+...
                    grad_tmp2(:,:, ones(EstimOpt.NVarstr,1)).*LV_der(:,:,:,LVindx)*bpoiss(2);
            end
            gmea(:,:,l+1) = grad_tmp1; % constant
            gmea(:,:,l+2:l+1+sum(EstimOpt.MeaMatrix(:,i)'== 1)) = grad_tmp1(:,:, ones(sum(EstimOpt.MeaMatrix(:,i)'== 1),1)).*LV_expand(:,:,EstimOpt.MeaMatrix(:,i)'== 1); % parameters for LV
            gmea(:,:,l+sum(EstimOpt.MeaMatrix(:,i)'== 1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp+2) = grad_tmp2; % constant
            gmea(:,:,l+sum(EstimOpt.MeaMatrix(:,i)'== 1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp+3:l+2*sum(EstimOpt.MeaMatrix(:,i)'== 1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp+2) = ...
                grad_tmp2(:,:, ones(sum(EstimOpt.MeaMatrix(:,i)'== 1),1)).*LV_expand(:,:,EstimOpt.MeaMatrix(:,i)'== 1); % parameters for LV
            
            if EstimOpt.MeaExpMatrix(i) ~= 0
                gmea(:,:,l+2+sum(EstimOpt.MeaMatrix(:,i)'== 1):l+1+sum(EstimOpt.MeaMatrix(:,i)'== 1)+EstimOpt.NVarmea_exp) = ...
                    grad_tmp1(:,:,ones(EstimOpt.NVarmea_exp,1)).*Xmea_exp_expand ;
                gmea(:,:,l+2*sum(EstimOpt.MeaMatrix(:,i)'== 1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp+3:l+2*sum(EstimOpt.MeaMatrix(:,i)'== 1)+2*EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp+2) = ...
                    grad_tmp2(:,:,ones(EstimOpt.NVarmea_exp,1)).*Xmea_exp_expand ;
            end
            l = l+2*size(X,2); ...
        end
    end
    
    probs = probs.*L_mea;
    p = max(realmin,mean(probs,2));
    f = -log(p);
    

    gstr = reshape(gstr, EstimOpt.NP, EstimOpt.NRep, EstimOpt.NVarstr*EstimOpt.NLatent);
    gmnl = reshape(gmnl, EstimOpt.NP, EstimOpt.NRep, EstimOpt.NVarA*(1+EstimOpt.NLatent));
%     size(probs(:,:, ones(EstimOpt.NVarA*(1+EstimOpt.NLatent),1)))
%     size(gmnl)
    g = squeeze(mean(probs(:,:, ones(EstimOpt.NVarA*(1+EstimOpt.NLatent),1)).*gmnl,2)); % 
    g3 = squeeze(mean(probs(:,:, ones(length(bmea),1)).*gmea,2)); % NP x NVarmea
    g2 = squeeze(mean(probs(:,:, ones(EstimOpt.NVarstr*EstimOpt.NLatent,1)).*gstr,2)); % NP x NVarstr*NLatent
    
    g = [g,g2,g3];
    g = -g./p(:, ones(1, length(B)));
end
% f = -log(mean(probs.*L_mea,2));

