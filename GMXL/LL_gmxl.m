function [f,g] = LL_gmxl(YY,XXa,XXm,XXs,XXt,err,EstimOpt,b0)

% save tmp_LL_gmxl
% return

NAlt = EstimOpt.NAlt;
NCT = EstimOpt.NCT;
NP = EstimOpt.NP;
NRep = EstimOpt.NRep;
NVarA = EstimOpt.NVarA;
NVarM = EstimOpt.NVarM;
NVarS = EstimOpt.NVarS;
NVarT = EstimOpt.NVarT;
Dist = EstimOpt.Dist;
DistS = EstimOpt.DistS;
WTP_space = EstimOpt.WTP_space;
WTP_matrix = EstimOpt.WTP_matrix;
NCTMiss = EstimOpt.NCTMiss;
NAltMiss = EstimOpt.NAltMiss;
FullCov = EstimOpt.FullCov;
Gamma0 = EstimOpt.Gamma0;
indx3 = EstimOpt.indx3;
if FullCov == 1 && nargout > 1
    indx1 = EstimOpt.indx1;
    indx2 = EstimOpt.indx2;
else
    indx1 = [];
    indx2 = [];
end
RealMin = EstimOpt.RealMin;

b0a = b0(1:NVarA);

if FullCov == 0
    b0v = (b0(NVarA+1:NVarA*2));
    %     b0v = b0(NVarA+1:NVarA*2);
    VC = diag(b0v);
    b0m = b0(NVarA*2+1:NVarA*(NVarM+2));
    b0m = reshape(b0m,NVarA, NVarM);
    b0s = b0(NVarA*(NVarM+2)+1:NVarA*(NVarM+2)+NVarS);
    b0t = b0(NVarA*(NVarM+2)+NVarS+1:NVarA*(NVarM+2)+NVarS+NVarT);
    tau = exp((b0(NVarA*(NVarM+2)+NVarS+NVarT+1)));
    %     tau = exp(b0(NVarA*(NVarM+2)+NVarS+NVarT+1));
    if Gamma0 == 0 || Gamma0 == 1
        gamma = Gamma0;
    else
        gamma = exp(b0(NVarA*(NVarM+2)+NVarS+NVarT+2)) ./ (1 + exp(b0(NVarA*(NVarM+2)+NVarS+NVarT+2)));
    end
else
    b0v = b0(NVarA+1:NVarA+sum(1:NVarA));
    VC = tril(ones(NVarA));
    VC(VC==1) = b0v;
    b0m = b0(NVarA*(NVarA/2+1.5)+1:NVarA*(NVarA/2+1.5+NVarM));
    b0m = reshape(b0m,NVarA, NVarM);
    b0s = b0(NVarA*(NVarA/2+1.5+NVarM)+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS);
    b0t = b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarT);
    tau = exp(b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarT+1));
    if Gamma0 == 0 || Gamma0 == 1
        gamma = Gamma0;
    else
        gamma = exp(b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarT+2)) ./ (1 + exp(b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarT+2)));
    end
end
b0n = b0a + b0m*XXm;
b0n = reshape(b0n([1:size(b0n,1)]'*ones(1,NRep),[1:size(b0n,2)]'),NVarA,NRep*NP);
if size(XXt,2) == 0
    if DistS == 1
        errx = reshape(err(end,:), NRep, NP)'; % NP x NRep
        sigma_bar = - log(mean(exp(tau*errx),1)); %1 x NRep
        %sigma_bar = -log(mean(exp(tau*(squeeze(err_sliced(1,:,:))')),1));
        sigma = exp(sigma_bar + tau*errx); % NP x NRep
    elseif DistS == 5
        shape_n = 1./tau;
        errx = reshape(err(end,:), NRep, NP)';
        sigma_bar = mean(errx.^shape_n,1); %1 x NRep
        sigma = (errx.^shape_n) ./ sigma_bar;
        %sigma_bar = mean(squeeze(err_sliced(1,:,:)).^shape_n(ones(EstimOpt.NRep,1),ones(1,EstimOpt.NP)),2)';
    end
else
    if DistS == 1
        cov_tau = exp(XXt*b0t);  % NP x 1
        %         cov_tau = cov_tau(:,ones(1,NRep));
        errx = reshape(err(end,:), NRep, NP)';
        sigma_bar = -log(mean(exp(tau*cov_tau.*errx),1)); % 1 x NRep
        %sigma_bar = -log(mean(exp(tau*cov_tau(:,ones(1,EstimOpt.NRep)).*(squeeze(err_sliced(1,:,:))')),1));
        sigma = exp(sigma_bar + tau*cov_tau.*errx);
    elseif DistS == 5
        shape_n = 1./(tau*exp(XXt*b0t));  % NP x 1
        %         shape_n = shape_n(:,  ones(1, NRep));
        errx = reshape(err(end,:), NRep, NP)';
        sigma_bar = mean(errx.^shape_n,1); % 1 x NRep
        %sigma_bar = mean(squeeze(err_sliced(1,:,:)).^shape_n(ones(EstimOpt.NRep,1),:),2)';
        sigma = (errx.^shape_n) ./ sigma_bar; % NP x NRep
    end
end

sigma = reshape(sigma',[1,NRep*NP]);

if sum(Dist > 0) == 0  % Normal
    if Gamma0 == 0
        b_mtx = sigma .* (b0n + VC*err(1:end-1,:));  % NVarA x NRep*NP
    elseif Gamma0 == 1
        b_mtx = sigma .* b0n + VC*err(1:end-1,:);
    else
        b_mtx = sigma .* (b0n + (1-gamma)*VC*err(1:end-1,:)) + gamma*VC*err(1:end-1,:);
    end
    if nargout == 2
        b_mtx_grad = zeros(NVarA,0, NP);
    end
elseif sum(Dist==1) > 0  % Log - normal
    if Gamma0 == 0
        b_mtx  = b0n + VC*err(1:end-1,:);
        b_mtx(Dist==1,:) = exp(b_mtx(Dist==1,:));
        if nargout == 2
            b_mtx_grad = reshape(b_mtx,[NVarA,NRep,NP]);
        else
            b_mtx_grad = zeros(NVarA,0,NP);
        end
        b_mtx = sigma .* b_mtx;  % NVarA x NRep*NP
    elseif Gamma0 == 1
        b_mtx = sigma .* b0n + VC*err(1:end-1,:);
        b_mtx(Dist==1,:) = exp(b_mtx(Dist==1,:));
    else
        b_mtx = sigma .* (b0n + (1-gamma)*VC*err(1:end-1,:)) + gamma*VC*err(1:end-1,:);
        b_mtx(Dist==1,:) = exp(b_mtx(Dist==1,:));
    end
elseif sum(Dist==2) > 0  % Spike
    if Gamma0 == 0
        b_mtx  = b0n + VC*err(1:end-1,:);
        b_mtx(Dist==1,:) = max(b_mtx(Dist==1,:),0);
        b_mtx = sigma .* b_mtx;  % NVarA x NRep*NP
    elseif Gamma0 == 1
        b_mtx = sigma .* b0n + VC*err(1:end-1,:);
        b_mtx(Dist==1,:) = max(b_mtx(Dist==1,:),0);
    else
        b_mtx = sigma .* (b0n + (1-gamma)*VC*err(1:end-1,:)) + gamma*VC*err(1:end-1,:);
        b_mtx(Dist==1,:) = max(b_mtx(Dist==1,:),0);
    end
    if nargout == 2
        b_mtx_grad = zeros(NVarA,0, NP);
    end
elseif sum(Dist==5) > 0
    if sum(sum(VC.*(1-eye(size(b0n,1)))~=0))~=0; error ('Weibull distribution can only be used with non-correlated parameters'); end
    if gma ~= 0
        error ('Weibull distributed attriute parameters possible only with G-MNL Type II');
    else
        b_mtx = zeros(NVarA,NP*NRep);
        err2 = err(1:end-1,:);
        b_mtx(Dist==0,:) = b0n(Dist==0,:) + VC(Dist==0,Dist==0)*err2(Dist==0,:);
        Wexp = (1./diag(VC(Dist==5,Dist==5)));
        b_mtx(Dist==5,:) = b0n(Dist==5).*err2(Dist==5,:).^Wexp;
        b_mtx = b_mtx.*sigma;
    end
    if nargout == 2
        b_mtx_grad = zeros(NVarA,0, NP);
    end
end

if WTP_space > 0 % this is wrong
    b_mtx(1:end-WTP_space,:) = b_mtx(1:end-WTP_space,:).*b_mtx(WTP_matrix,:);
end

cs = reshape(exp(XXs*b0s),NAlt*NCT,1,NP);
XXa = XXa .* cs;

b_mtx = reshape(b_mtx,NVarA,NRep,NP);

p0 = zeros(NP,1);

if nargout == 1 % function value only
    if any(isnan(XXa(:))) == 0 % faster version for complete dataset
        YYy = YY==1;
        parfor n = 1:NP
            U = reshape(XXa(:,:,n)*b_mtx(:,:,n),NAlt,NCT,NRep);
            U_max = max(U);
            U = exp(U - U_max); % rescale utility to avoid exploding
            U_sum = reshape(sum(U,1),NCT,NRep);
            YYy_n = YYy(:,n);
            U_selected = reshape(U(YYy_n(:,ones(NRep,1))),NCT,NRep);
            p0(n) = mean(prod(U_selected ./ U_sum,1),2);
        end
    else  % this works only if NAlt is constant for each respondent and if missing ALT is not the first in NCT
        parfor n = 1:NP
            YnanInd = ~isnan(YY(:,n));
            XXa_n = XXa(:,:,n);
            U = reshape(XXa_n(YnanInd,:,:)*b_mtx(:,:,n),NAltMiss(n),NCTMiss(n),NRep);
            U_max = max(U);
            U = exp(U - U_max);
            U_sum = reshape(sum(U,1),NCTMiss(n),NRep);
            YYy_n = YY(:,n)==1;
            U_selected = reshape(U(YYy_n(YnanInd,ones(NRep,1))),[NCTMiss(n),NRep]);
            p0(n) = mean(prod(U_selected ./ U_sum,1),2);
        end
    end
else
    if NVarS > 0
        Xs_sliced = reshape(XXs, [NAlt*NCT,NP,NVarS]);
    else
        Xs_sliced = reshape(XXs, [NAlt*NCT,NP,0]);
    end
    if FullCov == 0
        g = zeros(NP,2*NVarA + NVarS+NVarT + 1 + (Gamma0 ~= 0 && Gamma0 ~= 1));
        %         VC2 = reshape(2*diag(b0(NVarA+1:NVarA*2))*err, NVarA, NRep, NP);
    else
        g = zeros(NP, 2*NVarA+NVarA*(NVarA-1)/2 +  NVarS + NVarT+1 +(Gamma0 ~= 0 && Gamma0 ~= 1));
    end
    VC2 = reshape(err(1:end-1,:),[NVarA,NRep,NP]);
    Eta = reshape(VC*err(1:end-1,:),[NVarA,NRep,NP]);
    b0n = reshape(b0n,[NVarA,NRep,NP]);
    sigma = reshape(sigma,[1,NRep,NP]);
    if NVarT == 0
        mSig = -mean(exp(tau*errx).*(tau*errx),1)./mean(exp(tau*errx),1); %1 x NRep
        cov_tau = zeros(NP,NRep,0);
        mSigCov = zeros(NVarT,NRep,0);
    else
        mSig = -mean(exp(tau*errx.*cov_tau).*(tau*errx.*cov_tau),1)./mean(exp(tau*errx.*cov_tau),1); %1 x NRep
        mSigTmp = reshape(exp(tau*errx.*cov_tau).*(tau*errx.*cov_tau),1, NP, NRep);
        mSigTmp = bsxfun(@times, mSigTmp, XXt'); % NVarT x NP x NRep
        mSigCov = -reshape(mean(mSigTmp,2), NVarT, NRep);
        mSigTmp = mean(exp(tau*errx.*cov_tau),1);
        mSigCov = mSigCov ./ mSigTmp; % NVarT x NRep
    end
    %Distx = Dist(2:end);
    if any(isnan(XXa(:))) == 0  % faster version for complete dataset
        YYy = (YY==1);
        parfor n = 1:NP
            b_mtx_n = b_mtx(:,:,n);
            XXa_n = XXa(:,:,n);
            U = reshape(XXa_n*b_mtx_n,NAlt,NCT,NRep);  % NAlt x NCT x NRep
            U = exp(bsxfun(@minus,U,max(U)));  % rescale utility to avoid exploding
            U_prob = bsxfun(@rdivide,U,sum(U,1)); % NAlt x NCT x NRep
            YYy_n = YYy(:,n);
            U_prod = prod(reshape(U_prob(YYy_n(:,ones(NRep,1))),NCT,NRep),1);  % 1 x NRep
            p0(n) = mean(U_prod);
            
            % Calculations for gradient
            U_prob = reshape(U_prob, NAlt*NCT,1,NRep);  % NAlt*NCT x NVarA x NRep
            X_hat = sum(reshape(bsxfun(@times,U_prob,XXa_n), NAlt, NCT, NVarA, NRep),1);
            if NCT ~= 1
                F = bsxfun(@minus,XXa_n(YYy_n,:), reshape(X_hat,[NCT,NVarA,NRep]));  %NCT x NVarA x NRep
                sumFsqueezed = reshape(sum(F,1),[NVarA,NRep]);  %NVarA x NRep
            else
                sumFsqueezed = reshape(bsxfun(@minus,XXa_n(YYy_n,:,n),squeeze(X_hat)),[NCT,NVarA,NRep]); %NVarA x NRep
            end
            b_mtx_grad_n = b_mtx_grad(:,:,n)
            if Gamma0 == 0 && sum(Dist==1) > 0
                sumFsqueezed2 = sumFsqueezed;
                sumFsqueezed(Dist==1, :) = sumFsqueezed(Dist==1, :).*b_mtx_grad_n(Dist==1,:);
            else
                sumFsqueezed2 = sumFsqueezed;
            end
            sigma_n = sigma(:,:,n);
            errx_n = errx(n,:);
            mSig_n = mSig;
            if Gamma0 == 0 && sum(Dist==1)
                DerTau = b_mtx_grad_n;
            else
                DerTau = bsxfun(@plus,b0n(:,:,n),(1-gamma)*Eta(:,:,n));
            end
            DerTau = bsxfun(@times,DerTau,sigma_n);
            if NVarT == 0
                DerTau = bsxfun(@times,bsxfun(@plus,tau*errx_n,mSig_n), DerTau);
            else
                TauTmp = bsxfun(@times, tau*errx_n, cov_tau(n,:)); %1 x Nrep
                DerCovTau = bsxfun(@plus,bsxfun(@times,TauTmp,XXt(n,:)'),mSigCov)
                DerCovTau = bsxfun(@times, reshape(DerCovTau, 1,NVarT, NRep), reshape(DerTau, NVarA,1, NRep)); %  NVarA x NVarT x NRep
                DerCovTau = reshape(DerCovTau, NVarA*NVarT, NRep);
                DerTau = bsxfun(@times,bsxfun(@plus,TauTmp,mSig_n), DerTau);
            end
            if Gamma0 ~= 0 && Gamma0 ~= 1% gradient for gamma
                DerGam = bsxfun(@times,gamma*(1-gamma)*Eta(:,:,n),1-sigma_n);
            end
            
            if NVarS >0
                FScale = sum(sumFsqueezed2.*b_mtx_n,1); % 1 x NRep
                Xs_tmp = squeeze(Xs_sliced(1,n,:));
                FScale = FScale .* Xs_tmp; % NVarS x NRep
            end
            
            sumBeta = bsxfun(@times,sumFsqueezed, sigma_n);
            VC2tmp = (1-gamma)*bsxfun(@times,VC2(:,:,n), sigma_n) + gamma*VC2(:,:,n);
            if FullCov == 0
                sumVC2tmp = sumFsqueezed.*VC2tmp;  % NVarA x NRep
                gtmp = -mean([sumBeta.*U_prod; sumVC2tmp.*U_prod],2)./p0(n);
            else % FullCov = 1
                sumVC2tmp = sumFsqueezed(indx1,:).*VC2tmp(indx2,:);
                gtmp =  -mean([sumBeta.*U_prod; sumVC2tmp.*U_prod],2)./p0(n);
            end
            if NVarS > 0
                gtmp = [gtmp;-mean(FScale.*U_prod,2)./p0(n)];
            end
            
            sumTau = sum(sumFsqueezed2.*DerTau,1); %1 x NRep
            if Gamma0 ~= 0 && Gamma0 ~= 1% gradient for gamma
                sumGam = sum(sumFsqueezed.*DerGam,1);
                gtmp2 = -mean([sumTau.*U_prod; sumGam.*U_prod],2)./p0(n);
            else
                gtmp2 = -mean(sumTau.*U_prod,2)./p0(n);
            end
            if NVarT > 0
                sumCovTau = bsxfun(@times,sumFsqueezed2(indx3,:), DerCovTau);
                sumCovTau = reshape(sum(reshape(sumCovTau, NVarA, NVarT, NRep),1),NVarT, NRep);
                gtmp2 = [-mean(bsxfun(@times,sumCovTau, U_prod),2)./p0(n);gtmp2];
            end
            gtmp = [gtmp; gtmp2];
            g(n,:) = gtmp';
        end
    else
        YYy = (YY==1);
        parfor n = 1:NP
            YnanInd = ~isnan(YY(:,n));
            b_mtx_n = b_mtx(:,:,n);
            XXa_n = XXa(:,:,n);
            U = reshape(XXa_n(YnanInd,:)*b_mtx_n,NAltMiss(n),NCTMiss(n),NRep);  % NAlt x NCT x NRep
            U = exp(bsxfun(@minus,U,max(U)));  % rescale utility to avoid exploding
            U_prob = bsxfun(@rdivide,U,sum(U,1)); % NAlt x NCT x NRep
            YYy_n = YYy(:,n);
            U_prod = prod(reshape(U_prob(YYy_n(YnanInd,ones(NRep,1))),NCTMiss(n),NRep),1);  % 1 x NRep
            p0(n) = mean(U_prod);
            
            % Calculations for gradient
            U_prob = reshape(U_prob, NAltMiss(n)*NCTMiss(n),1,NRep);  % NAlt*NCT x NVarA x NRep
            X_hat = sum(reshape(bsxfun(@times,U_prob,XXa_n(YnanInd,:)), NAltMiss(n), NCTMiss(n), NVarA, NRep),1);
            %             if NCTMiss(n) ~= 1
            F = bsxfun(@minus,XXa_n(YYy_n,:), reshape(X_hat,[NCTMiss(n),NVarA,NRep]));  %NCT x NVarA x NRep
            sumFsqueezed = reshape(sum(F,1),[NVarA,NRep]);  %NVarA x NRep
            %             else
            %                 sumFsqueezed = reshape(bsxfun(@minus,XXa_n(YYy_n,:,n),squeeze(X_hat)),[NCTMiss(n),NVarA,NRep]); %NVarA x NRep
            %                 sumFsqueezed = reshape(bsxfun(@minus,XXa_n(YYy_n,:),reshape(X_hat,[NCTMiss(n),NVarA,NRep])),[NCTMiss(n),NVarA,NRep]); %NVarA x NRep
            %                 sumFsqueezed = reshape(sum(reshape(bsxfun(@minus,XXa_n(YYy_n,:),reshape(X_hat,[NCTMiss(n),NVarA,NRep])),[NCTMiss(n),NVarA,NRep]),1),[NVarA,NRep]);  %NVarA x NRep
            %             end
            b_mtx_grad_n = b_mtx_grad(:,:,n);
            if Gamma0 == 0 && sum(Dist==1) > 0
                sumFsqueezed2 = sumFsqueezed;
                sumFsqueezed(Dist==1, :) = sumFsqueezed(Dist==1, :).*b_mtx_grad_n(Dist==1,:);
            else
                sumFsqueezed2 = sumFsqueezed;
            end
            sigma_n = sigma(:,:,n);
            errx_n = errx(n,:);
            mSig_n = mSig;
            if Gamma0 == 0 && sum(Dist==1)
                DerTau = b_mtx_grad_n;
            else
                DerTau = bsxfun(@plus,b0n(:,:,n),(1-gamma)*Eta(:,:,n));
            end
            DerTau = bsxfun(@times,DerTau,sigma_n);
            if NVarT == 0
                DerTau = bsxfun(@times,bsxfun(@plus,tau*errx_n,mSig_n), DerTau);
            else
                TauTmp = bsxfun(@times, tau*errx_n, cov_tau(n,:)); %1 x Nrep
                DerCovTau = bsxfun(@plus,bsxfun(@times,TauTmp,XXt(n,:)'),mSigCov);
                DerCovTau = bsxfun(@times, reshape(DerCovTau, 1,NVarT, NRep), reshape(DerTau, NVarA,1, NRep)); %  NVarA x NVarT x NRep
                DerCovTau = reshape(DerCovTau, NVarA*NVarT, NRep);
                DerTau = bsxfun(@times,bsxfun(@plus,TauTmp,mSig_n), DerTau);
            end
            if Gamma0 ~= 0 && Gamma0 ~= 1% gradient for gamma
                DerGam = bsxfun(@times,gamma*(1-gamma)*Eta(:,:,n),1-sigma_n);
            end
            
            if NVarS >0
                FScale = sum(sumFsqueezed2.*b_mtx_n,1); % 1 x NRep
                Xs_tmp = squeeze(Xs_sliced(1,n,:));
                FScale = FScale .* Xs_tmp; % NVarS x NRep
            end
            
            sumBeta = bsxfun(@times,sumFsqueezed, sigma_n);
            VC2tmp = (1-gamma)*bsxfun(@times,VC2(:,:,n), sigma_n) + gamma*VC2(:,:,n);
            if FullCov == 0
                sumVC2tmp = sumFsqueezed.*VC2tmp;  % NVarA x NRep
                gtmp = -mean([sumBeta.*U_prod; sumVC2tmp.*U_prod],2)./p0(n);
            else % FullCov = 1
                sumVC2tmp = sumFsqueezed(indx1,:).*VC2tmp(indx2,:);
                gtmp =  -mean([sumBeta.*U_prod; sumVC2tmp.*U_prod],2)./p0(n);
            end
            if NVarS > 0
                gtmp = [gtmp;-mean(FScale.*U_prod,2)./p0(n)];
            end
            
            sumTau = sum(sumFsqueezed2.*DerTau,1); %1 x NRep
            if Gamma0 ~= 0 && Gamma0 ~= 1% gradient for gamma
                sumGam = sum(sumFsqueezed.*DerGam,1);
                gtmp2 = -mean([sumTau.*U_prod; sumGam.*U_prod],2)./p0(n);
            else
                gtmp2 = -mean(sumTau.*U_prod,2)./p0(n);
            end
            if NVarT > 0
                sumCovTau = bsxfun(@times,sumFsqueezed2(indx3,:), DerCovTau);
                sumCovTau = reshape(sum(reshape(sumCovTau, NVarA, NVarT, NRep),1),NVarT, NRep);
                gtmp2 = [-mean(bsxfun(@times,sumCovTau, U_prod),2)./p0(n);gtmp2];
            end
            gtmp = [gtmp; gtmp2];
            g(n,:) = gtmp';
        end
    end
    if NVarM > 0
        gm =  g(:,repmat(1:NVarA, 1, NVarM)).*(XXm(kron(1:NVarM,ones(1,NVarA)),:)');
        if EstimOpt.FullCov == 0
            g = [g(:,1:2*NVarA),gm, g(:,2*NVarA+1:end)];
        else
            g = [g(:,1:NVarA*(NVarA/2+1.5)),gm, g(:,NVarA*(NVarA/2+1.5)+1:end)];
        end
    end
end

% f = -log(p0);
% f = -log(max(p0,realmin));
if RealMin == 1
    f = -log(max(p0,realmin));
else
    f = -log(p0);
end
