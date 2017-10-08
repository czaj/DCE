function [f,g] = LL_lcmxl(YY,XXa,Xc,Xs,err,EstimOpt,B)

% save tmp_LL_lcmxl
% return

NVarA = EstimOpt.NVarA;
NClass = EstimOpt.NClass;
NVarC = EstimOpt.NVarC;
NP = EstimOpt.NP;
NAlt = EstimOpt.NAlt;
NCT = EstimOpt.NCT;
NRep = EstimOpt.NRep;
FullCov = EstimOpt.FullCov;
Dist = EstimOpt.Dist;
Dist1 = Dist(1:NVarA,1) == 1;
WTP_space = EstimOpt.WTP_space;
WTP_matrix = EstimOpt.WTP_matrix;
indx1 = EstimOpt.indx1;
indx2 = EstimOpt.indx2;
NCTMiss = EstimOpt.NCTMiss;
NAltMiss = EstimOpt.NAltMiss;
NVarS = EstimOpt.NVarS;

% b_mtx = zeros(NVarA, NP*NRep, NClass);


b_mtx = zeros(NVarA,NRep*NP,NClass);
VC = zeros(NVarA);
if EstimOpt.NumGrad == 0
    b_mtx_grad = zeros(NVarA,NRep,NClass,NP);% needed for gradient calculation in WTP_space
end

if FullCov == 0
    if NVarS > 0
        Scale = exp(Xs*reshape(B(2*NVarA*NClass+1:(2*NVarA+NVarS)*NClass),[NVarS,NClass])); % NP*NCT*NAlt x NClass
    end
    for c = 1:NClass
        b_mtx(:,:,c) = B((c-1)*NVarA+1:c*NVarA,:) + diag(B(NVarA*NClass+(c-1)*NVarA+1:NVarA*(NClass+c)))*err((c-1)*NVarA+1:c*NVarA,:);
        if sum(Dist(:,c) == 1) > 0  % lognormal
            b_mtx(Dist(:,c) == 1,:,c) = exp(b_mtx(Dist(:,c) == 1,:,c));
        elseif sum(Dist(:,c) == 2) > 0  % Spike
            b_mtx(Dist(:,c) == 2,:,c) = max(b_mtx(Dist(:,c) == 1,:,c),0);
        end
        if WTP_space > 0
            if EstimOpt.NumGrad == 0
                b_mtx_grad(:,:,c,:) = reshape(b_mtx(:,:,c),[NVarA,NRep,NP]);% needed for gradient calculation in WTP_space
            end
            b_mtx(1:end-WTP_space,:,c) = b_mtx(1:end-WTP_space,:,c).*b_mtx(WTP_matrix,:,c);
        end
    end
else
    if NVarS > 0
        Scale = exp(Xs*reshape(B(NClass*(NVarA+sum(1:NVarA))+1:NClass*(NVarA+NVarS+sum(1:NVarA))),[NVarS,NClass])); % NP*NCT*NAlt x NClass
    end
    VC_tmp = tril(ones(NVarA));
    for c = 1:NClass
        VC(VC_tmp == 1) = B(NVarA*NClass+(c-1)*sum(1:NVarA)+1:NVarA*NClass+c*sum(1:NVarA));
        b_mtx(:,:,c) = B((c-1)*NVarA+1:c*NVarA,ones(NP*NRep,1)) + VC*err((c-1)*NVarA+1:c*NVarA,:);
        if sum(Dist(:,c) == 1) > 0  % lognormal
            b_mtx(Dist(:,c) == 1,:,c) = exp(b_mtx(Dist(:,c) == 1,:,c));
        elseif sum(Dist(:,c) == 2) > 0  % Spike
            b_mtx(Dist(:,c) == 2,:,c) = max(b_mtx(Dist(:,c) == 1,:,c),0);
        end
        if WTP_space > 0
            if EstimOpt.NumGrad == 0
                b_mtx_grad(:,:,c,:) = reshape(b_mtx(:,:,c),[NVarA,NRep,NP]);% needed for gradient calculation in WTP_space
            end
            b_mtx(1:end-WTP_space,:,c) = b_mtx(1:end-WTP_space,:,c).*b_mtx(WTP_matrix,:,c);
        end
    end
end
if NVarS > 0
   Scale = reshape(Scale,[NAlt,NCT,NP,NClass]); 
else
   Scale = zeros(0,0,NP,0); 
end
if EstimOpt.NumGrad == 0
    b_mtx_grad = reshape(b_mtx_grad,[NVarA,NRep*NClass,NP]);% needed for gradient calculation in WTP_space
end

%     Dist = reshape(EstimOpt.Dist,[NVarA,1,NClass]);
%     b_mtx(Dist(:,ones(1,NP*NRep,1),:)==1) = exp(b_mtx(Dist(:,ones(1,NP*NRep,1),:)==1)); % this is slower

b_mtx = reshape(b_mtx,[NVarA,NRep,NP,NClass]);
b_mtx = permute(b_mtx,[1,2,4,3]);
b_mtx = reshape(b_mtx,[NVarA,NRep*NClass,NP]);

p = zeros(NP,NClass);

if nargout == 1 % function value only
    if any(isnan(XXa(:))) == 0  % faster version for complete dataset
        YYy = YY == 1;
        parfor n = 1:NP
            U = reshape(XXa(:,:,n)*b_mtx(:,:,n),[NAlt,NCT,NRep,NClass]);
            if NVarS > 0
               Scale_n = Scale(:,:,n,:)
               U = U.*Scale_n; 
            end
            U = exp(U-max(U));
            U_sum = reshape(sum(U,1),[NCT,NRep,NClass]);
            YYy_n = YYy(:,n);
            %             U_selected = reshape(U(YYn==1),NCT,NRep,NClass);
            U_selected = reshape(U(YYy_n(:,ones(1,NRep),ones(1,1,NClass))),[NCT,NRep,NClass]);
            p(n,:) = mean(prod(U_selected./U_sum,1),2);
        end
    else
        
        parfor n = 1:NP
            YnanInd = ~isnan(YY(:,n));
            XXa_n = XXa(:,:,n);
            %             U = reshape(XXa_n(YnanInd,:,:)*reshape(b_mtx(:,:,n,:),[NVarA,NRep*NClass]),[NAltMiss(n),NCTMiss(n),NRep,NClass]);
            U = reshape(XXa_n(YnanInd,:,:)*b_mtx(:,:,n),[NAltMiss(n),NCTMiss(n),NRep,NClass]);
            if NVarS > 0
               Scale_n = reshape(Scale(:,:,n,:),[NAlt*NCT,1,NClass]);
               Scale_n = reshape(Scale_n(YnanInd,:,:),[NAltMiss(n),NCTMiss(n),1,NClass]);
               U = U.*Scale_n; 
            end
            U = exp(U-max(U));
            U_sum = reshape(sum(U,1),[NCTMiss(n),NRep,NClass]);
            YYy_n = YY(:,n) == 1;
            U_selected = reshape(U(YYy_n(YnanInd,ones(1,NRep),ones(1,1,NClass))),[NCTMiss(n),NRep,NClass]);
            p(n,:) = mean(prod(U_selected./U_sum,1),2);
        end
    end
    
    if FullCov == 0
        Pclass = exp(Xc*reshape([B(NClass*(2*NVarA+NVarS)+1:end);zeros(NVarC,1)],[NVarC,NClass]));
    else
        Pclass = exp(Xc*reshape([B(NClass*(NVarA+NVarS+sum(1:NVarA))+1:end);zeros(NVarC,1)],[NVarC,NClass]));
    end
    Pclass = Pclass./sum(Pclass,2); % NP x NClass
    
    f = -log(sum(p.*Pclass,2));
    
elseif nargout == 2 % function value + gradient
    
    if FullCov == 0
        g = zeros(NP,(2*NVarA+NVarS)*NClass);
        VC2 = zeros(NVarA,NClass,NRep,NP);
        for c = 1:NClass
            VC2(:,c,:,:) = reshape(err((c-1)*NVarA+1:c*NVarA,:),[NVarA,1,NRep,NP]);
        end
        VC2 = reshape(VC2,[NVarA*NClass,NRep,NP]);
        VC2f = zeros(0,0,NP);
    else
        g = zeros(NP,(2*NVarA+NVarA*(NVarA-1)/2+NVarS)*NClass);
        VC2 = zeros(0,0,NP);
        VC2f = reshape(err,[NVarA*NClass,NRep,NP]);
    end
    if NVarS > 0
       Xs_sliced = permute(reshape(Xs(1:NAlt:end,:),[NCT,NP,NVarS]),[1 3 2]);
    else
        Xs_sliced = zeros(0,0,NP);
    end
    if FullCov == 0
        Pclass = exp(Xc*reshape([B(NClass*(2*NVarA+NVarS)+1:end);zeros(NVarC,1)],[NVarC,NClass]));
    else
        Pclass = exp(Xc*reshape([B(NClass*(NVarA+NVarS+sum(1:NVarA))+1:end);zeros(NVarC,1)],[NVarC,NClass]));
    end
    Pclass = Pclass./sum(Pclass,2); % NP x NClass
    f = zeros(NP,1);
    
    if any(isnan(XXa(:))) == 0  % faster version for complete dataset
        YYy = YY == 1;
        parfor n = 1:NP
            b_mtx_n = b_mtx(:,:,n);
            XXa_n = XXa(:,:,n);
            U = reshape(XXa_n*b_mtx_n,[NAlt,NCT,NRep,NClass]);
            if NVarS > 0
               Scale_n = Scale(:,:,n,:);
               U = U.*Scale_n; 
            end
            U = exp(U-max(U));
            U_prob = U./sum(U,1); % NAlt x NCT x NRep
            YYy_n = YYy(:,n);
            U_prod = prod(reshape(U_prob(YYy_n(:,ones(1,NRep),ones(1,1,NClass))),[NCT,NRep,NClass]),1); % 1 x NRep x NClass
            p(n,:) = mean(U_prod,2);
            f(n) = sum(p(n,:).*Pclass(n,:),2);
            U_prob = reshape(U_prob,[NAlt*NCT,1,NRep,NClass]); % NAlt*NCT x NVarA x NRep x NClass
            if WTP_space == 0
                if NVarS > 0
                   Scale_n = reshape(Scale_n,[NAlt*NCT,1,1,NClass]);
                   XXa_nn = XXa_n(:,:,ones(NRep,1),:).*Scale_n(:,:,ones(1,NRep),:);
                else
                   XXa_nn = XXa_n(:,:,ones(NRep,1),ones(NClass,1));
                end
%                 X_hat = sum(reshape(U_prob(:,ones(1,NVarA,1),:,:).*XXa_nn,[NAlt,NCT,NVarA,NRep,NClass]),1);
                X_hat = sum(reshape(U_prob.*XXa_nn,[NAlt,NCT,NVarA,NRep,NClass]),1);
                if NCT ~= 1
                    F = XXa_nn(YYy_n,:,:,:) - squeeze(X_hat); %NCT x NVarA x NRep x NClass
                    sumFsqueezed = squeeze(sum(F,1)); %NVarA x NRep x NClass
                else
                    sumFsqueezed = squeeze(XXa_nn(YYy_n,:,:,:)) - squeeze(X_hat); %NVarA x NRep x NClass
                end
                if any(Dist1)
                    b_mtx_n = reshape(b_mtx_n,[NVarA,NRep,NClass]);
                    sumFsqueezed(Dist1,:,:) = sumFsqueezed(Dist1,:,:).*b_mtx_n(Dist1,:,:);
                end
                if NVarS > 0
                    bss = reshape(b_mtx_n,[1,NVarA,NRep,NClass]);
%                     Fs = sum(F.*bss(ones(NCT,1),:,:,:),2);
                    Fs = sum(F.*bss,2);
                    Xs_n = Xs_sliced(:,:,n);
%                     Fs = Fs(:,ones(1,NVarS),:,:).*Xs_n(:,:,ones(1,NRep),ones(1,NClass)); % NCT x NVarS x NRep x NClass
                    Fs = Fs.*Xs_n; % NCT x NVarS x NRep x NClass
                    Fs = reshape(sum(Fs,1),[NVarS,NRep,NClass]); % NVarS x NRep x NClass
                end
            elseif WTP_space == 1
                % for non cost variables
                % NVarA x NRep x NClass
                b_mtx_wtp = reshape(b_mtx_n,[1,NVarA,NRep,NClass]);
                XXalpha = XXa_n(:,1:end-WTP_space,:,:).*b_mtx_wtp(:,WTP_matrix,:,:);
                if NVarS > 0
                   Scale_n = reshape(Scale_n,[NAlt*NCT,1,1,NClass]);
%                    XXalpha = XXalpha.*Scale_n(:,ones(1,NVarA-WTP_space),ones(1,NRep),:);
                   XXalpha = XXalpha.*Scale_n(:,ones(1,NVarA-WTP_space),:,:);
                end
%                 X_hat1 = sum(reshape(U_prob(:,ones(1,NVarA-WTP_space),:,:).*XXalpha,[NAlt,NCT,NVarA-WTP_space,NRep,NClass]),1);
                X_hat1 = sum(reshape(U_prob.*XXalpha,[NAlt,NCT,NVarA-WTP_space,NRep,NClass]),1);
                F1 = XXalpha(YY(:,n) == 1,:,:,:) - squeeze(X_hat1); %NCT x NVarA-WTP_space x NRep x NClass
                
                % for cost variable
                b_mtx_grad_n = b_mtx_grad(:,:,n);
%                 pX = squeeze(XXa_n(:,NVarA,ones(NClass*NRep,1))) + XXa_n(:,1:end-WTP_space)*b_mtx_grad_n(1:end-WTP_space,:); %
                pX = squeeze(XXa_n(:,NVarA,:)) + XXa_n(:,1:end-WTP_space)*b_mtx_grad_n(1:end-WTP_space,:); %
                pX = reshape(pX,[NAlt*NCT,NRep,NClass]);
                if NVarS > 0
                   Scale_n = reshape(Scale_n,[NAlt*NCT,1,NClass]);
%                    pX = pX.*Scale_n(:,ones(1,NRep),:);
                   pX = pX.*Scale_n;
                end
                X_hat2 = sum(reshape(squeeze(U_prob).*pX,[NAlt,NCT,NRep,NClass]),1);
                F2 = pX(YY(:,n) == 1,:,:) - squeeze(X_hat2); %NCT x NRep x NClass
                sumFsqueezed = zeros(NVarA,NRep,NClass);
                sumFsqueezed(1:NVarA-1,:,:) = squeeze(sum(F1,1));
                sumFsqueezed(NVarA,:,:) = squeeze(sum(F2,1));
                b_mtx_tmp = reshape(b_mtx_grad_n,[NVarA,NRep,NClass]);
                sumFsqueezed(Dist1,:,:) = sumFsqueezed(Dist1,:,:).* b_mtx_tmp(Dist1,:,:);
                if NVarS > 0
%                     Fs = reshape(F2.*b_mtx_tmp(NVarA*ones(NCT,1),:,:),[NCT,1,NRep,NClass]); % NCT x 1 x NRep x NClass
                    Fs = reshape(F2.*b_mtx_tmp(NVarA,:,:),[NCT,1,NRep,NClass]); % NCT x 1 x NRep x NClass
                    Xs_n = Xs_sliced(:,:,n);
%                     Fs = Fs(:,ones(1,NVarS),:,:).*Xs_n(:,:,ones(1,NRep),ones(1,NClass)); % NCT x NVarS x NRep x NClass
                    Fs = Fs.*Xs_n; % NCT x NVarS x NRep x NClass
                    Fs = reshape(sum(Fs,1),[NVarS,NRep,NClass]); % NVarS x NRep x NClass
                end
            end
            sumFsqueezed = reshape(permute(sumFsqueezed,[1 3 2]),[NVarA*NClass,NRep]);
            Pclass_n = Pclass(n,:); %reshape to 1x... ?
            if NVarS > 0
               Fs = reshape(permute(Fs,[1 3 2]),[NVarS*NClass,NRep]);
               U_prod2 = permute(U_prod(ones(NVarS,1),:,:),[1 3 2]);
               U_prod2 = reshape(U_prod2,[NVarS*NClass,NRep]);
               gtmps = -mean(Fs.*U_prod2,2).*reshape(Pclass_n(ones(NVarS,1),:),[NVarS*NClass,1])/f(n);
            end
            if FullCov == 0
                U_prod = permute(U_prod(ones(NVarA,1),:,:),[1 3 2]);
                U_prod = reshape(U_prod,[NVarA*NClass,NRep]);
                sumVC2tmp = sumFsqueezed.*VC2(:,:,n); % NVarA*NClass x NRep
                gtmp = mean([sumFsqueezed.*U_prod;sumVC2tmp.*U_prod],2);
                if NVarS > 0
                    g(n,:) = [-gtmp.*reshape(Pclass_n(ones(NVarA,1),:,ones(2,1)),[2*NVarA*NClass,1])/f(n);gtmps];
                else
                    g(n,:) = -gtmp.*reshape(Pclass_n(ones(NVarA,1),:,ones(2,1)),[2*NVarA*NClass,1])/f(n);
                end  
                
            else % FullCov = 1
                sumVC2tmp = sumFsqueezed(indx1,:).*VC2f(indx2,:,n);
                U_prod1 = permute(U_prod(ones(NVarA,1),:,:),[1 3 2]);
                U_prod1 = reshape(U_prod1,[NVarA*NClass,NRep]);
                U_prod2 = permute(U_prod(ones(NVarA*(NVarA-1)/2+NVarA,1),:,:),[1 3 2]);
                U_prod2 = reshape(U_prod2,[(NVarA*(NVarA-1)/2+NVarA)*NClass,NRep]);
                gtmp = mean([sumFsqueezed.*U_prod1;sumVC2tmp.*U_prod2],2);
                Pclass_tmp = [reshape(Pclass_n(ones(NVarA,1),:),[NVarA*NClass,1]);reshape(Pclass_n(ones(NVarA*(NVarA-1)/2+NVarA,1),:),[(NVarA*(NVarA-1)/2+NVarA)*NClass,1])];
                if NVarS > 0
                    g(n,:) = [-gtmp.*Pclass_tmp/f(n);gtmps];
                else
                    g(n,:) = -gtmp.*Pclass_tmp/f(n);
                end  
            end
        end
    else % Missing NCT
% save tmp1
% return
        if NVarS > 0
           YCT = reshape(sum(reshape(~isnan(YY),[NAlt,NCT,NP]),1) ~= 0,[NCT,NP]); 
        else
           YCT = zeros(0,NP);
        end
        parfor n = 1:NP
            YnanInd = ~isnan(YY(:,n));
            XXa_n = XXa(:,:,n);
            b_mtx_n = b_mtx(:,:,n);
            U = reshape(XXa_n(YnanInd,:,:)*b_mtx_n,[NAltMiss(n),NCTMiss(n),NRep,NClass]);
            if NVarS > 0
               Scale_n = reshape(Scale(:,:,n,:),[NAlt*NCT,1,NClass]);
               Scale_n = reshape(Scale_n(YnanInd,:,:),[NAltMiss(n),NCTMiss(n),1,NClass]);
               U = U.*Scale_n; 
            end
            U = exp(U - max(U));
            U_prob = U./sum(U,1); % NAlt x NCT x NRep
            YYy_n = YY(:,n) == 1;
            U_prod = prod(reshape(U_prob(YYy_n(YnanInd,ones(1,NRep),ones(1,1,NClass))),[NCTMiss(n),NRep,NClass]),1); % 1 x NRep x NClass
            p(n,:) = mean(U_prod,2);
            f(n) = sum(p(n,:).*Pclass(n,:),2);
            U_prob = reshape(U_prob,[NAltMiss(n)*NCTMiss(n),1,NRep,NClass]); % NAlt*NCT x NVarA x NRep x NClass
            if WTP_space == 0
                if NVarS > 0
                   Scale_n = reshape(Scale_n,[NAltMiss(n)*NCTMiss(n),1,1,NClass]);
                   XXa_n = XXa_n(YnanInd,:).*Scale_n;
                else
                   XXa_n = XXa_n(YnanInd,:);
                end
                X_hat = reshape(sum(reshape(U_prob.*XXa_n,[NAltMiss(n),NCTMiss(n),NVarA,NRep,NClass]),1),[NCTMiss(n),NVarA,NRep,NClass]);
                F = XXa_n(YYy_n(YnanInd),:,:,:) - X_hat;
                sumFsqueezed = reshape((sum(F,1)),[NVarA,NRep*NClass]);
                sumFsqueezed(Dist1,:) = sumFsqueezed(Dist1,:).*b_mtx_n(Dist1,:);
                sumFsqueezed = reshape(sumFsqueezed,[NVarA,NRep,NClass]);
                if NVarS > 0
                    bss = reshape(b_mtx_n,1,[NVarA,NRep,NClass]);
%                     Fs = sum(F.*bss(ones(NCTMiss(n),1),:,:,:),2);
                    Fs = sum(F.*bss,2);
                    Xs_n = Xs_sliced(:,:,n); 
%                     Fs = Fs(:,ones(1,NVarS),:,:).*Xs_n(YCT(:,n),:,ones(1,NRep),ones(1,NClass)); % NCT x NVarS x NRep x NClass
                    Fs = Fs.*Xs_n(YCT(:,n),:,:,:); % NCT x NVarS x NRep x NClass
                    Fs = reshape(sum(Fs,1),[NVarS,NRep,NClass]); % NVarS x NRep x NClass
                end
            elseif WTP_space == 1
                % for non cost variables
                % NVarA x NRep x NClass
                b_mtx_wtp = reshape(b_mtx_n,[1,NVarA,NRep,NClass]);
                XXalpha = XXa_n(YnanInd,1:end-WTP_space,:,:).*b_mtx_wtp(:,WTP_matrix,:,:);
                if NVarS > 0
                   Scale_n = reshape(Scale_n,[NAltMiss(n)*NCTMiss(n),1,1,NClass]);
                   XXalpha = XXalpha.*Scale_n(:,ones(1,NVarA-WTP_space),:,:);
                end
                X_hat1 = sum(reshape(U_prob.*XXalpha,[NAltMiss(n),NCTMiss(n),NVarA-WTP_space,NRep,NClass]),1);
                F1 = XXalpha(YYy_n(YnanInd) == 1,:,:,:) - reshape(X_hat1,[NCTMiss(n),NVarA-WTP_space,NRep,NClass]);
                
                % for cost variable
                b_mtx_grad_n = b_mtx_grad(:,:,n);
                pX = XXa_n(YnanInd,end) + XXa_n(YnanInd,1:end-WTP_space)*b_mtx_grad_n(1:end-WTP_space,:);
                pX = reshape(pX,[NAltMiss(n)*NCTMiss(n),NRep,NClass]);
                if NVarS > 0
                   Scale_n = reshape(Scale_n,[NAltMiss(n)*NCTMiss(n),1,NClass]);
%                    pX = pX.*Scale_n(:,ones(1,NRep),:);
                   pX = pX.*Scale_n;
                end
                X_hat2 = sum(reshape(reshape(U_prob,[NAltMiss(n)*NCTMiss(n),NRep,NClass]).*pX,[NAltMiss(n),NCTMiss(n),NRep,NClass]),1);
                F2 = pX(YYy_n(YnanInd) == 1,:,:) - reshape(X_hat2,[NCTMiss(n),NRep,NClass]); % NCT x NRep x NClass
                sumFsqueezed = zeros(NVarA,NRep,NClass);
                sumFsqueezed(1:NVarA-1,:,:) = reshape(sum(F1,1),[NVarA-1,NRep,NClass]);
                sumFsqueezed(NVarA,:,:) = reshape(sum(F2,1),[1,NRep,NClass]);
                
                b_mtx_tmp = reshape(b_mtx_grad_n,[NVarA,NRep,NClass]);
                sumFsqueezed(Dist1,:,:) = sumFsqueezed(Dist1,:,:).* b_mtx_tmp(Dist1,:,:);
                
                if NVarS > 0
%                     Fs = reshape(F2.*b_mtx_tmp(NVarA*ones(NCTMiss(n),1),:,:),[NCTMiss(n),1,NRep,NClass]); % NCT x 1 x NRep x NClass
                    Fs = reshape(F2.*b_mtx_tmp,[NCTMiss(n),1,NRep,NClass]); % NCT x 1 x NRep x NClass
                    Xs_n = Xs_sliced(:,:,n);
                    Fs = Fs.*Xs_n(YCT(:,n),:,:,:); % NCT x NVarS x NRep x NClass
                    Fs = reshape(sum(Fs,1),[NVarS,NRep,NClass]); % NVarS x NRep x NClass
                end
            end

            sumFsqueezed = reshape(permute(sumFsqueezed,[1 3 2]),[NVarA*NClass,NRep]);
            Pclass_n = Pclass(n,:); 
            if NVarS > 0
               Fs = reshape(permute(Fs,[1 3 2]),[NVarS*NClass,NRep]);
               U_prod2 = permute(U_prod(ones(NVarS,1),:,:),[1 3 2]);
               U_prod2 = reshape(U_prod2,[NVarS*NClass,NRep]);
               gtmps = -mean(Fs.*U_prod2,2).*reshape(Pclass_n(ones(NVarS,1),:),[NVarS*NClass,1])/f(n);
            end
            if FullCov == 0
                U_prod = permute(U_prod(ones(NVarA,1),:,:),[1 3 2]);
                U_prod = reshape(U_prod,[NVarA*NClass,NRep]);
                sumVC2tmp = sumFsqueezed.*VC2(:,:,n); % NVarA*NClass x NRep
                gtmp = mean([sumFsqueezed.*U_prod;sumVC2tmp.*U_prod],2);
                if NVarS > 0
                    g(n,:) = [-gtmp.*reshape(Pclass_n(ones(NVarA,1),:,ones(2,1)),[2*NVarA*NClass,1])/f(n);gtmps];
                else
                    g(n,:) = -gtmp.*reshape(Pclass_n(ones(NVarA,1),:,ones(2,1)),[2*NVarA*NClass,1])/f(n);
                end  
            else % FullCov = 1
                sumVC2tmp = sumFsqueezed(indx1,:).*VC2f(indx2,:,n);
                U_prod1 = permute(U_prod(ones(NVarA,1),:,:),[1 3 2]);
                U_prod1 = reshape(U_prod1,[NVarA*NClass,NRep]);
                U_prod2 = permute(U_prod(ones(NVarA*(NVarA-1)/2+NVarA,1),:,:),[1 3 2]);
                U_prod2 = reshape(U_prod2,[(NVarA*(NVarA-1)/2+NVarA)*NClass,NRep]);
                
                gtmp = mean([sumFsqueezed.*U_prod1;sumVC2tmp.*U_prod2],2);
                Pclass_tmp = [reshape(Pclass_n(ones(NVarA,1),:),[NVarA*NClass,1]);reshape(Pclass_n(ones(NVarA*(NVarA-1)/2+NVarA,1),:),[(NVarA*(NVarA-1)/2+NVarA)*NClass,1])];
                if NVarS > 0
                    g(n,:) = [-gtmp.*Pclass_tmp/f(n);gtmps];
                else
                    g(n,:) = -gtmp.*Pclass_tmp/f(n);
                end  
            end
        end
    end
    % class probabilities
    Pclass_tmp = reshape(Pclass(:,1:NClass-1),[NP,1,NClass-1]);
    Pclass_tmp = reshape(Pclass_tmp(:,ones(NClass-1,1),:),[NP,(NClass-1)^2]);
    Pclass = reshape(Pclass(:,:,ones(NClass-1,1)),[NP,NClass*(NClass-1)]);
    indx1 = 1:NClass-1;
    indx1 = indx1 + (indx1-1)*NClass;
    Pclass(:,indx1) = Pclass(:,indx1).*(1-Pclass(:,indx1));
    indx2 = 1:NClass*(NClass-1);
    indx2(indx1) = []; %
    
    Pclass(:,indx2) = -Pclass(:,indx2).*Pclass_tmp;
    Pclass = reshape(Pclass,[NP,NClass,NClass-1]);
    
%     g2 = sum(Pclass.*p(:,:,ones(NClass-1,1)),2); % NP x 1 xNClass-1
    g2 = sum(Pclass.*p,2); % NP x 1 xNClass-1
%     g2 = g2(:,ones(NVarC,1),:).*Xc(:,:, ones(NClass-1,1));
    g2 = g2.*Xc;
    g2 = reshape(g2,[NP,NVarC*(NClass-1)]); % gradient for class probability parameters
    g2 = -g2./f(:,ones(NVarC*(NClass-1),1));
    g = [g,g2];
    f = -log(f);
end