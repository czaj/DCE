function [f,g] = LL_lcmxl(YY,Xa,Xc,err_sliced, EstimOpt,B)

NVarA = EstimOpt.NVarA; 
NClass = EstimOpt.NClass; 
NVarC = EstimOpt.NVarC; 
NP = EstimOpt.NP; 
NAlt = EstimOpt.NAlt; 
NCT = EstimOpt.NCT; 
NRep = EstimOpt.NRep; 
FullCov = EstimOpt.FullCov; 
Dist = EstimOpt.Dist; 
WTP_space = EstimOpt.WTP_space; 
WTP_matrix = EstimOpt.WTP_matrix; 
indx1 = EstimOpt.indx1;
indx2 = EstimOpt.indx2;

b_mtx = zeros(NVarA, NP*NRep, NClass); ... 
VC = zeros(NVarA); ...
if EstimOpt.NumGrad == 0
    b_mtx_grad = zeros(NVarA,NRep,NClass,NP); ...% needed for gradient calculation in WTP_space
end
if FullCov == 0; 
    for c = 1:NClass; 
        b_mtx(:,:,c) = B((c-1)*NVarA+1:c*NVarA,ones(NP*NRep,1)) + diag(B(NVarA*NClass + (c - 1)*NVarA + 1 : NVarA*(NClass + c)).^2)*err_sliced((c-1)*NVarA+1:c*NVarA,:); 
        if sum(Dist(:,c)==1) > 0;  % lognormal
            b_mtx(Dist(:,c)==1,:,c) = exp(b_mtx(Dist(:,c)==1,:,c)); 
        elseif sum(Dist(:,c)==2) > 0;  % Spike
            b_mtx(Dist(:,c)==2,:,c) = max(b_mtx(Dist(:,c)==1,:,c),0); 
        end; ...
        if WTP_space > 0; 
            if EstimOpt.NumGrad == 0
                %size(b_mtx_grad(:,:,c,:))
                b_mtx_grad(:,:,c,:) = reshape(b_mtx(:,:,c) ,NVarA,NRep,NP); ...% needed for gradient calculation in WTP_space
            end
            b_mtx(1:end-WTP_space,:,c) = b_mtx(1:end-WTP_space,:,c).*b_mtx(WTP_matrix,:,c);      
        end; 
    end; 
else 
    VC_tmp = tril(ones(NVarA)); ...
    for c = 1:NClass; 
        VC(VC_tmp == 1) = B(NVarA*NClass + (c-1)*sum(1:NVarA) + 1 : NVarA*NClass + c*sum(1:NVarA)); ...
        b_mtx(:,:,c) = B((c-1)*NVarA+1:c*NVarA,ones(NP*NRep,1)) + VC*err_sliced((c-1)*NVarA+1:c*NVarA,:); 
        if sum(Dist(:,c)==1) > 0;  % lognormal
            b_mtx(Dist(:,c)==1,:,c) = exp(b_mtx(Dist(:,c)==1,:,c)); 
        elseif sum(Dist(:,c)==2) > 0;  % Spike
            b_mtx(Dist(:,c)==2,:,c) = max(b_mtx(Dist(:,c)==1,:,c),0); 
        end; 
        if WTP_space > 0; 
            if EstimOpt.NumGrad == 0
                b_mtx_grad(:,:,c,:) = reshape(b_mtx(:,:,c) ,NVarA,NRep,NP); ...% needed for gradient calculation in WTP_space
            end
            b_mtx(1:end - WTP_space,:,c) = b_mtx(1:end-WTP_space,:,c).*b_mtx(WTP_matrix,:,c); 
        end; 
    end; 
end;
if EstimOpt.NumGrad == 0
    b_mtx_grad = reshape(b_mtx_grad, NVarA, NRep*NClass, NP); ...% needed for gradient calculation in WTP_space
end

%     Dist = reshape(EstimOpt.Dist,[NVarA,1,NClass]);
%     b_mtx(Dist(:,ones(1,NP*NRep,1),:)==1) = exp(b_mtx(Dist(:,ones(1,NP*NRep,1),:)==1)); % this is slower  

p = zeros(NP,NClass); 
if EstimOpt.NumGrad == 1 % numerical gradient
    if any(isnan(Xa(:))) == 0;  % faster version for complete dataset      
        for n = 1:NP; 
            YYn = YY(:,n*ones(NRep,1),ones(1,1,NClass)); ...
            U = reshape(Xa(:,:,n)*reshape(b_mtx(:,((n-1)*NRep+1):n*NRep,:),[NVarA,NRep*NClass]),[NAlt,NCT,NRep,NClass]); ...
            U_max = max(U); ...
            U = exp(U - U_max(ones(NAlt,1),:,:,:)); ...
            U_sum = reshape(sum(U,1),NCT,NRep,NClass); ...
            U_selected = reshape(U(YYn==1),NCT,NRep,NClass); ...   
            p(n,:) = mean(prod(U_selected ./ U_sum,1),2); 
        end; 
    else 
        for n = 1:NP;        
            YYn = YY(:,n*ones(NRep,1),ones(1,1,NClass)); ... %it might be done once, not for every class    
            U = reshape(Xa(~isnan(YY(:,n)),:,n)*reshape(b_mtx(:,((n-1)*NRep+1):n*NRep,:),[NVarA,NRep*NClass]),[NAlt,NCT-sum(isnan(YY(1:NAlt:end,n))),NRep,NClass]); ...
            U_max = max(U); ...
            U = exp(U - U_max(ones(NAlt,1),:,:,:)); ... % NAlt x NCT - NaNs x NRep x NClass
            U_sum = reshape(nansum(U,1),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep,NClass); ...
            U_selected = reshape(U(YYn(~isnan(YY(:,n)),:)==1),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep, NClass); ...
            p(n,:) = mean(prod(U_selected ./ U_sum,1),2); 
        end; 
    end 

    if FullCov == 0; 
        Pclass = exp(Xc*reshape([B(2*NClass*NVarA+1:end);zeros(NVarC,1)],NVarC,NClass)); ...
    else 
        Pclass = exp(Xc*reshape([B(NClass*(NVarA+sum(1:NVarA))+1:end);zeros(NVarC,1)],NVarC,NClass)); 
    end; 
    Pclass_sum = sum(Pclass,2); ...
    Pclass = Pclass./Pclass_sum(:,ones(NClass,1)); ... % NP x NClass

    f = -log(sum(p.*Pclass,2)); 

else % analitycal gradient
    
    if FullCov == 0    
        g = zeros(NP, 2*NVarA*NClass);
        VC2 = zeros(NVarA, NClass,NRep, NP);
        for c = 1:NClass
            VC2(:,c,:,:) = reshape(2*diag(B(NVarA*NClass + (c - 1)*NVarA + 1 : NVarA*(NClass + c)))*err_sliced((c-1)*NVarA+1:c*NVarA,:), NVarA,1, NRep, NP);
           
        end
        VC2 = reshape(VC2, NVarA*NClass, NRep, NP);
    else
        g = zeros(NP, (2*NVarA+NVarA*(NVarA-1)/2)*NClass);
        VC2 = reshape(err_sliced, NVarA*NClass, NRep, NP);
    end
    
    if FullCov == 0; 
        Pclass = exp(Xc*reshape([B(2*NClass*NVarA+1:end);zeros(NVarC,1)],NVarC,NClass)); ...
    else 
        Pclass = exp(Xc*reshape([B(NClass*(NVarA+sum(1:NVarA))+1:end);zeros(NVarC,1)],NVarC,NClass)); 
    end; 
    Pclass_sum = sum(Pclass,2); ...
    Pclass = Pclass./Pclass_sum(:,ones(NClass,1)); ... % NP x NClass
    f= zeros(NP,1);
    if any(isnan(Xa(:))) == 0;  % faster version for complete dataset      
        for n = 1:NP; 
            YYn = YY(:,n*ones(NRep,1),ones(1,1,NClass)); ...
            U = reshape(Xa(:,:,n)*reshape(b_mtx(:,((n-1)*NRep+1):n*NRep,:),[NVarA,NRep*NClass]),[NAlt,NCT,NRep,NClass]); ...
            U_max = max(U); ...
            U = exp(U - U_max(ones(NAlt,1),:,:,:)); ...
            U_sum = sum(U,1); ...
            U_prob = U./U_sum(ones(NAlt,1,1),:,:,:); ... % NAlt x NCT x NRep
            U_prod = prod(reshape(U_prob(YYn==1),NCT,NRep, NClass),1); ... % 1 x NRep x NClass
            %U_selected = reshape(U(YYn==1),NCT,NRep,NClass); ...   
            p(n,:) = mean(U_prod,2); 
            f(n) = max(sum(p(n,:).*Pclass(n,:),2), realmin);
            
            U_prob = reshape(U_prob, NAlt*NCT,1, NRep, NClass); ... % NAlt*NCT x NVarA x NRep x NClass
            if WTP_space == 0
                X_hat = sum(reshape(U_prob(:,ones(1,NVarA,1),:,:).* Xa(:,:, n*ones(NRep,1), ones(NClass,1)), NAlt, NCT, NVarA, NRep, NClass),1); ...
                if NCT ~= 1
                    F = Xa(YY(:,n) == 1,:,n*ones(NRep,1), ones(NClass,1)) - squeeze(X_hat); ... %NCT x NVarA x NRep x NClass
                    sumFsqueezed = squeeze(sum(F,1)); ... %NVarA x NRep x NClass
                else
                    sumFsqueezed = squeeze(Xa(YY(:,n) == 1,:,n*ones(NRep,1), ones(NClass,1))) - squeeze(X_hat); %NVarA x NRep x NClass
                end
                
                sumFsqueezed(Dist(1:NVarA)==1, :,:) = sumFsqueezed(Dist(1:NVarA)==1, :,:).*b_mtx(Dist(1:NVarA)==1,((n-1)*NRep+1):n*NRep,:);
                
            elseif WTP_space == 1
                % for non cost variables
                 % NVarA x NRep x NClass
                b_mtx_wtp = reshape(b_mtx(:,((n-1)*NRep+1):n*NRep,:), 1, NVarA, NRep, NClass);
                Xalpha = Xa(:,1:end-WTP_space, n*ones(NRep,1), ones(1, NClass)).*b_mtx_wtp(ones(NAlt*NCT,1),WTP_matrix,:,:);
                % for non-cost variables
                X_hat1 = sum(reshape(U_prob(:,ones(1,NVarA-WTP_space),:,:).* Xalpha, NAlt, NCT, NVarA-WTP_space, NRep, NClass),1); ...
                F1 = Xalpha(YY(:,n) == 1,:,:,:) - squeeze(X_hat1); ... %NCT x NVarA-WTP_space x NRep x NClass
                
                % for cost variable
                pX = squeeze(Xa(:,NVarA, n*ones(NClass*NRep,1))) + Xa(:,1:end-WTP_space,n)*b_mtx_grad(1:end-WTP_space,:,n); % 
                pX = reshape(pX, NAlt*NCT, NRep, NClass);
                X_hat2 = sum(reshape(squeeze(U_prob).*pX, NAlt, NCT, NRep, NClass),1); ...
                F2 = pX(YY(:,n) == 1,:,:) - squeeze(X_hat2); ... %NCT x NRep x NClass
                sumFsqueezed = zeros(NVarA, NRep, NClass);
                sumFsqueezed(1:NVarA-1, :,:) = squeeze(sum(F1,1));
                sumFsqueezed(NVarA, :,:) = squeeze(sum(F2,1));
                
                b_mtx_tmp = reshape(b_mtx_grad(:,:,n), NVarA, NRep, NClass);
                sumFsqueezed(Dist(1:NVarA)==1, :,:) = sumFsqueezed(Dist(1:NVarA)==1, :,:).* b_mtx_tmp(Dist(1:NVarA)==1,:,:);
            end
            sumFsqueezed = reshape(permute(sumFsqueezed, [1 3 2]), NVarA*NClass, NRep);
            if FullCov == 0
                U_prod = permute(U_prod(ones(NVarA,1),:,:), [1 3 2]);
                U_prod = reshape(U_prod, NVarA*NClass, NRep);
                sumVC2tmp = sumFsqueezed.*VC2(:,:,n); ... % NVarA*NClass x NRep
                gtmp = mean([sumFsqueezed.*U_prod; sumVC2tmp.*U_prod],2);
                g(n,:) = -gtmp.*reshape(Pclass(n*ones(NVarA,1),:, ones(2,1)),2*NVarA*NClass,1)/f(n);
            else % FullCov = 1
                sumVC2tmp = sumFsqueezed(indx1,:).*VC2(indx2,:,n); ...   
                U_prod1 = permute(U_prod(ones(NVarA,1),:,:), [1 3 2]);
                U_prod1 = reshape(U_prod1, NVarA*NClass, NRep);  
                U_prod2 = permute(U_prod(ones(NVarA*(NVarA-1)/2+NVarA,1),:,:), [1 3 2]);
                U_prod2 = reshape(U_prod2, (NVarA*(NVarA-1)/2+NVarA)*NClass, NRep);

                gtmp = mean([sumFsqueezed.*U_prod1; sumVC2tmp.*U_prod2],2);    
                Pclass_tmp = [reshape(Pclass(n*ones(NVarA,1),:), NVarA*NClass,1); reshape(Pclass(n*ones(NVarA*(NVarA-1)/2+NVarA,1),:), (NVarA*(NVarA-1)/2+NVarA)*NClass,1)]; 
                g(n,:) = -gtmp.*Pclass_tmp/f(n);            
            end
        end; 
    else % Missing NCT
        for n = 1:NP; 
            
            
%             YYn = YY(:,n*ones(NRep,1),ones(1,1,NClass)); ... %it might be done once, not for every class    
%             U = reshape(Xa(~isnan(YY(:,n)),:,n)*reshape(b_mtx(:,((n-1)*NRep+1):n*NRep,:),[NVarA,NRep*NClass]),[NAlt,NCT-sum(isnan(YY(1:NAlt:end,n))),NRep,NClass]); ...
%             U_max = max(U); ...
%             U = exp(U - U_max(ones(NAlt,1),:,:,:)); ... % NAlt x NCT - NaNs x NRep x NClass
%             U_sum = reshape(nansum(U,1),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep,NClass); ...
%             U_selected = reshape(U(YYn(~isnan(YY(:,n)),:)==1),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep, NClass); ...
%             p(n,:) = mean(prod(U_selected ./ U_sum,1),2); 
%         
        
            YnanInd = ~isnan(YY(:,n));
            NCT_n =EstimOpt.NCT-sum(isnan(YY(1:EstimOpt.NAlt:end,n)));
            YYn = YY(:,n*ones(NRep,1),ones(1,1,NClass)); ...
            U = reshape(Xa(YnanInd,:,n)*reshape(b_mtx(:,((n-1)*NRep+1):n*NRep,:),[NVarA,NRep*NClass]),[NAlt,NCT_n,NRep,NClass]); ...
            U_max = max(U); ...
            U = exp(U - U_max(ones(NAlt,1),:,:,:)); ...
            U_sum = sum(U,1); ...
            U_prob = U./U_sum(ones(NAlt,1,1),:,:,:); ... % NAlt x NCT x NRep
            U_prod = prod(reshape(U_prob(YYn(YnanInd,:,:)==1),NCT_n,NRep, NClass),1); ... % 1 x NRep x NClass
            %U_selected = reshape(U(YYn==1),NCT,NRep,NClass); ...   
            p(n,:) = mean(U_prod,2); 
            f(n) = max(sum(p(n,:).*Pclass(n,:),2), realmin);
            
            U_prob = reshape(U_prob, NAlt*NCT_n,1, NRep, NClass); ... % NAlt*NCT x NVarA x NRep x NClass
            if WTP_space == 0
                X_hat = sum(reshape(U_prob(:,ones(1,NVarA,1),:,:).* Xa(YnanInd,:, n*ones(NRep,1), ones(NClass,1)), NAlt, NCT_n, NVarA, NRep, NClass),1); ...
                if NCT_n ~= 1
                    F = Xa(YY(:,n) == 1,:,n*ones(NRep,1), ones(NClass,1)) - squeeze(X_hat); ... %NCT x NVarA x NRep x NClass
                    sumFsqueezed = squeeze(sum(F,1)); ... %NVarA x NRep x NClass
                else
                    sumFsqueezed = squeeze(Xa(YY(:,n) == 1,:,n*ones(NRep,1), ones(NClass,1))) - squeeze(X_hat); %NVarA x NRep x NClass
                end
                
                sumFsqueezed(Dist(1:NVarA)==1, :,:) = sumFsqueezed(Dist(1:NVarA)==1, :,:).*b_mtx(Dist(1:NVarA)==1,((n-1)*NRep+1):n*NRep,:);
                
            elseif WTP_space == 1
                % for non cost variables
                 % NVarA x NRep x NClass
                b_mtx_wtp = reshape(b_mtx(:,((n-1)*NRep+1):n*NRep,:), 1, NVarA, NRep, NClass);
                Xalpha = Xa(YnanInd,1:end-WTP_space, n*ones(NRep,1), ones(1, NClass)).*b_mtx_wtp(ones(NAlt*NCT_n,1),WTP_matrix,:,:);
                X_hat1 = sum(reshape(U_prob(:,ones(1,NVarA-WTP_space),:,:).* Xalpha, NAlt, NCT_n, NVarA-WTP_space, NRep, NClass),1); ...
                F1 = Xalpha(YY(YnanInd,n) == 1,:,:,:) - squeeze(X_hat1); ... %NCT x NVarA-WTP_space x NRep x NClass
                
                % for cost variable
                pX = squeeze(Xa(YnanInd,NVarA, n*ones(NClass*NRep,1))) + Xa(YnanInd,1:end-WTP_space,n)*b_mtx_grad(1:end-WTP_space,:,n); % 
                pX = reshape(pX, NAlt*NCT_n, NRep, NClass);
                X_hat2 = sum(reshape(squeeze(U_prob).*pX, NAlt, NCT_n, NRep, NClass),1); ...
                F2 = pX(YY(YnanInd,n) == 1,:,:) - squeeze(X_hat2); ... %NCT x NRep x NClass
                sumFsqueezed = zeros(NVarA, NRep, NClass);
                sumFsqueezed(1:NVarA-1, :,:) = squeeze(sum(F1,1));
                sumFsqueezed(NVarA, :,:) = squeeze(sum(F2,1));
                
                b_mtx_tmp = reshape(b_mtx_grad(:,:,n), NVarA, NRep, NClass);
                sumFsqueezed(Dist(1:NVarA)==1, :,:) = sumFsqueezed(Dist(1:NVarA)==1, :,:).* b_mtx_tmp(Dist(1:NVarA)==1,:,:);
            end

            sumFsqueezed = reshape(permute(sumFsqueezed, [1 3 2]), NVarA*NClass, NRep);
            if FullCov == 0
                U_prod = permute(U_prod(ones(NVarA,1),:,:), [1 3 2]);
                U_prod = reshape(U_prod, NVarA*NClass, NRep);
                sumVC2tmp = sumFsqueezed.*VC2(:,:,n); ... % NVarA*NClass x NRep
                gtmp = mean([sumFsqueezed.*U_prod; sumVC2tmp.*U_prod],2);
                g(n,:) = -gtmp.*reshape(Pclass(n*ones(NVarA,1),:, ones(2,1)),2*NVarA*NClass,1)/f(n);
            else % FullCov = 1
                sumVC2tmp = sumFsqueezed(indx1,:).*VC2(indx2,:,n); ...   
                U_prod1 = permute(U_prod(ones(NVarA,1),:,:), [1 3 2]);
                U_prod1 = reshape(U_prod1, NVarA*NClass, NRep);  
                U_prod2 = permute(U_prod(ones(NVarA*(NVarA-1)/2+NVarA,1),:,:), [1 3 2]);
                U_prod2 = reshape(U_prod2, (NVarA*(NVarA-1)/2+NVarA)*NClass, NRep);

                gtmp = mean([sumFsqueezed.*U_prod1; sumVC2tmp.*U_prod2],2);    
                Pclass_tmp = [reshape(Pclass(n*ones(NVarA,1),:), NVarA*NClass,1); reshape(Pclass(n*ones(NVarA*(NVarA-1)/2+NVarA,1),:), (NVarA*(NVarA-1)/2+NVarA)*NClass,1)]; 
                g(n,:) = -gtmp.*Pclass_tmp/f(n);            
            end
        end; 
    end 
    % class probabilities 
   Pclass_tmp = reshape(Pclass(:,1:EstimOpt.NClass-1), EstimOpt.NP, 1, EstimOpt.NClass-1);
   Pclass_tmp = reshape(Pclass_tmp(:, ones(EstimOpt.NClass-1,1),:),EstimOpt.NP,(EstimOpt.NClass-1)^2);
   Pclass = reshape(Pclass(:,:, ones(EstimOpt.NClass-1,1)), EstimOpt.NP, EstimOpt.NClass*(EstimOpt.NClass-1));
   indx1 = 1:EstimOpt.NClass-1;
   indx1 = indx1 + (indx1-1)*EstimOpt.NClass;
   Pclass(:, indx1) = Pclass(:, indx1).*(1-Pclass(:, indx1));
   indx2 = 1:EstimOpt.NClass*(EstimOpt.NClass-1);
   indx2(indx1) = []; % 

   Pclass(:,indx2) =-Pclass(:,indx2).*Pclass_tmp;
   Pclass = reshape(Pclass, EstimOpt.NP, EstimOpt.NClass, EstimOpt.NClass-1);
   
   g2 = sum(Pclass.*p(:,:, ones(EstimOpt.NClass-1,1)),2); % NP x 1 xNClass-1
   g2 = g2(:, ones(EstimOpt.NVarC,1),:).*Xc(:,:, ones(EstimOpt.NClass-1,1));
   g2 = reshape(g2, EstimOpt.NP, EstimOpt.NVarC*(EstimOpt.NClass-1)); % gradient for class probability parameters
   g2 = -g2./f(:, ones(EstimOpt.NVarC*(EstimOpt.NClass-1),1));
   g = [g, g2];
    f = -log(f);

end