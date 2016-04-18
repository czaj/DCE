function [f,g,h] = LL_mxl(YY,XXa,XXm,Xs,err,EstimOpt,b0)

% save LL_mxl
% return

NAlt = EstimOpt.NAlt; 
NCT = EstimOpt.NCT; 
NP = EstimOpt.NP; 
NRep = EstimOpt.NRep; 
NVarA = EstimOpt.NVarA; 
NVarM = EstimOpt.NVarM; 
NVarS = EstimOpt.NVarS; 
Dist = EstimOpt.Dist(2:end); 
WTP_space = EstimOpt.WTP_space; 
WTP_matrix = EstimOpt.WTP_matrix; 
FullCov = EstimOpt.FullCov; 
DiagIndex = EstimOpt.DiagIndex; 
Triang = EstimOpt.Triang; 
% NVarNLT = EstimOpt.NVarNLT; 
% NLTVariables = EstimOpt.NLTVariables; 
% NLTType = EstimOpt.NLTType; 
Johnson = EstimOpt.Johnson; 
NCTMiss = EstimOpt.NCTMiss; 
NAltMiss = EstimOpt.NAltMiss; 

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
%     b0v = b0v.^2;
    VC = diag(b0v); 
    b0m = b0(NVarA*2+1:NVarA*(NVarM+2));     
    b0m = reshape(b0m,NVarA, NVarM); 
    b0s = b0(NVarA*(NVarM+2)+1:NVarA*(NVarM+2)+NVarS); 
    b0j = b0(NVarA*(NVarM+2)+NVarS+1:NVarA*(NVarM+2)+NVarS+2*Johnson); 
else 
	b0v = b0(NVarA+1:NVarA+sum(1:NVarA)); 
    tmp = b0v(DiagIndex);
    b0v(DiagIndex(Dist >=3 & Dist <=5)) = 1;
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
        VC(Dist >= 3 & Dist <= 5,:) = VC(Dist >= 3 & Dist <= 5,:)./tmp(:, ones(1,NVarA));
    end
    b0m = b0(NVarA*(NVarA/2+1.5)+1:NVarA*(NVarA/2+1.5+NVarM)); 
    b0m = reshape(b0m,NVarA, NVarM); 
    b0s = b0(NVarA*(NVarA/2+1.5+NVarM)+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS); 
    b0j = b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS+2*Johnson); 
end 

% b0n = b0a(:,ones(NP,1)) + b0m*XXm;
% b0n = reshape(b0n((1:size(b0n,1))'*ones(1,NRep),(1:size(b0n,2))'),NVarA,NRep*NP);  % NVarA x NRep*NP
% b_mtx = b0n + VC*err;  % NVarA x NRep*NP

b_mtx = b0a(:,ones(NP,1)) + b0m*XXm;
b_mtx = reshape(b_mtx((1:size(b_mtx,1))'*ones(1,NRep),(1:size(b_mtx,2))'),NVarA,NRep*NP) + VC*err;  % NVarA x NRep*NP

if sum(Dist==1) > 0 % Log - normal
    b_mtx(Dist==1,:) = exp(b_mtx(Dist==1,:)); 
end
if sum(Dist==2) > 0 % Spike
    b_mtx(Dist==2,:) = max(b_mtx(Dist==2,:),0); 
end
if sum(Dist==3) > 0 % Triangular
    tmp = normcdf(b_mtx(Dist==3,:)); 
    Triang = Triang(ones(NRep*NP,1),:)';
    b0triag_c = b0triag_c(:, ones(NRep*NP,1));
    b0triag_b = b0triag_b(:, ones(NRep*NP,1));
    Ftriang =  (b0triag_c - Triang)./(b0triag_b- Triang);
    bmtx_triang = zeros(size(tmp));
    tmp2 = (b0triag_b- Triang).*(b0triag_c - Triang);
    bmtx_triang(tmp < Ftriang) = Triang(tmp < Ftriang)+ sqrt(tmp(tmp < Ftriang).*tmp2(tmp < Ftriang));
    tmp2 = (b0triag_b- Triang).*(b0triag_b-b0triag_c);
    bmtx_triang(tmp >= Ftriang) = b0triag_b(tmp >= Ftriang)- sqrt((1-tmp(tmp >= Ftriang)).*tmp2(tmp >= Ftriang));  
    b_mtx(Dist==3,:) = bmtx_triang;
end
if sum(Dist==4) > 0 % Weibull
    tmp = -log(1-normcdf(b_mtx(Dist==4,:))); 
    b_mtx(Dist==4,:) = b0weibA(:, ones(1,NP*NRep,1)).*(tmp.^b0weibB(:, ones(1,NP*NRep,1)));
end
if sum(Dist>=5) > 0 % Johnson
    if sum(Dist==5) > 0 % Sinh-Arcsinh
        b_mtx(Dist==5,:) = b0sinhA(:,ones(NRep*NP,1))+ b0sinhB(:,ones(NRep*NP,1)).*asinh(b_mtx(Dist==5,:));

        b_mtx(Dist==5,:) = b0j(1:Johnson, ones(NRep*NP,1)) + exp(b0j(Johnson+1:end, ones(NRep*NP,1))).*sinh(b_mtx(Dist==5,:));  
    end
    if sum(Dist==6) > 0 % Johnson Sb
        tmp = exp(b_mtx(Dist==6,:));
        b_mtx(Dist==6,:) = tmp./(1+tmp); 
        b_mtx(Dist==6,:) = b0j(1:Johnson, ones(NRep*NP,1)) + exp(b0j(Johnson+1:end, ones(NRep*NP,1))).*b_mtx(Dist==6,:);
    end
    if sum(Dist==7) > 0 % Johnson Su
        b_mtx(Dist==7,:) = sinh(b_mtx(Dist==7,:)); 
        b_mtx(Dist==7,:) = b0j(1:Johnson, ones(NRep*NP,1)) + exp(b0j(Johnson+1:end, ones(NRep*NP,1))).*b_mtx(Dist==7,:);
    end

end

if WTP_space > 0
    if EstimOpt.NumGrad == 0
       b_mtx_grad = reshape(b_mtx,NVarA,NRep,NP); % needed for gradient calculation in WTP_space
    end
    b_mtx(1:end-WTP_space,:) = b_mtx(1:end-WTP_space,:).*b_mtx(WTP_matrix,:); 
else
    b_mtx_grad = zeros(0,0,NP);
end

if NVarS > 0
    cs = reshape(exp(Xs*b0s),NAlt*NCT,1,NP);
    XXa = XXa .* cs(:,ones(1,NVarA,1),:);
end

b_mtx = reshape(b_mtx,NVarA,NRep,NP); 
%max(max(max(b_mtx)))
p0 = zeros(NP,1); 

if nargout == 1 % function value only    

    
    if any(isnan(XXa(:))) == 0 % faster version for complete dataset
        parfor n = 1:NP
            U = reshape(XXa(:,:,n)*b_mtx(:,:,n),NAlt-1,NCT,NRep);            
            U_sum = reshape(sum(exp(U),1),NCT,NRep);
            p0(n) = mean(prod(1./(1+U_sum),1));
        end;
    else  % this works only if NAlt is constant for each respondent
        parfor n = 1:NP
            YnanInd = ~isnan(YY(:,n));
            XXa_n = XXa(:,:,n);
            U = reshape(XXa_n(YnanInd,:)*b_mtx(:,:,n),NAltMiss(n)-1,NCTMiss(n),NRep);
            U_sum = reshape(sum(exp(U),1),NCTMiss(n),NRep);
            p0(n) = mean(prod(1./(1+U_sum),1));
        end; 
    end 

%     na razie sobie to zostawi³em, jako notatkê

    % szybsze ale podatne na eksplozje:
    % for n = 1:NP
    %     U = exp(XXa_n(:,:,n)*b_mtx(:,:,n)); 
    %     U_selected = reshape(U(YY2(:,n*ones(NRep,1)) == 1),NCT,NRep); 
    %     U_sum = reshape((sum(reshape(U,NAlt,NCT,NRep),1)),NCT,NRep); 
    %     p0(n) = mean(prod(U_selected ./ U_sum),2);
    % %     U_selected = U(YY2(:,n*ones(NRep,1)) == 1);
    % %     U_sum = sum(reshape(U,NAlt,NCT*NRep),1)';
    % %     p0(n) = mean(prod(reshape(U_selected ./ U_sum,NCT,NRep)),2);
    % end;

    % wyci¹ganie rzeczy z pêtli niekoniecznie pomaga:
    % YY3 = reshape(YY,NAlt*NCT,1,NP);
    % tic ; for n = 1:NP
    %     U(:,:,n) = exp(XXa_n(:,:,n)*b_mtx(:,:,n));
    % end;
    % U_selected = reshape(U(YY3(:,ones(NRep,1),:) == 1),NCT,NRep,NP); 
    % U_sum = reshape((sum(reshape(U,NAlt,NCT,NRep,NP),1)),NCT,NRep,NP); 
    % p0 = reshape(mean(prod(U_selected ./ U_sum),2),NP,1); toc
    

elseif nargout == 2 %  function value + gradient
    
    if NVarS > 0
        Xs_sliced = reshape(Xs, NAlt*NCT, NP, NVarS);
    else
        Xs_sliced = reshape(Xs, NAlt*NCT, NP, 0);
    end
	if FullCov == 0    
        g = zeros(NP, 2*NVarA + NVarS);
%         VC2 = reshape(2*diag(b0(NVarA+1:NVarA*2))*err, NVarA, NRep, NP);
        VC2 = reshape(err, NVarA, NRep, NP);
        VC2f = zeros(0,0,NP);
	else
        g = zeros(NP, 2*NVarA+NVarA*(NVarA-1)/2 + NVarS);
        VC2 = zeros(0,0,NP);
        VC2f = reshape(err, NVarA, NRep, NP);
	end
    
    if any(isnan(XXa(:))) == 0 % faster version for complete dataset
        parfor n = 1:NP
            b_mtx_n = b_mtx(:,:,n);
            XXa_n = XXa(:,:,n);
            % With normalization
%             U = reshape(XXa_n*b_mtx_n,NAlt-1,NCT,NRep);
%             maxU = max(U);
%             eU = exp(bsxfun(@minus, U, maxU));
%             U_sum = reshape(exp(-maxU)+sum(eU,1),NCT,NRep);
%             U_prod = prod(bsxfun(@rdivide, reshape(exp(-maxU), [NCT,NRep]), U_sum),1);
%             p0(n) = mean(U_prod);
%             U_prob = reshape(eU, (NAlt-1)*NCT,1,NRep);  % (NAlt-1)*NCT x 1 x NRep
%             %p0(n) = max(mean(U_prod),realmin);
            
            % Without normalization
            U = reshape(exp(XXa_n*b_mtx_n),NAlt-1,NCT,NRep);  % NAlt -1 x NCT x NRep
            U_sum = reshape(1+sum(U,1),NCT,NRep);
            U_prod = prod(1./U_sum,1);
            p0(n) = mean(U_prod);
            %p0(n) = max(mean(U_prod),realmin);
            U_prob = reshape(U, (NAlt-1)*NCT,1,NRep);  % (NAlt-1)*NCT x 1 x NRep
            
            % calculations for gradient
            U_prob = min(U_prob, realmax);
            U_sum = min(U_sum, realmax);
            if WTP_space == 0   
                X_hat = sum(reshape(bsxfun(@times,U_prob,XXa_n), NAlt-1, NCT, NVarA, NRep),1); % 1 x NCT x NVarA x NRep
                if NCT ~= 1
                    X_hat = min(X_hat, realmax);
                    F = -bsxfun(@rdivide,reshape(X_hat, [NCT,NVarA,NRep]), reshape(U_sum, [NCT,1,NRep]));
                    sumFsqueezed = reshape(sum(F,1),[NVarA,NRep]);  %NVarA x NRep
                else
                    sumFsqueezed = -bsxfun(@rdivide,reshape(X_hat, [NVarA,NRep]), U_sum);
                end                
                sumFsqueezed(Dist==1, :) = sumFsqueezed(Dist==1, :).*b_mtx_n(Dist==1,:);                 
            else
                Xalpha = bsxfun(@times,XXa_n(:,1:end-WTP_space),reshape(b_mtx_n(WTP_matrix,:),[1,NVarA-WTP_space,NRep])); % (NAlt-1)*NCT x NVarA - WTP x NRep
                % for non-cost variables
                X_hat1 = sum(reshape(bsxfun(@times,U_prob,Xalpha), NAlt-1, NCT, NVarA-WTP_space, NRep),1); % 1 x NCT x NVarA-WTP_space x NRep
                F1 = -bsxfun(@rdivide,reshape(X_hat1, [NCT,NVarA-WTP_space,NRep]), reshape(U_sum, [NCT,1,NRep]));
                
                % for cost variables
                b_mtx_grad_n = b_mtx_grad(:,:,n);
                if WTP_space == 1 % without starting the loop
                    pX = bsxfun(@plus,XXa_n(:,NVarA),  XXa_n(:,1:end-WTP_space)*b_mtx_grad_n(1:end-WTP_space,:)); % (NAlt-1)*NCT x NRep
                    %pX = squeeze(XXa_n(:,NVarA,ones(NRep,1))) + XXa_n(:,1:end-WTP_space)*b_mtx_grad_n(1:end-WTP_space,:);
                    X_hat2 = sum(reshape(bsxfun(@times,squeeze(U_prob),pX), NAlt-1, NCT, WTP_space, NRep),1); 
                else
                    pX = zeros(NCT*(NAlt-1), WTP_space, NRep);
                    for i = 1:WTP_space
                        pX(:,i,:) = bsxfun(@plus,XXa_n(:,NVarA-WTP_space+i),  XXa_n(:,WTP_matrix == NVarA-WTP_space+i)*b_mtx_grad_n(WTP_matrix == NVarA-WTP_space+i,:)); % (NAlt-1)*NCT x NRep
                    end
                    X_hat2 = sum(reshape(bsxfun(@times,U_prob,pX), NAlt-1, NCT, WTP_space, NRep),1); 
                end
                F2 = -bsxfun(@rdivide,reshape(X_hat2, [NCT,WTP_space,NRep]), reshape(U_sum, [NCT,1,NRep]));
                 
                sumFsqueezed = [squeeze(sum(F1,1));reshape(sum(F2,1),WTP_space, NRep)];  %NVarA x NRep
                sumFsqueezed(Dist==1, :) = sumFsqueezed(Dist==1, :).*b_mtx_grad_n(Dist==1,:);   
                
            end
            sumFsqueezed  = min(sumFsqueezed , realmax);
            if NVarS >0
                if WTP_space == 0
                    FScale = sum(sumFsqueezed.*b_mtx_n,1); % 1 x NRep
                else
                    if WTP_space == 1
                        FScale = sum(squeeze(sum(F2,1)).*b_mtx_n(NVarA-WTP_space+1:end,:),1); % 1 x NRep
                    else
                        FScale = sum(squeeze(sum(F2,1))'.*b_mtx_n(NVarA-WTP_space+1:end,:),1); % 1 x NRep
                    end
                end
                Xs_tmp = squeeze(Xs_sliced(1,n,:));
                FScale = FScale(ones(NVarS,1),:).*Xs_tmp(:, ones(NRep,1)); % NVarS x NRep
            end
            
            if FullCov == 0
                sumVC2tmp = bsxfun(@times, sumFsqueezed, VC2(:,:,n));
                sumVC2tmp  = min(sumVC2tmp , realmax);
                gtmp = -mean([bsxfun(@times, sumFsqueezed, U_prod); bsxfun(@times, sumVC2tmp, U_prod)],2)./p0(n);
            else % FullCov = 1
                %sumVC2tmp = sumFsqueezed(indx1,:).*VC2f(indx2,:,n);   
                sumVC2tmp = bsxfun(@times, sumFsqueezed(indx1,:), VC2f(indx2,:,n));
                sumVC2tmp  = min(sumVC2tmp , realmax);
                gtmp = -mean([bsxfun(@times, sumFsqueezed, U_prod); bsxfun(@times, sumVC2tmp, U_prod)],2)./p0(n);
            end
            
            if NVarS > 0
                gtmp = [gtmp;-mean(bsxfun(@times, FScale, U_prod),2)./p0(n)];
            end
            g(n,:) = gtmp';
            
        end
        
    else
%         save tmp1
        parfor n = 1:NP 
            YnanInd = ~isnan(YY(:,n));
            b_mtx_n = b_mtx(:,:,n);
            XXa_n = XXa(:,:,n);
%             XXa_n = XXa(:,:,n);
%             U = reshape(XXa_n(YnanInd,:)*b_mtx(:,:,n),NAltMiss(n)-1,NCTMiss(n),NRep);
%             U_sum = reshape(sum(exp(U),1),NCTMiss(n),NRep);
%             p0(n) = mean(prod(1./(1+U_sum),1));
            U = reshape(exp(XXa_n(YnanInd,:)*b_mtx_n),NAltMiss(n)-1,NCTMiss(n),NRep);  % NAlt -1 x NCT x NRep
            U_sum = reshape(1+sum(U,1),NCTMiss(n),NRep);
            U_prod = prod(1./U_sum,1);
            p0(n) = mean(U_prod);
%             p0(n) = max(mean(U_prod),realmin); 
            % calculations for gradient
            U_prob = reshape(U, (NAltMiss(n)-1)*NCTMiss(n),1,NRep);  % (NAlt-1)*NCT x 1 x NRep
            if WTP_space == 0   
                X_hat = sum(reshape(bsxfun(@times,U_prob,XXa_n(YnanInd,:)), NAltMiss(n)-1, NCTMiss(n), NVarA, NRep),1); % 1 x NCT x NVarA x NRep
                if NCTMiss(n) ~= 1
                    F = -bsxfun(@rdivide,reshape(X_hat, [NCTMiss(n),NVarA,NRep]), reshape(U_sum, [NCTMiss(n),1,NRep]));
                    sumFsqueezed = reshape(sum(F,1),[NVarA,NRep]);  %NVarA x NRep
                else
                    sumFsqueezed = -bsxfun(@rdivide,reshape(X_hat, [NVarA,NRep]), U_sum);
                end                
                sumFsqueezed(Dist==1, :) = sumFsqueezed(Dist==1, :).*b_mtx_n(Dist==1,:);                 
            else
                Xalpha = bsxfun(@times,XXa_n(YnanInd,1:end-WTP_space),reshape(b_mtx_n(WTP_matrix,:),[1,NVarA-WTP_space,NRep])); % (NAlt-1)*NCT x NVarA - WTP x NRep
                % for non-cost variables
                X_hat1 = sum(reshape(bsxfun(@times,U_prob,Xalpha), NAltMiss(n)-1, NCTMiss(n), NVarA-WTP_space, NRep),1); % 1 x NCT x NVarA-WTP_space x NRep
                F1 = -bsxfun(@rdivide,reshape(X_hat1, [NCTMiss(n),NVarA-WTP_space,NRep]), reshape(U_sum, [NCTMiss(n),1,NRep]));
                % for cost variables
                b_mtx_grad_n = b_mtx_grad(:,:,n);
                if WTP_space == 1 % without starting the loop
                    pX = bsxfun(@plus,XXa_n(YnanInd,NVarA),  XXa_n(YnanInd,1:end-WTP_space)*b_mtx_grad_n(1:end-WTP_space,:)); % (NAlt-1)*NCT x NRep
                    %pX = squeeze(XXa_n(:,NVarA,ones(NRep,1))) + XXa_n(:,1:end-WTP_space)*b_mtx_grad_n(1:end-WTP_space,:);
                    X_hat2 = sum(reshape(bsxfun(@times,squeeze(U_prob),pX), NAltMiss(n)-1, NCTMiss(n), WTP_space, NRep),1); 
                else
                    pX = zeros(NCTMiss(n)*(NAltMiss(n)-1), WTP_space, NRep);
                    for i = 1:WTP_space
                        pX(:,i,:) = bsxfun(@plus,XXa_n(YnanInd,NVarA-WTP_space+i),  XXa_n(YnanInd,WTP_matrix == NVarA-WTP_space+i)*b_mtx_grad_n(WTP_matrix == NVarA-WTP_space+i,:)); % (NAlt-1)*NCT x NRep
                    end
                    X_hat2 = sum(reshape(bsxfun(@times,U_prob,pX), NAltMiss(n)-1, NCTMiss(n), WTP_space, NRep),1); 
                end
                F2 = -bsxfun(@rdivide,reshape(X_hat2, [NCTMiss(n),WTP_space,NRep]), reshape(U_sum, [NCTMiss(n),1,NRep]));
                 
                sumFsqueezed = [squeeze(sum(F1,1));reshape(sum(F2,1),WTP_space, NRep)];  %NVarA x NRep
                sumFsqueezed(Dist==1, :) = sumFsqueezed(Dist==1, :).*b_mtx_grad_n(Dist==1,:);   
               
            end
            if NVarS >0
                if WTP_space == 0
                    FScale = sum(sumFsqueezed.*b_mtx_n,1); % 1 x NRep
                else
                    if WTP_space == 1
                        FScale = sum(squeeze(sum(F2,1)).*b_mtx_n(NVarA-WTP_space+1:end,:),1); % 1 x NRep
                    else
                        FScale = sum(squeeze(sum(F2,1))'.*b_mtx_n(NVarA-WTP_space+1:end,:),1); % 1 x NRep
                    end
                end
                Xs_tmp = squeeze(Xs_sliced(1,n,:));
                FScale = FScale(ones(NVarS,1),:).*Xs_tmp(:, ones(NRep,1)); % NVarS x NRep
            end
            
            if FullCov == 0
                sumVC2tmp = bsxfun(@times, sumFsqueezed, VC2(:,:,n));
                gtmp = -mean([bsxfun(@times, sumFsqueezed, U_prod); bsxfun(@times, sumVC2tmp, U_prod)],2)./p0(n);
            else % FullCov = 1
                %sumVC2tmp = sumFsqueezed(indx1,:).*VC2f(indx2,:,n);   
                sumVC2tmp = bsxfun(@times, sumFsqueezed(indx1,:), VC2f(indx2,:,n));
                gtmp = -mean([bsxfun(@times, sumFsqueezed, U_prod); bsxfun(@times, sumVC2tmp, U_prod)],2)./p0(n);
            end
            if NVarS > 0
                gtmp = [gtmp;-mean(bsxfun(@times, FScale, U_prod),2)./p0(n)];
            end
            g(n,:) = gtmp';
            
        end
    end;
    if NVarM > 0
       gm =  g(:,repmat(1:NVarA, 1, NVarM)).*(XXm(kron(1:NVarM, ones(1,NVarA)),:)');
       if EstimOpt.FullCov == 0
            g = [g(:,1:2*NVarA),gm, g(:,2*NVarA+1:end)];
       else
           g = [g(:,1:NVarA*(NVarA/2+1.5)),gm, g(:,NVarA*(NVarA/2+1.5)+1:end)];    
       end
    end
    
end

f = -log(p0);

% save LL_mxl
% return

% f = -log(max(p0,realmin));
