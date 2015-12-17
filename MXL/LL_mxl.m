function [f,g,h] = LL_mxl(YY,XXa,XXm,Xs,err,EstimOpt,b0)

% save LL_mxl
% return

NAlt = EstimOpt.NAlt; ...
NCT = EstimOpt.NCT; ...
NP = EstimOpt.NP; ...
NRep = EstimOpt.NRep; ...
NVarA = EstimOpt.NVarA; ...
NVarM = EstimOpt.NVarM; ...
NVarS = EstimOpt.NVarS; ...
Dist = EstimOpt.Dist; ...
WTP_space = EstimOpt.WTP_space; ...
WTP_matrix = EstimOpt.WTP_matrix; ...
FullCov = EstimOpt.FullCov; ...
DiagIndex = EstimOpt.DiagIndex; ...
Triang = EstimOpt.Triang; ...
NVarNLT = EstimOpt.NVarNLT; ...
NLTVariables = EstimOpt.NLTVariables; ...
NLTType = EstimOpt.NLTType; ...
Johnson = EstimOpt.Johnson; ...

% if FullCov == 1 && EstimOpt.NumGrad == 0
    indx1 = EstimOpt.indx1;
    indx2 = EstimOpt.indx2;
% end
if EstimOpt.ApproxHess == 0
   XXX = EstimOpt.XXX;
   VCx = EstimOpt.VCx;
end

b0a = b0(1:NVarA); ...
if any(Dist(2:end) == 3)
    b0triag_c = exp(b0a(Dist(2:end) == 3)) + Triang';
    b0a(Dist(2:end) == 3) = 0;
end
if any(Dist(2:end) == 4)
    b0weibA = exp(b0a(Dist(2:end) == 4));
    b0a(Dist(2:end) == 4) = 0;
end
if any(Dist(2:end) == 5)
    b0sinhA = b0a(Dist(2:end) == 5);
    b0a(Dist(2:end) == 5) = 0;
end

if EstimOpt.FullCov == 0 
    b0v = (b0(NVarA+1:NVarA*2)); ...
    if any(Dist(2:end) == 3)
        b0triag_b = exp(b0v(Dist(2:end) == 3)) + b0triag_c; ...
        b0v(Dist(2:end) == 3) = 1;
    end
    if any(Dist(2:end) == 4)
        b0weibB = exp(-b0v(Dist(2:end) == 4));
        b0v(Dist(2:end) == 4) = 1;
    end
    if any(Dist(2:end) == 5)
        b0sinhB = b0v(Dist(2:end) == 5).^2;
        b0v(Dist(2:end) == 5) = 1;
    end
    b0v = b0v.^2;
    VC = diag(b0v); ...
    b0m = b0(NVarA*2+1:NVarA*(NVarM+2)); ...    
    b0m = reshape(b0m,NVarA, NVarM); ...
    b0s = b0(NVarA*(NVarM+2)+1:NVarA*(NVarM+2)+NVarS); ...
    b0t = b0(NVarA*(NVarM+2)+NVarS+1:NVarA*(NVarM+2)+NVarS+NVarNLT);
    b0j = b0(NVarA*(NVarM+2)+NVarS+NVarNLT+1:NVarA*(NVarM+2)+NVarS+NVarNLT+2*Johnson); ...
else ...
	b0v = b0(NVarA+1:NVarA+sum(1:NVarA)); ...
    tmp = b0v(DiagIndex);
    b0v(DiagIndex(Dist(2:end) >=3 & Dist(2:end) <=5)) = 1;
    if any(Dist(2:end) == 3)
        b0triag_b = exp(tmp(Dist(2:end) == 3))+b0triag_c;
    end    
    if any(Dist(2:end) == 4)
        b0weibB = exp(-tmp(Dist(2:end) == 4));
    end    
    if any(Dist(2:end) == 5)
        b0sinhB = tmp(Dist(2:end) == 5).^2;
    end  
    VC = tril(ones(NVarA)); ...
    VC(VC==1) = b0v; ...
    if any(Dist(2:end) >= 3 & Dist(2:end) <= 5)
        tmp = sqrt(sum(VC(Dist(2:end) >= 3 & Dist(2:end) <= 5,:).^2,2));
        VC(Dist(2:end) >= 3 & Dist(2:end) <= 5,:) = VC(Dist(2:end) >= 3 & Dist(2:end) <= 5,:)./tmp(:, ones(1,NVarA));
    end
    b0m = b0(NVarA*(NVarA/2+1.5)+1:NVarA*(NVarA/2+1.5+NVarM)); ...
    b0m = reshape(b0m,NVarA, NVarM); ...
    b0s = b0(NVarA*(NVarA/2+1.5+NVarM)+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS); ...
    b0t = b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarNLT); ...
    b0j = b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarNLT+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarNLT+2*Johnson); ...
end 

if NVarNLT > 0
    % IndTransNon0 = (abs(bt) > 0.00001)';
    IndTransNon0 = (abs(b0t) > eps)';

    Xt = XXa(:, EstimOpt.NLTVariables,:);
    
%     bt_tmp = permute(b0t(:, ones(NAlt*NCT,1), ones(NP,1)), [2 1 3]);
    bt_tmp = b0t(:,ones(NAlt*NCT,1))';...
	bt_tmp = bt_tmp(:,:,ones(1,1,NP));

    
    
    if NLTType == 1 % BC
        Xt(:,IndTransNon0,:) = -(Xt(:,IndTransNon0,:).^bt_tmp(:,IndTransNon0,:) - 1)./bt_tmp(:,IndTransNon0,:);
        Xt(:,~IndTransNon0,:) = -log(Xt(:,~IndTransNon0,:));
    elseif NLTType == 2 % YJ
        IndXtNon0 = (Xt >= 0);
        IndXtCase1 = IndXtNon0 & IndTransNon0(ones(size(IndXtNon0,1),1),:,ones(1,1,NP)); % X >= 0, lam ~= 0
        IndXtCase2 = IndXtNon0 & ~IndTransNon0(ones(size(IndXtNon0,1),1),:,ones(1,1,NP)); % X >= 0, lam = 0
    %     IndTransNon2 = (abs(bt - 2) < 0.00001)';
        IndTransNon2 = (abs(b0t - 2) > eps)';
        IndXtCase3 = ~IndXtNon0 & IndTransNon2(ones(size(IndXtNon0,1),1),:,ones(1,1,NP));  % X < 0, lam ~= 2
        IndXtCase4 = ~IndXtNon0 & ~IndTransNon2(ones(size(IndXtNon0,1),1),:,ones(1,1,NP)); % X < 0, lam = 2
%         bt_tmp = b0t(:,ones(size(XXa,1),1))';
        Xt(IndXtCase1) = ((Xt(IndXtCase1) + 1).^bt_tmp(IndXtCase1) - 1)./bt_tmp(IndXtCase1);
        Xt(IndXtCase2) = log(Xt(IndXtCase2) + 1);
        Xt(IndXtCase3) = -((-Xt(IndXtCase3) + 1).^(2 - bt_tmp(IndXtCase3)) - 1)./(2 - bt_tmp(IndXtCase3));
        Xt(IndXtCase4) = -log(-Xt(IndXtCase4) + 1);
    end
    
    if EstimOpt.NumGrad == 0 % 
       if NLTType == 1 % BC
           XXt = XXa(:,EstimOpt.NLTVariables,:);
           XXt(:,IndTransNon0,:) = -(XXt(:,IndTransNon0,:).^bt_tmp(:,IndTransNon0,:).*(bt_tmp(:,IndTransNon0,:).*log(XXt(:, IndTransNon0,:))-1)+1)./(bt_tmp(:,IndTransNon0,:).^2);
           XXt(:,IndTransNon0 == 0,:) = -0.5*log(XXt(:, IndTransNon0 == 0)).^2;
       elseif NLTType == 2 % YJ
           XXt = XXa(:,EstimOpt.NLTVariables,:);
           XXt(IndXtCase1) = ((XXt(IndXtCase1)+1).^bt_tmp(IndXtCase1).*(bt_tmp(IndXtCase1).*log(XXt(IndXtCase1)+1)-1)+1)./(bt_tmp(IndXtCase1).^2);% X >= 0, lam ~= 0
           XXt(IndXtCase2) = 0.5*log(XXt(IndXtCase2)+1).^2;% X >= 0, lam == 0
           XXt(IndXtCase3) = -((-XXt(IndXtCase3)+1).^(2-bt_tmp(IndXtCase3)).*(1-(2-bt_tmp(IndXtCase3)).*log(-XXt(IndXtCase3)+1))-1)./((2-bt_tmp(IndXtCase3)).^2);% X < 0, lam ~= 2
           XXt(IndXtCase4) = -0.5*log(-XXt(IndXtCase4)+1).^2;% X < 0, lam == 2
       end
    end
    XXa(:, NLTVariables,:) = Xt;
else
    if EstimOpt.NumGrad == 0 % 
        XXt = [];
    end
end

b0n = b0a(:,ones(NP,1)) + b0m*XXm; ...
b0n = reshape(b0n((1:size(b0n,1))'*ones(1,NRep),(1:size(b0n,2))'),NVarA,NRep*NP); ... % NVarA x NRep*NP

b_mtx_n = b0n + VC*err; ... % NVarA x NRep*NP
    
if sum(Dist(2:end)==1) > 0 % Log - normal
    b_mtx_n(Dist(2:end)==1,:) = exp(b_mtx_n(Dist(2:end)==1,:)); ...
end
if sum(Dist(2:end)==2) > 0 % Spike
    b_mtx_n(Dist(2:end)==2,:) = max(b_mtx_n(Dist(2:end)==2,:),0); ...
end
if sum(Dist(2:end) ==3) > 0 % Triangular
    tmp = normcdf(b_mtx_n(Dist(2:end)==3,:)); ...
    Triang = Triang(ones(NRep*NP,1),:)';
    b0triag_c = b0triag_c(:, ones(NRep*NP,1));
    b0triag_b = b0triag_b(:, ones(NRep*NP,1));
    Ftriang =  (b0triag_c - Triang)./(b0triag_b- Triang);
    bmtx_triang = zeros(size(tmp));
    tmp2 = (b0triag_b- Triang).*(b0triag_c - Triang);
    bmtx_triang(tmp < Ftriang) = Triang(tmp < Ftriang)+ sqrt(tmp(tmp < Ftriang).*tmp2(tmp < Ftriang));
    tmp2 = (b0triag_b- Triang).*(b0triag_b-b0triag_c);
    bmtx_triang(tmp >= Ftriang) = b0triag_b(tmp >= Ftriang)- sqrt((1-tmp(tmp >= Ftriang)).*tmp2(tmp >= Ftriang));  
    b_mtx_n(Dist(2:end)==3,:) = bmtx_triang;
end
if sum(Dist(2:end) ==4) > 0 % Weibull
    tmp = -log(1-normcdf(b_mtx_n(Dist(2:end)==4,:))); ...
    b_mtx_n(Dist(2:end)==4,:) = b0weibA(:, ones(1,NP*NRep,1)).*(tmp.^b0weibB(:, ones(1,NP*NRep,1)));
end
if sum(Dist(2:end)>=5) > 0 % Johnson
    if sum(Dist(2:end)==5) > 0 % Sinh-Arcsinh
        b_mtx_n(Dist(2:end)==5,:) = b0sinhA(:,ones(NRep*NP,1))+ b0sinhB(:,ones(NRep*NP,1)).*asinh(b_mtx_n(Dist(2:end)==5,:));

        b_mtx_n(Dist(2:end)==5,:) = b0j(1:Johnson, ones(NRep*NP,1)) + exp(b0j(Johnson+1:end, ones(NRep*NP,1))).*sinh(b_mtx_n(Dist(2:end)==5,:));  
    end
    if sum(Dist(2:end)==6) > 0 % Johnson Sb
        tmp = exp(b_mtx_n(Dist(2:end) ==6,:));
        b_mtx_n(Dist(2:end)==6,:) = tmp./(1+tmp); ...
        b_mtx_n(Dist(2:end) == 6,:) = b0j(1:Johnson, ones(NRep*NP,1)) + exp(b0j(Johnson+1:end, ones(NRep*NP,1))).*b_mtx_n(Dist(2:end)==6,:);
    end
    if sum(Dist(2:end)==7) > 0 % Johnson Su
        b_mtx_n(Dist(2:end)==7,:) = sinh(b_mtx_n(Dist(2:end)==7,:)); ...
        b_mtx_n(Dist(2:end) ==7,:) = b0j(1:Johnson, ones(NRep*NP,1)) + exp(b0j(Johnson+1:end, ones(NRep*NP,1))).*b_mtx_n(Dist(2:end)==7,:);
    end

end

%sum(isnan(b_mtx_n),2)
if WTP_space > 0
    if EstimOpt.NumGrad == 0
       b_mtx_grad = reshape(b_mtx_n ,NVarA,NRep,NP); ...% needed for gradient calculation in WTP_space
    end
    b_mtx_n(1:end-WTP_space,:) = b_mtx_n(1:end-WTP_space,:).*b_mtx_n(WTP_matrix,:); ...
else
    b_mtx_grad = [];
end

cs = reshape(exp(Xs*b0s),NAlt*NCT,1,NP); ...
XXa_n = XXa .* cs(:,ones(1,NVarA,1),:); ...

b_mtx_n = reshape(b_mtx_n,NVarA,NRep,NP); ...

p0 = zeros(NP,1); ...

if EstimOpt.NumGrad == 1 % with numerical gradient
    if any(isnan(XXa(:))) == 0 % faster version for complete dataset
        for n = 1:NP             
            U = reshape(XXa_n(:,:,n)*b_mtx_n(:,:,n),NAlt,NCT,NRep); ...
            U_max = max(U); ...
            U = exp(U - U_max(ones(NAlt,1),:,:)); ... % rescale utility to avoid exploding
            U_sum = reshape(sum(U,1),NCT,NRep); ...
            U_selected = reshape(U(YY(:,n*ones(NRep,1))==1),NCT,NRep); ...   
            p0(n) = mean(prod(U_selected ./ U_sum,1));
        end; ...
    else ... % this works only if NAlt is constant for each respondent and if missing ALT is not the first in NCT
        for n = 1:NP             
%             U = reshape(XXa_n(~isnan(YY(:,n)),:,n)*b_mtx_n(:,:,n),NAlt,NCT-sum(isnan(YY(1:NAlt:end,n))),NRep);... % this would be faster if there are no ALT missing
            U = reshape(XXa_n(~isnan(YY(:,n)),:,n)*b_mtx_n(:,:,n),numel(YY(~isnan(YY(:,n)),n))./(NCT-sum(isnan(YY(1:NAlt:end,n)))),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep); ...
            U_max = max(U); ...
%             U = exp(U - U_max(ones(NAlt,1),:,:)); ...
            U = exp(U - U_max(ones(numel(YY(~isnan(YY(:,n)),n))./(NCT-sum(isnan(YY(1:NAlt:end,n)))),1),:,:)); ...
            U_sum = reshape(nansum(U,1),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep); ...
            U_selected = reshape(U(YY(~isnan(YY(:,n)),n*ones(NRep,1))==1),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep); ...   
            p0(n) = mean(prod(U_selected ./ U_sum,1));
        end; ...
    end ...

% na razie sobie to zostawi³em, jako notatkê

    % szybsze ale podatne na eksplozje:
    % for n = 1:NP
    %     U = exp(XXa_n(:,:,n)*b_mtx_n(:,:,n)); ...
    %     U_selected = reshape(U(YY2(:,n*ones(NRep,1)) == 1),NCT,NRep); ...
    %     U_sum = reshape((sum(reshape(U,NAlt,NCT,NRep),1)),NCT,NRep); ...
    %     p0(n) = mean(prod(U_selected ./ U_sum),2);
    % %     U_selected = U(YY2(:,n*ones(NRep,1)) == 1);...
    % %     U_sum = sum(reshape(U,NAlt,NCT*NRep),1)';...
    % %     p0(n) = mean(prod(reshape(U_selected ./ U_sum,NCT,NRep)),2);...
    % end;

    % wyci¹ganie rzeczy z pêtli niekoniecznie pomaga:
    % YY3 = reshape(YY,NAlt*NCT,1,NP);
    % tic ; for n = 1:NP
    %     U(:,:,n) = exp(XXa_n(:,:,n)*b_mtx_n(:,:,n));
    % end;
    % U_selected = reshape(U(YY3(:,ones(NRep,1),:) == 1),NCT,NRep,NP); ...
    % U_sum = reshape((sum(reshape(U,NAlt,NCT,NRep,NP),1)),NCT,NRep,NP); ...
    % p0 = reshape(mean(prod(U_selected ./ U_sum),2),NP,1); toc

    
elseif EstimOpt.NumGrad == 0 && EstimOpt.ApproxHess == 1 % analitical gradient + approximated hessian

%     save tmp2
%     return
    
	if FullCov == 0    
%         g = zeros(NP, 2*NVarA + 1);
        g = zeros(NP, 2*NVarA + NVarNLT);
        VC2 = reshape(2*diag(b0(NVarA+1:NVarA*2))*err, NVarA, NRep, NP);
    else
%         g = zeros(NP, 2*NVarA+NVarA*(NVarA-1)/2+1);
        g = zeros(NP, 2*NVarA+NVarA*(NVarA-1)/2+NVarNLT);
        VC2 = reshape(err, NVarA, NRep, NP);       
    end
    
    if any(isnan(XXa(:))) == 0 % faster version for complete dataset

%         Xalpha0 = bsxfun(@times,reshape(XXa(:,1:end-WTP_space,:),[NCT*NAlt,NVarA-WTP_space,1,NP]),reshape(b_mtx_n(WTP_matrix,:,:),[1,NVarA-WTP_space,NRep,NP]));
        parfor n = 1:NP
            U = reshape(XXa_n(:,:,n)*b_mtx_n(:,:,n),NAlt,NCT,NRep); ... % NAlt x NCT x NRep
            U_max = max(U); ...
            U = exp(U - U_max(ones(NAlt,1),:,:)); ... % rescale utility to avoid exploding 
            U_sum = reshape(sum(U,1),1,NCT,NRep); ... 
            U_prob = U./U_sum(ones(NAlt,1,1),:,:); ... % NAlt x NCT x NRep
            U_prod = prod(reshape(U_prob(YY(:,n*ones(NRep,1))==1),NCT,NRep),1); ... % 1 x NRep
            %p0(n) = mean(U_prod);
            p0(n) = max(mean(U_prod),realmin); ...
                
            % calculations for gradient
            U_prob = reshape(U_prob, NAlt*NCT,1, NRep); ... % NAlt*NCT x NVarA x NRep
            if WTP_space == 0   
                X_hat = sum(reshape(U_prob(:,ones(1,NVarA,1),:).* XXa(:,:, n*ones(NRep,1)), NAlt, NCT, NVarA, NRep),1); ...
                if NCT ~= 1
                    F = XXa(YY(:,n) == 1,:,n*ones(NRep,1)) - squeeze(X_hat); ... %NCT x NVarA x NRep 
                    sumFsqueezed = squeeze(sum(F,1)); ... %NVarA x NRep
                else
                    sumFsqueezed = squeeze(XXa(YY(:,n) == 1,:,n*ones(NRep,1))) - squeeze(X_hat); %NVarA x NRep 
                end                
                sumFsqueezed(Dist(2:end)==1, :) = sumFsqueezed(Dist(2:end)==1, :).*b_mtx_n(Dist(2:end)==1,:,n);
                
                if NVarNLT == 1
                    XXtt = XXt(:,:,n)*b_mtx_n(NLTVariables,:,n); %NAlt*NCT x NRep
                    X_hat_lam = sum(reshape(squeeze(U_prob(:,1,:)).* XXtt, NAlt, NCT, NRep),1);
                    F3 =   XXtt(YY(:,n) == 1,:)- squeeze(X_hat_lam); % CT x NRep
                    F3sum = sum(F3,1); % 1  x NRep
                elseif NVarNLT > 1

                    XXtt = XXt(:,:,n*ones(NRep,1)).*permute(b_mtx_n(NLTVariables,:,n*ones(NCT*NAlt,1)), [3 1 2]); %NAlt*NCT x NVarNLT x NRep
%                     X_hat_lam = sum(reshape(U_prob(:,ones(1,NVarNLT,1),:).* XXtt, NAlt, NCT, NRep),1);
                    X_hat_lam = sum(reshape(U_prob(:,ones(1,NVarNLT,1),:).* XXtt, NAlt, NCT,NVarNLT,NRep),1);
                    F3 =  XXtt(YY(:,n) == 1,:,:)- squeeze(X_hat_lam); % CT x NVarNLT x NRep
                    F3sum = squeeze(sum(F3,1)); % NVarNLT  x NRep

                end
            else
%                 b_mtx_wtp = reshape(b_mtx_n(:,:,n), 1, NVarA, NRep);...
%                 Xalpha = XXa(:,1:end-WTP_space, n*ones(NRep,1)).*b_mtx_wtp(ones(NAlt*NCT,1),WTP_matrix,:); 
                Xalpha = bsxfun(@times,XXa(:,1:end-WTP_space, n),reshape(b_mtx_n(WTP_matrix,:,n),[1,NVarA-WTP_space,NRep]));
                % for non-cost variables
                X_hat1 = sum(reshape(U_prob(:,ones(1,NVarA-WTP_space),:).* Xalpha, NAlt, NCT, NVarA-WTP_space, NRep),1); ...
                F1 = Xalpha(YY(:,n) == 1,:,:) - squeeze(X_hat1); ... %NCT x NVarA-WTP_space x NRep   
                % for cost variables
                if WTP_space == 1 % without starting the loop
                    pX = squeeze(XXa(:,NVarA, n*ones(NRep,1))) + XXa(:,1:end-WTP_space,n)*b_mtx_grad(1:end-WTP_space,:,n);
                    X_hat2 = sum(reshape(squeeze(U_prob).*pX, NAlt, NCT, WTP_space, NRep),1); ...
                else
                    pX = zeros(NCT*NAlt, WTP_space, NRep);
                    for i = 1:WTP_space
                        pX(:,i,:) = squeeze(XXa(:,NVarA-WTP_space+i, n*ones(NRep,1))) + XXa(:,WTP_matrix == NVarA-WTP_space+i,n)*b_mtx_grad(WTP_matrix == NVarA-WTP_space+i,:,n);
                    end
                    X_hat2 = sum(reshape(U_prob(:,ones(1,WTP_space),:).* pX, NAlt, NCT, WTP_space, NRep),1); ...
                end
                F2 = pX(YY(:,n) == 1,:,:) - squeeze(X_hat2); ... %NCT x WTP_space x NRep       
                sumFsqueezed = [squeeze(sum(F1,1));squeeze(sum(F2,1)) ]; ... %NVarA x NRep
                sumFsqueezed(Dist(2:end)==1, :) = sumFsqueezed(Dist(2:end)==1, :).*b_mtx_grad(Dist(2:end)==1,:,n);
                if NVarNLT == 1
                    XXtt = XXt(:,:,n)*b_mtx_n(NLTVariables,:,n); %NAlt*NCT x NRep
                    X_hat_lam = sum(reshape(squeeze(U_prob(:,1,:)).* XXtt, NAlt, NCT, NRep),1);
                    F3 = XXtt(YY(:,n) == 1,:)- squeeze(X_hat_lam) ; % CT x NRep
                    F3sum = sum(F3,1); % 1  x NRep

                elseif NVarNLT > 1
                    XXtt = XXt(:,:,n*ones(NRep,1)).*permute(b_mtx_n(NLTVariables,:,n*ones(NCT*NAlt,1)), [3 1 2]); %NAlt*NCT x NVarNLT x NRep
                    X_hat_lam = sum(reshape(U_prob(:,ones(NVarNLT,1),:).* XXtt, NAlt, NCT, NVarNLT, NRep),1);
                    F3 =   XXtt(YY(:,n) == 1,:,:) - squeeze(X_hat_lam); % CT x NVarNLT x NRep
                    F3sum = squeeze(sum(F3,1)); % NVarNLT  x NRep
                end
            end
            
            
            if FullCov == 0
                sumVC2tmp = sumFsqueezed.*VC2(:,:,n); ... % NVarA x NRep
                if NVarNLT == 0
                    g(n,:) = -mean([sumFsqueezed.*U_prod(ones(NVarA,1),:); sumVC2tmp.*U_prod(ones(NVarA,1),:)],2)./p0(n);
                else
                    g(n,:) = (-mean([sumFsqueezed.*U_prod(ones(NVarA,1),:); sumVC2tmp.*U_prod(ones(NVarA,1),:); F3sum.*U_prod(ones(NVarNLT,1),:)],2)./p0(n))';
                end
            else % FullCov = 1
                sumVC2tmp = sumFsqueezed(indx1,:).*VC2(indx2,:,n); ...  
                if NVarNLT == 0
                    g(n,:) = -mean([sumFsqueezed.*U_prod(ones(NVarA,1),:); sumVC2tmp.*U_prod(ones(NVarA*(NVarA-1)/2+NVarA,1),:)],2)./p0(n);    
                else
                    g(n,:) = -mean([sumFsqueezed.*U_prod(ones(NVarA,1),:); sumVC2tmp.*U_prod(ones(NVarA*(NVarA-1)/2+NVarA,1),:); F3sum.*U_prod(ones(NVarNLT,1),:) ],2)./p0(n);    
                end
                
            end     
            
            
        end; ...
    else ...
        parfor n = 1:NP ...
            U = reshape(XXa_n(~isnan(YY(:,n)),:,n)*b_mtx_n(:,:,n),NAlt,NCT-sum(isnan(YY(1:NAlt:end,n))),NRep); ...
            U_max = max(U); ...
            U = exp(U - U_max(ones(NAlt,1),:,:)); ... % rescale utility to avoid exploding 
            U_sum = reshape(nansum(U,1),1,NCT-sum(isnan(YY(1:NAlt:end,n))),NRep); ...
            U_prob = U./U_sum(ones(NAlt,1,1),:,:); ... % NAlt x NCT x NRep
            U_prod = prod(reshape(U_prob(YY(~isnan(YY(:,n)),n*ones(NRep,1))==1),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep),1); ... % 1 x NRep
            %p0(n) = mean(U_prod);
            p0(n) = max(mean(U_prod),realmin); ...
                
            % calculations for gradient
            
            U_prob = reshape(U_prob, NAlt*(NCT-sum(isnan(YY(1:NAlt:end,n)))),1, NRep); ... % NAlt*NCT x NVarA x NRep   
            if WTP_space == 0
                X_hat = sum(reshape(U_prob(:,ones(1,NVarA),:).* XXa(~isnan(YY(:,n)),:, n*ones(NRep,1)), NAlt, NCT-sum(isnan(YY(1:NAlt:end,n))), NVarA, NRep),1); ...
                F = XXa(YY(:,n) == 1,:,n*ones(NRep,1)) - squeeze(X_hat); ... %NCT x NVarA x NRep    
                sumFsqueezed = squeeze(sum(F,1)); ... %NVarA x NRep

                sumFsqueezed(Dist(2:end)==1, :) = sumFsqueezed(Dist(2:end)==1, :).*b_mtx_n(Dist(2:end)==1,:,n);
            else
                b_mtx_wtp = reshape(b_mtx_n(:,:,n), 1, NVarA, NRep);
                Xalpha = XXa(~isnan(YY(:,n)),1:end-WTP_space, n*ones(NRep,1)).*b_mtx_wtp(ones(NAlt*(NCT- sum(isnan(YY(1:NAlt:end,n)))),1),WTP_matrix,:);
                % for non-cost variables
                X_hat1 = sum(reshape(U_prob(:,ones(1,NVarA-WTP_space),:).* Xalpha, NAlt, NCT- sum(isnan(YY(1:NAlt:end,n))), NVarA-WTP_space, NRep),1); ...
                F1 = Xalpha(YY(~isnan(YY(:,n)),n) == 1,:,:) - squeeze(X_hat1); ... %NCT x NVarA-WTP_space x NRep   
                % for cost variables
                if WTP_space == 1 % without starting the loop
                    pX = squeeze(XXa(~isnan(YY(:,n)),NVarA, n*ones(NRep,1))) + XXa(~isnan(YY(:,n)),1:end-WTP_space,n)*b_mtx_grad(1:end-WTP_space,:,n);
                    X_hat2 = sum(reshape(squeeze(U_prob).*pX, NAlt, NCT-sum(isnan(YY(1:NAlt:end,n))), WTP_space, NRep),1); ...
                else
                    pX = zeros((NCT- sum(isnan(YY(1:NAlt:end,n))))*NAlt, WTP_space, NRep);
                    for i = 1:WTP_space
                        pX(:,i,:) = squeeze(XXa(~isnan(YY(:,n)),NVarA-WTP_space+i, n*ones(NRep,1))) + XXa(~isnan(YY(:,n)),WTP_matrix == NVarA-WTP_space+i,n)*b_mtx_grad(WTP_matrix == NVarA-WTP_space+i,:,n);
                    end
                    X_hat2 = sum(reshape(U_prob(:,ones(1,WTP_space),:).* pX, NAlt, sum(isnan(YY(1:NAlt:end,n))), WTP_space, NRep),1); ...
                end
                F2 = pX(YY(~isnan(YY(:,n)),n) == 1,:,:) - squeeze(X_hat2); ... %NCT x WTP_space x NRep       
                sumFsqueezed = [squeeze(sum(F1,1));squeeze(sum(F2,1)) ]; ... %NVarA x NRep
                sumFsqueezed(Dist(2:end)==1, :) = sumFsqueezed(Dist(2:end)==1, :).*b_mtx_grad(Dist(2:end)==1,:,n);
            end
            if FullCov == 0
                sumVC2tmp = sumFsqueezed.*VC2(:,:,n); ... % NVarA x NRep
                g(n,:) = -mean([sumFsqueezed.*U_prod(ones(NVarA,1),:); sumVC2tmp.*U_prod(ones(NVarA,1),:)],2)./p0(n);
            else % FullCov = 1
                sumVC2tmp = sumFsqueezed(indx1,:).*VC2(indx2,:,n); ...   
                g(n,:) = -mean([sumFsqueezed.*U_prod(ones(NVarA,1),:); sumVC2tmp.*U_prod(ones(NVarA*(NVarA-1)/2+NVarA,1),:)],2)./p0(n);    
            end
        end; 
    end;
    if NVarM > 0
       gm =  g(:,repmat(1:NVarA, 1, NVarM)).*(XXm(kron(1:NVarM, ones(1,NVarA)),:)');
       g = [g,gm];
    end

elseif EstimOpt.NumGrad == 0 && EstimOpt.ApproxHess == 0 % Analitical gradient and hesjan    
    
    if FullCov == 0    
        g = zeros(NP, 2*NVarA);
        VC2 = reshape(2*diag(b0(NVarA+1:NVarA*2))*err, NVarA, NRep, NP);
        SIG = 4*b0(NVarA+1:NVarA*2)*b0(NVarA+1:NVarA*2)';
        SIG = SIG(:,:, ones(NRep,1));
        VC3 = reshape(2*err, NVarA, NRep, NP);
        NVarG = 2*NVarA; % general no. of parameters
    else
        g = zeros(NP, 2*NVarA+NVarA*(NVarA-1)/2);
        VC2 = reshape(err, NVarA, NRep, NP);
        VC3 = []; % for parfor
        SIG = []; % for parfor
        NVarG = 2*NVarA+NVarA*(NVarA-1)/2;
    end
    hx1 = zeros(NVarG, NVarG, NP);
    hx2 = zeros(NVarG, NVarG, NP);
    if any(isnan(XXa(:))) == 0 % faster version for complete dataset
        parfor n = 1:NP
            U = reshape(XXa_n(:,:,n)*b_mtx_n(:,:,n),NAlt,NCT,NRep); ... % NAlt x NCT x NRep
            U_max = max(U); ...
            U = exp(U - U_max(ones(NAlt,1),:,:)); ... % rescale utility to avoid exploding 
            U_sum = reshape(sum(U,1),1,NCT,NRep); ... 
            U_prob = U./U_sum(ones(NAlt,1,1),:,:); ... % NAlt x NCT x NRep
            U_prod = prod(reshape(U_prob(YY(:,n*ones(NRep,1))==1),NCT,NRep),1); ... % 1 x NRep
            %p0(n) = mean(U_prod);
%             dlaczego max(x,0)? czy to nie powoduje problemów póŸniej przy dzieleniu przez 0? 
%             chodzi³o o realmin - mój b³¹d
            p0(n) = max(mean(U_prod),realmin); ...
                
            % calculations for gradient
            U_prob = reshape(U_prob, NAlt*NCT,1, NRep); ... % NAlt*NCT x 1 x NRep
            X_hat = sum(reshape(U_prob(:,ones(1,NVarA,1),:).* XXa(:,:, n*ones(NRep,1)), NAlt, NCT, NVarA, NRep),1); ... % 1 x NCT x NVarA x NRep
            F = XXa(YY(:,n) == 1,:,n*ones(NRep,1)) - squeeze(X_hat); ... %NCT x NVarA x NRep    
            sumFsqueezed = squeeze(sum(F,1)); ... %NVarA x NRep
            if FullCov == 0
                sumVC2tmp = sumFsqueezed.*VC2(:,:,n); ... % NVarA x NRep
                g(n,:) = -mean([sumFsqueezed.*U_prod(ones(NVarA,1),:); sumVC2tmp.*U_prod(ones(NVarA,1),:)],2)./p0(n);
            else % FullCov = 1
                sumVC2tmp = sumFsqueezed(indx1,:).*VC2(indx2,:,n); ...   
                g(n,:) = -mean([sumFsqueezed.*U_prod(ones(NVarA,1),:); sumVC2tmp.*U_prod(ones(NVarA*(NVarA-1)/2+NVarA,1),:)],2)./p0(n);    
            end
            % calculations for hessian
            %second term of hessian
            gtmp = [sumFsqueezed; sumVC2tmp];   
            hx1(:,:,n) = ((gtmp.*U_prod(ones(NVarG,1),:))*gtmp')./(p0(n)*NRep); %NVarG x NVarG
            
            %third term of hessian
            U_prod = reshape(U_prod,1,1,NRep);
            U_prob = reshape(U_prob,NAlt*NCT,1,1,NRep);
            X_hat_tmp = reshape(permute(squeeze(X_hat(1,:,:,:)).*U_prod(ones(NCT,1), ones(NVarA,1),:), [1 3 2]), NCT*NRep, NVarA);
            X_hat = reshape(permute(X_hat, [2 4 3 1]), NCT*NRep, NVarA);
            
            hbbx = sum(permute(-sum(reshape(U_prob(:, ones(NVarA,1), ones(NVarA,1), :).*XXX(:, :,:, n*ones(NRep,1)), NAlt,NCT, NVarA, NVarA, NRep),1), [3 4 5 2 1]),4);
            hbb = (mean(hbbx.*U_prod(ones(NVarA,1), ones(NVarA,1),:),3)+(X_hat_tmp')*X_hat/NRep)./p0(n);
            % partial derivative of standard deviations and b
            
            if FullCov == 0
                VC2_n = reshape(VC2(:,:,n), 1, NVarA, NRep);
                VC2_n_tmp = reshape(permute(VC2_n(ones(NCT,1),:,:), [1 3 2]), NCT*NRep, NVarA);

                hsb = (mean(hbbx.*U_prod(ones(NVarA,1), ones(NVarA,1),:).*VC2_n(ones(NVarA,1),:,:),3)+ (X_hat_tmp)'*(X_hat.*VC2_n_tmp)/NRep)./p0(n);
                VCxx  = diag(ones(NVarA,1));
                VCxx = VCxx(:,:, ones(NRep,1));
                VCxx(VCxx == 1) = sumFsqueezed.*VC3(:,:,n); 

                VCx_n = VCx(:,:,:,n).*SIG;
                hss = (mean(U_prod(ones(NVarA,1), ones(NVarA,1),:).*(hbbx.*VCx_n +VCxx),3)+(X_hat_tmp.*VC2_n_tmp)'*(X_hat.*VC2_n_tmp)/NRep)./p0(n);
            else
                VC2_n = reshape(VC2(indx2,:,n), 1, NVarA*(NVarA-1)/2+NVarA, NRep);
                VC2_n_tmp = reshape(permute(VC2_n(ones(NCT,1),:,:), [1 3 2]), NCT*NRep, NVarA*(NVarA-1)/2+NVarA);
                hsb = (mean(hbbx(:, indx1, :).*U_prod(ones(NVarA,1), ones(NVarA*(NVarA-1)/2+NVarA,1),:).*VC2_n(ones(NVarA,1),:,:),3)+(X_hat_tmp)'*(X_hat(:, indx1).*VC2_n_tmp)/NRep)./p0(n);

                hss = (mean(U_prod(ones(NVarA*(NVarA-1)/2+NVarA,1), ones(NVarA*(NVarA-1)/2+NVarA,1),:).*hbbx(indx1, indx1, :).*VCx(:,:,:,n),3)+(X_hat_tmp(:, indx1).*VC2_n_tmp)'*(X_hat(:, indx1).*VC2_n_tmp)/NRep)./p0(n);
            end

            hx2(:,:,n) = [hbb, hsb; hsb', hss];
        end; ...

    else ...
        parfor n = 1:NP ...
            U = reshape(XXa_n(~isnan(YY(:,n)),:,n)*b_mtx_n(:,:,n),NAlt,NCT-sum(isnan(YY(1:NAlt:end,n))),NRep); ...
            U_max = max(U); ...
            U = exp(U - U_max(ones(NAlt,1),:,:)); ... % rescale utility to avoid exploding 
            U_sum = reshape(nansum(U,1),1,NCT-sum(isnan(YY(1:NAlt:end,n))),NRep); ...
            U_prob = U./U_sum(ones(NAlt,1,1),:,:); ... % NAlt x NCT x NRep
            U_prod = prod(reshape(U_prob(YY(~isnan(YY(:,n)),n*ones(NRep,1))==1),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep),1); ... % 1 x NRep
            %p0(n) = mean(U_prod);
            p0(n) = max(mean(U_prod),realmin); ...
                
            % calculations for gradient
            U_prob = reshape(U_prob, NAlt*(NCT-sum(isnan(YY(1:NAlt:end,n)))),1, NRep); ... % NAlt*NCT x NVarA x NRep   
            X_hat = sum(reshape(U_prob(:,ones(1,NVarA,1),:).* XXa(:,:, n*ones(NRep,1)), NAlt, NCT-sum(isnan(YY(1:NAlt:end,n))), NVarA, NRep),1); ...
            F = XXa(YY(:,n) == 1,:,n*ones(NRep,1)) - squeeze(X_hat); ... %NCT x NVarA x NRep    
            sumFsqueezed = squeeze(sum(F,1)); ... %NVarA x NRep
            if FullCov == 0
                sumVC2tmp = sumFsqueezed.*VC2(:,:,n); ... % NVarA x NRep
                g(n,:) = -mean([sumFsqueezed.*U_prod(ones(NVarA,1),:); sumVC2tmp.*U_prod(ones(NVarA,1),:)],2)./p0(n);
            else % FullCov = 1
                sumVC2tmp = sumFsqueezed(indx1,:).*VC2(indx2,:,n); ...   
                g(n,:) = -mean([sumFsqueezed.*U_prod(ones(NVarA,1),:); sumVC2tmp.*U_prod(ones(NVarA*(NVarA-1)/2+NVarA,1),:)],2)./p0(n);    
            end
        end; 
    end;
%     if NVarM > 0
%        gm =  g(:,repmat(1:NVarA, 1, NVarM)).*(XXm(kron(1:NVarM, ones(1,NVarA)),:)');
%        g = [g,gm];
%     end
    h = g'*g - sum(hx1 + hx2,3);
end

f = -log(p0);
% f = -log(max(p0,realmin));
