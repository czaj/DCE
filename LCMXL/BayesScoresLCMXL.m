function [PClass_i,Betas_i] = BayesScoresLCMXL(YY,Xa,Xc,varargin)

% save scores1_tmp

if numel(varargin) == 3
    Xs = [];
    err_sliced = varargin{1};
    EstimOpt = varargin{2};
    B = varargin{3};
elseif numel(varargin) == 4
    Xs = varargin{1};
    err_sliced = varargin{2};
    EstimOpt = varargin{3};
    B = varargin{4};
else
    error('BayesScoresLCMXL: invalid number of inputs');
end

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
NVarS = EstimOpt.NVarS;

b_mtx = zeros(NVarA,NP*NRep,NClass);  
VC = zeros(NVarA); 
if FullCov == 0 
    for c = 1:NClass
        b_mtx(:,:,c) = B((c-1)*NVarA+1:c*NVarA,ones(NP*NRep,1)) + diag(B(NVarA*NClass+(c-1)*NVarA+1:NVarA*(NClass+c)))*err_sliced((c-1)*NVarA+1:c*NVarA,:); 
        if sum(Dist(:,c) == 1) > 0  % lognormal
            b_mtx(Dist(:,c) == 1,:,c) = exp(b_mtx(Dist(:,c) == 1,:,c)); 
        end
        if sum(Dist(:,c) == 2) > 0  % Spike
            b_mtx(Dist(:,c) == 2,:,c) = max(b_mtx(Dist(:,c) == 2,:,c),0); 
        end
        if WTP_space > 0 
            b_mtx(1:end-WTP_space,:,c) = b_mtx(1:end-WTP_space,:,c).*b_mtx(WTP_matrix,:,c); 
        end
    end
else 
    VC_tmp = tril(ones(NVarA));
    for c = 1:NClass
        VC(VC_tmp == 1) = B(NVarA*NClass+(c-1)*sum(1:NVarA)+1:NVarA*NClass+c*sum(1:NVarA));
        b_mtx(:,:,c) = B((c-1)*NVarA+1:c*NVarA,ones(NP*NRep,1)) + VC*err_sliced((c-1)*NVarA+1:c*NVarA,:); 
        if sum(Dist(:,c) == 1) > 0 % lognormal
            b_mtx(Dist(:,c) == 1,:,c) = exp(b_mtx(Dist(:,c) == 1,:,c)); 
        end
        if sum(Dist(:,c) == 2) > 0 % Spike
            b_mtx(Dist(:,c) == 2,:,c) = max(b_mtx(Dist(:,c) == 2,:,c),0); 
        end
        if WTP_space > 0
            b_mtx(1:end-WTP_space,:,c) = b_mtx(1:end-WTP_space,:,c).*b_mtx(WTP_matrix,:,c); 
        end
    end
end

%     Dist = reshape(EstimOpt.Dist,[NVarA,1,NClass]);
%     b_mtx(Dist(:,ones(1,NP*NRep,1),:)==1) = exp(b_mtx(Dist(:,ones(1,NP*NRep,1),:)==1)); % this is slower  

if NVarS > 0
    if isempty(Xs)
        error('BayesScoresLCMXL: Xs is required when EstimOpt.NVarS > 0.');
    end
    if FullCov == 0
        bs = reshape(B(2*NClass*NVarA+1:(2*NVarA+NVarS)*NClass),[NVarS,NClass]);
    else
        bs = reshape(B(NClass*(NVarA+sum(1:NVarA))+1:NClass*(NVarA+NVarS+sum(1:NVarA))),[NVarS,NClass]);
    end
    Scale = reshape(exp(Xs*bs),[NAlt,NCT,NP,NClass]);
else
    Scale = [];
end

p = zeros(NP,NClass); 
p_tmp = zeros(NP,NRep,NClass);

if any(isnan(Xa(:))) == 0  % faster version for complete dataset      
	for n = 1:NP 
        YYn = YY(:,n*ones(NRep,1),ones(1,1,NClass));
        U = reshape(Xa(:,:,n)*reshape(b_mtx(:,((n-1)*NRep+1):n*NRep,:),[NVarA,NRep*NClass]),[NAlt,NCT,NRep,NClass]);
        if NVarS > 0
            U = U.*Scale(:,:,n,:);
        end
        U = exp(U - max(U));
        U_sum = reshape(sum(U,1),[NCT,NRep,NClass]);
        U_selected = reshape(U(YYn == 1),[NCT,NRep,NClass]);
%         p(n,:) = mean(prod(U_selected ./ U_sum,1),2); 
        p_tmp(n,:,:) = prod(U_selected./U_sum,1); % NP x NRep x NClass
        p(n,:) = mean(p_tmp(n,:,:),2); % likelihoods, NP x NClass
    end
else 
    for n = 1:NP
        YYn = YY(:,n*ones(NRep,1),ones(1,1,NClass));  %it might be done once, not for every class    
        U = reshape(Xa(~isnan(YY(:,n)),:,n)*reshape(b_mtx(:,((n-1)*NRep+1):n*NRep,:),[NVarA,NRep*NClass]),[NAlt,NCT-sum(isnan(YY(1:NAlt:end,n))),NRep,NClass]);
        if NVarS > 0
            Scale_n = reshape(Scale(:,:,n,:),[NAlt*NCT,1,NClass]);
            Scale_n = reshape(Scale_n(~isnan(YY(:,n)),:,:),[NAlt,NCT-sum(isnan(YY(1:NAlt:end,n))),1,NClass]);
            U = U.*Scale_n;
        end
        U = exp(U - max(U)); ... % NAlt x NCT - NaNs x NRep x NClass
        U_sum = reshape(nansum(U,1),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep,NClass);
        U_selected = reshape(U(YYn(~isnan(YY(:,n)),:) == 1),[NCT-sum(isnan(YY(1:NAlt:end,n))),NRep,NClass]);
%         p(n,:) = mean(prod(U_selected ./ U_sum,1),2); 
        p_tmp(n,:,:) = prod(U_selected./U_sum,1);
        p(n,:) = mean(p_tmp(n,:,:),2); 
    end
end 

if FullCov == 0
    PClass = Xc*reshape([B(NClass*(2*NVarA+NVarS)+1:end);zeros(NVarC,1)],[NVarC,NClass]); 
else 
    PClass = Xc*reshape([B(NClass*(NVarA+NVarS+sum(1:NVarA))+1:end);zeros(NVarC,1)],[NVarC,NClass]); 
end
PClass = exp(PClass - max(PClass,[],2));
PClass = PClass./sum(PClass,2); ... % class probabilities (prior)
PClass_i = reshape(p.*PClass./max(sum(p.*PClass,2),realmin),[1,NP,NClass]); %class probabilities (posterior)

% b_mtx_sliced = squeeze(mean(reshape(b_mtx,NVarA,NRep,NP,NClass),2)); ...
% Betas_i = sum(b_mtx_sliced .* PClass_i(ones(NVarA,1,1),:,:),3)'; ...

% why was it here? seems scores are not right with it. 
% % get WTP-space parameters (clean): 
% if WTP_space > 0
%     b_mtx = zeros(NVarA,NP*NRep,NClass);  
%     VC = zeros(NVarA); 
%     if FullCov == 0
%         for c = 1:NClass
%             b_mtx(:,:,c) = B((c-1)*NVarA+1:c*NVarA,:) + diag(B(NVarA*NClass+(c-1)*NVarA+1:NVarA*(NClass+c)).^2)*err_sliced((c-1)*NVarA+1:c*NVarA,:); 
%             if sum(Dist(:,c) == 1) > 0  % lognormal
%                 b_mtx(Dist(:,c) == 1,:,c) = exp(b_mtx(Dist(:,c) == 1,:,c)); 
%             elseif sum(Dist(:,c) == 2) > 0  % Spike
%                 b_mtx(Dist(:,c) == 2,:,c) = max(b_mtx(Dist(:,c) == 1,:,c),0);   
%             end
% %             if WTP_space > 0; 
% %                 b_mtx(1:end-WTP_space,:,c) = b_mtx(1:end-WTP_space,:,c).*b_mtx(WTP_matrix,:,c); 
% %             end; 
%         end
%     else 
%         VC_tmp = tril(ones(NVarA));
%         for c = 1:NClass
%             VC(VC_tmp == 1) = B(NVarA*NClass+(c-1)*sum(1:NVarA)+1:NVarA*NClass+c*sum(1:NVarA));
%             b_mtx(:,:,c) = B((c-1)*NVarA+1:c*NVarA,ones(NP*NRep,1)) + VC*err_sliced((c-1)*NVarA+1:c*NVarA,:); 
%             if sum(Dist(:,c) == 1) > 0  % lognormal
%                 b_mtx(Dist(:,c) == 1,:,c) = exp(b_mtx(Dist(:,c) == 1,:,c)); 
%             elseif sum(Dist(:,c) == 2) > 0  % Spike
%                 b_mtx(Dist(:,c) == 2,:,c) = max(b_mtx(Dist(:,c) == 1,:,c),0); 
%             end 
% %             if WTP_space > 0; 
% %                 b_mtx(1:end-WTP_space,:,c) = b_mtx(1:end-WTP_space,:,c).*b_mtx(WTP_matrix,:,c); 
% %             end; 
%         end
%     end
% end


b_mtx_sliced = reshape(b_mtx,[NVarA,NRep,NP,NClass]);
% p_tmp = reshape(p_tmp, [1, NP, NRep, NClass]); ... 
% p_tmp = permute(p_tmp, [1 3 2 4]); ... % 1 x NRep x NP x NClass
% Betas_i = squeeze(mean(b_mtx_sliced.*p_tmp(ones(NVarA,1),:,:,:),2)); ...% NVarA x NP x NClass
p_tmp = reshape(permute(p_tmp,[2 1 3]),[1,NRep,NP,NClass]); ... % 1 x NRep x NP x NClass 
% Betas_i = squeeze(mean(bsxfun(@times,p_tmp,b_mtx_sliced),2)); ...% NVarA x NP x NClass
% p = reshape(p, [1, NP, NClass]); ...
% Betas_i = sum(Betas_i./p(ones(NVarA,1),:,:).* PClass_i(ones(NVarA,1,1),:,:),3)'; ...

% using posterior probabilities for weights:
% Betas_i = sum(squeeze(sum(bsxfun(@times,bsxfun(@rdivide,p_tmp,sum(p_tmp,2)),b_mtx_sliced),2)).*PClass_i(ones(NVarA,1,1),:,:),3)'; ...

PClass2 = PClass_i;
ClassBetas_i = reshape(sum((p_tmp./max(sum(p_tmp,2),realmin)).*b_mtx_sliced,2),[NVarA,NP,NClass]);
Betas_i = sum(ClassBetas_i.*PClass2,3)';

PClass_i = squeeze(PClass_i);

