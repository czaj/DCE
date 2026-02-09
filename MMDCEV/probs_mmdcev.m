function [f, grad] = probs_mmdcev(data, EstimOpt, b0)
% Function calculates probabilities of the MDCEV model.
% It returns log(probabilities)
% See:
% Bhat, Chandra R. 2008. 
% ‘The Multiple Discrete-Continuous Extreme Value (MDCEV) Model: Role of 
% Utility Function Parameters, Identification Considerations, and Model 
% Extensions’. Transportation Research Part B: Methodological 42 (3): 
% 274–303. https://doi.org/10.1016/j.trb.2007.06.002.
%
% Estimating probabilities as defined in eq. (19) or (20) in the article.
%
%
% EstimOpt.Profile -- sets the utility profile
% 1 -- "alpha" profile
% 2 -- "gamma" profile
%
% variables -- includes variables to be estimated:
%   1. betas    -- alternative attributes parameters (NVarA)
%   2. alphas   -- satiation parameters (NAlt)
%      or
%   2. gammas   -- translation parameters (enable corner solutions) (NAlt)
%   3. sigma    -- scale parameter (1)
%
% data contains:
%   Y           -- dependent variables (quantities demanded); (NAltxN)
%   Xa          -- covariates (alternatives attributes) (NxNVar)
%   priceMat    -- matrix of prices of each alternative (NAltxN)


% save LL_mnl
% return

y = data.Y; % dependent variables (quantities demands); size: NAltxN
Xu = data.Xu; % covariates (socio-demographic) (NxNVarU)
Xa = data.Xa; % covariates (alternatives attributes) (NxNVar)
Xm = data.Xm; % covariates (socio-demographic) (NxNVarM)

priceMat = data.priceMat; % matrix of prices of each alternative (NAltxN)
err = data.err;

FullCov = EstimOpt.FullCov;
NVarA = EstimOpt.NVarA; % Number of attributes
NVarP = EstimOpt.NVarP; % Number of parameters for alpha/gamma profile
NVarU = EstimOpt.NVarU;
NVarM = EstimOpt.NVarM;

NRep = EstimOpt.NRep;
Profile = EstimOpt.Profile; % Utility function version
SpecProfile = EstimOpt.SpecProfile;
Dist = EstimOpt.Dist;
NAlt = EstimOpt.NAlt; % Number of alternatives
NP = EstimOpt.NP;
NCT = EstimOpt.NCT;
RealMin = EstimOpt.RealMin;
indx1 = EstimOpt.indx1;
indx2 = EstimOpt.indx2;
N = size(y,2); % Number of decisions

% define variables to be optimized
betas = b0(1:NVarA); % betas
if FullCov == 0
    b0v = b0(NVarA+1:2*NVarA);
    VC = diag(b0v);
    l = 2*NVarA;
else
    b0v = b0(NVarA+1:NVarA+sum(1:NVarA));
    VC = tril(ones(NVarA));
    VC(VC == 1) = b0v;
    l = NVarA+sum(1:NVarA);
end
b_mtx = betas + VC*err;
if NVarM > 0
    b0m = reshape(b0(l+1:l+NVarA*NVarM), [NVarA,NVarM]);
    XmFit = b0m*Xm; % NVarA x NP
    XmFit = reshape(permute(XmFit(:,:, ones(1, NRep)), [1 3 2]), [NVarA, NRep*NP]);
    b_mtx = b_mtx + XmFit;
    l = l+NVarA*NVarM;
end
if sum(Dist == 1) > 0 % Log-normal
    b_mtx(Dist == 1,:) = exp(b_mtx(Dist == 1,:));
end
b_mtx = reshape(b_mtx,[NVarA,NRep,NP]);

b_profile = b0(l+1:l+NVarP*(1+NVarU));
scale = exp(b0(l+NVarP*(1+NVarU)+1)); % scale parameter (sigma); one for all dataset
% scale = exp(variables(end)); % scale parameter (sigma); one for all dataset
% exp() for ensuring > 0

% alphas and gammas
if Profile == 1 % alpha profile
    if NVarU == 0
        alphas = b_profile;
        alphas = 1 - exp(-alphas); % between 0 and 1
    else
        alpha = reshape(b_profile, [NVarU+1, NAlt]);
        alphas = 1 - exp(-Xu*alpha)'; % NAlt x N
    end
    gammas = ones(size(alphas));
elseif Profile == 2 % gamma profile
    if NVarU == 0
        gammas = b_profile;
        gammas = exp(gammas); % greater than 0
    else
        gamma = reshape(b_profile, [NVarU+1, NAlt]);
        gammas= exp(Xu*gamma)'; % NAlt x N
    end
    alphas = zeros(size(gammas));
elseif Profile == 3 % SpecProfile
    indx = length(unique(SpecProfile(1,SpecProfile(1,:) ~= 0))); % unique alphas
    if NVarU == 0
        a = b_profile(1:indx);
        a = 1 - exp(-a);
        a = a(SpecProfile(1,SpecProfile(1,:) ~= 0)); 
        alphas = zeros(NAlt,1);
        alphas(SpecProfile(1,:) ~= 0) = a;
    
        g = b_profile(indx+1:end);
        g = exp(g); % greater than 0
        g = g(SpecProfile(2,SpecProfile(2,:) ~= 0)); 
        gammas = ones(NAlt,1);
        gammas(SpecProfile(2,:) ~= 0) = g;
    else

        a = reshape(b_profile(1:indx*(NVarU+1)), [NVarU+1, indx]);
        alpha = 1 - exp(-Xu*a)';
        alpha = alpha(SpecProfile(1,SpecProfile(1,:) ~= 0),:); 
        alphas = zeros(NAlt,N);
        alphas(SpecProfile(1,:) ~= 0,:) = alpha;
    
        g = reshape(b_profile(indx*(NVarU+1)+1:end), [NVarU+1, NVarP - indx]);
        gamma = exp(Xu*g)'; % greater than 0
        gamma = gamma(SpecProfile(2,SpecProfile(2,:) ~= 0),:); 
        gammas = ones(NAlt,N);
        gammas(SpecProfile(2,:) ~= 0,:) = gamma;
    end
end
 
isChosen = (y ~= 0); % if the alternative was chosen or not
M = sum(isChosen, 1); % number of consumed goods in each decision

% computing baseline utility levels
betasZ = zeros(NAlt*NCT, NRep, NP);
for i = 1:NP
    betasZ(:,:,i) = Xa(:,:,i)*b_mtx(:,:,i);
end
betasZ = reshape(permute(betasZ, [1, 3, 2]), [NAlt, N, NRep]);

%% logarithms
% logf_i = log(1 - alphas) - log(y + gammas); % (NAltxN)
logc_i = log(1 - alphas) - log(y + gammas) - log(priceMat); % (NAltxN)
AlphaDerV = log(y ./ gammas + 1);
V = betasZ + (alphas - 1) .* AlphaDerV - log(priceMat);

logV = V / scale; % NAlt x N x NRep

priceRatio = exp(-logc_i);
logV = logV - max(logV,[], 1);
sumV = sum(exp(logV), 1);

logprobs = (1 - M) .* log(scale) + ...
    sum(isChosen .* logc_i, 1) + ...
    log(sum(isChosen .* priceRatio, 1)) + ...
    (sum(isChosen .* logV, 1) - M .* log(sumV)) + ...
    gammaln(M); % log(factorial(M-1)) = gammaln(M)

probs = reshape(exp(logprobs), [NCT, NP, NRep]);
ProbsProd = prod(probs,1); % 1 x NP x NRep
f = mean(ProbsProd,3)'; % NP x 1

if nargout == 2 % gradient
    if FullCov == 0
        grad = zeros(NP,(2+NVarM)*NVarA+NVarP*(1+NVarU)+1);
        VC2 = permute(reshape(err,[NVarA,NRep,NP]),[3, 2, 1]);
    else        
        grad = zeros([NP,(1+NVarM)*NVarA+sum(1:NVarA)+NVarP*(1+NVarU)+1]);
        VC2 = permute(reshape(err,[NVarA,NRep,NP]),[3, 2, 1]);
    end
    if NVarU > 0
        Xu = permute(Xu, [3, 1, 4, 2]); % 1 x N x 1 x 1+NVarU
    end
    % Gradient for random parameters
    Xa = reshape(permute(reshape(Xa, [NAlt, NCT, NVarA, NP]), [1 2 4 3]), [NAlt, NCT*NP, 1, NVarA])/scale;
    % New version
    GradV = Xa - sum((exp(logV)./sumV).*Xa,1);
    GradV = sum(isChosen .* GradV,1);
    GradV = squeeze(sum(reshape(GradV, [NCT, NP, NRep, NVarA]),1).*ProbsProd); % NP x NRep x NVarA
    
    GradVX = GradV;
    if sum(Dist == 1) > 0 % Log-normal
        GradV(:,:,Dist == 1) = GradV(:,:,Dist == 1).*permute(b_mtx(Dist == 1,:,:), [3, 2, 1]);
        GradVX(:,:,Dist == 1) = GradVX(:,:,Dist == 1).*permute(b_mtx(Dist == 1,:,:), [3, 2, 1]);
    end

    if FullCov == 0
        GradV_cov = GradVX.*VC2;
        grad(:, (NVarA + 1):2*NVarA) = -squeeze(mean(GradV_cov,2))./f;
        l = 2*NVarA;
    else % FullCov = 1
        GradV_cov = GradVX(:,:, indx1).*VC2(:,:,indx2);
        grad(:, (NVarA + 1):(2*NVarA + NVarA*(NVarA-1)/2)) = -squeeze(mean(GradV_cov,2))./f;
        l = (2*NVarA + NVarA*(NVarA-1)/2);
    end
    
    if NVarM > 0
        gradA = -squeeze(mean(GradV,2))./f;
        grad(:, 1:NVarA) = gradA;
%         grad(:, l+1:l+NVarA*NVarM) = repmat(gradA,[1,NVarM]).*Xm(reshape(repmat((1:NVarM), [NVarA, 1]), [NVarA*NVarM,1]),:)';
        grad(:, l+1:l+NVarA*NVarM) = reshape(permute(Xm, [2 3 1]).*reshape(gradA, [NP, NVarA, 1]), [NP, NVarA*NVarM]);
        l = l+NVarA*NVarM;
    else
        grad(:, 1:NVarA) = -squeeze(mean(GradV,2))./f;
    end
    % Gradient for profile parameters
      % Alphas
      indx = length(unique(SpecProfile(1,SpecProfile(1,:) ~= 0))); % unique alphas
      if indx > 0
          Tmp = (y + gammas).*priceMat;
          GradAlpha = -isChosen./(1-alphas) + isChosen .* (exp(-2*logc_i)./Tmp)./sum(isChosen .* priceRatio,1); % NAlt x N
    %       GradAlpha2 = (AlphaDerV/scale - sum((exp(logV)./sumV).*AlphaDerV/scale,1)).*isChosen;
          GradAlpha2 = (AlphaDerV.*isChosen/scale - M.*(exp(logV)./sumV).*AlphaDerV/scale);
          GradAlpha = -(GradAlpha + GradAlpha2).*(alphas-1);
          GradAlphaSum = zeros(indx, NCT*NP, NRep);
          for i = 1:indx
               GradAlphaSum(i,:,:) = sum(GradAlpha(SpecProfile(1,:) == i,:,:),1);
          end
          if NVarU > 0
            GradAlphaSum = reshape(permute(GradAlphaSum.*Xu, [4 1 2 3]), [indx*(1+NVarU), NCT*NP, NRep]); 
            indx = indx*(1+NVarU);
          end
          % 1 x NP x NRep
          GradAlphaSum = sum(permute(reshape(GradAlphaSum, [indx, NCT, NP, NRep]), [2, 3, 4, 1]),1); % 1 x NP x NRep x indx
          GradAlphaSum = reshape(GradAlphaSum.*ProbsProd, [NP, NRep, indx]);
          grad(:, (l+1):(l+indx)) = -reshape(mean(GradAlphaSum,2), [NP, indx])./f;
          l = (l+indx);
      end

      % Gammas
      indx = length(unique(SpecProfile(2,SpecProfile(2,:) ~= 0))); % unique gammas
      if indx > 0
          Tmp = 1./(y + gammas);
          Tmp2 = (Tmp.^2).*(1-alphas)./priceMat;
          GradGamma = -isChosen.*Tmp + isChosen .* (exp(-2*logc_i).*Tmp2)./sum(isChosen .* priceRatio,1); % NAlt x N
          GammaDerV = -Tmp.*y.*(alphas-1)./gammas;
          GradGamma2 = (GammaDerV.*isChosen/scale - M.*(exp(logV)./sumV).*GammaDerV/scale);
          GradGamma = (GradGamma + GradGamma2).*gammas;
          GradGammaSum = zeros(indx, NCT*NP, NRep);
          for i = 1:indx
               GradGammaSum(i,:,:) = sum(GradGamma(SpecProfile(2,:) == i,:,:),1);
          end
          if NVarU > 0
            GradGammaSum = reshape(permute(GradGammaSum.*Xu, [4 1 2 3]), [indx*(1+NVarU), NCT*NP, NRep]); 
            indx = indx*(1+NVarU);
          end

          % 1 x NP x NRep
          GradGammaSum = sum(permute(reshape(GradGammaSum, [indx, NCT, NP, NRep]), [2, 3, 4, 1]),1); % 1 x NP x NRep x indx
          GradGammaSum = reshape(GradGammaSum.*ProbsProd, [NP, NRep, indx]);
          grad(:, (l+1):(l+indx)) = -reshape(mean(GradGammaSum,2), [NP, indx])./f;
      end

    % Gradient for scale
    ScaleGrad = (1-M) - sum(isChosen .* logV, 1) + M.*sum(exp(logV).*logV,1)./sumV;
    ScaleGrad = reshape(sum(reshape(ScaleGrad, [NCT, NP, NRep]),1).*ProbsProd, [NP, NRep]);
    grad(:, end) = -mean(ScaleGrad,2)./f;

end

if RealMin == 1
    f = -log(max(f,realmin));
else
    f = -log(f);
end
end
