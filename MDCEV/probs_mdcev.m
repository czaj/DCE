function f = probs_mdcev(data, EstimOpt, variables)
% Per-task negative log-likelihood contributions for the MDCEV model.
%
% INPUT.Xm shifts attribute coefficients, consistently with MXL/MMDCEV.
% INPUT.Xu shifts alpha/gamma utility-profile parameters.

y = data.Y;
Xa = data.Xa;
Xm = data.Xm;
Xu = data.Xu;
priceMat = data.priceMat;

NVarA = EstimOpt.NVarA;
NVarP = EstimOpt.NVarP;
NVarM = EstimOpt.NVarM;
NVarU = EstimOpt.NVarU;
NAlt = EstimOpt.NAlt;
N = size(y,2);

betas = variables(1:NVarA);
l = NVarA;
if NVarM > 0
    meanShifters = reshape(variables(l+1:l+NVarA*NVarM),NVarA,NVarM);
    l = l+NVarA*NVarM;
end
bProfile = variables(l+1:l+NVarP*(1+NVarU));
scale = exp(variables(l+NVarP*(1+NVarU)+1));

SpecProfile = EstimOpt.SpecProfile;
if EstimOpt.Profile == 1
    coef = reshape(bProfile,1+NVarU,NVarP);
    fit = profileFit(Xu,coef,NVarU);
    alphas = (1-exp(-fit(:,SpecProfile(1,:))))';
    gammas = ones(size(alphas));
elseif EstimOpt.Profile == 2
    coef = reshape(bProfile,1+NVarU,NVarP);
    fit = profileFit(Xu,coef,NVarU);
    gammas = exp(fit(:,SpecProfile(2,:)))';
    alphas = zeros(size(gammas));
else
    nAlpha = numel(unique(SpecProfile(1,SpecProfile(1,:) ~= 0)));
    nRows = size(profileFit(Xu,zeros(1+NVarU,1),NVarU),1);
    alphas = zeros(NAlt,nRows);
    gammas = ones(NAlt,nRows);
    if nAlpha > 0
        coef = reshape(bProfile(1:nAlpha*(1+NVarU)),1+NVarU,nAlpha);
        fit = profileFit(Xu,coef,NVarU);
        mask = SpecProfile(1,:) ~= 0;
        alphas(mask,:) = (1-exp(-fit(:,SpecProfile(1,mask))))';
    end
    nGamma = NVarP-nAlpha;
    if nGamma > 0
        coef = reshape(bProfile(nAlpha*(1+NVarU)+1:end),1+NVarU,nGamma);
        fit = profileFit(Xu,coef,NVarU);
        mask = SpecProfile(2,:) ~= 0;
        gammas(mask,:) = exp(fit(:,SpecProfile(2,mask)))';
    end
end

betasZ = reshape(Xa*betas,[NAlt,N]);
if NVarM > 0
    shift = Xm*meanShifters';
    betasZ = betasZ + sum(reshape(Xa,[NAlt,N,NVarA]).*reshape(shift,[1,N,NVarA]),3);
end

isChosen = y ~= 0;
M = sum(isChosen,1);
if EstimOpt.RealMin == 1
    logc = log(max(1-alphas,realmin))-log(max(y+gammas,realmin))-log(max(priceMat,realmin));
    V = betasZ+(alphas-1).*log(max(y./gammas+1,realmin))-log(max(priceMat,realmin));
else
    logc = log(1-alphas)-log(y+gammas)-log(priceMat);
    V = betasZ+(alphas-1).*log(y./gammas+1)-log(priceMat);
end

logV = V/scale;
logV = logV-max(logV,[],1);
sumV = sum(exp(logV),1);
priceRatio = exp(-logc);
if EstimOpt.RealMin == 1
    logprobs = (1-M).*log(max(scale,realmin)) + sum(isChosen.*logc,1) + ...
        log(max(sum(isChosen.*priceRatio,1),realmin)) + ...
        sum(isChosen.*logV,1)-M.*log(max(sumV,realmin)) + gammaln(M);
else
    logprobs = (1-M).*log(scale) + sum(isChosen.*logc,1) + ...
        log(sum(isChosen.*priceRatio,1)) + ...
        sum(isChosen.*logV,1)-M.*log(sumV) + gammaln(M);
end
f = -logprobs';
end

function fit = profileFit(Xu,coef,NVarU)
if NVarU == 0
    fit = coef;
else
    fit = Xu*coef;
end
end
