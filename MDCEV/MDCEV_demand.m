function [D,Err] = MDCEV_demand(Results,num)
% Predict demand from an estimated MDCEV model.
% num is empty or the index of an alternative that must be consumed.

EstimOpt = Results.EstimOpt;
if ~isfield(EstimOpt,'NSim')
    EstimOpt.NSim = 500;
end
validateattributes(EstimOpt.NSim,{'numeric'},{'scalar','integer','positive','finite'});
NAlt = EstimOpt.NAlt;
N = EstimOpt.NCT*EstimOpt.NP;
if ~isempty(num)
    validateattributes(num,{'numeric'},{'scalar','integer','>=',1,'<=',NAlt})
end

Xa = Results.INPUT.Xa;
Y = Results.INPUT.Y;
priceMat = Results.INPUT.priceMat;
income = Results.INPUT.I(:)';
Xm = Results.INPUT.Xm;
Xu = Results.INPUT.Xu;
b = Results.bhat;

NVarA = EstimOpt.NVarA;
NVarM = EstimOpt.NVarM;
NVarU = EstimOpt.NVarU;
NVarP = EstimOpt.NVarP;
betas = b(1:NVarA);
l = NVarA;
if NVarM > 0
    meanShifters = reshape(b(l+1:l+NVarA*NVarM),NVarA,NVarM);
    l = l+NVarA*NVarM;
end
bProfile = b(l+1:l+NVarP*(1+NVarU));
scale = exp(b(l+NVarP*(1+NVarU)+1));
[alphas,gammas] = profileValues(bProfile,Xu,EstimOpt,N);

betasZ = reshape(Xa*betas,[NAlt,N]);
if NVarM > 0
    shift = Xm*meanShifters';
    betasZ = betasZ+sum(reshape(Xa,[NAlt,N,NVarA]).*reshape(shift,[1,N,NVarA]),3);
end

sequence = scramble(sobolset(NAlt,'Skip',1),'MatousekAffineOwen');
block = min(N,max(1,floor(64*1024^2/(8*8*NAlt*EstimOpt.NSim))));
if isfield(EstimOpt,'VCEBlockSize') && ~isempty(EstimOpt.VCEBlockSize), block = EstimOpt.VCEBlockSize; end
validateattributes(block,{'numeric'},{'scalar','integer','positive','finite'});
D = zeros(NAlt,N);
for first = 1:block:N
    ids = first:min(first+block-1,N);
    epsDraw = sequence((first-1)*EstimOpt.NSim+1:ids(end)*EstimOpt.NSim,:);
    epsDraw = permute(reshape(-log(-log(epsDraw')),[NAlt,EstimOpt.NSim,numel(ids)]),[1 3 2]);
    logMU = betasZ(:,ids)+scale*epsDraw-log(priceMat(:,ids));
    aa = alphas(:,min(ids,size(alphas,2))); gg = gammas(:,min(ids,size(gammas,2)));
    D(:,ids) = mean(mdcevDemand(logMU,priceMat(:,ids),aa,gg,income(ids),num),3);
end
Err = Y-D;
end

function [alphas,gammas] = profileValues(bProfile,Xu,opt,N)
spec = opt.SpecProfile;
nRows = 1+(opt.NVarU > 0)*(N-1);
alphas = zeros(opt.NAlt,nRows);
gammas = ones(opt.NAlt,nRows);
nAlpha = numel(unique(spec(1,spec(1,:) ~= 0)));
if nAlpha > 0
    coef = reshape(bProfile(1:nAlpha*(1+opt.NVarU)),1+opt.NVarU,nAlpha);
    fit = profileFit(Xu,coef,opt.NVarU);
    mask = spec(1,:) ~= 0;
    alphas(mask,:) = (1-exp(-fit(:,spec(1,mask))))';
end
nGamma = opt.NVarP-nAlpha;
if nGamma > 0
    coef = reshape(bProfile(nAlpha*(1+opt.NVarU)+1:end),1+opt.NVarU,nGamma);
    fit = profileFit(Xu,coef,opt.NVarU);
    mask = spec(2,:) ~= 0;
    gammas(mask,:) = exp(fit(:,spec(2,mask)))';
end
end

function fit = profileFit(Xu,coef,NVarU)
if NVarU == 0
    fit = coef;
else
    fit = Xu*coef;
end
end
