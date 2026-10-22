function [D,Err] = MMDCEV_demand(Results,num)
% Predict demand from an estimated diagonal or full-covariance MMDCEV model.
% num is empty or the index of an alternative that must be consumed.

EstimOpt = Results.EstimOpt;
if ~isfield(EstimOpt,'NSim')
    EstimOpt.NSim = min(500,EstimOpt.NRep);
end
validateattributes(EstimOpt.NSim,{'numeric'},{'scalar','integer','positive','finite'});
assert(EstimOpt.NSim <= EstimOpt.NRep,'NSim must not exceed NRep.')
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
if isfield(Results.INPUT,'err')
    draws = Results.INPUT.err;
else
    draws = generateRandomDraws(EstimOpt)';
end
b = Results.bhat;

NVarA = EstimOpt.NVarA;
NVarM = EstimOpt.NVarM;
NVarU = EstimOpt.NVarU;
NVarP = EstimOpt.NVarP;
betas = b(1:NVarA);
if EstimOpt.FullCov == 0
    VC = diag(b(NVarA+1:2*NVarA));
    l = 2*NVarA;
else
    VC = tril(ones(NVarA));
    VC(VC == 1) = b(NVarA+1:NVarA+sum(1:NVarA));
    l = NVarA+sum(1:NVarA);
end
draws = reshape(draws,NVarA,EstimOpt.NRep,EstimOpt.NP);
draws = draws(:,1:EstimOpt.NSim,:);
bMatrix = betas+pagemtimes(VC,draws);
if NVarM > 0
    meanShifters = reshape(b(l+1:l+NVarA*NVarM),NVarA,NVarM);
    meanFit = meanShifters*Xm;
    meanFit = reshape(meanFit,NVarA,1,EstimOpt.NP);
    bMatrix = bMatrix+meanFit;
    l = l+NVarA*NVarM;
end
if any(EstimOpt.Dist == 1)
    bMatrix(EstimOpt.Dist == 1,:,:) = exp(bMatrix(EstimOpt.Dist == 1,:,:));
end

bProfile = b(l+1:l+NVarP*(1+NVarU));
scale = exp(b(l+NVarP*(1+NVarU)+1));
[alphas,gammas] = profileValues(bProfile,Xu,EstimOpt,N);

sequence = scramble(sobolset(NAlt,'Skip',1),'MatousekAffineOwen');
block = min(EstimOpt.NP,max(1,floor(64*1024^2/(8*8*NAlt*EstimOpt.NSim*EstimOpt.NCT))));
if isfield(EstimOpt,'VCEBlockSize') && ~isempty(EstimOpt.VCEBlockSize), block = EstimOpt.VCEBlockSize; end
validateattributes(block,{'numeric'},{'scalar','integer','positive','finite'});
D = zeros(NAlt,N);
for first = 1:block:EstimOpt.NP
    ids = first:min(first+block-1,EstimOpt.NP);
    tasks = (first-1)*EstimOpt.NCT+1:ids(end)*EstimOpt.NCT;
    betasZ = pagemtimes(Xa(:,:,ids),bMatrix(:,:,ids));
    betasZ = reshape(permute(reshape(betasZ,[NAlt,EstimOpt.NCT,EstimOpt.NSim,numel(ids)]),[1 2 4 3]),NAlt,numel(tasks),EstimOpt.NSim);
    epsDraw = sequence((tasks(1)-1)*EstimOpt.NSim+1:tasks(end)*EstimOpt.NSim,:);
    epsDraw = permute(reshape(-log(-log(epsDraw')),[NAlt,EstimOpt.NSim,numel(tasks)]),[1 3 2]);
    logMU = betasZ+scale*epsDraw-log(priceMat(:,tasks));
    aa = alphas(:,min(tasks,size(alphas,2))); gg = gammas(:,min(tasks,size(gammas,2)));
    D(:,tasks) = mean(mdcevDemand(logMU,priceMat(:,tasks),aa,gg,income(tasks),num),3);
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
