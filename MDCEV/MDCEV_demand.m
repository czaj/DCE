function [D,Err] = MDCEV_demand(Results,num)
% Predict demand from an estimated MDCEV model.
% num is empty or the index of an alternative that must be consumed.

EstimOpt = Results.EstimOpt;
if ~isfield(EstimOpt,'NSim')
    EstimOpt.NSim = 500;
end
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
epsDraw = net(sequence,N*EstimOpt.NSim);
epsDraw = permute(reshape(-log(-log(epsDraw')),[NAlt,EstimOpt.NSim,N]),[1 3 2]);
logMU = betasZ+scale*epsDraw-log(priceMat);
D = mean(solveDemand(logMU,priceMat,alphas,gammas,income,num),3);
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

function D = solveDemand(logMU,prices,alphas,gammas,income,num)
[nAlt,n,nSim] = size(logMU);
D = zeros(nAlt,n,nSim);
for i=1:n
    col = min(i,size(alphas,2));
    inva = 1./(1-alphas(:,col));
    gamma = gammas(:,min(i,size(gammas,2)));
    for j=1:nSim
        mu = logMU(:,i,j);
        high = max(mu)+1;
        if ~isempty(num)
            high = min(high,mu(num));
        end
        spend = @(lambda) sum(prices(:,i).*gamma.*max(exp((mu-lambda).*inva)-1,0));
        if spend(high) > income(i)
            error('MDCEV_demand:Numeraire','The requested numeraire cannot be positive for task %d, draw %d.',i,j)
        end
        low = high-10;
        while spend(low) < income(i)
            low = low-10;
        end
        for k=1:80
            mid = (low+high)/2;
            if spend(mid) > income(i)
                low = mid;
            else
                high = mid;
            end
        end
        D(:,i,j) = gamma.*max(exp((mu-(low+high)/2).*inva)-1,0);
    end
end
end
