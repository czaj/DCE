function [f,g] = vce_block(d,b,ids,need)
% Contract draw scores inside each block, without an alternatives-by-parameters tensor.
o = d.opt; A = o.NAlt; K = o.NVarA; T = d.T; R = d.R; B = numel(ids);
tasks = reshape((ids-1)*T+(1:T)',1,[]); N = T*B;
X = d.X(:,:,ids); E = d.E(:,:,ids); Xm = d.M(:,:,ids);
y = d.Y(:,tasks); price = d.price(:,tasks); U = d.U(tasks,:);
pos = K; beta = reshape(b(1:K),K,1);
covIndex = [];
if d.mixed
    if o.FullCov == 0
        covIndex = (1:K+1:K*K)';
    else
        covIndex = find(tril(true(K)));
    end
    V = zeros(K,'like',b); V(covIndex) = b(pos+(1:numel(covIndex))); pos = pos+numel(covIndex);
    beta = beta+pagemtimes(V,E);
end
if o.NVarM > 0
    bm = reshape(b(pos+(1:K*o.NVarM)),K,o.NVarM); pos = pos+K*o.NVarM;
    beta = beta+pagemtimes(bm,Xm);
end
if size(beta,3) == 1, beta = repmat(beta,1,1,B); end
if d.mixed, beta(o.Dist == 1,:,:) = exp(beta(o.Dist == 1,:,:)); end
spec = o.SpecProfile; nAlpha = max(spec(1,:)); nGamma = max(spec(2,:));
assert(nAlpha+nGamma == o.NVarP && numel(b) == pos+o.NVarP*(1+o.NVarU)+1,'DCE:VCE:Parameters','Incorrect parameter count.');
coef = reshape(b(pos+(1:o.NVarP*(1+o.NVarU))),1+o.NVarU,o.NVarP);
fit = U*coef;
fa = zeros(A,N,'like',b); gamma = ones(A,N,'like',b);
mask = spec(1,:) > 0; fa(mask,:) = fit(:,spec(1,mask))';
mask = spec(2,:) > 0; gamma(mask,:) = exp(fit(:,nAlpha+spec(2,mask)))';
oneMinusAlpha = exp(-fa); scale = exp(b(end));
logY = log(y+gamma); logGamma = log(gamma); logPrice = log(price);
logAlpha = -fa; alphaActive = ones(A,N,'like',b);
if o.RealMin == 1
    logY = log(max(y+gamma,realmin)); logGamma = log(max(gamma,realmin));
    logAlpha = max(logAlpha,log(realmin)); alphaActive = oneMinusAlpha > realmin;
end
L = logY-logGamma;
logc = logAlpha-logY-logPrice;
chosen = y > 0; M = sum(chosen,1);
% Normalize the Jacobian sum in log space as well as the choice probabilities.
logRatio = -logc; logRatio(~chosen) = -Inf;
shift = max(logRatio,[],1); ratio = exp(logRatio-shift);
q = ratio./sum(ratio,1);
logJacobian = sum(chosen.*logc,1)+shift+log(sum(ratio,1));
utility = pagemtimes(X,beta);
utility = reshape(permute(reshape(utility,A,T,R,B),[1 2 4 3]),A,N,R);
utility = (utility-oneMinusAlpha.*L-logPrice)/scale;
utility = utility-max(utility,[],1);
ev = exp(utility); P = ev./sum(ev,1);
logp = (1-M)*b(end)+logJacobian+sum(chosen.*utility,1)-M.*log(sum(ev,1))+gammaln(M);
panel = reshape(sum(reshape(logp,T,B,R),1),B,R);
shift = max(panel,[],2); w = exp(panel-shift); total = sum(w,2);
f = -shift-log(total/R); g = zeros(B,0,'like',b);
if ~need, return; end
w = w./total; wp = reshape(w',1,R,B);
residual = chosen-M.*P;
resPages = reshape(permute(reshape(residual,A,T,B,R),[1 2 4 3]),A*T,R,B);
score = pagemtimes(X,'transpose',resPages,'none')/scale;
if d.mixed, score(o.Dist == 1,:,:) = score(o.Dist == 1,:,:).*beta(o.Dist == 1,:,:); end
weighted = score.*wp; meanScore = sum(weighted,2);
scores = meanScore;
if d.mixed
    covariance = reshape(pagemtimes(weighted,'none',E,'transpose'),K*K,B);
    scores = [scores;reshape(covariance(covIndex,:),[],1,B)];
end
for m = 1:o.NVarM, scores = [scores;meanScore.*Xm(m,1,:)]; end %#ok<AGROW>
taskWeights = reshape(permute(repmat(wp,T,1,1),[1 3 2]),1,N,R);
for kind = 1:2
    if kind == 1
        profileScore = (-chosen+q).*alphaActive+residual.*oneMinusAlpha.*L/scale;
    else
        profileScore = (-chosen+q).*gamma./(y+gamma)+residual.*oneMinusAlpha.*y./(y+gamma)/scale;
    end
    for group = 1:max(spec(kind,:))
        s = sum(sum(profileScore(spec(kind,:) == group,:,:),1).*taskWeights,3)';
        s = reshape(sum(reshape(U.*s,T,B,1+o.NVarU),1),B,1+o.NVarU);
        scores = [scores;reshape(s',1+o.NVarU,1,B)]; %#ok<AGROW>
    end
end
s = (1-M)-sum(residual.*utility,1);
s = reshape(sum(reshape(s,T,B,R),1),B,R);
scores = [scores;reshape(sum(s.*w,2),1,1,B)];
g = -reshape(scores,numel(b),B)';
end
