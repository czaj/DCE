function [f,g] = LL_mxl_gpu(data,B,blockSize)
% Blocked double-precision likelihood for the supported MXL GPU path.
opt = data.opt;
K = opt.NVarA;
R = opt.NRep;
NP = opt.NP;
T = opt.NAlt*opt.NCT;
MCount = opt.NVarM;
SCount = opt.NVarS;
QCount = K + opt.FullCov*K*(K-1)/2;
assert(isscalar(blockSize) && blockSize >= 1 && blockSize == fix(blockSize));
assert(numel(B) == K+QCount+K*MCount+SCount,'Unexpected parameter layout.');
B = gpuArray(double(B(:)));
mu = B(1:K);
covIndex = find(tril(ones(K)));
VC = zeros(K,K,'like',data.X);
if opt.FullCov
    VC(covIndex) = B(K+(1:QCount));
else
    VC(1:K+1:end) = B(K+(1:K));
end
bm = reshape(B(K+QCount+(1:K*MCount)),[K,MCount]);
bs = B(K+QCount+K*MCount+(1:SCount));
f = zeros(NP,1,'like',data.X);
g = zeros(NP,numel(B)*(nargout > 1),'like',data.X);
for first = 1:blockSize:NP
    ids = first:min(first+blockSize-1,NP);
    count = numel(ids);
    X = data.X(:,:,ids);
    E = data.E(:,:,ids);
    if MCount > 0, M = data.M(:,:,ids); else, M = zeros(0,1,count,'like',X); end
    if SCount > 0, Xs = data.Xs(:,:,ids); else, Xs = zeros(T,0,count,'like',X); end
    chosen = data.Y(:,:,ids);
    available = data.A(:,:,ids);
    if nargout == 1
        panel = mxl_choice_pages(X,E,M,Xs,chosen,available,opt,mu,bm,VC,bs);
    else
        [panel,score,Sm,Ss] = mxl_choice_pages(X,E,M,Xs,chosen,available,opt,mu,bm,VC,bs,0,1,data.chosenX(:,:,ids));
    end
    p = mean(panel,2);
    if opt.RealMin == 1, p = max(p,realmin); end
    f(ids) = reshape(-log(p),[count,1]);
    if nargout == 1, continue; end
    gm = -mean(score.*panel,2)./p;
    gmm = -mean(Sm.*panel,2)./p;
    gs = -mean(Ss.*panel,2)./p;
    if opt.FullCov
        weighted = score.*panel;
        covarianceScore = -pagemtimes(weighted,'none',E,'transpose')./(R*p);
        covarianceScore = reshape(covarianceScore,[K*K,count]);
        gv = reshape(covarianceScore(covIndex,:),[QCount,1,count]);
        if opt.RealMin ~= 1 && any(gather(p > 0 & p < realmin),'all')
            [ci,cj] = ind2sub([K,K],covIndex);
            gv = -mean((score(ci,:,:).*E(cj,:,:)).*panel,2)./p;
        end
    else
        gv = -mean((score.*E).*panel,2)./p;
    end
    g(ids,:) = reshape([gm;gv;gmm;gs],[numel(B),count])';
end
end
