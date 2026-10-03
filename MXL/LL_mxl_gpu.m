function [f,g] = LL_mxl_gpu(data,B,blockSize)
% Blocked double-precision likelihood for the supported MXL GPU path.
opt = data.opt;
K = opt.NVarA;
R = opt.NRep;
NP = opt.NP;
T = opt.NAlt*opt.NCT;
W = opt.WTP_space;
mapping = opt.WTP_matrix;
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
lognormal = find(opt.Dist == 1);
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
    if SCount > 0, X = X.*exp(pagemtimes(Xs,bs)); end
    base = mu + pagemtimes(VC,E);
    varying = opt.mCT ~= 0 && MCount > 0;
    if ~varying
        z = base;
        if MCount > 0, z = z + pagemtimes(bm,M); end
        z(lognormal,:,:) = exp(z(lognormal,:,:));
        beta = z;
        if W > 0, beta(1:K-W,:,:) = z(1:K-W,:,:).*z(mapping,:,:); end
        V = pagemtimes(X,beta);
    else
        fit = pagemtimes(bm,M);
        V = zeros(T,R,count,'like',X);
        costZ = cell(W,1);
        costX = cell(W,1);
        for c = 1:W
            k = K-W+c;
            zc = base(k,:,:) + permute(fit(k,:,:),[2 1 3]);
            if opt.Dist(k) == 1, zc = exp(zc); end
            costZ{c} = zc;
            costX{c} = X(:,k,:);
        end
        for k = 1:K-W
            zk = base(k,:,:) + permute(fit(k,:,:),[2 1 3]);
            if opt.Dist(k) == 1, zk = exp(zk); end
            term = X(:,k,:).*zk;
            if W > 0
                c = mapping(k)-(K-W);
                costX{c} = costX{c} + term;
            else
                V = V + term;
            end
        end
        for c = 1:W, V = V + costX{c}.*costZ{c}; end
    end
    unavailable = ~repmat(available,[1 R 1]);
    V(unavailable) = 0;
    U = V;
    U(unavailable) = -Inf;
    U = reshape(U,[opt.NAlt,opt.NCT,R,count]);
    missing = repmat(data.missing(:,:,:,ids),[1,1,R,1]);
    shift = max(U,[],1);
    shift(missing) = 0;
    U = exp(U-shift);
    denominator = sum(U,1);
    denominator(missing) = 1;
    P4 = U./denominator;
    selected = reshape(P4,[T,R,count]);
    selected(~repmat(chosen,[1,R,1])) = 1;
    panel = prod(selected,1);
    p = mean(panel,2);
    if opt.RealMin == 1, p = max(p,realmin); end
    f(ids) = reshape(-log(p),[count,1]);
    if nargout == 1, continue; end
    P = reshape(P4,[T,R,count]);
    if SCount > 0
        scaleScore = pagemtimes(Xs,'transpose',(chosen-P).*V,'none');
        gs = -mean(scaleScore.*panel,2)./p;
    else
        gs = zeros(0,1,count,'like',X);
    end
    gmm = zeros(K*MCount,1,count,'like',X);
    if ~varying
        if SCount > 0
            chosenX = permute(sum(X.*chosen,1),[2 1 3]);
        else
            chosenX = data.chosenX(:,:,ids);
        end
        D = chosenX - pagemtimes(X,'transpose',P,'none');
        score = D;
        if W > 0
            score(1:K-W,:,:) = D(1:K-W,:,:).*z(mapping,:,:);
            for k = K-W+1:K
                linked = mapping == k;
                score(k,:,:) = D(k,:,:) + sum(D(linked,:,:).*z(linked,:,:),1);
            end
        end
        score(lognormal,:,:) = score(lognormal,:,:).*z(lognormal,:,:);
        gm = -mean(score.*panel,2)./p;
        for m = 1:MCount
            gmm((m-1)*K+(1:K),:,:) = -mean((score.*M(m,1,:)).*panel,2)./p;
        end
    else
        score = zeros(K,R,count,'like',X);
        residual = chosen-P;
        for k = 1:K
            if W > 0 && k > K-W
                c = k-(K-W);
                J = costX{c};
                zk = costZ{c};
            else
                J = X(:,k,:);
                if W > 0, J = J.*costZ{mapping(k)-(K-W)}; end
                zk = base(k,:,:) + permute(fit(k,:,:),[2 1 3]);
                if opt.Dist(k) == 1, zk = exp(zk); end
            end
            if opt.Dist(k) == 1, J = J.*zk; end
            A = residual.*J;
            A(unavailable) = 0;
            score(k,:,:) = sum(A,1);
            meanScore = pagemtimes(M,A);
            gmm(k:K:end,:,:) = -mean(meanScore.*panel,2)./p;
        end
        gm = -mean(score.*panel,2)./p;
    end
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
