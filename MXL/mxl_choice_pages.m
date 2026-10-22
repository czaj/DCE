function [panel,score,Sm,Ss,Su] = mxl_choice_pages(X,E,M,Xs,chosen,available,opt,mu,bm,VC,bs,meanDraw,scaleDraw,chosenX)
% Draw-level choice probabilities and scores; shared CPU/GPU page kernel.
K = opt.NVarA; R = opt.NRep; T = opt.NAlt*opt.NCT;
W = opt.WTP_space; mapping = opt.WTP_matrix;
MCount = size(M,1); SCount = size(Xs,2); count = size(X,3);
lognormal = find(opt.Dist == 1);
if nargin < 12, meanDraw = 0; end
if nargin < 13, scaleDraw = 1; end
if SCount > 0, X = X.*exp(pagemtimes(Xs,bs)); end
base = mu + meanDraw + pagemtimes(VC,E);
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
V = V.*scaleDraw;
unavailable = ~repmat(available,[1 R 1]);
V(unavailable) = 0;
U = V;
U(unavailable) = -Inf;
U = reshape(U,[opt.NAlt,opt.NCT,R,count]);
missing = repmat(~any(reshape(available,opt.NAlt,opt.NCT,1,count),1),[1,1,R,1]);
shift = max(U,[],1);
shift(missing) = 0;
U = exp(U-shift);
denominator = sum(U,1);
denominator(missing) = 1;
P4 = U./denominator;
selected = reshape(P4,[T,R,count]);
selected(~repmat(chosen,[1,R,1])) = 1;
panel = prod(selected,1);
if nargout == 1, return; end
P = reshape(P4,[T,R,count]);
Su = sum((chosen-P).*V,1);
Ss = pagemtimes(Xs,'transpose',(chosen-P).*V,'none');
Sm = zeros(K*MCount,R,count,'like',X);
if ~varying
    if nargin < 14 || SCount > 0, chosenX = permute(sum(X.*chosen,1),[2 1 3]); end
    D = (chosenX - pagemtimes(X,'transpose',P,'none')).*scaleDraw;
    score = D;
    if W > 0
        score(1:K-W,:,:) = D(1:K-W,:,:).*z(mapping,:,:);
        for k = K-W+1:K
            linked = mapping == k;
            score(k,:,:) = D(k,:,:) + sum(D(linked,:,:).*z(linked,:,:),1);
        end
    end
    score(lognormal,:,:) = score(lognormal,:,:).*z(lognormal,:,:);
    for m = 1:MCount
        Sm((m-1)*K+(1:K),:,:) = score.*M(m,1,:);
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
        J = J.*scaleDraw;
        if opt.Dist(k) == 1, J = J.*zk; end
        A = residual.*J;
        A(unavailable) = 0;
        score(k,:,:) = sum(A,1);
        meanScore = pagemtimes(M,A);
        Sm(k:K:end,:,:) = meanScore;
    end
end
end
