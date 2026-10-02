function [panel,S,Sm,Ss,Su] = mxl_choice(Y,X,M,Xs,E,opt,mu,bm,VC,bs,meanDraw,scaleDraw)
% Conditional panel probabilities and draw-level log-probability derivatives.
K = opt.NVarA;
R = opt.NRep;
W = opt.WTP_space;
mapping = opt.WTP_matrix;
lognormal = opt.Dist(:) == 1;
if nargin < 11 || isempty(meanDraw), meanDraw = 0; end
if nargin < 12 || isempty(scaleDraw), scaleDraw = 1; end
available = ~isnan(Y);
chosen = Y == 1;
X(~available,:) = 0;
if size(M,2) > 1, M(:,~available) = 0; end
if isempty(Xs), Xs = zeros(numel(Y),0); end
Xs(~available,:) = 0;
if isempty(bs), scale = 1; else, scale = exp(Xs*bs); end
X = X.*scale;
base = mu + meanDraw + VC*E;
varying = size(M,2) > 1;
if ~varying
    z = base + bm*M;
    z(lognormal,:) = exp(z(lognormal,:));
    beta = z;
    if W > 0, beta(1:K-W,:) = z(1:K-W,:).*z(mapping,:); end
    V = (X*beta).*scaleDraw;
else
    % One alternative-by-draw buffer per attribute, not rows-by-K-by-R.
    fit = bm*M;
    V = zeros(size(X,1),R);
    costZ = cell(W,1);
    costX = cell(W,1);
    for c = 1:W
        k = K-W+c;
        zc = base(k,:) + fit(k,:)';
        if lognormal(k), zc = exp(zc); end
        costZ{c} = zc;
        costX{c} = X(:,k);
    end
    for k = 1:K-W
        zk = base(k,:) + fit(k,:)';
        if lognormal(k), zk = exp(zk); end
        term = X(:,k).*zk;
        if W > 0
            c = mapping(k)-(K-W);
            costX{c} = costX{c} + term;
        else
            V = V + term;
        end
    end
    for c = 1:W, V = V + costX{c}.*costZ{c}; end
    V = V.*scaleDraw;
end
V(~available,:) = 0;
U = V;
U(~available,:) = -Inf;
U = reshape(U,[opt.NAlt,opt.NCT,R]);
missingTask = ~any(reshape(available,[opt.NAlt,opt.NCT]),1);
shift = max(U,[],1);
shift(:,missingTask,:) = 0;
U = exp(U-shift);
denominator = sum(U,1);
denominator(:,missingTask,:) = 1;
P = reshape(U./denominator,[opt.NAlt*opt.NCT,R]);
panel = prod(P(chosen,:),1);
if nargout == 1, return; end
Ss = zeros(size(Xs,2),R);
Su = [];
if size(Xs,2) > 0 || nargout > 4
    weightedV = (chosen-P).*V;
    Su = sum(weightedV,1);
    Ss = Xs'*weightedV;
end
S = zeros(K,R);
Sm = zeros(K*size(M,1),R);
if ~varying
    D = (sum(X(chosen,:),1)' - X'*P).*scaleDraw;
    S = D;
    if W > 0
        S(1:K-W,:) = D(1:K-W,:).*z(mapping,:);
        for k = K-W+1:K
            linked = mapping == k;
            S(k,:) = D(k,:) + sum(D(linked,:).*z(linked,:),1);
        end
    end
    S(lognormal,:) = S(lognormal,:).*z(lognormal,:);
    for m = 1:size(M,1), Sm((m-1)*K+(1:K),:) = S.*M(m); end
else
    Q = chosen-P;
    for k = 1:K
        if W > 0 && k > K-W
            c = k-(K-W);
            J = costX{c}.*scaleDraw;
            zk = costZ{c};
        else
            J = X(:,k).*scaleDraw;
            if W > 0, J = J.*costZ{mapping(k)-(K-W)}; end
            zk = base(k,:) + fit(k,:)';
            if lognormal(k), zk = exp(zk); end
        end
        if lognormal(k), J = J.*zk; end
        A = Q.*J;
        A(~available,:) = 0;
        S(k,:) = sum(A,1);
        Sm(k:K:end,:) = M*A;
    end
end
end
