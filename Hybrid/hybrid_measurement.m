function [likelihood,dLV,g] = hybrid_measurement(d,ids,LV,b)
% Measurement likelihood and log scores, limited to the current respondent block.
o = d.opt; R = o.NRep; N = numel(ids); L = o.NLatent;
likelihood = ones(1,R,N,'like',LV);
need = nargout > 1;
dLV = zeros(L,R,N*need,'like',LV);
g = zeros(numel(b),R,N*need,'like',LV);
pos = 0;
for i = 1:size(d.Ymea,2)
    type = o.MeaSpecMatrix(i);
    latent = find(o.MeaMatrix(:,i) == 1);
    X = LV(latent,:,:);
    if o.MeaExpMatrix(i)
        ex = reshape(d.Xexp(ids,:)',o.NVarMeaExp,1,N);
        X = [X;repmat(ex,1,R,1)];
    end
    if type ~= 2, X = [ones(1,R,N,'like',LV);X]; end
    P = size(X,1);
    y = reshape(d.Ymea(ids,i),1,1,N);
    observed = reshape(d.observed(ids,i),1,1,N);
    if type == 1
        categories = d.categories{i}; C = numel(categories);
        nb = P*(C-1); beta = [zeros(P,1,'like',LV),reshape(b(pos+(1:nb)),P,C-1)];
        u = pagemtimes(beta',X); u = exp(u-max(u,[],1)); pr = u./sum(u,1);
        selected = reshape(y == reshape(categories,[],1),C,1,N);
        prob = sum(pr.*selected,1);
        if need
            residual = selected-pr;
            score = zeros(nb,R,N,'like',LV);
            for c = 2:C, score((c-2)*P+(1:P),:,:) = X.*residual(c,:,:); end
            dl = pagemtimes(beta(2:1+numel(latent),:),residual);
        end
    elseif type == 2
        categories = d.categories{i}; C = numel(categories); k = C-1;
        nb = P+k; beta = b(pos+(1:P));
        threshold = reshape(b(pos+P+(1:k)),k,1);
        tail = reshape(threshold(2:end),k-1,1);
        increments = [ones(1,1,'like',LV);exp(min(max(tail,-o.MaxExp),o.MaxExp))];
        alpha = cumsum([threshold(1);reshape(increments(2:end),k-1,1)]);
        fit = sum(X.*beta,1);
        hi = inf(1,1,N,'like',LV); lo = -hi;
        category = zeros(1,1,N);
        for c = 1:C
            ix = d.Ymea(ids,i) == categories(c);
            category(1,1,ix) = c;
            if c <= k, hi(1,1,ix) = alpha(c); end
            if c > 1, lo(1,1,ix) = alpha(c-1); end
        end
        hi = hi-fit; lo = lo-fit;
        prob = normcdf(hi)-normcdf(lo);
        % Use the upper tail when subtracting two CDFs close to one.
        upper = lo > 0;
        prob(upper) = normcdf(-lo(upper))-normcdf(-hi(upper));
        if need
            ph = normpdf(hi); pl = normpdf(lo);
            denom = max(prob,realmin);
            df = (pl-ph)./denom;
            score = zeros(nb,R,N,'like',LV); score(1:P,:,:) = X.*df;
            for j = 1:k
                score(P+j,:,:) = (ph.*(category >= j & category <= k)-...
                    pl.*(category > j))./denom*increments(j);
                if j > 1 && abs(threshold(j)) >= o.MaxExp, score(P+j,:,:) = 0; end
            end
            dl = beta(1:numel(latent)).*df;
        end
    else
        inflated = ismember(type,[5 6]); negativeBinomial = ismember(type,[4 6]);
        nb = P*(1+inflated)+double(type == 0 || negativeBinomial);
        beta = b(pos+inflated*P+(1:P)); fit = sum(X.*beta,1);
        if type == 0
            logSigma = b(pos+nb);
            sigma = max(exp(min(max(logSigma,-o.MaxExp),o.MaxExp)),o.SigmaMin);
            residual = y-fit;
            prob = normpdf(residual,0,sigma);
            if need
                df = residual/(sigma^2);
                extra = -1+(residual/sigma).^2;
                if abs(logSigma) >= o.MaxExp || exp(logSigma) <= o.SigmaMin, extra(:) = 0; end
            end
        else
            lambda = exp(fit);
            if negativeBinomial
                theta = exp(b(pos+nb)); u = theta./(theta+lambda);
                prob = exp(gammaln(theta+y)-gammaln(theta)-gammaln(y+1)).*u.^theta.*(1-u).^y;
                if need
                    df = y.*u-theta*(1-u);
                    extra = theta*(psi(theta+y)-psi(theta)+log(u)+1-u)-y.*u;
                end
            else
                prob = exp(y.*fit-lambda-gammaln(y+1));
                if need, df = y-lambda; end
            end
        end
        if inflated
            zi = b(pos+(1:P)); z = sum(X.*zi,1);
            pzero = 1./(1+exp(-z)); iszero = y == 0;
            mix = (1-pzero).*prob+pzero.*iszero;
            if need
                responsibility = (1-pzero).*prob./max(mix,realmin);
                dz = pzero.*(iszero./max(mix,realmin)-1);
                df = df.*responsibility;
                if negativeBinomial, extra = extra.*responsibility; end
            end
            prob = mix;
        end
        if need
            score = zeros(nb,R,N,'like',LV);
            score(inflated*P+(1:P),:,:) = X.*df;
            dl = beta(2:1+numel(latent)).*df;
            if inflated
                score(1:P,:,:) = X.*dz;
                dl = dl+zi(2:1+numel(latent)).*dz;
            end
            if type == 0 || negativeBinomial, score(end,:,:) = extra; end
        end
    end
    if o.RealMin >= 2 && ismember(type,[0 2])
        clamped = prob <= o.ProbMin | ~isfinite(prob);
        prob(clamped) = o.ProbMin;
        if need
            score(repmat(clamped,size(score,1),1,1)) = 0;
            dl(repmat(clamped,numel(latent),1,1)) = 0;
        end
    end
    prob(repmat(~observed,1,R,1)) = 1;
    likelihood = likelihood.*prob;
    if need
        score(repmat(~observed,size(score,1),R,1)) = 0;
        dl(repmat(~observed,numel(latent),R,1)) = 0;
        g(pos+(1:nb),:,:) = score;
        dLV(latent,:,:) = dLV(latent,:,:)+dl;
    end
    pos = pos+nb;
end
assert(pos == numel(b),'DCE:Hybrid:Parameters','Unexpected measurement parameter count.');
end
