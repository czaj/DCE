function p = hybrid_parameters(d,b)
o = d.opt; L = o.NLatent; D = o.NVarStr; J = size(d.E,1);
b = b(:); pos = 0; K = 0;
if ~strcmp(d.model,'MIMIC')
    K = o.NVarA;
    if strcmp(d.model,'HLC')
        p.beta = reshape(b(1:K*o.NClass),K,o.NClass); pos = K*o.NClass;
        nc = (o.NVarC+L)*(o.NClass-1);
        p.class = [reshape(b(pos+(1:nc)),o.NVarC+L,[]),zeros(o.NVarC+L,1)];
        pos = pos+nc;
    else
        p.mu = b(1:K); pos = K;
        p.V = zeros(K,J); p.covIndex = [];
        if strcmp(d.model,'HMXL')
            if o.FullCov == 0
                p.covIndex = (1:K+1:K*K)'; q = K;
                p.V(:,1:K) = diag(b(pos+(1:K)));
            elseif o.FullCov == 1
                p.covIndex = find(tril(ones(K))); q = numel(p.covIndex);
                v = zeros(K); v(p.covIndex) = b(pos+(1:q)); p.V(:,1:K) = v;
            else
                mask = tril(true(K+L));
                mask(sub2ind([K+L,K+L],K+(1:L),K+(1:L))) = false;
                p.covIndex = find(mask); q = numel(p.covIndex);
                p.rawV = zeros(K+L); p.rawV(p.covIndex) = b(pos+(1:q));
                p.rawV(sub2ind([K+L,K+L],K+(1:L),K+(1:L))) = 1;
                p.normV = sqrt(sum(p.rawV(K+1:end,:).^2,2));
                p.V = p.rawV(1:K,:);
            end
            pos = pos+q;
        end
        p.bm = reshape(b(pos+(1:K*o.NVarM)),K,o.NVarM); pos = pos+K*o.NVarM;
        p.bl = reshape(b(pos+(1:K*L)),K,L); pos = pos+K*L;
        ns = o.NVarS-o.ScaleLV*L;
        p.bs = b(pos+(1:ns)); pos = pos+ns;
        p.bsLV = zeros(L,1);
        if o.ScaleLV, p.bsLV = b(pos+(1:L)); pos = pos+L; end
    end
end
p.bstr = reshape(b(pos+(1:D*L)),D,L); pos = pos+D*L;
p.bmea = b(pos+1:end);
if strcmp(d.model,'HMXL') && o.FullCov == 2
    p.latentV = p.rawV(K+1:end,:)./p.normV;
else
    p.latentV = [zeros(L,J-L),eye(L)];
end
p.A = [p.bstr',p.latentV];
p.mean = p.A*d.moments.mean;
p.s = sqrt(sum((p.A*d.moments.covariance).*p.A,2));
assert(all(p.s > 0 & isfinite(p.s)),'DCE:Hybrid:Normalization','Latent variance must be positive.');
p.correction = (p.A*d.moments.covariance)./(p.s.^2);
p.npar = numel(b);
end
