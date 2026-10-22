function [f,g] = hybrid_block(d,p,ids,need)
o = d.opt; R = o.NRep; N = numel(ids); L = o.NLatent; D = o.NVarStr;
e = d.E(:,:,ids);
x = reshape(d.Xstr(ids,:)',D,1,N);
LV = (pagemtimes(p.bstr',x)+pagemtimes(p.latentV,e)-p.mean)./p.s;
if need
    [measurement,dl,gm] = hybrid_measurement(d,ids,LV,p.bmea);
else
    measurement = hybrid_measurement(d,ids,LV,p.bmea);
end
if strcmp(d.model,'MIMIC')
    panel = ones(1,R,N,'like',LV);
elseif strcmp(d.model,'HLC')
    beta = p.beta;
    if o.WTP_space, beta(1:end-o.WTP_space,:) = beta(1:end-o.WTP_space,:).*beta(o.WTP_matrix,:); end
    u = pagemtimes(d.X(:,:,ids),beta);
    u(repmat(~d.A(:,:,ids),1,o.NClass,1)) = -Inf;
    u = reshape(u,o.NAlt,o.NCT,o.NClass,N);
    shift = max(u,[],1); shift(~isfinite(shift)) = 0;
    u = exp(u-shift); den = sum(u,1); den(den == 0) = 1;
    pr = reshape(u./den,o.NAlt*o.NCT,o.NClass,N);
    pr(repmat(~d.Y(:,:,ids),1,o.NClass,1)) = 1;
    choice = prod(pr,1);
    xc = [repmat(reshape(d.Xc(ids,:)',o.NVarC,1,N),1,R,1);LV];
    membership = pagemtimes(p.class',xc);
    membership = exp(membership-max(membership,[],1)); membership = membership./sum(membership,1);
    panel = sum(membership.*permute(choice,[2 1 3]),1);
else
    scale = exp(pagemtimes(p.bsLV',LV));
    if need
        [panel,S,Sm,Ss,Su] = mxl_choice_pages(d.X(:,:,ids),e,d.M(:,:,ids),d.Xs(:,:,ids),...
            d.Y(:,:,ids),d.A(:,:,ids),o,p.mu,p.bm,p.V,p.bs,pagemtimes(p.bl,LV),scale);
        dl = dl+pagemtimes(p.bl',S)+p.bsLV.*Su;
    else
        panel = mxl_choice_pages(d.X(:,:,ids),e,d.M(:,:,ids),d.Xs(:,:,ids),...
            d.Y(:,:,ids),d.A(:,:,ids),o,p.mu,p.bm,p.V,p.bs,pagemtimes(p.bl,LV),scale);
    end
end
joint = panel.*measurement;
prob = mean(joint,2);
if o.RealMin == 1, prob = max(prob,realmin); end
if o.RealMin >= 2
    prob(~isfinite(prob)) = 0; prob = max(prob,o.ProbMin);
end
f = reshape(-log(prob),N,1); g = zeros(N,0,'like',LV);
if ~need, return; end
w = joint./(R*prob);
if o.RealMin >= 2, w(repmat(prob <= o.ProbMin,1,R,1)) = 0; end
dw = dl.*w;
z = sum(dw.*LV,2);
xc = x-reshape(d.moments.mean(1:D),D,1,1);
ec = e-reshape(d.moments.mean(D+1:end),[],1,1);
gs = pagemtimes(sum(dw,2)./p.s,'none',xc,'transpose')-z.*reshape(p.correction(:,1:D),L,D,1);
gs = reshape(permute(gs,[2 1 3]),D*L,1,N);
gmea = sum(gm.*w,2);
if strcmp(d.model,'MIMIC')
    scores = [gs;gmea];
else
    K = o.NVarA;
    gv = zeros(0,1,N,'like',LV);
    if strcmp(d.model,'HMXL')
        v = pagemtimes(S.*w,'none',e,'transpose');
        if o.FullCov == 0
            v = reshape(v(:,1:K,:),K*K,N); gv = reshape(v(p.covIndex,:),K,1,N);
        elseif o.FullCov == 1
            v = reshape(v(:,1:K,:),K*K,N); gv = reshape(v(p.covIndex,:),[],1,N);
        else
            vl = pagemtimes(dw,'none',ec,'transpose')./p.s-...
                z.*reshape(p.correction(:,D+1:end),L,[],1);
            vl = (vl-sum(vl.*p.latentV,2).*p.latentV)./p.normV;
            v = reshape([v;vl],(K+L)^2,N); gv = reshape(v(p.covIndex,:),[],1,N);
        end
    end
    gl = pagemtimes(S.*w,'none',LV,'transpose');
    gl = reshape(gl,K*L,1,N);
    scaleScore = sum(Ss.*w,2);
    if o.ScaleLV, scaleScore = [scaleScore;sum(Su.*LV.*w,2)]; end
    scores = [sum(S.*w,2);gv;sum(Sm.*w,2);gl;scaleScore;gs;gmea];
end
g = -reshape(scores,p.npar,N)';
end
