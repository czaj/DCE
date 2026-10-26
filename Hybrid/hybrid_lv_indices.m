function indices = hybrid_lv_indices(model,o)
% Choice means, covariance parameters and observed mean interactions precede LV links.
K = o.NVarA; offset = K*(1+o.NVarM);
if strcmp(model,'HMXL')
    switch o.FullCov
        case 0, offset = offset+K;
        case 1, offset = offset+K*(K+1)/2;
        case 2, offset = offset+(K+o.NLatent)*(K+o.NLatent+1)/2-o.NLatent;
        otherwise, error('DCE:Hybrid:FullCov','FullCov must be 0, 1 or 2.');
    end
else
    assert(strcmp(model,'HMNL'),'DCE:Hybrid:Model','Expected HMNL or HMXL.');
end
indices = offset+(1:K*o.NLatent);
end
