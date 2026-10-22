function [f,g] = hybrid_evaluate(d,b,blockSize,need,nWorkers,workerData)
p = hybrid_parameters(d,b);
onGPU = isa(d.E,'gpuArray');
if onGPU
    for name = {'mu','V','bm','bl','bs','bsLV','bstr','bmea','latentV','mean','s','correction','normV','beta','class'}
        if isfield(p,name{1}), p.(name{1}) = gpuArray(p.(name{1})); end
    end
end
N = d.opt.NP; starts = 1:blockSize:N;
ff = cell(numel(starts),1); gg = ff;
if nWorkers > 0 && ~onGPU
    parfor (j = 1:numel(starts),nWorkers)
        ids = starts(j):min(starts(j)+blockSize-1,N);
        [ff{j},gg{j}] = hybrid_block(workerData.Value,p,ids,need);
    end
else
    for j = 1:numel(starts)
        ids = starts(j):min(starts(j)+blockSize-1,N);
        [ff{j},gg{j}] = hybrid_block(d,p,ids,need);
    end
end
f = vertcat(ff{:}); g = vertcat(gg{:});
if onGPU, f = gather(f); g = gather(g); end
end
