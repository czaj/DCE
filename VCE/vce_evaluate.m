function [f,g] = vce_evaluate(d,b,block,need,nWorkers,resident)
N = d.N; starts = 1:block:N; ff = cell(numel(starts),1); gg = ff;
onGPU = isa(d.E,'gpuArray');
if onGPU, b = gpuArray(b); end
if nWorkers > 0 && ~onGPU
    parfor (j = 1:numel(starts),nWorkers)
        ids = starts(j):min(starts(j)+block-1,N);
        [ff{j},gg{j}] = vce_block(resident.Value,b,ids,need);
    end
else
    for j = 1:numel(starts)
        ids = starts(j):min(starts(j)+block-1,N);
        [ff{j},gg{j}] = vce_block(d,b,ids,need);
    end
end
f = vertcat(ff{:}); g = vertcat(gg{:});
if onGPU, f = gather(f); g = gather(g); end
end
