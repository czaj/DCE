function test_vce_saved(file,model,baselineDir,outDir)
% Validate copies of saved results without rerunning or changing their estimation.
addpath(baselineDir); if ~isfolder(outDir), mkdir(outDir); end
s = load(file,'Results'); r = s.Results.(model); d = r.INPUT; o = r.EstimOpt;
o.GPU = 'cpu'; o.VCEBlockSize = 16;
if strcmp(model,'MMDCEV')
    o.NRep = 64; d.err = generateRandomDraws(o)'; d.err(o.Dist == -1,:) = 0;
end
fun = @(b) feval(['probs_',lower(model)],d,o,b);
[f,g] = fun(r.bhat);
old = @(b) feval(['baseline_probs_',lower(model)],d,o,b);
fo = old(r.bhat); eg = NaN;
if strcmp(model,'MMDCEV'), [~,go] = old(r.bhat); eg = relative(g,go); assert(eg < 1e-6); end
ev = relative(f,fo); assert(ev < 1e-8 && all(isfinite([f(:);g(:)])));
indices = unique(round(linspace(1,numel(r.bhat),min(16,numel(r.bhat)))));
gn = zeros(size(g,1),numel(indices));
for j = 1:numel(indices)
    k = indices(j); h = 1e-5*max(1,abs(r.bhat(k))); bp = r.bhat; bm = bp; bp(k) = bp(k)+h; bm(k) = bm(k)-h;
    gn(:,j) = (fun(bp)-fun(bm))/(2*h);
end
fd = relative(g(:,indices),gn); assert(fd < 1e-5);
report = struct('model',model,'NP',o.NP,'parameters',numel(r.bhat),'ValueError',ev,'GradientError',eg,'FD',fd);
save(fullfile(outDir,'saved.mat'),'report','f','g','d','o');
disp(report); fprintf('VCE_SAVED_PASSED %s\n',model);
end
function r = relative(a,b)
r = max(abs(a-b),[],'all')/max(1,max(abs(b),[],'all'));
end
