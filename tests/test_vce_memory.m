function report = test_vce_memory(outDir,baselineDir,workers)
if nargin < 3, workers = 0; end
if ~isfolder(outDir), mkdir(outDir); end
addpath(baselineDir);
settings = parallel.Settings; oldAuto = settings.Pool.AutoCreate;
restore = onCleanup(@() setAuto(settings,oldAuto)); %#ok<NASGU>
settings.Pool.AutoCreate = false; assert(isempty(gcp('nocreate')));
cases = {};
for model = {'MDCEV','MMDCEV'}
    for p = 1:3
        for cov = 0:strcmp(model{1},'MMDCEV')
            for u = 0:1
                for m = 0:1
                    cases{end+1} = vce_fixture(model{1},p,cov,u,m); %#ok<AGROW>
                    cases{end}.opt.RealMin = m;
                end
            end
        end
    end
end
for counts = {[1 7 2 3],[4 1 3 3],[1 1 1 1],[3 7 2 1],[1 7 200 3]}
    v = counts{1}; cases{end+1} = vce_fixture('MMDCEV',3,1,1,1,v(1),v(2),v(3),v(4)); %#ok<AGROW>
end
for cov = 0:1
    c = vce_fixture('MMDCEV',3,cov,1,1);
    c.name = [c.name,'_fixed']; c.opt.Dist(1) = -1; c.data.err(1,:) = 0;
    if cov, c.b(4:6) = 0; else, c.b(4) = 0; end
    cases{end+1} = c; %#ok<AGROW>
end
rows = struct('Case',{},'FD',{},'BaselineValue',{});
for model = {'MDCEV','MMDCEV'}
    c = vce_fixture(model{1},3,1,1,1); c.name = [c.name,'_logical_covariates'];
    c.data.Xm = c.data.Xm > 0; c.data.Xu = c.data.Xu > 0;
    cases{end+1} = c; %#ok<AGROW>
end
for model = {'MDCEV','MMDCEV'}
    c = vce_fixture(model{1},3,1,0,0); c.name = [c.name,'_extreme_alpha'];
    c.b(end-2) = 50; c.compareBaseline = false;
    cases{end+1} = c; %#ok<AGROW>
end
serial = cell(size(cases));
for k = 1:numel(cases)
    c = cases{k}; entry = ['probs_',lower(c.model)];
    [f,g] = feval(entry,c.data,c.opt,c.b); fv = feval(entry,c.data,c.opt,c.b);
    assert(all(isfinite([f(:);g(:)])) && relative(f,fv) < 1e-12,c.name);
    gn = zeros(size(g));
    for j = 1:numel(c.b)
        h = 1e-5*max(1,abs(c.b(j))); bp = c.b; bm = bp; bp(j) = bp(j)+h; bm(j) = bm(j)-h;
        gn(:,j) = (feval(entry,c.data,c.opt,bp)-feval(entry,c.data,c.opt,bm))/(2*h);
    end
    fd = relative(g,gn); assert(fd < 5e-6,'Gradient mismatch %s: %.3g',c.name,fd);
    fo = feval(['baseline_',entry],c.data,c.opt,c.b);
    delta = relative(f,fo);
    if ~isfield(c,'compareBaseline') || c.compareBaseline
        assert(delta < 1e-9,'Value mismatch %s: %.3g',c.name,delta);
    end
    o = c.opt; o.VCEBlockSize = 1; [f1,g1] = feval(entry,c.data,o,c.b);
    o.VCEBlockSize = 100; [f2,g2] = feval(entry,c.data,o,c.b);
    assert(relative(f1,f2) < 1e-10 && relative(g1,g2) < 1e-9,c.name);
    oo = struct('GradObj','on','Hessian','user-supplied','FinDiffType','central');
    o.BActive(end) = 0;
    [fw,gw,H] = feval(['LL_',lower(c.model),'_MATlike'],c.data,o,oo,c.b);
    W = c.data.W; if strcmp(c.model,'MMDCEV'), W = W(1:o.NCT:end); end
    j = g.*W; j(:,end) = 0;
    assert(relative(fw,sum(f.*W)) < 1e-10 && relative(gw(:),sum(j,1)') < 1e-5 && relative(H,j'*j) < 1e-5,c.name);
    rows(k) = struct('Case',c.name,'FD',fd,'BaselineValue',delta); %#ok<AGROW>
    serial{k} = struct('f',f,'g',g);
    fprintf('%d/%d %s FD %.3g baseline %.3g\n',k,numel(cases),c.name,fd,delta);
end
if workers > 0
    pool = parpool('Processes',workers); cleanup = onCleanup(@() delete(pool)); %#ok<NASGU>
    spmd, maxNumCompThreads(2); end
    for k = 1:numel(cases)
        c = cases{k}; [f,g] = feval(['probs_',lower(c.model)],c.data,c.opt,c.b);
        assert(relative(f,serial{k}.f) < 1e-10 && relative(g,serial{k}.g) < 1e-9,c.name);
    end
end
report = struct2table(rows); writetable(report,fullfile(outDir,'vce.csv'));
save(fullfile(outDir,'vce.mat'),'cases','serial','report');
fprintf('VCE_TESTS_PASSED %d cases, workers=%d\n',numel(cases),workers);
end

function r = relative(a,b)
r = max(abs(a-b),[],'all')/max(1,max(abs(b),[],'all'));
end
function setAuto(settings,value)
settings.Pool.AutoCreate = value;
end
