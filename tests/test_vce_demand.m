function test_vce_demand(baselineDir,outDir)
addpath(baselineDir); if ~isfolder(outDir), mkdir(outDir); end
rows = struct('Case',{},'Error',{},'BudgetError',{},'OldSeconds',{},'NewSeconds',{});
for model = {'MDCEV','MMDCEV'}
    for profile = 1:3
        for u = 0:1
            c = vce_fixture(model{1},profile,1,u,1,5,32,2,3);
            c.opt.NSim = 16; c.data.I = 100*ones(c.opt.NP*c.opt.NCT,1);
            r = struct('INPUT',c.data,'EstimOpt',c.opt,'bhat',c.b);
            fun = [c.model,'_demand'];
            rng(123); t = tic; expected = feval(['baseline_',fun],r,[]); oldSeconds = toc(t);
            rng(123); t = tic; actual = feval(fun,r,[]); newSeconds = toc(t);
            r.EstimOpt.VCEBlockSize = 1; rng(123); single = feval(fun,r,[]);
            err = max(abs(actual-expected),[],'all');
            budget = max(abs(sum(actual.*c.data.priceMat,1)-c.data.I'));
            assert(err < 1e-9 && budget < 1e-9 && max(abs(actual-single),[],'all') < 1e-9,c.name);
            rows(end+1) = struct('Case',c.name,'Error',err,'BudgetError',budget,'OldSeconds',oldSeconds,'NewSeconds',newSeconds); %#ok<AGROW>
        end
    end
end
% Feasible/infeasible numeraire requests and zero income keep their semantics.
mu = zeros(3,2,4); prices = ones(3,2); a = zeros(3,2); gamma = ones(3,2);
d = mdcevDemand(mu,prices,a,gamma,[100 0],[]);
assert(max(abs(sum(d,1)-reshape([100 0],1,2,1)),[],'all') < 1e-9);
d = mdcevDemand(mu,prices,a,gamma,[100 100],1); assert(all(d(1,:,:) > 0,'all'));
mu(1,:,:) = -100; failed = false;
try, mdcevDemand(mu,prices,a,gamma,[1 1],1); catch ex, failed = strcmp(ex.identifier,'Tools:MDCEV:Numeraire'); end
assert(failed);
writetable(struct2table(rows),fullfile(outDir,'demand.csv'));
fprintf('VCE_DEMAND_PASSED %d cases and budget/numeraire edge cases\n',numel(rows));
end
