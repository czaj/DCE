function report = test_hybrid_gpu(outDir,fixtureFile)
if ~exist(outDir,'dir'), mkdir(outDir); end
s = load(fixtureFile,'cases');
device = gpuDevice; fprintf('GPU: %s, capability %s\n',device.Name,device.ComputeCapability);
rows = struct([]);
for k = 1:numel(s.cases)
    c = s.cases{k}; a = c.args;
    a{end-1}.GPU = 'cpu'; [fc,gc] = hybrid_likelihood(c.model,a,true);
    for block = [1 3]
        a{end-1}.GPU = 'gpu'; a{end-1}.HybridBlockSize = block;
        [f,g,backend] = hybrid_likelihood(c.model,a,true);
        [fv,~,valueBackend] = hybrid_likelihood(c.model,a,false);
        assert(strcmp(backend,'gpu') && strcmp(valueBackend,'gpu'),'GPU not exercised: %s',c.name);
        ef = relative(f,fc); eg = relative(g,gc);
        assert(ef < 1e-8 && eg < 1e-6 && relative(f,fv) < 1e-8,c.name);
        j = numel(rows)+1; rows(j).Case = c.name; rows(j).Block = block;
        rows(j).ValueError = ef; rows(j).GradientError = eg;
        fprintf('%s block %d GPU value %.3g gradient %.3g\n',c.name,block,ef,eg);
    end
end
report = struct2table(rows); writetable(report,fullfile(outDir,'gpu.csv'));
save(fullfile(outDir,'gpu.mat'),'report');
fprintf('HYBRID_GPU_PASSED %d comparisons\n',height(report));
end

function r = relative(a,b)
r = max(abs(a-b),[],'all')/max(1,max(abs(b),[],'all'));
end
