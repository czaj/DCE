function test_vce_gpu(fixtureFile,outDir)
if ~isfolder(outDir), mkdir(outDir); end
s = load(fixtureFile,'cases','serial'); gpuDevice;
rows = struct('Case',{},'Block',{},'ValueError',{},'GradientError',{});
for k = 1:numel(s.cases)
    c = s.cases{k};
    for block = [1 3]
        c.opt.GPU = 'gpu'; c.opt.VCEBlockSize = block;
        [f,g,backend] = vce_likelihood(c.model,c.data,c.opt,c.b,true);
        [fv,~,bv] = vce_likelihood(c.model,c.data,c.opt,c.b,false);
        assert(strcmp(backend,'gpu') && strcmp(bv,'gpu'),'GPU not exercised: %s',c.name);
        ef = relative(f,s.serial{k}.f); eg = relative(g,s.serial{k}.g);
        assert(ef < 1e-8 && eg < 1e-6 && relative(f,fv) < 1e-10,c.name);
        rows(end+1) = struct('Case',c.name,'Block',block,'ValueError',ef,'GradientError',eg); %#ok<AGROW>
        fprintf('GPU %s block %d value %.3g gradient %.3g\n',c.name,block,ef,eg);
    end
end
writetable(struct2table(rows),fullfile(outDir,'gpu.csv'));
fprintf('VCE_GPU_PASSED %d comparisons\n',numel(rows));
end
function r = relative(a,b)
r = max(abs(a-b),[],'all')/max(1,max(abs(b),[],'all'));
end
