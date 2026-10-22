function result = bench_vce_memory(outDir,baselineDir,model,mode,NP,R)
% Run each version in a fresh process on the same idle node.
assert(~isfolder(outDir),'Use a new output directory.'); mkdir(outDir);
outDir = char(java.io.File(outDir).getCanonicalPath()); addpath(baselineDir);
c = vce_fixture(model,3,1,1,1,NP,R,6,3); c.opt = rmfield(c.opt,'VCEBlockSize'); c.opt.GPU = mode;
settings = parallel.Settings; previous = settings.Pool.AutoCreate;
restore = onCleanup(@() setAuto(settings,previous)); %#ok<NASGU>
settings.Pool.AutoCreate = false; assert(isempty(gcp('nocreate')));
for k = 1:2, [f,g,backend] = evaluate(c,mode); end
assert(all(isfinite([f(:);g(:)])));
if strcmp(mode,'gpu'), assert(strcmp(backend,'gpu'),'GPU benchmark fell back to CPU.'); end
pid = feature('getpid'); write_json(fullfile(outDir,'ready.json'),struct('pids',pid,'NP',NP));
sampler = System.Diagnostics.Process;
sampler.StartInfo.FileName = 'powershell.exe';
sampler.StartInfo.Arguments = sprintf('-NoProfile -ExecutionPolicy Bypass -File "%s" -OutDir "%s" -AttachPid %d -SampleMilliseconds 200',...
    fullfile(fileparts(mfilename('fullpath')),'measure_mxl_memory.ps1'),outDir,pid);
sampler.StartInfo.UseShellExecute = false; sampler.StartInfo.CreateNoWindow = true;
assert(sampler.Start()); wait_for(fullfile(outDir,'sample.start'));
seconds = []; timer = tic;
while numel(seconds) < 5 || toc(timer) < 20
    one = tic; [f,g] = evaluate(c,mode); seconds(end+1) = toc(one); %#ok<AGROW>
end
result = struct('model',model,'mode',mode,'backend',backend,'workers',0,'NP',NP,'NRep',R,...
    'NCT',c.opt.NCT,'parameters',numel(c.b),'evalSeconds',seconds,'medianEvalSeconds',median(seconds),...
    'evaluations',numel(seconds),'negativeLL',sum(f),'pids',pid,'MATLAB',version);
write_json(fullfile(outDir,'evaluation.json'),result); wait_for(fullfile(outDir,'sample.done'));
save(fullfile(outDir,'evaluation.mat'),'result','f','g');
fprintf('VCE_BENCHMARK_COMPLETE %s %s NP=%d R=%d median=%.6fs\n',model,mode,NP,R,result.medianEvalSeconds);
end

function [f,g,backend] = evaluate(c,mode)
if strcmp(mode,'legacy')
    fun = @(b) feval(['baseline_probs_',lower(c.model)],c.data,c.opt,b);
    if strcmp(c.model,'MDCEV')
        f = fun(c.b); g = numdiff(fun,f,c.b,true,c.opt.BActive);
    else
        [f,g] = fun(c.b);
    end
    backend = 'legacy';
else
    [f,g,backend] = vce_likelihood(c.model,c.data,c.opt,c.b,true);
end
end
function write_json(path,value)
fid = fopen(path,'w'); assert(fid >= 0); close = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(value));
end
function wait_for(path)
t = tic;
while ~isfile(path), assert(toc(t) < 120,'Memory sampler timeout.'); pause(.1); end
end
function setAuto(settings,value)
settings.Pool.AutoCreate = value;
end
