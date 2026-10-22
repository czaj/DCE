function result = bench_hybrid_memory(fitFile,outDir,mode,workers,nrep,repeats)
% Dedicated-session benchmark of a saved HMXL fit; native sampler includes workers.
if nargin < 6, repeats = 5; end
assert(~isfolder(outDir),'Use a new output directory.'); mkdir(outDir);
outDir = char(java.io.File(outDir).getCanonicalPath());
saved = load(fitFile,'fit'); fit = saved.fit; I = fit.INPUT; o = fit.EstimOpt;
if nargin >= 5 && ~isempty(nrep), o.NRep = nrep; end
rng(o.Seed1); h = sobolset(o.NVarA+o.NLatent,'Skip',o.HaltonSkip,'Leap',o.HaltonLeap);
h = scramble(h,'MatousekAffineOwen'); err = norminv(net(h,o.NP*o.NRep))';
o.GPU = mode;
if strcmp(mode,'legacy'), o.HybridBlockSize = 0; o.GPU = 'cpu'; end
args = {I.YY,I.XXa,I.Xm,I.Xs,I.Xstr,I.Xmea,I.Xmea_exp,err,o,fit.bhat};
clear saved fit I err;
settings = parallel.Settings; previousAuto = settings.Pool.AutoCreate;
restoreAuto = onCleanup(@() setAutoCreate(settings,previousAuto)); %#ok<NASGU>
settings.Pool.AutoCreate = false;
assert(isempty(gcp('nocreate')),'Use a dedicated session.');
pids = feature('getpid');
if workers > 0
    pool = parpool('Processes',workers);
    spmd, workerPid = feature('getpid'); maxNumCompThreads(2); end
    pids = [pids,workerPid{:}];
end
for i = 1:2
    [f,g,backend] = evaluate(args,mode);
    assert(all(isfinite([f(:);g(:)])),'Nonfinite benchmark evaluation.');
end
if strcmp(mode,'gpu'), assert(strcmp(backend,'gpu'),'The requested GPU benchmark fell back to CPU.'); end
write_json(fullfile(outDir,'ready.json'),struct('pids',pids,'NP',o.NP));
sampler = System.Diagnostics.Process;
sampler.StartInfo.FileName = 'powershell.exe';
sampler.StartInfo.Arguments = sprintf('-NoProfile -ExecutionPolicy Bypass -File "%s" -OutDir "%s" -AttachPid %d -SampleMilliseconds 200',...
    fullfile(fileparts(mfilename('fullpath')),'measure_mxl_memory.ps1'),outDir,pids(1));
sampler.StartInfo.UseShellExecute = false; sampler.StartInfo.CreateNoWindow = true;
assert(sampler.Start(),'Unable to start native memory sampler.');
wait_for(fullfile(outDir,'sample.start'),120);
seconds = []; timer = tic;
while numel(seconds) < repeats || toc(timer) < 20
    one = tic; [f,g] = evaluate(args,mode); seconds(end+1) = toc(one); %#ok<AGROW>
end
result = struct('mode',mode,'backend',backend,'workers',workers,'NP',o.NP,'NRep',o.NRep,...
    'parameters',numel(args{end}),'evalSeconds',seconds,'medianEvalSeconds',median(seconds),...
    'evaluations',numel(seconds),'negativeLL',sum(f),'pids',pids,'MATLAB',version);
write_json(fullfile(outDir,'evaluation.json'),result);
wait_for(fullfile(outDir,'sample.done'),120);
save(fullfile(outDir,'evaluation.mat'),'result','f','g');
if workers > 0, delete(pool); end
fprintf('BENCHMARK_COMPLETE %s NP=%d R=%d median=%.6fs\n',mode,o.NP,o.NRep,result.medianEvalSeconds);
end

function [f,g,backend] = evaluate(args,mode)
if strcmp(mode,'legacy')
    [f,g] = LL_hmxl(args{:}); backend = 'legacy';
else
    [f,g,backend] = hybrid_likelihood('HMXL',args,true);
end
end

function write_json(path,value)
fid = fopen(path,'w'); assert(fid >= 0); cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(value));
end

function wait_for(path,seconds)
t = tic;
while ~isfile(path), assert(toc(t) < seconds,'Memory sampler handshake timed out.'); pause(.1); end
end

function setAutoCreate(settings,value)
settings.Pool.AutoCreate = value;
end
