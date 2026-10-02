function result = bench_mxl_memory(caseName, functionName, outDir, replicationRoot, baselineDir, repeats, minSeconds)
% Called by measure_mxl_memory.ps1; warm-up and pool startup are not timed.
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'MXL'));
if ~isempty(baselineDir)
    addpath(baselineDir);
end
sourceImplementation = functionName;
if strcmp(functionName,'LL_mxl_baseline')
    sourceImplementation = 'dd6f704 (unmodified HEAD)';
elseif strcmp(functionName,'LL_mxl_baseline_parfor')
    source = fileread(fullfile(baselineDir,'LL_mxl_baseline.m'));
    header = 'function [f,g,h] = LL_mxl_baseline(';
    assert(startsWith(source,header),'Expected the renamed dd6f704 baseline.');
    first = strfind(source,'elseif nargout == 2 %% function value + gradient');
    last = strfind(source,'elseif nargout == 3 % function value + gradient + hessian');
    assert(isscalar(first) && isscalar(last) && first < last,'Unexpected baseline gradient section.');
    gradient = source(first:last-1);
    pattern = '(?m)^([ \t]*)for n = 1:NP';
    assert(numel(regexp(gradient,pattern)) == 2,'Expected exactly two respondent gradient loops.');
    gradient = regexprep(gradient,pattern,'$1parfor n = 1:NP');
    source = [source(1:first-1),gradient,source(last:end)];
    source = strrep(source,header,'function [f,g,h] = LL_mxl_baseline_parfor(');
    write_file(fullfile(outDir,'LL_mxl_baseline_parfor.m'),source);
    addpath(outDir);
    sourceImplementation = 'dd6f704 + parfor (only the two gradient respondent loops)';
end
inputFile = fullfile(fileparts(outDir),['input_',caseName,'.mat']);
if isfile(inputFile)
    saved = load(inputFile,'C');
    C = saved.C;
else
    C = mxl_memory_case(caseName,replicationRoot);
    save(inputFile,'C','-v7.3');
end
pool = gcp('nocreate');
assert(isempty(pool),'Run this benchmark in a fresh MATLAB batch session.');
pool = parpool('Processes',3);
spmd
    workerPid = feature('getpid');
end
pids = [feature('getpid'),workerPid{:}];
fun = str2func(functionName);
for k = 1:2
    [f,g] = fun(C.YY,C.XXa,C.XXm,C.Xs,C.err,C.EstimOpt,C.b);
    assert(all(isfinite(f),'all') && all(isfinite(g),'all'));
end
write_file(fullfile(outDir,'ready.json'),jsonencode(struct('pids',pids,'NP',C.EstimOpt.NP)));
wait_for(fullfile(outDir,'sample.start'),120);
seconds = [];
started = utc_now();
series = tic;
while numel(seconds) < repeats || toc(series) < minSeconds
    one = tic;
    [f,g] = fun(C.YY,C.XXa,C.XXm,C.Xs,C.err,C.EstimOpt,C.b);
    seconds(end+1) = toc(one); %#ok<AGROW>
end
ended = utc_now();
result = struct('caseName',caseName,'functionName',functionName, ...
    'sourceImplementation',sourceImplementation, ...
    'NP',C.EstimOpt.NP,'workers',pool.NumWorkers,'matlabVersion',version, ...
    'matlabRelease',version('-release'),'startedUTC',started,'endedUTC',ended, ...
    'evalSeconds',seconds,'medianEvalSeconds',median(seconds), ...
    'totalEvalSeconds',sum(seconds),'evaluations',numel(seconds), ...
    'weightedNegativeLL',sum(C.W(:).*f(:)), ...
    'weightedLL',-sum(C.W(:).*f(:)), ...
    'weightedGradient',sum(C.W(:).*g,1),'pids',pids);
write_file(fullfile(outDir,'evaluation.json'),jsonencode(result));
wait_for(fullfile(outDir,'sample.done'),120);
save(fullfile(outDir,'evaluation.mat'),'result','f','g');
delete(pool);
end

function write_file(path,value)
fid = fopen([path,'.tmp'],'w');
assert(fid ~= -1,'Cannot write benchmark metadata.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s',value);
clear cleanup;
movefile([path,'.tmp'],path);
end

function wait_for(path,seconds)
timer = tic;
while ~isfile(path)
    assert(toc(timer) < seconds,'Performance sampler handshake timed out.');
    pause(.1);
end
end

function value = utc_now()
value = char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
end
