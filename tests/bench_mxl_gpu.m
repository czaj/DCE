function result = bench_mxl_gpu(caseName,outDir,blockSizes,phase,inputFile)
% Prototype benchmark in a fresh session. Root must first verify local idleness.
% Split into cpu_serial/gpu/cpu_parallel phases and one GPU block size per
% invocation to keep local runs below the shared five-minute offload threshold.
if nargin < 3 || isempty(blockSizes), blockSizes = 128; end
if nargin < 4 || isempty(phase), phase = 'gpu'; end
phase = char(phase);
caseName = char(caseName);
assert(ismember(caseName,{'CH','pooled'}),'Use CH or pooled saved benchmark inputs.');
assert(ismember(phase,{'all','cpu_serial','gpu','cpu_parallel'}),'Unknown benchmark phase.');
assert(all(isfinite(blockSizes)) && all(blockSizes >= 1) &&...
    all(blockSizes == fix(blockSizes)),'Block sizes must be positive integers.');
repo = fileparts(fileparts(mfilename('fullpath')));
documents = fileparts(fileparts(fileparts(repo)));
addpath(fullfile(repo,'MXL'),fullfile(repo,'tests'));
if nargin < 5 || isempty(inputFile)
    inputFile = fullfile(documents,'_dce_memory','performance',['input_',caseName,'.mat']);
end
inputFile = char(inputFile);
assert(isfile(inputFile),'Saved performance input is required: %s',inputFile);
outDir = char(java.io.File(outDir).getCanonicalPath());
frozen = fullfile(documents,'lasy','replication_package');
assert(~strcmpi(outDir,frozen) && ~startsWith(lower(outDir),[lower(frozen),filesep]),...
    'GPU benchmark output must be outside replication_package.');
if isfolder(outDir)
    entries = dir(outDir);
    assert(isempty(entries(~ismember({entries.name},{'.','..'}))),...
        'Use a fresh output directory for every benchmark invocation.');
else
    mkdir(outDir);
end
assert(isempty(gcp('nocreate')),'Run the benchmark in a fresh session without a pool.');
settings = parallel.Settings;
oldAutoCreate = settings.Pool.AutoCreate;
settings.Pool.AutoCreate = false;
restoreSettings = onCleanup(@() set_auto_create(settings,oldAutoCreate));
saved = load(inputFile,'C');
C = saved.C;
cpu = @() LL_mxl(C.YY,C.XXa,C.XXm,C.Xs,C.err,C.EstimOpt,C.b);
[cpuF,cpuG] = cpu();
assert(all(isfinite(cpuF),'all') && all(isfinite(cpuG),'all'));
result = struct('caseName',caseName,'dataName',C.Name,'phase',phase,'inputFile',inputFile,...
    'matlabVersion',version,'matlabRelease',version('-release'),...
    'NP',C.EstimOpt.NP,'draws',C.EstimOpt.NRep,...
    'startedUTC',utc_now(),'endedUTC','','GPU',[],...
    'staticPrepareUploadSeconds',NaN,'staticNumericPayloadBytes',NaN,...
    'memorySnapshots',struct([]),'rows',struct([]),...
    'timingMethod','wait before tic and before toc; two f/g warmups; at least 5 calls and 20 seconds');
raw = struct([]);
if ismember(phase,{'all','cpu_serial'})
    [times,f,g,window] = measure(cpu,[],5,20);
    append('CPU_serial',NaN,0,times,f,g,NaN,window);
end
if ismember(phase,{'all','gpu'})
    initialization = tic;
    device = gpuDevice;
    wait(device);
    result.GPU = struct('Name',device.Name,'ComputeCapability',device.ComputeCapability,...
        'GraphicsDriverVersion',device.GraphicsDriverVersion,'DriverModel',device.DriverModel,...
        'CachePolicy',device.CachePolicy,'TotalMemory',double(device.TotalMemory),...
        'SingleDoubleRatio',double(device.SingleDoubleRatio),...
        'InitializationSeconds',toc(initialization));
    snapshot('before_static_upload',NaN);
    upload = tic;
    data = mxl_gpu_prepare(C);
    wait(device);
    result.staticPrepareUploadSeconds = toc(upload);
    result.staticNumericPayloadBytes = payload_bytes(data);
    snapshot('after_static_upload',NaN);
    for blockIndex = 1:numel(blockSizes)
        blockSize = blockSizes(blockIndex);
        gpu = @() LL_mxl_gpu(data,C.b,blockSize);
        [checkF,checkG] = gpu();
        checkF = gather(checkF);
        checkG = gather(checkG);
        onlyF = gather(LL_mxl_gpu(data,C.b,blockSize));
        validate(checkF,checkG,onlyF,cpuF,cpuG,C);
        [times,f,g,window] = measure(gpu,device,5,20);
        append('GPU_resident_device',blockSize,0,times,gather(f),gather(g),...
            relative(onlyF,cpuF),window);
        snapshot('after_resident_device',blockSize);
        practical = @() gather_evaluation(data,C.b,blockSize);
        [times,f,g,window] = measure(practical,device,5,20);
        append('GPU_resident_gather',blockSize,0,times,f,g,relative(onlyF,cpuF),window);
        snapshot('after_resident_gather',blockSize);
        % Three separate cold-data calls include static preparation and upload.
        clear data gpu practical;
        wait(device);
        inclusive = @() transfer_evaluation(C,blockSize);
        [times,f,g,window] = measure(inclusive,device,3,0,0);
        append('GPU_transfer_inclusive',blockSize,0,times,f,g,relative(onlyF,cpuF),window);
        snapshot('after_transfer_inclusive',blockSize);
        if blockIndex < numel(blockSizes)
            data = mxl_gpu_prepare(C);
            wait(device);
        end
    end
    clear data f g gpu practical inclusive;
end
if ismember(phase,{'all','cpu_parallel'})
    pool = parpool('Processes',3);
    closePool = onCleanup(@() delete(pool));
    [times,f,g,window] = measure(cpu,[],5,20);
    append('CPU_pool3',NaN,3,times,f,g,NaN,window);
    clear closePool;
end
result.endedUTC = utc_now();
report = struct2table(result.rows,'AsArray',true);
% Delay disk writes until all timed series finish, avoiding sync interference.
writetable(report,fullfile(outDir,'gpu_benchmark.csv'));
fid = fopen(fullfile(outDir,'gpu_benchmark.json'),'w');
assert(fid ~= -1,'Cannot write GPU benchmark metadata.');
closeFile = onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(result));
clear closeFile;
save(fullfile(outDir,'gpu_benchmark.mat'),'result','report','raw','cpuF','cpuG','blockSizes');
clear restoreSettings;

    function append(mode,blockSize,workers,times,f,g,valueError,window)
        errors = validate(f,g,f,cpuF,cpuG,C);
        row = struct('Case',C.Name,'Mode',mode,'BlockSize',blockSize,...
            'Workers',workers,'Calls',numel(times),'MedianSeconds',median(times),...
            'TotalSeconds',sum(times),'StartedUTC',window.StartedUTC,...
            'EndedUTC',window.EndedUTC,'WeightedLL',-sum(C.W(:).*f(:)),...
            'LLRelativeError',errors.LL,'ValueRelativeError',errors.Value,...
            'GradientRelativeError',errors.Gradient,...
            'WeightedGradientRelativeError',errors.WeightedGradient,...
            'ValueOnlyRelativeError',valueError);
        one = struct('Mode',mode,'BlockSize',blockSize,'Seconds',times,'f',f,'g',g);
        if isempty(result.rows)
            result.rows = row;
            raw = one;
        else
            result.rows(end+1) = row;
            raw(end+1) = one;
        end
    end

    function snapshot(label,blockSize)
        % Current allocator/device availability, NOT a sampled VRAM peak.
        one = struct('Label',label,'BlockSize',blockSize,...
            'AvailableBytes',double(device.AvailableMemory),'TimestampUTC',utc_now());
        if isempty(result.memorySnapshots)
            result.memorySnapshots = one;
        else
            result.memorySnapshots(end+1) = one;
        end
    end
end

function [seconds,f,g,window] = measure(fun,device,repeats,minSeconds,warmups)
if nargin < 5, warmups = 2; end
for k = 1:warmups
    [f,g] = fun();
    if ~isempty(device), wait(device); end
end
seconds = [];
window = struct('StartedUTC',utc_now(),'EndedUTC','');
series = tic;
while numel(seconds) < repeats || toc(series) < minSeconds
    if ~isempty(device), wait(device); end
    one = tic;
    [f,g] = fun();
    if ~isempty(device), wait(device); end
    seconds(end+1) = toc(one); %#ok<AGROW>
end
window.EndedUTC = utc_now();
end

function [f,g] = gather_evaluation(data,b,blockSize)
[f,g] = LL_mxl_gpu(data,b,blockSize);
f = gather(f);
g = gather(g);
end

function [f,g] = transfer_evaluation(C,blockSize)
data = mxl_gpu_prepare(C);
[f,g] = gather_evaluation(data,C.b,blockSize);
end

function errors = validate(f,g,onlyF,cpuF,cpuG,C)
assert(isequal(size(f),size(cpuF)) && isequal(size(g),size(cpuG)) &&...
    all(isfinite(f),'all') && all(isfinite(g),'all'),'Invalid GPU output dimensions or values.');
errors = struct('LL',relative(sum(C.W(:).*f(:)),sum(C.W(:).*cpuF(:))),...
    'Value',relative(f,cpuF),'Gradient',relative(g,cpuG),...
    'WeightedGradient',relative(sum(C.W(:).*g,1),sum(C.W(:).*cpuG,1)));
assert(errors.LL <= 1e-8 && errors.Value <= 1e-8 &&...
    relative(onlyF,cpuF) <= 1e-8 && relative(onlyF,f) <= 1e-8,...
    'GPU likelihood or value-only branch differs from the current CPU control.');
assert(errors.Gradient <= 1e-6 && errors.WeightedGradient <= 1e-6,...
    'GPU gradient differs from the current CPU control.');
if strcmp(C.Name,'CH')
    assert(abs(-sum(C.W(:).*f(:))+4793.8164) < 5e-5,...
        'The CH published likelihood check failed.');
end
end

function n = payload_bytes(data)
fields = {'X','Y','A','missing','E','M','Xs','chosenX'};
n = 0;
for k = 1:numel(fields)
    x = data.(fields{k});
    n = n + numel(x)*(1+7*~strcmp(classUnderlying(x),'logical'));
end
end

function r = relative(actual,expected)
r = max(abs(actual(:)-expected(:)))/max(1,max(abs(expected(:))));
end

function set_auto_create(settings,value)
settings.Pool.AutoCreate = value;
end

function value = utc_now()
value = char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
end
