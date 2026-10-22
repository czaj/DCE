function [f,g,backend] = vce_likelihood(model,input,o,b,need)
% Exact-input cache and warmed CPU/GPU calibration, as in the hybrid evaluator.
persistent previous d workers poolPrevious gpuData device decisions failed
mode = 'auto';
if isfield(o,'GPU') && ~isempty(o.GPU)
    if isequal(o.GPU,false) || isequal(o.GPU,0), mode = 'cpu';
    elseif isequal(o.GPU,true) || isequal(o.GPU,1), mode = 'gpu';
    elseif (ischar(o.GPU) && isrow(o.GPU)) || (isstring(o.GPU) && isscalar(o.GPU)), mode = lower(char(o.GPU));
    else, error('DCE:VCE:GPUOption','GPU must be auto, cpu, gpu, false or true.'); end
end
assert(ismember(mode,{'auto','cpu','gpu'}),'DCE:VCE:GPUOption','GPU must be auto, cpu or gpu.');
pool = [];
if license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()), pool = gcp('nocreate'); end
key = {model,input,o};
if isempty(previous) || ~isequaln(previous,key) || ~isequal(pool,poolPrevious)
    if ~isempty(workers) && isvalid(workers), delete(workers); end
    d = vce_prepare(model,input,o); previous = key; poolPrevious = pool;
    workers = []; gpuData = []; decisions = [0 0]; failed = false;
    if ~isempty(pool), workers = parallel.pool.Constant(d); end
end
b = b(:);
assert(isa(b,'double') && isreal(b) && all(isfinite(b)),'DCE:VCE:Parameters','Parameters must be finite real doubles.');
perPerson = 8*d.R*(24*o.NAlt*d.T+8*o.NVarA+16)+8*numel(b)*d.T;
block = min(d.N,max(1,floor(64*1024^2/perPerson)));
if d.mixed, block = min(64,block); end
if isfield(o,'VCEBlockSize') && ~isempty(o.VCEBlockSize), block = o.VCEBlockSize; end
assert(isscalar(block) && isfinite(block) && block >= 1 && block == fix(block),...
    'DCE:VCE:BlockSize','VCEBlockSize must be a positive integer.');
nWorkers = 0; if ~isempty(pool), nWorkers = pool.NumWorkers; end
cpu = @() vce_evaluate(d,b,block,need,nWorkers,workers);
backend = 'cpu'; slot = 1+need;
if strcmp(mode,'cpu') || failed || decisions(slot) < 0 || ...
        ~license('test','Distrib_Computing_Toolbox') || exist('gpuDeviceCount','file') == 0 || ~isempty(getCurrentTask())
    [f,g] = cpu(); return;
end
try
    if gpuDeviceCount('available') == 0, failed = true; [f,g] = cpu(); return; end
    current = gpuDevice;
    if isempty(device) || current.Index ~= device.Index, gpuData = []; decisions = [0 0]; end
    device = current;
    if ~device.SupportsDouble, failed = true; [f,g] = cpu(); return; end
    if isempty(gpuData)
        names = {'X','E','M','U','Y','price'}; bytes = 0;
        for name = names, bytes = bytes+8*numel(d.(name{1})); end
        budget = .75*double(device.AvailableMemory)-bytes-256*1024^2;
        gpuBlock = min(d.N,floor(budget/perPerson));
        if isfield(o,'VCEBlockSize') && ~isempty(o.VCEBlockSize), gpuBlock = min(gpuBlock,block); end
        if gpuBlock < 1, failed = true; [f,g] = cpu(); return; end
        gpuData = d;
        for name = names, gpuData.(name{1}) = gpuArray(d.(name{1})); end
        gpuData.block = gpuBlock;
    end
    gpu = @() vce_evaluate(gpuData,b,gpuData.block,need,0,[]);
    if decisions(slot) == 0
        [~,~] = cpu(); t = tic; [fc,gc] = cpu(); tc = toc(t);
        [~,~] = gpu(); wait(device); t = tic; [f,g] = gpu(); wait(device); tg = toc(t);
        agrees = all(isfinite([fc(:);gc(:);f(:);g(:)])) && relative(f,fc) < 1e-8 && relative(g,gc) < 1e-6;
        decisions(slot) = 2*(agrees && (strcmp(mode,'gpu') || tg < .9*tc))-1;
        if decisions(slot) < 0, f = fc; g = gc; return; end
    else
        [f,g] = gpu();
    end
    backend = 'gpu';
catch ex
    if ~startsWith(ex.identifier,'parallel:gpu:') && ~startsWith(ex.identifier,'gpuArray:') && ...
            ~ismember(ex.identifier,{'MATLAB:nomem','MATLAB:OutOfMemory'}), rethrow(ex); end
    gpuData = []; failed = true; [f,g] = cpu();
end
end

function r = relative(a,b)
if isempty(a) && isempty(b), r = 0; else, r = max(abs(a-b),[],'all')/max(1,max(abs(b),[],'all')); end
end
