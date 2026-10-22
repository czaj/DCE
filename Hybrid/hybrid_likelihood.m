function [f,g,backend] = hybrid_likelihood(model,args,need)
% One exact-input cache; pooled moments stay fixed across respondent blocks.
persistent previous data resident decisions device failed poolPrevious workers
o = args{end-1}; b = args{end}; inputs = args(1:end-1);
mode = 'auto';
if isfield(o,'GPU') && ~isempty(o.GPU)
    if isequal(o.GPU,false) || isequal(o.GPU,0), mode = 'cpu';
    elseif isequal(o.GPU,true) || isequal(o.GPU,1), mode = 'gpu';
    elseif (ischar(o.GPU) && isrow(o.GPU)) || (isstring(o.GPU) && isscalar(o.GPU)), mode = lower(char(o.GPU));
    else, error('DCE:Hybrid:GPUOption','GPU must be auto, cpu, gpu, false or true.'); end
end
assert(ismember(mode,{'auto','cpu','gpu'}),'DCE:Hybrid:GPUOption','GPU must be auto, cpu or gpu.');
pool = [];
if license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()), pool = gcp('nocreate'); end
key = {model,inputs};
if isempty(previous) || ~isequaln(previous,key) || ~isequal(pool,poolPrevious)
    if ~isempty(workers) && isvalid(workers), delete(workers); end
    data = hybrid_prepare(model,args);
    previous = key; poolPrevious = pool; resident = []; decisions = [0 0]; failed = false; workers = [];
    if ~isempty(pool), workers = parallel.pool.Constant(data); end
end
choiceRows = 0;
if isfield(o,'NAlt') && isfield(o,'NCT'), choiceRows = o.NAlt*o.NCT; end
perPerson = 8*o.NRep*(4*numel(b)+8*o.NLatent+16*choiceRows+32);
blockSize = min(32,max(1,floor(64*1024^2/perPerson)));
if isfield(o,'HybridBlockSize') && ~isempty(o.HybridBlockSize), blockSize = o.HybridBlockSize; end
assert(isscalar(blockSize) && isfinite(blockSize) && blockSize >= 1 && blockSize == fix(blockSize),...
    'DCE:Hybrid:BlockSize','HybridBlockSize must be a positive integer (or 0 for the legacy path).');
nWorkers = 0; if ~isempty(pool), nWorkers = pool.NumWorkers; end
cpu = @() hybrid_evaluate(data,b,blockSize,need,nWorkers,workers);
backend = 'cpu'; slot = 1+need;
if strcmp(mode,'cpu') || failed || decisions(slot) < 0 ||...
        ~license('test','Distrib_Computing_Toolbox') || exist('gpuDeviceCount','file') == 0 || ~isempty(getCurrentTask())
    [f,g] = cpu(); return;
end
try
    if gpuDeviceCount('available') == 0, failed = true; [f,g] = cpu(); return; end
    current = gpuDevice;
    if isempty(device) || current.Index ~= device.Index
        resident = []; decisions = [0 0]; failed = false;
    end
    device = current;
    if ~device.SupportsDouble, failed = true; [f,g] = cpu(); return; end
    if isempty(resident)
        names = {'E','Xstr','Xexp','X','M','Xs','Xc','Y','A'};
        bytes = 0;
        for name = names
            if isfield(data,name{1}), bytes = bytes+8*numel(data.(name{1})); end
        end
        budget = .75*double(device.AvailableMemory)-bytes-256*1024^2;
        gpuBlock = min([128,o.NP,floor(budget/perPerson)]);
        if isfield(o,'HybridBlockSize') && ~isempty(o.HybridBlockSize), gpuBlock = min(gpuBlock,blockSize); end
        if gpuBlock < 1, failed = true; [f,g] = cpu(); return; end
        resident = data;
        for name = names
            if isfield(resident,name{1}), resident.(name{1}) = gpuArray(resident.(name{1})); end
        end
        resident.blockSize = gpuBlock;
    end
    gpu = @() hybrid_evaluate(resident,b,resident.blockSize,need,0,[]);
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
catch exception
    if ~startsWith(exception.identifier,'parallel:gpu:') && ~startsWith(exception.identifier,'gpuArray:') &&...
            ~ismember(exception.identifier,{'MATLAB:nomem','MATLAB:OutOfMemory'})
        rethrow(exception);
    end
    resident = []; failed = true; [f,g] = cpu();
end
end

function r = relative(a,b)
if isempty(a) && isempty(b), r = 0; else, r = max(abs(a-b),[],'all')/max(1,max(abs(b),[],'all')); end
end
