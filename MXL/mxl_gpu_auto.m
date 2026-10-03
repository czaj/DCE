function [f,g,backend] = mxl_gpu_auto(C,b,needGradient,cpu)
% ponytail: cache one dataset/device; clear this function to retry a disabled GPU.
persistent previous poolPrevious resident device decisions blockSize available failed
backend = 'cpu';
mode = 'auto';
if isfield(C.EstimOpt,'GPU') && ~isempty(C.EstimOpt.GPU)
    setting = C.EstimOpt.GPU;
    if isequal(setting,false) || isequal(setting,0)
        mode = 'cpu';
    elseif isequal(setting,true) || isequal(setting,1)
        mode = 'gpu';
    elseif (ischar(setting) && isrow(setting)) || (isstring(setting) && isscalar(setting))
        mode = lower(char(setting));
    else
        error('DCE:MXL:GPUOption','EstimOpt.GPU must be auto, cpu, gpu, false or true.');
    end
end
assert(ismember(mode,{'auto','cpu','gpu'}),'DCE:MXL:GPUOption',...
    'EstimOpt.GPU must be auto, cpu or gpu.');
opt = C.EstimOpt;
supported = ismember(opt.FullCov,[0,1]) && all(ismember(opt.Dist,[-1,0,1])) &&...
    opt.NVarNLT == 0 && opt.Johnson == 0 && isempty(opt.ExpB) &&...
    (opt.WTP_space == 0 || all(ismember(opt.WTP_matrix,opt.NVarA-opt.WTP_space+1:opt.NVarA))) &&...
    isa(C.YY,'double') && isa(C.XXa,'double') && isa(C.err,'double') &&...
    (isempty(C.XXm) || isa(C.XXm,'double')) && (isempty(C.Xs) || isa(C.Xs,'double'));
if supported && opt.FullCov == 1 && needGradient
    [i,j] = find(tril(ones(opt.NVarA)));
    supported = isequal(opt.indx1(:),i) && isequal(opt.indx2(:),j);
end
if strcmp(mode,'cpu') || ~supported || ~license('test','Distrib_Computing_Toolbox') ||...
        exist('gpuDeviceCount','file') == 0 || exist('pagemtimes') == 0
    resident = []; previous = []; decisions = [];
    [f,g] = cpu();
    return
end
if ~isempty(getCurrentTask())
    [f,g] = cpu();
    return
end
pool = gcp('nocreate');
if isfield(C.EstimOpt,'GPU'), C.EstimOpt = rmfield(C.EstimOpt,'GPU'); end
if isempty(previous) || ~isequaln(C,previous) || ~isequal(pool,poolPrevious)
    resident = [];
    previous = C;
    poolPrevious = pool;
    decisions = [0,0];
    failed = false;
end
slot = 1+needGradient;
if isempty(available)
    try
        available = gpuDeviceCount('available') > 0;
    catch exception
        if ~gpu_failure(exception), rethrow(exception); end
        available = false;
    end
end
if ~available
    [f,g] = cpu();
    return
end
try
    current = gpuDevice;
    if ~isempty(device) && current.Index ~= device.Index
        resident = [];
        decisions = [0,0];
        failed = false;
    end
    device = current;
    if failed || ~device.DeviceAvailable || ~device.SupportsDouble ||...
            (strcmp(mode,'auto') && decisions(slot) < 0)
        [f,g] = cpu();
        return
    end
    if isempty(resident)
        K = opt.NVarA;
        Q = K+opt.FullCov*K*(K-1)/2;
        staticBytes = 8*(numel(C.XXa)+numel(C.err)+numel(C.XXm)+numel(C.Xs)+K*opt.NP)+...
            2*numel(C.YY)+opt.NCT*opt.NP;
        resultBytes = 8*opt.NP*(1+K+Q+K*opt.NVarM+opt.NVarS);
        perPersonBytes = 8*opt.NRep*((16+2*opt.WTP_space)*opt.NAlt*opt.NCT+6*K);
        budget = .8*double(device.AvailableMemory)-staticBytes-resultBytes-256*1024^2;
        blockSize = min([128,opt.NP,floor(budget/perPersonBytes)]);
        if blockSize < 1
            decisions(:) = -1;
            [f,g] = cpu();
            return
        end
        resident = mxl_gpu_prepare(C,needGradient);
    end
    if strcmp(mode,'auto') && decisions(slot) == 0
        [~,~] = cpu();
        timer = tic;
        [cpuF,cpuG] = cpu();
        cpuSeconds = toc(timer);
        [~,~] = evaluate(resident,b,blockSize,needGradient,device);
        wait(device);
        timer = tic;
        [f,g] = evaluate(resident,b,blockSize,needGradient,device);
        gpuSeconds = toc(timer);
        finite = all(isfinite([cpuF(:);cpuG(:);f(:);g(:)]));
        if ~finite
            f = cpuF; g = cpuG;
            return
        end
        agrees = relative(f,cpuF) <= 1e-8 && (~needGradient || relative(g,cpuG) <= 1e-6);
        decisions(slot) = 2*(agrees && gpuSeconds < .9*cpuSeconds)-1;
        if decisions(slot) < 0
            f = cpuF; g = cpuG;
            if all(decisions < 0), resident = []; end
            return
        end
    else
        [f,g] = evaluate(resident,b,blockSize,needGradient,device);
    end
    backend = 'gpu';
catch exception
    if ~gpu_failure(exception), rethrow(exception); end
    resident = [];
    decisions(:) = -1;
    failed = true;
    [f,g] = cpu();
end
end

function [f,g] = evaluate(data,b,blockSize,needGradient,device)
if needGradient
    [f,g] = LL_mxl_gpu(data,b,blockSize);
    g = gather(g);
else
    f = LL_mxl_gpu(data,b,blockSize);
    g = zeros(data.opt.NP,0);
end
f = gather(f);
wait(device);
end

function yes = gpu_failure(exception)
yes = startsWith(exception.identifier,'parallel:gpu:') ||...
    startsWith(exception.identifier,'gpuArray:') ||...
    ismember(exception.identifier,{'MATLAB:nomem','MATLAB:OutOfMemory'});
end

function value = relative(actual,reference)
value = max(abs(actual-reference),[],'all')/max(1,max(abs(reference),[],'all'));
end
