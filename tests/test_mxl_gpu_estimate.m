function [report,result] = test_mxl_gpu_estimate(outDir,blockSize,maxSeconds)
% Bounded CH optimizer check using only the experimental GPU likelihood.
% Root must verify local idleness before launching this dedicated session.
repo = fileparts(fileparts(mfilename('fullpath')));
documents = fileparts(fileparts(fileparts(repo)));
if nargin < 1 || isempty(outDir)
    outDir = fullfile(documents,'_dce_memory',...
        ['gpu_estimate_CH_',datestr(now,'yyyymmdd_HHMMSS')]);
end
if nargin < 2 || isempty(blockSize), blockSize = 128; end
if nargin < 3 || isempty(maxSeconds), maxSeconds = 120; end
assert(isscalar(blockSize) && isfinite(blockSize) && blockSize >= 1 &&...
    blockSize == fix(blockSize),'Block size must be a positive integer.');
assert(isscalar(maxSeconds) && isfinite(maxSeconds) &&...
    maxSeconds > 0 && maxSeconds <= 120,'The optimizer budget must be at most 120 seconds.');
outDir = char(java.io.File(outDir).getCanonicalPath());
frozen = char(java.io.File(fullfile(documents,'lasy','replication_package')).getCanonicalPath());
assert(~strcmpi(outDir,frozen) &&...
    ~startsWith(lower(outDir),[lower(frozen),filesep]),...
    'GPU estimation output must be outside replication_package.');
if isfolder(outDir)
    entries = dir(outDir);
    assert(isempty(entries(~ismember({entries.name},{'.','..'}))),...
        'Use a fresh output directory for every GPU estimation check.');
else
    mkdir(outDir);
end
assert(isempty(gcp('nocreate')),'Use a dedicated session with no existing pool.');
settings = parallel.Settings;
oldAutoCreate = settings.Pool.AutoCreate;
settings.Pool.AutoCreate = false;
restoreSettings = onCleanup(@() setAutoCreate(settings,oldAutoCreate));
addpath(fullfile(repo,'MXL'),fullfile(repo,'tests'));
inputFile = fullfile(documents,'_dce_memory','performance','input_CH.mat');
assert(isfile(inputFile),'The exact saved CH performance input is required.');
saved = load(inputFile,'C');
C = saved.C;
E = C.EstimOpt;
E.GPU = 'cpu';
b0 = C.b(:);
W = C.W(:);
assert(strcmp(C.Name,'CH') && E.NP == 644 && E.NRep == 1000 &&...
    E.NVarA == 17 && E.FullCov == 1 && E.WTP_space == 1,...
    'This optimizer check is restricted to the published CH specification.');
assert(isequal(b0,C.Published.bhat(:)) &&...
    abs(C.Published.LL+4793.8164) <= 5e-5,...
    'The saved input must start at the published CH estimates.');
assert(numel(W) == E.NP && all(isfinite(W)) && all(W >= 0),...
    'Invalid saved CH weights.');
assert(~isfield(E,'ConstVarActive') || isempty(E.ConstVarActive) ||...
    E.ConstVarActive == 0,'Equality-constrained CH is not supported by this check.');
assert(~isfield(E,'BActive') || isempty(E.BActive) ||...
    (numel(E.BActive) == numel(b0) && all(E.BActive(:) == 1)),...
    'Inactive parameters must not be silently optimized.');
assert(~any(E.Dist == -1) && ismember(E.RealMin,[0 1]),...
    'Fixed-coefficient masks and repaired objectives are outside this optimizer check.');
assert(isfield(C,'OptimOpt') && isa(C.OptimOpt,'optim.options.Fminunc'),...
    'Saved unconstrained fminunc options are required.');
options = C.OptimOpt;
options.Algorithm = 'quasi-newton';
options.GradObj = 'on';
options.Hessian = 'off';
options.Display = 'off';
options.OutputFcn = @stopForTime;
maxAbsB = 1000;
if isfield(E,'MaxAbsB') && ~isempty(E.MaxAbsB), maxAbsB = E.MaxAbsB; end
penalty = 1e50;
if isfield(E,'BadEvalPenalty') && ~isempty(E.BadEvalPenalty), penalty = E.BadEvalPenalty; end
penaltySlope = 1e6;
if isfield(E,'PenaltySlope') && ~isempty(E.PenaltySlope), penaltySlope = E.PenaltySlope; end
assert(all(isfinite(b0)) && max(abs(b0)) <= maxAbsB,...
    'Published CH parameters exceed the production parameter guard.');
totalTimer = tic;
device = gpuDevice;
uploadTimer = tic;
data = mxl_gpu_prepare(C);
wait(device);
prepareUploadSeconds = toc(uploadTimer);
evaluations = 0;
penaltyEvaluations = 0;
stoppedForTime = false;
[initialObjective,initialGradient] = objective(b0);
assert(abs(-initialObjective-C.Published.LL)/max(1,abs(C.Published.LL)) <= 1e-8,...
    'The initial GPU likelihood does not match the published CH result.');
optimizerTimer = tic;
[bhat,negativeLL,exitflag,output,optimizerGradient] = fminunc(@objective,b0,options);
optimizerSeconds = toc(optimizerTimer);
[gpuF,gpuG] = LL_mxl_gpu(data,bhat,blockSize);
gpuF = gather(gpuF);
gpuG = gather(gpuG);
[cpuF,cpuG] = LL_mxl(C.YY,C.XXa,C.XXm,C.Xs,C.err,E,bhat);
assert(all(isfinite([bhat(:);gpuF(:);gpuG(:);cpuF(:);cpuG(:)])),...
    'The final estimate or likelihood evaluation is non-finite.');
gpuLL = -sum(W.*gpuF);
cpuLL = -sum(W.*cpuF);
gpuGradient = sum(W.*gpuG,1)';
cpuGradient = sum(W.*cpuG,1)';
LLRelativeError = abs(gpuLL-C.Published.LL)/max(1,abs(C.Published.LL));
parameterRelativeError = max(abs(bhat-C.Published.bhat(:)))/...
    max(1,max(abs(C.Published.bhat(:))));
CPUValueRelativeError = relative(gpuF,cpuF);
CPUGradientRelativeError = relative(gpuG,cpuG);
CPUWeightedGradientRelativeError = relative(gpuGradient,cpuGradient);
objectiveRelativeError = abs(negativeLL+gpuLL)/max(1,abs(gpuLL));
passed = exitflag > 0 && ~stoppedForTime && LLRelativeError <= 1e-8 &&...
    parameterRelativeError <= 1e-6 && CPUValueRelativeError <= 1e-8 &&...
    CPUGradientRelativeError <= 1e-6 && CPUWeightedGradientRelativeError <= 1e-6 &&...
    objectiveRelativeError <= 1e-8;
report = table(blockSize,maxSeconds,optimizerSeconds,exitflag,stoppedForTime,...
    C.Published.LL,gpuLL,cpuLL,LLRelativeError,parameterRelativeError,...
    CPUValueRelativeError,CPUGradientRelativeError,CPUWeightedGradientRelativeError,passed,...
    'VariableNames',{'BlockSize','TimeBudgetSeconds','OptimizerSeconds','ExitFlag',...
    'StoppedForTime','PublishedLL','GPULL','CPUFinalLL','LLRelativeError',...
    'ParameterRelativeError','CPUValueRelativeError','CPUGradientRelativeError',...
    'CPUWeightedGradientRelativeError','Passed'});
recordedOptions = options;
recordedOptions.OutputFcn = [];
result = struct('inputFile',inputFile,'matlabVersion',version,...
    'GPU',device.Name,'ComputeCapability',device.ComputeCapability,...
    'Driver',device.GraphicsDriverVersion,'options',recordedOptions,'EstimOpt',E,...
    'outputFunctionPolicy','Stop at an optimizer callback after the time budget.',...
    'b0',b0,'bhat',bhat,'published',C.Published,'negativeLL',negativeLL,...
    'exitflag',exitflag,'output',output,'optimizerGradient',optimizerGradient,...
    'initialObjective',initialObjective,'initialGradient',initialGradient,...
    'gpuF',gpuF,'gpuG',gpuG,'cpuF',cpuF,'cpuG',cpuG,...
    'gpuGradient',gpuGradient,'cpuGradient',cpuGradient,...
    'maxPersonGradientDifference',max(abs(gpuG-cpuG),[],'all'),...
    'objectiveRelativeError',objectiveRelativeError,...
    'evaluations',evaluations,'penaltyEvaluations',penaltyEvaluations,...
    'prepareUploadSeconds',prepareUploadSeconds,'totalSeconds',toc(totalTimer),...
    'scope','GPU likelihood optimizer only; no MXL Hessian, standard errors or CPU reoptimization.');
save(fullfile(outDir,'gpu_CH_estimation.mat'),'report','result');
writetable(report,fullfile(outDir,'gpu_CH_estimation.csv'));
disp(report);
assert(passed,'GPU CH estimation failed: exit %d, LL error %.3g, parameter error %.3g.',...
    exitflag,LLRelativeError,parameterRelativeError);
clear restoreSettings;

    function [f,g] = objective(b)
        evaluations = evaluations+1;
        assert(all(isfinite(b)),'Non-finite optimizer candidate.');
        if max(abs(b)) > maxAbsB
            % Match LL_mxl_MATlike's parameter-limit penalty, without diagnostics.
            penaltyEvaluations = penaltyEvaluations+1;
            excess = max(0,abs(b)-maxAbsB);
            f = penalty+penaltySlope*sum(excess.^2);
            g = 2*penaltySlope*excess.*sign(b);
            return
        end
        [fv,j] = LL_mxl_gpu(data,b,blockSize);
        fv = gather(fv);
        j = gather(j);
        assert(all(isfinite([fv(:);j(:)])),'Non-finite GPU objective.');
        f = sum(W.*fv);
        g = sum(W.*j,1)';
    end

    function stop = stopForTime(~,~,state)
        stop = ~strcmp(state,'done') && toc(optimizerTimer) >= maxSeconds;
        stoppedForTime = stoppedForTime || stop;
    end
end

function value = relative(actual,reference)
value = max(abs(actual-reference),[],'all')/max(1,max(abs(reference),[],'all'));
end

function setAutoCreate(settings,value)
settings.Pool.AutoCreate = value;
end
