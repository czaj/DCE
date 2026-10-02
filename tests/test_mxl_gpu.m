function [report,raw] = test_mxl_gpu(outDir,fixtureFile,blockSizes)
% Compare the experimental GPU likelihood with the serial CPU implementation.
% fixtureFile must contain cases saved by test_mxl_extended. Regenerate with
% test_mxl_extended(outDir,[0 3],'6bf56b8',baselineDir) when needed.
repo = fileparts(fileparts(mfilename('fullpath')));
documents = fileparts(fileparts(fileparts(repo)));
if nargin < 1 || isempty(outDir), outDir = fullfile(tempdir,'DCE_mxl_gpu'); end
if nargin < 2 || isempty(fixtureFile)
    fixtureFile = fullfile(documents,'_dce_memory','extended_cluster_20261002',...
        'output','extended_final','fixtures.mat');
end
if nargin < 3 || isempty(blockSizes), blockSizes = [1 3 5]; end
assert(all(isfinite(blockSizes)) && all(blockSizes >= 1) &&...
    all(blockSizes == fix(blockSizes)),'Block sizes must be positive integers.');
outDir = char(java.io.File(outDir).getCanonicalPath());
frozen = fullfile(documents,'lasy','replication_package');
assert(~strcmpi(outDir,frozen) && ~startsWith(lower(outDir),[lower(frozen),filesep]),...
    'GPU test output must be outside replication_package.');
assert(exist(fixtureFile,'file') == 2,...
    'Missing fixtures. Run test_mxl_extended first or supply fixtureFile.');
assert(isempty(gcp('nocreate')),'Use a dedicated session with no existing pool.');
settings = parallel.Settings;
oldAutoCreate = settings.Pool.AutoCreate;
settings.Pool.AutoCreate = false;
restoreSettings = onCleanup(@() setAutoCreate(settings,oldAutoCreate));
addpath(fullfile(repo,'MXL'),fullfile(repo,'tests'));
assert(exist('mxl_gpu_prepare','file') == 2 && exist('LL_mxl_gpu','file') == 2,...
    'The experimental mxl_gpu_prepare and LL_mxl_gpu helpers are required.');
if ~exist(outDir,'dir'), mkdir(outDir); end
device = gpuDevice;
environment = struct('MATLAB',version,'GPU',device.Name,...
    'ComputeCapability',device.ComputeCapability,'TotalMemory',double(device.TotalMemory));
loaded = load(fixtureFile,'cases');
assert(isfield(loaded,'cases'),'Fixture MAT file must contain cases.');
cases = loaded.cases(strcmp({loaded.cases.Model},'MXL'));
assert(~isempty(cases),'No MXL fixtures were found.');
base = [];
singleton = [];
fixed = [];
for k = 1:numel(cases)
    opt = cases(k).Args{end-1};
    if opt.NP == 1 && opt.NRep == 1, singleton = cases(k); end
    if opt.WTP_space == 0 && opt.FullCov == 0 && opt.NVarM == 0 &&...
            opt.NVarS == 0 && opt.Dist(1) == 0
        fixed = cases(k);
    end
    if opt.WTP_space == 1 && opt.FullCov == 1 && opt.mCT == 1 &&...
            opt.NVarM == 1 && opt.NVarS == 1 && any(isnan(cases(k).Args{1}),'all')
        base = cases(k);
    end
end
assert(~isempty(base),'A row-varying WTP/full-covariance fixture is required.');
extra = twoCovariates(base);
cases(end+1) = extra;
cases(end+1) = missingRespondent(extra);
assert(~isempty(singleton),'A singleton MXL fixture is required.');
cases(end+1) = subnormalMean(singleton);
extra = cases(end);
extra.Name = [extra.Name,'_full_covariance'];
extra.Args{6}.FullCov = 1;
extra.Args{5} = .7;
cases(end+1) = extra;
extra.Name = [extra.Name,'_realmin'];
extra.Args{6}.RealMin = 1;
cases(end+1) = extra;
assert(~isempty(fixed),'A simple normal MXL fixture is required.');
opt = fixed.Args{6};
fixed.Name = [fixed.Name,'_gpu_fixed_zero_draws'];
fixed.Args{6}.Dist(1) = -1;
draws = reshape(fixed.Args{5},[opt.NVarA,opt.NRep,opt.NP]);
draws(1,:,:) = 0;
fixed.Args{5} = reshape(draws,[opt.NVarA,opt.NRep*opt.NP]);
cases(end+1) = fixed;
rows = struct([]);
raw = struct([]);
for k = 1:numel(cases)
    fixture = cases(k);
    C = asInput(fixture.Args);
    [cpuF,cpuG] = LL_mxl(fixture.Args{:});
    cpuValue = LL_mxl(fixture.Args{:});
    checkOutputs(cpuF,cpuG,cpuValue,C,fixture.Name);
    prepared = mxl_gpu_prepare(C);
    for blockSize = blockSizes(:)'
        [gpuF,gpuG] = LL_mxl_gpu(prepared,C.b,blockSize);
        gpuValue = LL_mxl_gpu(prepared,C.b,blockSize);
        gpuF = gather(gpuF);
        gpuG = gather(gpuG);
        gpuValue = gather(gpuValue);
        checkOutputs(gpuF,gpuG,gpuValue,C,fixture.Name);
        r = struct('Case',fixture.Name,'BlockSize',blockSize,...
            'LLRelativeError',relative(sum(gpuF),sum(cpuF)),...
            'ValueRelativeError',relative(gpuF,cpuF),...
            'GradientRelativeError',relative(gpuG,cpuG),...
            'GradientSumRelativeError',relative(sum(gpuG,1),sum(cpuG,1)),...
            'ValueOnlyRelativeError',relative(gpuValue,cpuValue),...
            'GPUBranchRelativeError',relative(gpuF,gpuValue),'Passed',true);
        assert(r.LLRelativeError <= 1e-8 && r.ValueRelativeError <= 1e-8 &&...
            r.ValueOnlyRelativeError <= 1e-8 && r.GPUBranchRelativeError <= 1e-8,...
            'GPU likelihood mismatch: %s, block %d.',fixture.Name,blockSize);
        assert(r.GradientRelativeError <= 1e-6 && r.GradientSumRelativeError <= 1e-6,...
            'GPU gradient mismatch: %s, block %d.',fixture.Name,blockSize);
        if contains(fixture.Name,'_all_missing_respondent')
            assert(abs(gpuF(2)) <= 1e-14 && all(abs(gpuG(2,:)) <= 1e-14),...
                'An all-missing respondent must have neutral likelihood and gradient.');
        end
        j = numel(rows)+1;
        if isempty(rows), rows = r; else, rows(j) = r; end
        result = struct('Case',fixture.Name,'BlockSize',blockSize,'Input',C,...
            'CPUValue',cpuValue,'CPUf',cpuF,'CPUg',cpuG,...
            'GPUValue',gpuValue,'GPUf',gpuF,'GPUg',gpuG);
        if isempty(raw), raw = result; else, raw(j) = result; end
        report = struct2table(rows,'AsArray',true);
        writetable(report,fullfile(outDir,'gpu_correctness.csv'));
        save(fullfile(outDir,'gpu_correctness.mat'),'report','raw','fixtureFile',...
            'blockSizes','environment');
    end
end
guards = struct('ExpBRejected',rejectsOption(asInput(base.Args),'ExpB',1),...
    'NLTRejected',rejectsOption(asInput(base.Args),'NVarNLT',1));
assert(guards.ExpBRejected && guards.NLTRejected,...
    'The GPU prototype must reject unsupported ExpB and NLT options.');
save(fullfile(outDir,'gpu_correctness.mat'),'report','raw','fixtureFile',...
    'blockSizes','environment','guards');
fprintf('GPU correctness passed: %d fixtures, %d block comparisons.\n',...
    numel(cases),height(report));
clear restoreSettings;
end

function C = asInput(args)
C = struct('YY',args{1},'XXa',args{2},'XXm',args{3},'Xs',args{4},...
    'err',args{5},'EstimOpt',args{6},'b',args{7});
end

function C = twoCovariates(C)
opt = C.Args{6};
K = opt.NVarA;
rows = opt.NAlt*opt.NCT;
indices = 1:rows*opt.NP;
unavailable = isnan(C.Args{1}(:));
secondMean = .11*cos(indices/3)+.013*mod(indices,rows);
secondMean(unavailable) = NaN;
secondScale = .09*sin(indices'/4)+.005*indices';
secondScale(unavailable) = NaN;
C.Args{3} = [C.Args{3};secondMean];
C.Args{4} = [C.Args{4},secondScale];
head = K+K+opt.FullCov*K*(K-1)/2;
b = C.Args{7};
C.Args{7} = [b(1:head);b(head+(1:K));.03;-.02;.04;b(end);.05];
C.Args{6}.NVarM = 2;
C.Args{6}.NVarS = 2;
C.Name = [C.Name,'_gpu_two_mean_two_scale'];
end

function C = missingRespondent(C)
opt = C.Args{6};
rows = opt.NAlt*opt.NCT;
indices = rows+(1:rows);
C.Args{1}(:,2) = NaN;
C.Args{2}(:,:,2) = NaN;
C.Args{3}(:,indices) = NaN;
C.Args{4}(indices,:) = NaN;
C.Args{6}.MissingCT(:,2) = true;
C.Args{6}.NCTMiss(2) = 0;
C.Args{6}.NAltMiss(2) = 0;
C.Args{6}.NAltMissInd(:,2) = opt.NAlt;
C.Args{6}.NAltMissIndExp(:,2) = opt.NAlt;
C.Name = [C.Name,'_all_missing_respondent'];
end

function C = subnormalMean(C)
opt = C.Args{6};
opt.NVarA = 1;
opt.NAlt = 2;
opt.NCT = 1;
opt.NP = 1;
opt.NRep = 1;
opt.NVarM = 1;
opt.NVarS = 0;
opt.FullCov = 0;
opt.Dist = 0;
opt.WTP_space = 0;
opt.WTP_matrix = [];
opt.mCT = 0;
opt.RealMin = 0;
opt.indx1 = 1;
opt.indx2 = 1;
opt.DiagIndex = 1;
opt.MissingCT = false;
opt.NCTMiss = 1;
opt.NAltMiss = 2;
opt.NAltMissInd = 2;
opt.NAltMissIndExp = [2;2];
C.Args = {[1;0],[.7;0],2,zeros(2,0),0,opt,...
    [log(realmin*eps)/.7;0;0]};
C.Name = 'MXL_subnormal_respondent_mean';
end

function checkOutputs(f,g,value,C,name)
assert(isequal(size(f),[C.EstimOpt.NP,1]) &&...
    isequal(size(value),[C.EstimOpt.NP,1]) &&...
    isequal(size(g),[C.EstimOpt.NP,numel(C.b)]),...
    'Unexpected output dimensions: %s.',name);
assert(all(isfinite([f(:);g(:);value(:)])),...
    'Non-finite output: %s.',name);
assert(relative(f,value) <= 1e-8,'Value/gradient branch mismatch: %s.',name);
end

function value = relative(actual,reference)
assert(isequal(size(actual),size(reference)),'Comparison dimensions differ.');
value = max(abs(actual-reference),[],'all')/max(1,max(abs(reference),[],'all'));
end

function rejected = rejectsOption(C,field,value)
C.EstimOpt.(field) = value;
if strcmp(field,'ExpB')
    expectedIdentifier = 'DCE:GPU:UnsupportedExpB';
else
    expectedIdentifier = 'DCE:GPU:UnsupportedNLT';
end
rejected = false;
try
    mxl_gpu_prepare(C);
catch exception
    if strcmp(exception.identifier,expectedIdentifier)
        rejected = true;
    else
        rethrow(exception);
    end
end
end

function setAutoCreate(settings,value)
settings.Pool.AutoCreate = value;
end
