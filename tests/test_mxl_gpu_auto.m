function [report,raw] = test_mxl_gpu_auto(outDir,fixtureFile,nWorkers)
% Public LL_mxl routing and automatic GPU/CPU fallback regressions.
% fixtureFile contains cases written by test_mxl_extended; GPU-free hosts work.
repo = fileparts(fileparts(mfilename('fullpath')));
documents = fileparts(fileparts(fileparts(repo)));
if nargin < 1 || isempty(outDir), outDir = fullfile(tempdir,'DCE_mxl_gpu_auto'); end
if nargin < 2 || isempty(fixtureFile)
    fixtureFile = fullfile(documents,'_dce_memory','extended_cluster_20261002',...
        'output','extended_final','fixtures.mat');
end
if nargin < 3 || isempty(nWorkers), nWorkers = 3; end
assert(isscalar(nWorkers) && isfinite(nWorkers) && nWorkers >= 0 &&...
    nWorkers == fix(nWorkers),'nWorkers must be a nonnegative integer.');
outDir = char(java.io.File(outDir).getCanonicalPath());
frozen = fullfile(documents,'lasy','replication_package');
assert(~strcmpi(outDir,frozen) && ~startsWith(lower(outDir),[lower(frozen),filesep]),...
    'Test output must be outside replication_package.');
assert(exist(fixtureFile,'file') == 2,'Supply fixtures saved by test_mxl_extended.');
assert(isempty(gcp('nocreate')),'Use a dedicated session with no existing pool.');
settings = parallel.Settings;
oldAutoCreate = settings.Pool.AutoCreate;
settings.Pool.AutoCreate = false;
restoreSettings = onCleanup(@() setAutoCreate(settings,oldAutoCreate));
addpath(fullfile(repo,'MXL'),fullfile(repo,'tests'));
assert(exist('mxl_gpu_auto','file') == 2,'The production GPU dispatcher is required.');
if ~exist(outDir,'dir'), mkdir(outDir); end
loaded = load(fixtureFile,'cases');
assert(isfield(loaded,'cases'),'The fixture MAT file must contain cases.');
normal = [];
varying = [];
for fixture = loaded.cases
    if ~strcmp(fixture.Model,'MXL'), continue; end
    opt = fixture.Args{6};
    if opt.WTP_space == 0 && opt.FullCov == 0 && opt.NVarM == 0 &&...
            opt.NVarS == 0 && opt.Dist(1) == 0
        normal = asInput(fixture.Args);
        normal.EstimOpt.Dist(:) = 0;
    end
    if opt.WTP_space == 2 && opt.FullCov == 1 && opt.mCT == 1 &&...
            opt.NVarM == 1 && opt.NVarS == 1
        varying = asInput(fixture.Args);
    end
end
assert(~isempty(normal) && ~isempty(varying),'Required MXL fixtures were not found.');
if isfield(normal.EstimOpt,'GPU'), normal.EstimOpt = rmfield(normal.EstimOpt,'GPU'); end
varying.EstimOpt.GPU = 'gpu';
cases = testCases(normal,varying);
references = cell(numel(cases),1);
% Compute references before the candidate sequence, preserving resident-cache tests.
for k = 1:numel(cases)
    references{k} = reference(cases(k).Input);
end
try
    gpuAvailable = gpuDeviceCount('available') > 0;
catch
    gpuAvailable = false;
end
environment = struct('MATLAB',version,'GPUAvailable',gpuAvailable,...
    'Workers',nWorkers,'FixtureFile',fixtureFile);
rows = struct([]);
raw = struct([]);
for k = 1:numel(cases)
    for wantGradient = [false true false]
        expectedCPU = cases(k).ExpectedCPU ||...
            (wantGradient && cases(k).ExpectedCPUGradient);
        [r,s] = compareCase(cases(k).Name,cases(k).Input,wantGradient,...
            'serial',expectedCPU,references{k});
        append(r,s);
    end
end
[r,s] = compareHessian(normal);
append(r,s);
if nWorkers > 0
    pool = parpool('Processes',nWorkers);
    closePool = onCleanup(@() deletePool(pool));
    C = varying;
    future = parfeval(pool,@workerProbe,2,C);
    [r,s] = fetchOutputs(future);
    append(r,s);
    [r,s] = compareCase('client_with_pool',C,true,'pool',false,reference(C));
    append(r,s);
    clear closePool;
    [r,s] = compareCase('after_pool_deleted',C,true,'serial',false,reference(C));
    append(r,s);
end
if ~gpuAvailable
    assert(all(strcmp({rows.Backend},'cpu')),...
        'A GPU-free host must use graceful CPU fallback.');
else
    required = ismember({rows.Case},{'forced_gpu','row_varying_covariates'}) &...
        [rows.Gradient];
    assert(sum(required) == 2 && all(strcmp({rows(required).Backend},'gpu')),...
        'Forced supported gradients must use GPU on a GPU-available host.');
end
report = struct2table(rows,'AsArray',true);
writetable(report,fullfile(outDir,'gpu_auto.csv'));
save(fullfile(outDir,'gpu_auto.mat'),'report','raw','environment');
fprintf('Automatic GPU/CPU regressions passed: %d comparisons, %d GPU results.\n',...
    height(report),sum(strcmp(report.Backend,'gpu')));
clear restoreSettings;

    function append(r,s)
        j = numel(rows)+1;
        if isempty(rows), rows = r; else, rows(j) = r; end
        if isempty(raw), raw = s; else, raw(j) = s; end
        report = struct2table(rows,'AsArray',true);
        writetable(report,fullfile(outDir,'gpu_auto.csv'));
        save(fullfile(outDir,'gpu_auto.mat'),'report','raw','environment');
    end
end

function cases = testCases(normal,varying)
cases = entry('default_auto',normal,false);
C = normal; C.EstimOpt.GPU = 'auto';
cases(end+1) = entry('explicit_auto',C,false);
C = normal; C.EstimOpt.GPU = 'cpu';
cases(end+1) = entry('explicit_cpu',C,true);
C = normal; C.EstimOpt.GPU = false;
cases(end+1) = entry('logical_cpu_alias',C,true);
C = normal; C.EstimOpt.GPU = 'gpu';
cases(end+1) = entry('forced_gpu',C,false);
C = normal; C.EstimOpt.GPU = true;
cases(end+1) = entry('logical_gpu_alias',C,false);
cases(end+1) = entry('row_varying_covariates',varying,false);
C = varying;
C.EstimOpt.indx1 = flip(C.EstimOpt.indx1);
C.EstimOpt.indx2 = flip(C.EstimOpt.indx2);
cases(end+1) = entry('custom_Cholesky_order',C,false,true);
C = varying; C.b(1) = C.b(1)+.17;
cases(end+1) = entry('changed_b',C,false);
C = varying;
task = 1:C.EstimOpt.NAlt;
old = task(find(C.YY(task,1) == 1,1));
next = task(find(C.YY(task,1) == 0,1));
C.YY(old,1) = 0; C.YY(next,1) = 1;
cases(end+1) = entry('changed_Y_same_shape',C,false);
C = varying; C.XXa(1,1,1) = C.XXa(1,1,1)+.137;
cases(end+1) = entry('changed_X_same_shape',C,false);
C = varying; C.err(1,1) = C.err(1,1)+.21;
cases(end+1) = entry('changed_draws_same_shape',C,false);
C = varying; C.XXm(1,1) = C.XXm(1,1)+.11;
cases(end+1) = entry('changed_Xm_same_shape',C,false);
C = varying; C.Xs(1,1) = C.Xs(1,1)+.12;
cases(end+1) = entry('changed_Xs_same_shape',C,false);
C = varying; C.EstimOpt.WTP_matrix = C.EstimOpt.NVarA;
cases(end+1) = entry('changed_WTP_mapping',C,false);
C = varying; C.EstimOpt.Dist(1) = 1;
cases(end+1) = entry('changed_distribution',C,false);
C = varying; C.EstimOpt.RealMin = 1;
cases(end+1) = entry('changed_RealMin',C,false);
cases(end+1) = entry('return_to_original_dataset',varying,false);
C = varying; C.EstimOpt.ExpB = 1;
cases(end+1) = entry('ExpB_CPU_fallback',C,true);
C = normal; C.EstimOpt.GPU = 'gpu'; C.EstimOpt.Dist(1) = 2;
cases(end+1) = entry('spike_CPU_fallback',C,true);
C = normal; C.EstimOpt.GPU = 'gpu'; C.YY(:) = 0;
cases(end+1) = entry('legacy_zero_chosen',C,false);
C = normal; C.EstimOpt.GPU = 'gpu'; C.YY(2,1) = 1;
cases(end+1) = entry('legacy_two_chosen',C,false);
C = tinyInput(normal); C.XXa(1) = Inf;
cases(end+1) = entry('infinite_available_utility',C,false);
C = tinyInput(normal); C.XXa(1) = NaN;
cases(end+1) = entry('NaN_available_utility',C,false);
end

function C = tinyInput(C)
opt = C.EstimOpt;
opt.NVarA = 1; opt.NAlt = 2; opt.NCT = 1; opt.NP = 1; opt.NRep = 1;
opt.NVarM = 0; opt.NVarS = 0; opt.FullCov = 1; opt.Dist = 0;
opt.WTP_space = 0; opt.WTP_matrix = []; opt.mCT = 0; opt.RealMin = 0;
opt.indx1 = 1; opt.indx2 = 1; opt.DiagIndex = 1;
opt.MissingCT = false; opt.NCTMiss = 1; opt.NAltMiss = 2;
opt.NAltMissInd = 2; opt.NAltMissIndExp = [2;2]; opt.GPU = 'gpu';
C = struct('YY',[1;0],'XXa',[.7;0],'XXm',zeros(0,1),...
    'Xs',zeros(2,0),'err',.3,'EstimOpt',opt,'b',[-.2;.4]);
end

function value = entry(name,C,expectedCPU,expectedCPUGradient)
if nargin < 4, expectedCPUGradient = expectedCPU; end
value = struct('Name',name,'Input',C,'ExpectedCPU',expectedCPU,...
    'ExpectedCPUGradient',expectedCPUGradient);
end

function C = asInput(args)
C = struct('YY',args{1},'XXa',args{2},'XXm',args{3},'Xs',args{4},...
    'err',args{5},'EstimOpt',args{6},'b',args{7});
end

function result = reference(C)
[result.f,result.g] = cpuReference(C,true);
result.value = cpuReference(C,false);
end

function [f,g] = cpuReference(C,wantGradient)
C.EstimOpt.GPU = 'cpu';
args = {C.YY,C.XXa,C.XXm,C.Xs,C.err,C.EstimOpt,C.b};
if wantGradient, [f,g] = LL_mxl(args{:}); else, f = LL_mxl(args{:}); g = []; end
end

function [r,s] = compareCase(name,C,wantGradient,stage,expectedCPU,ref)
args = {C.YY,C.XXa,C.XXm,C.Xs,C.err,C.EstimOpt,C.b};
if wantGradient
    [publicF,publicG] = LL_mxl(args{:});
    referenceF = ref.f; referenceG = ref.g;
else
    publicF = LL_mxl(args{:}); publicG = [];
    referenceF = ref.value; referenceG = [];
end
% The public call calibrates with its real CPU kernel before this numerical probe.
input = rmfield(C,'b');
[f,g,backend] = mxl_gpu_auto(input,C.b,wantGradient,@() referencePair(ref,wantGradient));
assert(ismember(backend,{'cpu','gpu'}),'Unexpected backend.');
if expectedCPU, assert(strcmp(backend,'cpu'),'This case must use CPU: %s.',name); end
assert(isa(f,'double') && ~isa(f,'gpuArray') &&...
    isa(g,'double') && ~isa(g,'gpuArray'),'Public results must be CPU double.');
r = struct('Case',name,'Stage',stage,'Gradient',wantGradient,'Backend',backend,...
    'DispatcherValueRelativeError',difference(f,referenceF),...
    'DispatcherGradientRelativeError',difference(g,referenceG),...
    'PublicValueRelativeError',difference(publicF,referenceF),...
    'PublicGradientRelativeError',difference(publicG,referenceG),...
    'HessianRelativeError',NaN,'Passed',true);
assert(r.DispatcherValueRelativeError <= 1e-8 && r.PublicValueRelativeError <= 1e-8,...
    'Likelihood regression: %s.',name);
assert(r.DispatcherGradientRelativeError <= 1e-6 &&...
    r.PublicGradientRelativeError <= 1e-6,'Gradient regression: %s.',name);
s = struct('Case',name,'Stage',stage,'Gradient',wantGradient,'Backend',backend,...
    'Input',C,'Referencef',referenceF,'Referenceg',referenceG,...
    'Dispatcherf',f,'Dispatcherg',g,'Publicf',publicF,'Publicg',publicG,...
    'ReferenceHessian',[],'PublicHessian',[]);
end

function [f,g] = referencePair(ref,wantGradient)
if wantGradient, f = ref.f; g = ref.g; else, f = ref.value; g = []; end
end

function [r,s] = compareHessian(C)
opt = C.EstimOpt; opt.GPU = 'cpu';
args = {C.YY,C.XXa,C.XXm,C.Xs,C.err,opt,C.b};
[refF,refG,refH] = LL_mxl(args{:});
opt = rmfield(opt,'GPU'); args{6} = opt;
[f,g,h] = LL_mxl(args{:});
r = struct('Case','Hessian_CPU_bypass','Stage','serial','Gradient',true,'Backend','cpu',...
    'DispatcherValueRelativeError',0,'DispatcherGradientRelativeError',0,...
    'PublicValueRelativeError',difference(f,refF),...
    'PublicGradientRelativeError',difference(g,refG),...
    'HessianRelativeError',difference(h,refH),'Passed',true);
assert(r.PublicValueRelativeError <= 1e-8 && r.PublicGradientRelativeError <= 1e-6 &&...
    r.HessianRelativeError <= 1e-6,'Hessian CPU bypass regression.');
s = struct('Case',r.Case,'Stage',r.Stage,'Gradient',true,'Backend','cpu',...
    'Input',C,'Referencef',refF,'Referenceg',refG,'Dispatcherf',refF,...
    'Dispatcherg',refG,'Publicf',f,'Publicg',g,...
    'ReferenceHessian',refH,'PublicHessian',h);
end

function [r,s] = workerProbe(C)
assert(~isempty(getCurrentTask()),'Worker probe must run on a worker.');
[r,s] = compareCase('GPU_disabled_on_worker',C,true,'worker',true,reference(C));
end

function value = difference(actual,expected)
if isempty(actual) && isempty(expected), value = 0; return; end
assert(isequal(size(actual),size(expected)),'Output dimensions differ.');
nonfinite = ~isfinite(actual) | ~isfinite(expected);
assert(isequaln(actual(nonfinite),expected(nonfinite)),...
    'Inf/NaN output semantics changed.');
actual = actual(~nonfinite); expected = expected(~nonfinite);
if isempty(actual), value = 0; return; end
value = max(abs(actual-expected),[],'all')/max(1,max(abs(expected),[],'all'));
end

function deletePool(pool)
if isvalid(pool), delete(pool); end
end

function setAutoCreate(settings,value)
settings.Pool.AutoCreate = value;
end
