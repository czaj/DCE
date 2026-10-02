function [report,raw] = test_mxl_memory(outDir,workers,mode,replicationRoot,caseNames,baselineDir)
% Compare dd6f704 with working-tree LL_mxl on identical inputs and draws.
% test_mxl_memory(outDir,[0 3],'compare') checks six cases, no optimizer.
% test_mxl_memory(outDir,3,'estimate') re-estimates CH from published bhat.
% Optional replicationRoot and baselineDir keep all output outside the paper.
global B_backup
repo = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(outDir), outDir = fullfile(tempdir,'DCE_mxl_memory'); end
if nargin < 2 || isempty(workers), workers = [0 3]; end
if nargin < 3 || isempty(mode), mode = 'compare'; end
if nargin < 4, replicationRoot = []; end
if nargin < 5 || isempty(caseNames)
    caseNames = {'demo_pref_diag','demo_pref_full','demo_wtp_diag',...
        'demo_wtp_full','CH','pooled'};
end
if nargin < 6 || isempty(baselineDir), baselineDir = fullfile(outDir,'baseline'); end
if ischar(caseNames) || isstring(caseNames), caseNames = cellstr(caseNames); end
if isempty(replicationRoot)
    replicationRoot = fullfile(fileparts(fileparts(fileparts(repo))),...
        'lasy','replication_package');
end
outDir = char(java.io.File(outDir).getCanonicalPath());
baselineDir = char(java.io.File(baselineDir).getCanonicalPath());
protectedDir = char(java.io.File(replicationRoot).getCanonicalPath());
assert(~strcmpi(outDir,protectedDir) &&...
    ~startsWith(lower(outDir),[lower(protectedDir),filesep]),...
    'Test output must be outside replication_package.');
assert(~strcmpi(baselineDir,protectedDir) &&...
    ~startsWith(lower(baselineDir),[lower(protectedDir),filesep]),...
    'Baseline export must be outside replication_package.');
if ~exist(outDir,'dir'), mkdir(outDir); end
addpath(fullfile(repo,'tests'),fullfile(repo,'MXL'));
assert(ismember(mode,{'compare','estimate'}),'Unknown test mode: %s',mode);
assert(isempty(gcp('nocreate')),...
    'Run in a dedicated MATLAB session without an existing pool.');
settings = parallel.Settings;
autoCreate = settings.Pool.AutoCreate;
settings.Pool.AutoCreate = false;
restore = onCleanup(@() setAutoCreate(settings,autoCreate));

if strcmp(mode,'estimate')
    assert(isscalar(workers),'Select one worker count for CH estimation.');
    if workers > 0
        pool = parpool('Processes',workers);
        closePool = onCleanup(@() deletePool(pool));
    end
    C = mxl_memory_case('CH',replicationRoot);
    oldDir = cd(outDir);
    restoreDir = onCleanup(@() cd(oldDir));
    E = C.EstimOpt;
    E.ProjectName = 'mxl_memory_CH';
    E.Display = 0;
    OptimOpt = C.OptimOpt;
    oldBackup = B_backup;
    restoreBackup = onCleanup(@() setBackup(oldBackup));
    B_backup = C.b;
    C.Results.MXL.b0 = C.b;
    OptimOpt.Algorithm = 'trust-region';
    OptimOpt.Hessian = 'user-supplied';
    timer = tic;
    trustRegion = MXL(C.INPUT,C.Results,E,OptimOpt);
    C.Results.MXL.b0 = trustRegion.bhat;
    OptimOpt.Algorithm = 'quasi-newton';
    OptimOpt.Hessian = 'off';
    estimate = MXL(C.INPUT,C.Results,E,OptimOpt);
    elapsed = toc(timer);
    LLerror = abs(estimate.LL-C.Published.LL)/max(1,abs(C.Published.LL));
    parameterError = max(abs(estimate.bhat(:)-C.b))/max(1,max(abs(C.b)));
    report = table(workers,elapsed,C.Published.LL,estimate.LL,LLerror,parameterError,...
        'VariableNames',{'Workers','Seconds','PublishedLL','EstimatedLL',...
        'LLRelativeError','ParameterRelativeError'});
    raw = struct('trustRegion',trustRegion,'estimate',estimate,'published',C.Published);
    save(fullfile(outDir,'CH_estimation.mat'),'report','raw','-v7.3');
    writetable(report,fullfile(outDir,'CH_estimation.csv'));
    disp(report);
    assert(LLerror <= 1e-8,'CH estimate differs from published LL.');
    assert(parameterError <= 1e-6,'CH estimate differs from published parameters.');
    return
end

baselineFile = fullfile(baselineDir,'LL_mxl_baseline.m');
if ~exist(baselineFile,'file')
    [status,source] = system(sprintf('git -C "%s" show dd6f704:MXL/LL_mxl.m',repo));
    assert(status == 0,'Cannot extract dd6f704 baseline.');
    source = regexprep(source,'^(function[^\r\n]*=\s*)LL_mxl\(',...
        '$1LL_mxl_baseline(','once');
    if ~exist(baselineDir,'dir'), mkdir(baselineDir); end
    fid = fopen(baselineFile,'w');
    assert(fid >= 0,'Cannot write baseline: %s',baselineFile);
    closeBaseline = onCleanup(@() fclose(fid));
    fwrite(fid,source,'char');
    clear closeBaseline
end
addpath(baselineDir);
rows = struct([]);
raw = struct([]);
for nWorkers = workers(:)'
    if nWorkers > 0
        pool = parpool('Processes',nWorkers);
        closePool = onCleanup(@() deletePool(pool));
    end
    for k = 1:numel(caseNames)
        C = mxl_memory_case(caseNames{k},replicationRoot);
        args = {C.YY,C.XXa,C.XXm,C.Xs,C.err,C.EstimOpt,C.b};
        [~,~] = LL_mxl_baseline(args{:});
        [~,~] = LL_mxl(args{:});
        % Warm-up includes JIT, worker setup and any immutable-data transfer.
        timer = tic;
        [f0,g0] = LL_mxl_baseline(args{:});
        baselineSeconds = toc(timer);
        timer = tic;
        [f1,g1] = LL_mxl(args{:});
        newSeconds = toc(timer);
        timer = tic;
        v0 = LL_mxl_baseline(args{:});
        baselineValueSeconds = toc(timer);
        timer = tic;
        v1 = LL_mxl(args{:});
        newValueSeconds = toc(timer);
        grad0 = sum(C.W.*g0,1);
        grad1 = sum(C.W.*g1,1);
        LL0 = -sum(C.W.*f0);
        LL1 = -sum(C.W.*f1);
        r = struct('Case',C.Name,'Workers',nWorkers,'NP',C.EstimOpt.NP,...
            'BaselineLL',LL0,'NewLL',LL1,...
            'LLRelativeError',abs(LL1-LL0)/max(1,abs(LL0)),...
            'GradientRelativeError',max(abs(grad1-grad0))/max(1,max(abs(grad0))),...
            'MaxPersonValueError',max(abs(f1-f0)),...
            'MaxPersonGradientError',max(abs(g1-g0),[],'all'),...
            'PersonGradientRelativeError',max(abs(g1-g0),[],'all')/max(1,max(abs(g0),[],'all')),...
            'ValueOnlyRelativeError',max(abs(v1-v0))/max(1,max(abs(v0))),...
            'BaselineGradientSeconds',baselineSeconds,'NewGradientSeconds',newSeconds,...
            'BaselineValueSeconds',baselineValueSeconds,'NewValueSeconds',newValueSeconds);
        j = numel(rows)+1;
        if isempty(rows), rows = r; else, rows(j) = r; end
        raw(j).Case = C.Name;
        raw(j).Workers = nWorkers;
        raw(j).b = C.b;
        raw(j).W = C.W;
        raw(j).EstimOpt = C.EstimOpt;
        raw(j).fBaseline = f0;
        raw(j).gBaseline = g0;
        raw(j).fNew = f1;
        raw(j).gNew = g1;
        raw(j).valueBaseline = v0;
        raw(j).valueNew = v1;
        report = struct2table(rows);
        save(fullfile(outDir,'correctness.mat'),'report','raw','-v7.3');
        writetable(report,fullfile(outDir,'correctness.csv'));
        disp(report(end,:));
        assert(all(isfinite([f0(:);f1(:);g0(:);g1(:);v0(:);v1(:)])),...
            'Non-finite likelihood or gradient: %s',C.Name);
        assert(r.LLRelativeError <= 1e-8 && r.ValueOnlyRelativeError <= 1e-8,...
            'Likelihood mismatch: %s (%d workers)',C.Name,nWorkers);
        assert(r.GradientRelativeError <= 1e-6 && r.PersonGradientRelativeError <= 1e-6,...
            'Gradient mismatch: %s (%d workers)',C.Name,nWorkers);
        assert(max(abs(f0-v0)) <= 1e-8*max(1,max(abs(f0))) &&...
            max(abs(f1-v1)) <= 1e-8*max(1,max(abs(f1))),...
            'Value-only and gradient branches differ: %s',C.Name);
        if strcmp(C.Name,'CH')
            assert(abs(LL0-(-4793.8164)) <= 5e-5,...
                'CH likelihood does not match published -4793.8164.');
        end
    end
    testSupplemental(outDir,nWorkers);
    if nWorkers > 0, clear closePool; end
end
end

function testSupplemental(outDir,nWorkers)
% Small cases exercise derivative contraction and immutable-data cache refresh.
E = struct('NAlt',3,'NCT',2,'NP',3,'NRep',8,'NVarA',4,'NVarM',1,...
    'NVarS',0,'Dist',[0 0 1 1],'WTP_space',2,'WTP_matrix',[3 4],...
    'FullCov',1,'Triang',[],'NVarNLT',0,'NLTVariables',[],'NLTType',[],...
    'Johnson',0,'NCTMiss',2*ones(3,1),'NAltMiss',3*ones(3,1),...
    'NAltMissInd',3*ones(2,3),'NAltMissIndExp',3*ones(6,3),...
    'MissingCT',false(2,3),'RealMin',0,'ExpB',[],'mCT',0,...
    'DiagIndex',[1;5;8;10],'indx1',[1 2 3 4 2 3 4 3 4 4],...
    'indx2',[1 1 1 1 2 2 2 3 3 4]);
YY = repmat([1;0;0;0;1;0],1,3);
X = reshape(sin(1:6*4*3),6,4,3);
X(:,3:4,:) = -abs(X(:,3:4,:))-.2;
Xm = [-.5,0,.5];
err = reshape(linspace(-1.2,1.1,4*8*3),4,[]);
L = [.3,0,0,0;.1,.2,0,0;-.05,.03,.15,0;.02,-.04,.05,.2];
b = [.2;-.1;-.2;.1;L(tril(true(4)));.1;-.2;.05;-.08];
Xs = zeros(6*3,0);
names = {'two_costs_means','same_shape_Xa','same_shape_err',...
    'same_shape_YY','same_shape_Xm','diagonal','fixed','realmin_underflow'};
rows = struct([]);
raw = struct([]);
for k = 1:numel(names)
    opt = E; y = YY; x = X; xm = Xm; draws = err; beta = b;
    switch names{k}
        case 'same_shape_Xa'
            x(1,1,1) = x(1,1,1)+.3;
        case 'same_shape_err'
            draws(2,3) = draws(2,3)+.4;
        case 'same_shape_YY'
            y(1:3,1) = [0;1;0];
        case 'same_shape_Xm'
            xm(1) = xm(1)+.25;
        case 'diagonal'
            opt.FullCov = 0;
            beta = [b(1:4);diag(L);b(15:end)];
        case 'fixed'
            opt.Dist(1) = -1;
            fixedL = L;
            fixedL(1,:) = 0;
            fixedL(:,1) = 0;
            beta(5:14) = fixedL(tril(true(4)));
            draws(1,:) = 0;
        case 'realmin_underflow'
            opt.RealMin = 1;
            beta(3:4) = 10;
            x(:) = 0;
            x(:,3:4,:) = -1;
            for n = 1:3
                x(y(:,n) == 1,3:4,n) = -100;
            end
    end
    args = {y,x,xm,Xs,draws,opt,beta};
    [f0,g0] = LL_mxl_baseline(args{:});
    [f1,g1] = LL_mxl(args{:});
    v0 = LL_mxl_baseline(args{:});
    v1 = LL_mxl(args{:});
    rows(k).Case = names{k};
    rows(k).Workers = nWorkers;
    rows(k).ValueRelativeError = max(abs([f1-f0;v1-v0]))/max(1,max(abs([f0;v0])));
    rows(k).GradientRelativeError = max(abs(g1-g0),[],'all')/max(1,max(abs(g0),[],'all'));
    result = struct('fBaseline',f0,'gBaseline',g0,'fNew',f1,'gNew',g1,...
        'valueBaseline',v0,'valueNew',v1);
    if isempty(raw), raw = result; else, raw(k) = result; end
    report = struct2table(rows);
    file = sprintf('supplemental_%d',nWorkers);
    save(fullfile(outDir,[file,'.mat']),'report','raw');
    writetable(report,fullfile(outDir,[file,'.csv']));
    assert(all(isfinite([f0(:);f1(:);g0(:);g1(:);v0(:);v1(:)])),...
        'Non-finite supplemental result: %s',names{k});
    assert(rows(k).ValueRelativeError <= 1e-8 &&...
        rows(k).GradientRelativeError <= 1e-6,...
        'Supplemental baseline mismatch: %s',names{k});
end
disp(report);
end

function deletePool(pool)
if isvalid(pool), delete(pool); end
end

function setAutoCreate(settings,value)
settings.Pool.AutoCreate = value;
end

function setBackup(value)
global B_backup
B_backup = value;
end
