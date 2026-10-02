function [report,raw] = test_mxl_covariates(baselineDir,outDir,workers)
% Small baseline/new checks for Xm, Xs and subnormal simulated probabilities.
% test_mxl_covariates(baselineDir,outDir,[0 3]); no external data or optimizer.
repo = fileparts(fileparts(mfilename('fullpath')));
if nargin < 3 || isempty(workers), workers = 0; end
outDir = char(java.io.File(outDir).getCanonicalPath());
protectedDir = fullfile(fileparts(fileparts(fileparts(repo))),'lasy','replication_package');
assert(~strcmpi(outDir,protectedDir) &&...
    ~startsWith(lower(outDir),[lower(protectedDir),filesep]),...
    'Test output must be outside replication_package.');
assert(exist(fullfile(baselineDir,'LL_mxl_baseline.m'),'file') == 2,...
    'Provide the dd6f704 baseline exported by test_mxl_memory.');
assert(isempty(gcp('nocreate')),'Run in a dedicated MATLAB session without a pool.');
if ~exist(outDir,'dir'), mkdir(outDir); end
addpath(baselineDir,fullfile(repo,'MXL'));
settings = parallel.Settings;
autoCreate = settings.Pool.AutoCreate;
settings.Pool.AutoCreate = false;
restore = onCleanup(@() restoreAutoCreate(settings,autoCreate));
rows = struct([]);
raw = struct([]);
for nWorkers = workers(:)'
    if nWorkers > 0
        pool = parpool('Processes',nWorkers);
        closePool = onCleanup(@() delete(pool));
    end
    for wtp = 0:1
        for fullCov = 0:1
            for meanMode = 0:3
                for scaleMode = 0:2
                    C = covariateCase(wtp,fullCov,meanMode,scaleMode);
                    [r,result] = compareCase(C,nWorkers);
                    if isempty(rows), rows = r; raw = result;
                    else, rows(end+1) = r; raw(end+1) = result; end
                    saveResults(outDir,rows,raw);
                    assert(r.Passed,'Covariate regression: %s',C.Name);
                end
            end
        end
    end
    % The pre-existing special path for CT-specific means with all-normal WTP.
    C = covariateCase(1,0,2,0);
    C.Name = 'wtp_diag_Xm_CT_all_normal';
    C.EstimOpt.Dist(:) = 0;
    [r,result] = compareCase(C,nWorkers);
    rows(end+1) = r; raw(end+1) = result;
    saveResults(outDir,rows,raw);
    assert(r.Passed,'All-normal mCT WTP regression.');
    for fullCov = 0:1
        C = subnormalCase(fullCov);
        [r,result] = compareCase(C,nWorkers);
        rows(end+1) = r; raw(end+1) = result;
        saveResults(outDir,rows,raw);
        assert(r.Passed,'Subnormal probability regression: %s',C.Name);
        expected = [1 0 2 zeros(1,1+fullCov)];
        assert(isempty(result.new.ErrorID) &&...
            max(abs(result.new.g-expected)) <= 1e-12,...
            'Subnormal gradient must preserve the baseline rounding order.');
    end
    if nWorkers > 0, clear closePool; end
end
report = struct2table(rows,'AsArray',true);
disp(report);
end

function C = covariateCase(wtp,fullCov,meanMode,scaleMode)
E = options(3,2,3,8,3,fullCov);
E.Dist = [0 0 1];
E.WTP_space = wtp;
E.WTP_matrix = [];
if wtp, E.WTP_matrix = [3 3]; end
E.NVarM = double(meanMode > 0);
E.NVarS = double(scaleMode > 0);
E.mCT = double(meanMode > 1);
C.YY = repmat([1;0;0;0;1;0],1,3);
C.XXa = reshape(sin(1:6*3*3),6,3,3);
C.XXa(:,3,:) = -abs(C.XXa(:,3,:))-.3;
C.err = reshape(linspace(-1.1,1.2,3*8*3),3,[]);
switch meanMode
    case 0, C.XXm = zeros(0,3);
    case 1, C.XXm = [-.5,0,.5];
    case 2, C.XXm = repelem([-.5,.2,.1,.6,-.2,.3],3);
    case 3, C.XXm = reshape(linspace(-.5,.5,18),1,[]);
end
switch scaleMode
    case 0, C.Xs = zeros(18,0);
    case 1, C.Xs = repelem([-.3;0;.3],6);
    case 2, C.Xs = repelem([-.4;.1;0;.3;-.2;.2],3);
end
L = [.3 0 0;.04 .2 0;-.03 .05 .15];
if fullCov, covariance = L(tril(true(3))); else, covariance = diag(L); end
C.b = [.2;-.1;-.3;covariance];
if E.NVarM, C.b = [C.b;.1;-.2;.05]; end
if E.NVarS, C.b = [C.b;.08]; end
C.EstimOpt = E;
C.Name = sprintf('WTP%d_FullCov%d_Xm%d_Xs%d',wtp,fullCov,meanMode,scaleMode);
end

function C = subnormalCase(fullCov)
E = options(2,1,1,1,2,fullCov);
C.YY = [1;0];
C.XXa = [-.7,0;0,0];
C.XXm = zeros(0,1);
C.Xs = zeros(2,0);
C.err = [3;0];
C.b = [-log(realmin*eps)/.7;0;zeros(2+fullCov,1)];
C.EstimOpt = E;
C.Name = sprintf('subnormal_FullCov%d',fullCov);
end

function E = options(nAlt,nCT,nP,nRep,k,fullCov)
E = struct('NAlt',nAlt,'NCT',nCT,'NP',nP,'NRep',nRep,'NVarA',k,...
    'NVarM',0,'NVarS',0,'Dist',zeros(1,k),'WTP_space',0,'WTP_matrix',[],...
    'FullCov',fullCov,'Triang',[],'NVarNLT',0,'NLTVariables',[],'NLTType',[],...
    'Johnson',0,'NCTMiss',nCT*ones(nP,1),'NAltMiss',nAlt*ones(nP,1),...
    'NAltMissInd',nAlt*ones(nCT,nP),'NAltMissIndExp',nAlt*ones(nAlt*nCT,nP),...
    'MissingCT',false(nCT,nP),'RealMin',0,'ExpB',[],'mCT',0);
index = tril(ones(k));
index(index == 1) = 1:sum(1:k);
E.DiagIndex = diag(index);
[E.indx1,E.indx2] = find(tril(true(k)));
E.indx1 = E.indx1';
E.indx2 = E.indx2';
end

function [r,result] = compareCase(C,nWorkers)
result = struct('baseline',evaluate(@LL_mxl_baseline,C),'new',evaluate(@LL_mxl,C),...
    'numericalGradient',[]);
old = result.baseline;
new = result.new;
r = struct('Case',C.Name,'Workers',nWorkers,'BaselineError',old.ErrorID,...
    'NewError',new.ErrorID,'BaselineValueError',old.ValueErrorID,...
    'NewValueError',new.ValueErrorID,'ValueRelativeError',NaN,'GradientRelativeError',NaN,...
    'FiniteDifferenceError',NaN,...
    'Passed',false);
if ~isempty(old.ValueErrorID) || ~isempty(new.ValueErrorID)
    valuePassed = ~isempty(old.ValueErrorID) &&...
        strcmp(old.ValueErrorID,new.ValueErrorID) &&...
        strcmp(old.ValueErrorMessage,new.ValueErrorMessage);
else
    r.ValueRelativeError = max(abs(new.v-old.v))/max(1,max(abs(old.v)));
    valuePassed = all(isfinite([old.v(:);new.v(:)])) && r.ValueRelativeError <= 1e-8;
end
if ~isempty(old.ErrorID) || ~isempty(new.ErrorID)
    if ~isempty(old.ErrorID) && isempty(new.ErrorID)
        result.numericalGradient = finiteDifference(C);
        r.FiniteDifferenceError = max(abs(new.g-result.numericalGradient),[],'all')/...
            max(1,max(abs(result.numericalGradient),[],'all'));
        r.Passed = valuePassed && all(isfinite([new.f(:);new.g(:)])) &&...
            max(abs(new.f-new.v)) <= 1e-8*max(1,max(abs(new.f))) &&...
            r.FiniteDifferenceError <= 5e-6;
    end
else
    r.GradientRelativeError = max(abs(new.g-old.g),[],'all')/max(1,max(abs(old.g),[],'all'));
    r.Passed = valuePassed && all(isfinite([old.f(:);new.f(:);old.g(:);new.g(:)])) &&...
        max(abs(new.f-old.f)) <= 1e-8*max(1,max(abs(old.f))) &&...
        r.GradientRelativeError <= 1e-6;
    if isempty(new.ValueErrorID) && isempty(old.ValueErrorID)
        r.Passed = r.Passed && max(abs(new.f-new.v)) <= 1e-8*max(1,max(abs(new.f)));
    end
end
end

function J = finiteDifference(C)
args = {C.YY,C.XXa,C.XXm,C.Xs,C.err,C.EstimOpt,C.b};
J = zeros(C.EstimOpt.NP,numel(C.b));
for j = 1:numel(C.b)
    h = 1e-5*max(1,abs(C.b(j)));
    plus = args;
    minus = args;
    plus{end}(j) = C.b(j)+h;
    minus{end}(j) = C.b(j)-h;
    J(:,j) = (LL_mxl(plus{:})-LL_mxl(minus{:}))/(2*h);
end
end

function result = evaluate(fun,C)
result = struct('f',[],'g',[],'v',[],'ErrorID','','ErrorMessage',...
    '','ValueErrorID','','ValueErrorMessage','');
args = {C.YY,C.XXa,C.XXm,C.Xs,C.err,C.EstimOpt,C.b};
try
    [result.f,result.g] = fun(args{:});
catch ME
    result.ErrorID = ME.identifier;
    if isempty(result.ErrorID), result.ErrorID = 'error_without_identifier'; end
    result.ErrorMessage = ME.message;
end
try
    result.v = fun(args{:});
catch ME
    result.ValueErrorID = ME.identifier;
    if isempty(result.ValueErrorID), result.ValueErrorID = 'error_without_identifier'; end
    result.ValueErrorMessage = ME.message;
end
end

function saveResults(outDir,rows,raw)
report = struct2table(rows,'AsArray',true);
save(fullfile(outDir,'covariates.mat'),'report','raw');
writetable(report,fullfile(outDir,'covariates.csv'));
end

function restoreAutoCreate(settings,value)
settings.Pool.AutoCreate = value;
end
