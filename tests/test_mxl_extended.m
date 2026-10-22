function [report,raw] = test_mxl_extended(outDir,workers,baselineRef,baselineDir)
% Small deterministic regressions for MXL, HMXL and latent-class MXL.
if nargin < 1 || isempty(outDir), outDir = fullfile(tempdir,'DCE_mxl_extended'); end
if nargin < 2 || isempty(workers), workers = [0 3]; end
if nargin < 3 || isempty(baselineRef), baselineRef = '6bf56b8'; end
repo = fileparts(fileparts(mfilename('fullpath')));
outDir = char(java.io.File(outDir).getCanonicalPath());
frozen = fullfile(fileparts(fileparts(fileparts(repo))),'lasy','replication_package');
assert(~strcmpi(outDir,frozen) && ~startsWith(lower(outDir),[lower(frozen),filesep]),...
    'Test output must be outside replication_package.');
assert(isempty(gcp('nocreate')),'Use a dedicated session with no existing pool.');
settings = parallel.Settings;
oldAutoCreate = settings.Pool.AutoCreate;
settings.Pool.AutoCreate = false;
restoreSettings = onCleanup(@() setAutoCreate(settings,oldAutoCreate)); %#ok<NASGU>
addpath(fullfile(repo,'MXL'),fullfile(repo,'HMXL'),fullfile(repo,'LCMXL'));
if ~exist(outDir,'dir'), mkdir(outDir); end
if nargin < 4 || isempty(baselineDir), baselineDir = fullfile(outDir,'baseline'); end
baselineDir = char(java.io.File(baselineDir).getCanonicalPath());
assert(~strcmpi(baselineDir,frozen) &&...
    ~startsWith(lower(baselineDir),[lower(frozen),filesep]),...
    'Baseline export must be outside replication_package.');
if ~exist(baselineDir,'dir'), mkdir(baselineDir); end
models = {'MXL','HMXL','LCMXL'};
entries = {'LL_mxl','LL_hmxl','LL_lcmxl'};
for k = 1:numel(models)
    exportBaseline(repo,baselineRef,models{k},entries{k},baselineDir);
end
addpath(baselineDir);

cases = struct([]);
for model = models
    for wtp = 0:2
        for fullCov = 0:1
            for config = 0:3
                meanMode = config;
                if strcmp(model{1},'LCMXL'), meanMode = 0; end
                missing = max(0,config-1);
                C = hybrid_fixture(model{1},wtp,fullCov,meanMode,config,missing,0);
                if isempty(cases), cases = C; else, cases(end+1) = C; end %#ok<AGROW>
            end
            if strcmp(model{1},'HMXL')
                cases(end+1) = hybrid_fixture(model{1},wtp,fullCov,3,3,2,1); %#ok<AGROW>
            end
        end
    end
    C = hybrid_fixture(model{1},1,1,0,0,0,0);
    C.Name = [C.Name,'_empty_Xs'];
    C.Args{4} = [];
    cases(end+1) = C; %#ok<AGROW>
    if ~strcmp(model{1},'LCMXL')
        C = hybrid_fixture(model{1},1,1,0,0,0,0);
        C.Name = [C.Name,'_empty_Xm'];
        C.Args{3} = [];
        cases(end+1) = C; %#ok<AGROW>
    end
    C = hybrid_fixture(model{1},0,0,0,0,0,0);
    C.Name = [C.Name,'_fixed'];
    if strcmp(model{1},'LCMXL')
        C.Args{6}.Dist(1,:) = -1;
        C.Args{5}([1 4],:) = 0;
        C.Args{7}([7 10]) = 0;
    else
        optionIndex = numel(C.Args)-1;
        drawIndex = optionIndex-1;
        C.Args{optionIndex}.Dist(1) = -1;
        C.Args{drawIndex}(1,:) = 0;
        C.Args{end}(4) = 0;
    end
    cases(end+1) = C; %#ok<AGROW>
end
C = hybrid_fixture('HMXL',1,1,3,3,2,0);
C.Name = [C.Name,'_flat_Xm'];
C.Args{3} = reshape(C.Args{3},9,[]);
cases(end+1) = C;
C = hybrid_fixture('LCMXL',0,0,0,0,0,0);
C.Name = [C.Name,'_zero_class_probability'];
C.Args{1}(:) = 0;
C.Args{1}(1:3:end,:) = 1;
C.Args{2}(:,1,:) = 0;
C.Args{2}(1:3:end,1,:) = -1;
C.Args{end}(1) = 1000;
cases(end+1) = C;
C = hybrid_fixture('HMXL',1,1,3,3,2,1,[2 2 4 7]);
C.Name = [C.Name,'_two_latent_variables'];
cases(end+1) = C;
C = hybrid_fixture('LCMXL',1,1,0,3,2,0,[1 3 4 7]);
C.Name = [C.Name,'_three_classes'];
cases(end+1) = C;
for model = {'MXL','LCMXL'}
    C = hybrid_fixture(model{1},1,1,0,1,0,0,[1 2 1 1]);
    C.Name = [C.Name,'_NP1_R1'];
    cases(end+1) = C; %#ok<AGROW>
end
save(fullfile(outDir,'fixtures.mat'),'cases','baselineRef');

% Finite differences run locally, even when the later comparisons use a pool.
serial = cell(numel(cases),1);
numerical = cell(numel(cases),1);
baseline = cell(numel(cases),1);
baselineFD = cell(numel(cases),1);
for k = 1:numel(cases)
    C = cases(k);
    serial{k} = evaluate(C.Entry,C.Args);
    assertGood(serial{k},C.Name);
    numerical{k} = finiteDifference(C.Entry,C.Args);
    FDerror = relative(serial{k}.g,numerical{k});
    fprintf('Finite differences %d/%d: %s, relative error %.3g\n',...
        k,numel(cases),C.Name,FDerror);
    assert(FDerror <= 5e-6,...
        'Finite-difference gradient mismatch: %s',C.Name);
    baseline{k} = evaluate(['baseline_',C.Entry],C.Args);
    if isempty(baseline{k}.ValueError)
        baselineFD{k} = finiteDifference(['baseline_',C.Entry],C.Args);
    end
end

rows = struct([]);
raw = struct([]);
for nWorkers = workers(:)'
    if nWorkers > 0
        pool = parpool('Processes',nWorkers);
        closePool = onCleanup(@() deletePool(pool));
    end
    for k = 1:numel(cases)
        C = cases(k);
        if nWorkers == 0, current = serial{k}; else, current = evaluate(C.Entry,C.Args); end
        assertGood(current,C.Name);
        old = baseline{k};
        rf = NaN;
        rg = NaN;
        oldFDerror = NaN;
        if isempty(old.ValueError)
            rf = relative(current.f,old.value);
            assert(rf <= 1e-8,'Baseline value mismatch: %s',C.Name);
        end
        if isempty(old.GradientError)
            rg = relative(current.g,old.g);
            if ~isempty(baselineFD{k})
                oldFDerror = relative(old.g,baselineFD{k});
                if oldFDerror <= 5e-6
                    assert(rg <= 1e-6,'Baseline gradient regression: %s',C.Name);
                end
            end
        end
        r = struct('Case',C.Name,'Model',C.Model,'Workers',nWorkers,...
            'BaselineValueError',old.ValueError,'BaselineGradientError',old.GradientError,...
            'BaselineValueRelativeError',rf,'BaselineGradientRelativeError',rg,...
            'BaselineFiniteDifferenceError',oldFDerror,...
            'FiniteDifferenceError',relative(current.g,numerical{k}),...
            'SerialValueRelativeError',relative(current.f,serial{k}.f),...
            'SerialGradientRelativeError',relative(current.g,serial{k}.g),...
            'ValueOnlyRelativeError',relative(current.f,current.value),'Passed',true);
        assert(r.SerialValueRelativeError <= 1e-8 && r.SerialGradientRelativeError <= 1e-6,...
            'Parallel/serial mismatch: %s',C.Name);
        j = numel(rows)+1;
        if isempty(rows), rows = r; else, rows(j) = r; end
        result = struct('Case',C.Name,'Workers',nWorkers,'Args',{C.Args},...
            'baseline',old,'current',current,'numericalGradient',numerical{k});
        if isempty(raw), raw = result; else, raw(j) = result; end
        report = struct2table(rows,'AsArray',true);
        writetable(report,fullfile(outDir,'extended.csv'));
        save(fullfile(outDir,'extended.mat'),'report','raw');
    end
    if nWorkers > 0, clear closePool; end
end
fprintf('Extended MXL/HMXL/LCMXL regressions passed: %d cases, %d comparisons.\n',...
    numel(cases),height(report));
end


function r = evaluate(entry,args)
r = struct('f',[],'g',[],'value',[],'ValueError','','GradientError','');
opt = args{end-1};
valueSize = [opt.NP,1];
gradientSize = [opt.NP,numel(args{end})];
try
    r.value = feval(entry,args{:});
    if ~isequal(size(r.value),valueSize)
        r.ValueError = shapeError('value',r.value,valueSize);
    elseif any(~isfinite(r.value),'all')
        r.ValueError = 'nonfinite_output';
    end
catch exception
    r.ValueError = [exception.identifier,': ',exception.message];
end
try
    [r.f,r.g] = feval(entry,args{:});
    if ~isequal(size(r.f),valueSize)
        r.GradientError = shapeError('gradient_branch_value',r.f,valueSize);
    elseif ~isequal(size(r.g),gradientSize)
        r.GradientError = shapeError('gradient',r.g,gradientSize);
    elseif any(~isfinite([r.f(:);r.g(:)]))
        r.GradientError = 'nonfinite_output';
    end
catch exception
    r.GradientError = [exception.identifier,': ',exception.message];
end
end

function message = shapeError(output,value,expected)
message = sprintf('unsupported_output_shape: %s is %s; expected %s',...
    output,mat2str(size(value)),mat2str(expected));
end

function assertGood(r,name)
assert(isempty(r.ValueError) && isempty(r.GradientError),...
    'Candidate raised an exception: %s\n%s\n%s',name,r.ValueError,r.GradientError);
assert(all(isfinite([r.f(:);r.g(:);r.value(:)])),...
    'Non-finite candidate output: %s',name);
assert(relative(r.f,r.value) <= 1e-8,'Value/gradient branch mismatch: %s',name);
end

function J = finiteDifference(entry,args)
b = args{end};
f = feval(entry,args{:});
J = zeros(numel(f),numel(b));
for j = 1:numel(b)
    h = 1e-5*max(1,abs(b(j)));
    plus = args;
    minus = args;
    plus{end}(j) = b(j)+h;
    minus{end}(j) = b(j)-h;
    J(:,j) = (feval(entry,plus{:})-feval(entry,minus{:}))/(2*h);
end
end

function value = relative(actual,reference)
if isempty(actual) && isempty(reference), value = 0; return; end
value = max(abs(actual-reference),[],'all')/max(1,max(abs(reference),[],'all'));
end

function exportBaseline(repo,ref,folder,entry,destination)
path = fullfile(destination,['baseline_',entry,'.m']);
if exist(path,'file'), return; end
[status,source] = system(sprintf('git -C "%s" show "%s:%s/%s.m"',...
    repo,ref,folder,entry));
assert(status == 0,'Cannot export baseline %s at %s.',entry,ref);
source = regexprep(source,['^(function[^\r\n]*=\s*)',entry,'\('],...
    ['$1baseline_',entry,'('],'once');
fid = fopen(path,'w');
assert(fid >= 0,'Cannot write baseline alias.');
closeFile = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid,source,'char');
end

function deletePool(pool)
if isvalid(pool), delete(pool); end
end

function setAutoCreate(settings,value)
settings.Pool.AutoCreate = value;
end
