function report = test_hybrid_memory(outDir,workers)
if nargin < 2, workers = 0; end
if ~exist(outDir,'dir'), mkdir(outDir); end
repo = fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(repo));
settings = parallel.Settings; previousAuto = settings.Pool.AutoCreate;
restoreAuto = onCleanup(@() setAutoCreate(settings,previousAuto)); %#ok<NASGU>
settings.Pool.AutoCreate = false;
assert(isempty(gcp('nocreate')),'Use a dedicated MATLAB session.');
test_hybridLatentMoments;
cases = {};
for model = {'HMXL','HMNL','MIMIC'}
    for type = 0:6
        for variant = 0:2
            C = hybrid_fixture('HMXL',mod(variant,3),mod(variant,2),variant,variant,0,mod(variant,2),[2 2 4 7]);
            a = C.Args; o = a{9}; o.GPU = 'cpu'; o.HybridBlockSize = 2;
            o.NVarMeaExp = double(variant > 0); o.MeaExpMatrix = double(variant > 0);
            o.MeaSpecMatrix = type; o.MeaMatrix = ones(2,1); o.RealMin = mod(variant,3);
            y = [-.2;.4;.7;-.1];
            if ismember(type,[1 2]), y = [-2;1;4;1]; elseif type >= 3, y = [0;1;3;2]; end
            if variant == 2, y(2) = NaN; end
            o.MissingIndMea = isnan(y);
            a{6} = y; a{7} = zeros(4,0);
            if o.NVarMeaExp, a{7} = [-.3;.2;.1;.5]; end
            q = 1+o.NLatent+o.NVarMeaExp;
            switch type
                case 0, bm = [.1;.15;-.2;.07*ones(o.NVarMeaExp,1);log(.8)];
                case 1, bm = .15*sin((1:2*q)');
                case 2, bm = [.15;-.2;.07*ones(o.NVarMeaExp,1);-.3;log(.8)];
                case 3, bm = .1*sin((1:q)');
                case 4, bm = [.1*sin((1:q)');log(1.7)];
                case 5, bm = .1*sin((1:2*q)');
                case 6, bm = [.1*sin((1:2*q)');log(1.7)];
            end
            a{10} = [a{10}(1:end-4);bm]; a{9} = o;
            if strcmp(model{1},'HMNL')
                nc = 3+3*o.FullCov;
                a{10}(4:3+nc) = []; a{8} = a{8}(4:end,:);
                if ~o.mCT, a{3} = a{3}'; end
            elseif strcmp(model{1},'MIMIC')
                a = {a{5},a{6},a{7},a{8}(4:end,:),o,a{10}(end-numel(bm)-3:end)};
            end
            cases{end+1} = struct('name',sprintf('%s_type%d_variant%d',model{1},type,variant),...
                'model',model{1},'args',{a}); %#ok<AGROW>
        end
    end
end
% Binary categories exercise empty threshold-increment arrays on the GPU.
for k = 1:numel(cases)
    c = cases{k}; a = c.args; o = a{end-1};
    if ~ismember(o.MeaSpecMatrix,[1 2]) || o.NVarMeaExp || any(o.MissingIndMea,'all'), continue; end
    if strcmp(c.model,'MIMIC'), iy = 2; else, iy = 6; end
    a{iy}(a{iy} == 4) = 1;
    if o.MeaSpecMatrix == 1, a{end}(end-o.NLatent:end) = []; else, a{end}(end) = []; end
    c.args = a; c.name = [c.name,'_binary']; cases{end+1} = c; %#ok<AGROW>
end
% The fully correlated latent-error layout is also checked with two latent variables.
C = hybrid_fixture('HMXL',1,1,3,3,2,1,[2 2 4 7]); a = C.Args;
o = a{9}; o.FullCov = 2; o.GPU = 'cpu'; o.HybridBlockSize = 1;
v = tril(.08*ones(5)); v(1:6:end) = .4;
mask = tril(true(5)); mask([19 25]) = false;
a{10} = [a{10}(1:3);v(mask);a{10}(10:end)]; a{9} = o;
cases{end+1} = struct('name','HMXL_FullCov2_two_latents','model','HMXL','args',{a});
for counts = {[2 2 1 7],[2 2 4 1]}
    C = hybrid_fixture('HMXL',1,1,0,0,0,0,counts{1}); C.Args{9}.HybridBlockSize = 1;
    cases{end+1} = struct('name',sprintf('HMXL_NP%d_R%d',counts{1}(3:4)),...
        'model','HMXL','args',{C.Args}); %#ok<AGROW>
end
rows = struct([]);
serial = cell(size(cases));
for k = 1:numel(cases)
    c = cases{k}; a = c.args; entry = ['LL_',lower(c.model)];
    [f,g] = feval(entry,a{:}); fv = feval(entry,a{:});
    assert(all(isfinite([f(:);g(:)])) && relative(f,fv) < 1e-10,c.name);
    J = zeros(size(g));
    for j = 1:numel(a{end})
        h = 1e-5*max(1,abs(a{end}(j))); ap = a; am = a;
        ap{end}(j) = ap{end}(j)+h; am{end}(j) = am{end}(j)-h;
        J(:,j) = (feval(entry,ap{:})-feval(entry,am{:}))/(2*h);
    end
    fd = relative(g,J); assert(fd < 5e-6,'Gradient mismatch %s: %.3g',c.name,fd);
    old = a; old{end-1}.HybridBlockSize = 0;
    baselineError = ''; delta = NaN;
    try
        fo = feval(entry,old{:}); delta = relative(f,fo);
    catch ex
        baselineError = ex.message;
    end
    % Missing categorical/count measurements were not handled by every legacy model.
    if isempty(baselineError) && all(isfinite(fo))
        observedRows = ~any(a{end-1}.MissingIndMea,2);
        observedDelta = relative(f(observedRows),fo(observedRows));
        assert(observedDelta < 1e-8,'Legacy observed value mismatch %s: %.3g',c.name,observedDelta);
    end
    a{end-1}.HybridBlockSize = 1; [f1,g1] = feval(entry,a{:});
    a{end-1}.HybridBlockSize = 4; [f4,g4] = feval(entry,a{:});
    assert(relative(f1,f4) < 1e-10 && relative(g1,g4) < 1e-9,'Block dependence %s',c.name);
    serial{k} = struct('f',f,'g',g);
    rows(k).Case = c.name; rows(k).FD = fd; rows(k).LegacyValue = delta;
    rows(k).LegacyError = baselineError; rows(k).BlockGradient = relative(g1,g4);
    fprintf('%d/%d %s FD %.3g legacy %.3g\n',k,numel(cases),c.name,fd,delta);
    report = struct2table(rows,'AsArray',true); writetable(report,fullfile(outDir,'hybrid.csv'));
end
if workers > 0
    pool = parpool('Processes',workers); cleanup = onCleanup(@() delete(pool)); %#ok<NASGU>
    for k = 1:numel(cases)
        c = cases{k}; [f,g] = feval(['LL_',lower(c.model)],c.args{:});
        assert(relative(f,serial{k}.f) < 1e-10 && relative(g,serial{k}.g) < 1e-9,c.name);
    end
end
save(fullfile(outDir,'hybrid.mat'),'cases','serial','report');
fprintf('HYBRID_TESTS_PASSED: %d cases, workers=%d\n',numel(cases),workers);
end

function r = relative(a,b)
r = max(abs(a-b),[],'all')/max(1,max(abs(b),[],'all'));
end

function setAutoCreate(settings,value)
settings.Pool.AutoCreate = value;
end
