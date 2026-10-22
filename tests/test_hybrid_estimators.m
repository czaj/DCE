function test_hybrid_estimators(outDir)
% Exercise HessEstFix=5 through public preprocessing, optimization and inference.
if ~isfolder(outDir), mkdir(outDir); end
settings = parallel.Settings; previousAuto = settings.Pool.AutoCreate;
restoreAuto = onCleanup(@() setAutoCreate(settings,previousAuto));
settings.Pool.AutoCreate = false;
assert(isempty(gcp('nocreate')),'Use a dedicated MATLAB session.');
C = hybrid_fixture('HMXL',0,1,0,0,0,0); a = C.Args;
T = 9; N = 4;
I = struct('Y',a{1}(:),'Xa',reshape(permute(a{2},[1 3 2]),T*N,3),...
    'Xstr',repelem(a{5}(:,2),T,1),'Xmea',repelem(a{6},T,1));
base = struct('NAlt',3,'NCT',3,'NP',N,'NLatent',1,'NRep',32,'Draws',1,...
    'Seed1',1234,'FullCov',1,'WTP_space',0,'StrNorm',0,'MeaSpecMatrix',0,...
    'MeaMatrix',1,'Display',0,'SaveTxtOutput',0,'Scores',0,'HessEstFix',5,...
    'RobustStd',0,'NumGrad',0,'CheckSeparation',0,'GPU','cpu','OutputDir',outDir);
models = {'HMXL','HMNL','MIMIC','MIMIC','MIMIC'}; draws = [1 1 1 2 6];
for k = 1:numel(models)
    name = models{k}; e = base; e.Draws = draws(k); e.Drawskeep = 1;
    e.ProjectName = sprintf('test_%s_draws%d',name,e.Draws);
    b = a{10}; b(end-4) = []; % Remove the normalized-away structural intercept.
    if strcmp(name,'HMNL'), b(4:9) = []; end
    if strcmp(name,'MIMIC'), b = b(end-3:end); end
    [input,old,e,oo] = DataCleanDCE(I,e);
    oo.Algorithm = 'quasi-newton'; oo.GradObj = 'on'; oo.Hessian = 'off';
    oo.MaxIterations = 1; oo.Display = 'off';
    old.(name) = struct('b0',b);
    r = feval(name,input,old,e,oo);
    assert(isfield(r,'hess') && all(isfinite(r.hess),'all') && any(r.hess(:) ~= 0),name);
    assert(r.EstimOpt.HessEstFix == 5,name);
    assert(size(r.err,1) == e.NLatent+strcmp(name,'HMXL')*r.EstimOpt.NVarA,name);
    save(fullfile(outDir,[e.ProjectName,'.mat']),'r');
    fprintf('PUBLIC_ESTIMATOR_HESSIAN_PASSED %s draws %d\n',name,e.Draws);
end
end

function setAutoCreate(settings,value)
settings.Pool.AutoCreate = value;
end
