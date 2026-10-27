function test_sobol_draws(outDir)
% New uniform grids pass through existing likelihood/score and public paths.
if ~isfolder(outDir),mkdir(outDir);end
settings=parallel.Settings;oldAuto=settings.Pool.AutoCreate;
restore=onCleanup(@() setAuto(settings,oldAuto)); %#ok<NASGU>
settings.Pool.AutoCreate=false;assert(isempty(gcp('nocreate')));
for model={'MXL','HMXL','LCMXL'}
    c=hybrid_fixture(model{1},1,1,0,0,0,0,[1 2 4 32]);a=c.Args;e=a{end-1};b=a{end};
    e=dceSimulationDefaults(e);e.GPU='cpu';e.Draws=6;dim=size(a{end-2},1);
    a{end-2}=norminv(dceQuasiDraws(e,dim))';a{end-1}=e;
    [f,J]=feval(c.Entry,a{:});assert(all(isfinite(f)) && all(isfinite(J),'all'));
    for j=1:numel(b)
        plus=a;minus=a;plus{end}(j)=b(j)+1e-5;minus{end}(j)=b(j)-1e-5;
        numeric=(feval(c.Entry,plus{:})-feval(c.Entry,minus{:}))/2e-5;
        assert(max(abs(numeric-J(:,j))./(1+abs(numeric)))<2e-5,model{1});
    end
    fprintf('SOBOL_LIKELIHOOD_GRADIENT_PASS %s %d parameters\n',model{1},numel(b));
end
c=hybrid_fixture('HMXL',0,1,0,0,0,0);a=c.Args;T=9;N=4;
I=struct('Y',a{1}(:),'Xa',reshape(permute(a{2},[1 3 2]),T*N,3),...
    'Xstr',repelem(a{5}(:,2),T,1),'Xmea',repelem(a{6},T,1));
base=struct('NAlt',3,'NCT',3,'NP',N,'NLatent',1,'Seed1',1234,...
    'FullCov',1,'WTP_space',0,'StrNorm',0,'MeaSpecMatrix',0,'MeaMatrix',1,...
    'Display',0,'SaveTxtOutput',0,'Scores',0,'RobustStd',0,'NumGrad',0,...
    'CheckSeparation',0,'GPU','cpu','OutputDir',outDir,'Drawskeep',1,'HybridFitOnly',1);
for model={'HMXL','HMNL','MIMIC'}
    name=model{1};e=base;e.ProjectName=['sobol_',name];b=a{10};b(end-4)=[];
    if strcmp(name,'HMNL'),b(4:9)=[];end
    if strcmp(name,'MIMIC'),b=b(end-3:end);end
    [input,old,e,oo]=DataCleanDCE(I,e);assert(e.NRep==1024 && e.HaltonSkip==0 && e.HaltonLeap==0);
    oo.Algorithm='quasi-newton';oo.GradObj='on';oo.Hessian='off';oo.MaxIterations=1;oo.Display='off';
    old.(name)=struct('b0',b);r=feval(name,input,old,e,oo);
    if strcmp(name,'HMXL'),dim=r.EstimOpt.NVarA+1;elseif strcmp(name,'HMNL'),dim=1;else,dim=r.EstimOpt.NVarA+2;end
    expected=norminv(dceQuasiDraws(r.EstimOpt,dim))';
    if strcmp(name,'MIMIC'),expected=expected(end,:);end
    assert(isequal(r.err,expected) && isfinite(r.LL),name);
    save(fullfile(outDir,[name,'.mat']),'r');
    fprintf('SOBOL_PUBLIC_DEFAULT_PASS %s\n',name);
end
% Demand keeps budgets and block invariance with independent error scrambles.
for model={'MDCEV','MMDCEV'}
    c=vce_fixture(model{1},3,1,1,1,5,32,2,3);c.opt.NSim=16;c.opt.Seed1=123;c.opt.Seed2=456;
    c.opt=dceSimulationDefaults(c.opt);c.data.I=100*ones(c.opt.NP*c.opt.NCT,1);
    r=struct('INPUT',c.data,'EstimOpt',c.opt,'bhat',c.b);fun=[model{1},'_demand'];
    x=feval(fun,r,[]);r.EstimOpt.VCEBlockSize=1;y=feval(fun,r,[]);
    assert(max(abs(x-y),[],'all')<1e-9);
    assert(max(abs(sum(x.*c.data.priceMat,1)-c.data.I'))<1e-9);
    fprintf('SOBOL_DEMAND_BUDGET_BLOCK_PASS %s\n',model{1});
end
end

function setAuto(settings,value)
settings.Pool.AutoCreate=value;
end
