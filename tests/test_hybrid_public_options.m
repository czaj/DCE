function test_hybrid_public_options(outDir)
% Public preprocessing -> optimization -> native tables, with mean interactions.
if ~isfolder(outDir), mkdir(outDir); end
settings=parallel.Settings; previous=settings.Pool.AutoCreate;
cleanup=onCleanup(@() restore(settings,previous)); %#ok<NASGU>
settings.Pool.AutoCreate=false;
C=hybrid_fixture('HMXL',0,1,0,0,0,0);a=C.Args;T=9;N=4;
I=struct('Y',a{1}(:),'Xa',reshape(permute(a{2},[1 3 2]),T*N,3),...
 'Xm',repelem([-.4;-.1;.3;.8],T,1),'Xstr',repelem(a{5}(:,2),T,1),...
 'Xmea',repelem([0 2;1 0;2 1;1 2],T,1));
base=struct('NAlt',3,'NCT',3,'NP',N,'NLatent',1,'NRep',32,'Draws',6,...
 'Seed1',1234,'FullCov',1,'WTP_space',0,'StrNorm',0,'MeaSpecMatrix',[2 2],...
 'MeaMatrix',[1 1],'Display',0,'SaveTxtOutput',0,'Scores',0,'HessEstFix',5,...
 'MaxExp',[],'RobustStd',0,'NumGrad',0,'CheckSeparation',0,'GPU','cpu','OutputDir',outDir);
b=a{10};b(end-4)=[];
b=[b(1:9);.05;.06;.07;b(10:end-3);.4;-.3;log(.8);.5;.2;log(1.1)];
for model={'HMXL','HMNL','MIMIC'}
 name=model{1}; e=base; bb=b; ii=I;
 if strcmp(name,'HMNL'), bb(4:9)=[]; end
 if strcmp(name,'MIMIC'), bb=bb(end-6:end);ii=rmfield(ii,'Xm');end
 e.BActive=ones(size(bb));e.BActive(end-2)=0; % nonzero fixed second loading
 if ~strcmp(name,'MIMIC'),e.LVMatrix=[1;0;1];end
 for method=1:3
  ee=e;ee.ProjectName=sprintf('public_%s_%d',name,method);
  ee.HybridOptimize=method~=2;ee.ConstVarActive=method==2;
  if method==3,ee.HybridFitOnly=1;ee.HessEstFix=0;ee.HybridStartHessian=2*eye(numel(bb));end
  [input,old,ee,oo]=DataCleanDCE(ii,ee);oo.MaxIterations=1;oo.Display='off';
  if method~=2,oo.Algorithm='quasi-newton';else,oo.Algorithm='sqp';end
  old.(name)=struct('b0',bb);
  % Explicit start fields differ across HMXL variants; B_backup is its public fallback.
  global B_backup;B_backup=bb;clearStart=onCleanup(@() resetStart); %#ok<NASGU>
  r=feval(name,input,old,ee,oo);B_backup=[];
  assert(r.exitf~=-999 && abs(r.bhat(end-2)-bb(end-2))<1e-12);
  if method~=2,assert(r.output.hybridActiveParameters==nnz(r.EstimOpt.BActive));end
  if ~strcmp(name,'MIMIC')
   ix=hybrid_lv_indices(name,r.EstimOpt);
   assert(r.EstimOpt.BActive(ix(2))==0 && abs(r.bhat(ix(2))-bb(ix(2)))<1e-12);
  end
  if method==3
   assert(~isfield(r,'hess') && ~isfield(r,'std') && ~isfield(r,'DetailsM') && isfield(r,'InferenceStatus'));
   assert(r.output.hybridPreconditioned);
   fprintf('HYBRID_SCREENING_PASS %s\n',name);
   continue
  end
  offset=numel(r.bhat)-size(r.DetailsM,1);
  for item=1:2
   idx=(item-1)*3+(2:3);raw=r.bhat(offset+idx);
   gradient=[1 exp(raw(2))];expected=sqrt(gradient*r.ihess(offset+idx,offset+idx)*gradient');
   assert(abs(r.DetailsM(idx(2),1)-(raw(1)+exp(raw(2))))<1e-10);
   if isreal(expected)&&isfinite(expected)&&expected>0
    assert(abs(r.DetailsM(idx(2),3)-expected)<1e-8);
   else
    assert(isnan(r.DetailsM(idx(2),3)));
   end
  end
  save(fullfile(outDir,[ee.ProjectName '.mat']),'r');
  fprintf('HYBRID_PUBLIC_OPTIONS_PASS %s %d\n',name,method);
 end
end
end
function restore(settings,value)
settings.Pool.AutoCreate=value;
end
function resetStart
global B_backup;B_backup=[];
end
