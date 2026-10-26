function test_hlc_cutoff_output(outDir)
% Exercise the public HLC cutoff offsets with fixed nonzero class/loadings.
if ~isfolder(outDir),mkdir(outDir);end
settings=parallel.Settings;previous=settings.Pool.AutoCreate;
cleanup=onCleanup(@() restore(settings,previous)); %#ok<NASGU>
settings.Pool.AutoCreate=false;
C=hybrid_fixture('HMXL',0,1,0,0,0,0);a=C.Args;
I=struct('Y',a{1}(:),'Xa',reshape(permute(a{2},[1 3 2]),36,3),...
 'Xmea',repelem([0 2;1 0;2 1;1 2],9,1));
E=struct('NAlt',3,'NCT',3,'NP',4,'NLatent',1,'NClass',2,'NRep',16,...
 'Draws',6,'Seed1',1234,'WTP_space',0,'StrNorm',0,'MeaSpecMatrix',[2 2],...
 'MeaMatrix',[1 1],'Display',0,'SaveTxtOutput',0,'Scores',0,'HessEstFix',1,...
 'RobustStd',0,'NumGrad',1,'CheckSeparation',0,'GPU','cpu','OutputDir',outDir,...
 'ProjectName','hlc_cutoffs','ConstVarActive',1);
b=[-.2;.1;-.3;.3;-.1;-.4;.2;.3;.4;-.3;log(.8);.5;.2;log(1.1)];
E.BActive=zeros(size(b));E.BActive([10 11 13 14])=1;
[I,old,E,O]=DataCleanDCE(I,E);old.HLC.b0=b;
O.Algorithm='sqp';O.SpecifyObjectiveGradient=true;O.MaxIterations=1;O.Display='off';
r=HLC(I,old,E,O);
assert(max(abs(r.bhat(E.BActive==0)-b(E.BActive==0)))<1e-12);
offset=numel(b)-size(r.DetailsM,1);
for item=1:2
 idx=(item-1)*3+(2:3);raw=r.bhat(offset+idx);grad=[1 exp(raw(2))];
 expected=sqrt(grad*r.ihess(offset+idx,offset+idx)*grad');
 assert(abs(r.DetailsM(idx(2),1)-sum([raw(1);exp(raw(2))]))<1e-10);
 if isreal(expected)&&isfinite(expected)&&expected>0
  assert(abs(r.DetailsM(idx(2),3)-expected)<1e-8);
 else
  assert(isnan(r.DetailsM(idx(2),3)));
 end
end
save(fullfile(outDir,'hlc_cutoffs.mat'),'r');fprintf('HLC_PUBLIC_CUTOFF_PASS\n');
end
function restore(settings,value)
settings.Pool.AutoCreate=value;
end
