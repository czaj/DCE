function test_hybrid_rebase
cases=0;
for model={'HMXL','HMNL','MIMIC'}
 for fullCov=0:2
  if fullCov==2 && ~strcmp(model{1},'HMXL'),continue;end
  for measurement=[0 2]
   C=hybrid_fixture('HMXL',1,min(fullCov,1),1,0,0,0,[2 2 4 7]);a=C.Args;
   if fullCov==2
    v=tril(.08*ones(5));v(1:6:end)=.4;mask=tril(true(5));mask([19 25])=false;
    a{10}=[a{10}(1:3);v(mask);a{10}(10:end)];a{9}.FullCov=2;
   end
   if measurement==2
    a{6}=[0;1;2;1];a{9}.MeaSpecMatrix=2;a{10}(end-3:end)=[.4;-.25;-.3;log(.8)];
   end
   if strcmp(model{1},'HMNL')
    nc=3+3*fullCov;a{10}(4:3+nc)=[];a{8}=a{8}(4:end,:);a{3}=a{3}';
   elseif strcmp(model{1},'MIMIC')
    a={a{5},a{6},a{7},a{8}(4:end,:),a{9},a{10}(end-7:end)};
   end
   o=a{end-1};di=numel(a)-2;old=a;
   draws=reshape(sin((1:size(a{di},1)*4*13)/7),size(a{di},1),[]);
   [b,change]=hybrid_rebase_start(model{1},a,draws);
   a{di}=draws;a{end-1}.NRep=13;a{end}=b;
   d0=hybrid_prepare(model{1},old);d1=hybrid_prepare(model{1},a);
   p0=hybrid_parameters(d0,old{end});p1=hybrid_parameters(d1,b);
   raw=reshape(cos(1:o.NLatent*13*4),o.NLatent,13,4);
   lv0=(raw-p0.mean)./p0.s;lv1=(raw-p1.mean)./p1.s;
   assert(max(abs(lv0-change.multiplier.*lv1-change.shift),[],'all')<1e-12);
   % Same raw latent values, so the transformed measurement indices/probabilities agree.
   d0.opt.NRep=13;
   v0=hybrid_measurement(d0,1:4,lv0,p0.bmea);
   v1=hybrid_measurement(d1,1:4,lv1,p1.bmea);
   assert(max(abs(v0-v1),[],'all')<1e-12);
   if ~strcmp(model{1},'MIMIC')
    u0=p0.mu+pagemtimes(p0.bl,lv0);u1=p1.mu+pagemtimes(p1.bl,lv1);
    assert(max(abs(u0-u1),[],'all')<1e-12);
   end
   cases=cases+1;
  end
 end
end
fprintf('HYBRID_REBASE_PASS %d cases\n',cases);
end
