function [b,change] = hybrid_rebase_start(model,args,newDraws)
% Re-express an existing start when changing draws, preserving latent indices.
% args are the prepared LL_hmxl/LL_hmnl/LL_mimic arguments, including old b.
% Supports unchanged model/sample, OLS/ordered measurements and no latent scale.
assert(ismember(model,{'HMXL','HMNL','MIMIC'}),'DCE:Hybrid:Rebase','Unsupported model.');
o=args{end-1};b=args{end}(:);oldb=b;
assert(all(ismember(o.MeaSpecMatrix,[0 2])) && (~isfield(o,'ScaleLV') || ~o.ScaleLV),...
 'DCE:Hybrid:Rebase','Rebasing currently supports OLS/ordered measurements without latent scale effects.');
drawIndex=numel(args)-2;
assert(isa(newDraws,'double') && isreal(newDraws) && all(isfinite(newDraws(:))) && ...
 size(newDraws,1)==size(args{drawIndex},1) && mod(size(newDraws,2),o.NP)==0,...
 'DCE:Hybrid:Rebase','New draws must preserve draw dimensions and respondent ordering.');
next=args;next{drawIndex}=newDraws;next{end-1}.NRep=size(newDraws,2)/o.NP;
d0=hybrid_prepare(model,args);d1=hybrid_prepare(model,next);
p0=hybrid_parameters(d0,b);p1=hybrid_parameters(d1,b);
a=p1.s./p0.s;shift=(p1.mean-p0.mean)./p0.s;
if ~strcmp(model,'MIMIC')
 indices=hybrid_lv_indices(model,o);gamma=reshape(b(indices),o.NVarA,o.NLatent);
 b(1:o.NVarA)=b(1:o.NVarA)+gamma*shift;
 b(indices)=reshape(gamma.*a',[],1);
end
pos=numel(b)-numel(p0.bmea);
for i=1:numel(o.MeaSpecMatrix)
 latent=find(o.MeaMatrix(:,i)==1);n=numel(latent);
 extra=o.NVarMeaExp*(o.MeaExpMatrix(i)~=0);
 if o.MeaSpecMatrix(i)==2
  idx=pos+(1:n);delta=b(idx)'*shift(latent);
  b(idx)=b(idx).*a(latent);b(pos+n+extra+1)=b(pos+n+extra+1)-delta;
  pos=pos+n+extra+numel(d0.categories{i})-1;
 else
  idx=pos+1+(1:n);delta=b(idx)'*shift(latent);
  b(idx)=b(idx).*a(latent);b(pos+1)=b(pos+1)+delta;
  pos=pos+2+n+extra;
 end
end
assert(pos==numel(b) && all(isfinite(b)),'DCE:Hybrid:Rebase','Parameter layout mismatch.');
if isfield(o,'BActive') && ~isempty(o.BActive)
 fixed=o.BActive(:)==0;
 assert(all(abs(b(fixed)-oldb(fixed))<1e-12),'DCE:Hybrid:Rebase',...
  'Rebasing would change fixed coefficients. Review identification/normalization before changing draws.');
end
change=struct('oldDraws',o.NRep,'newDraws',next{end-1}.NRep,...
 'oldMean',p0.mean,'newMean',p1.mean,'oldSD',p0.s,'newSD',p1.s,...
 'multiplier',a,'shift',shift);
end
