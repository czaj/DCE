function d = hybrid_prepare(model,args)
% Immutable double-precision inputs and pooled latent moments, cached by caller.
o = args{end-1};
N = o.NP; R = o.NRep; L = o.NLatent;
d = struct('model',model,'opt',o);
assert(all(ismember(o.MeaSpecMatrix,0:6)),'DCE:Hybrid:Measurement','Unsupported measurement specification.');
if strcmp(model,'MIMIC')
    d.Xstr = args{1}; d.Ymea = args{2}; d.Xexp = args{3}; e = args{4};
    K = 0;
else
    K = o.NVarA;
    T = o.NAlt*o.NCT;
    if strcmp(model,'HLC')
        y = reshape(args{1},o.NAlt,o.NCT,N,o.NClass);
        d.Y = reshape(y(:,:,:,1),T,1,N) == 1;
        d.A = ~isnan(reshape(y(:,:,:,1),T,1,N));
        d.X = permute(reshape(args{2},T,N,K),[1 3 2]);
        d.Xc = args{3}; d.Xstr = args{4}; d.Ymea = args{5};
        d.Xexp = args{6}; e = args{7};
    else
        d.Y = reshape(args{1},T,1,N) == 1;
        d.A = ~isnan(reshape(args{1},T,1,N));
        d.X = reshape(args{2},T,K,N);
        d.Xstr = args{5}; d.Ymea = args{6}; d.Xexp = args{7}; e = args{8};
        if ~isfield(o,'mCT'), o.mCT = 0; end
        if o.NVarM == 0
            d.M = zeros(0,1,N);
        elseif o.mCT
            d.M = permute(reshape(args{3},T,o.NVarM,N),[2 1 3]);
            d.M(repmat(permute(~d.A,[2 1 3]),o.NVarM,1,1)) = 0;
        elseif strcmp(model,'HMNL')
            d.M = reshape(args{3}',o.NVarM,1,N);
        else
            d.M = reshape(args{3},o.NVarM,1,N);
        end
        ns = o.NVarS-o.ScaleLV*L;
        d.Xs = permute(reshape(args{4},T,N,ns),[1 3 2]);
        d.Xs(repmat(~d.A,1,ns,1)) = 0;
    end
    d.X(repmat(~d.A,1,K,1)) = 0;
end
assert(isa(e,'double') && isa(d.Xstr,'double') && size(e,2) == N*R,...
    'DCE:Hybrid:Inputs','Hybrid calculations require double inputs with NP*NRep draw columns.');
assert(size(e,1) == L+K*strcmp(model,'HMXL'),'DCE:Hybrid:Inputs','Unexpected number of draw rows.');
assert(isequal(size(d.Xstr),[N,o.NVarStr]) && size(d.Ymea,1) == N,...
    'DCE:Hybrid:Inputs','Unexpected structural or measurement dimensions.');
if ~isfield(o,'RealMin'), o.RealMin = 1; end
if ~isfield(o,'ProbMin') || isempty(o.ProbMin), o.ProbMin = 1e-300; end
o.ProbMin = max(realmin,min(o.ProbMin,1e-6));
if ~isfield(o,'MaxExp') || isempty(o.MaxExp), o.MaxExp = 50; end
if ~isfield(o,'SigmaMin') || isempty(o.SigmaMin), o.SigmaMin = 1e-8; end
d.observed = isfinite(d.Ymea);
if isfield(o,'MissingIndMea') && ~isempty(o.MissingIndMea)
    assert(isequal(size(o.MissingIndMea),size(d.Ymea)),'Unexpected measurement missing mask.');
    d.observed = d.observed & o.MissingIndMea == 0;
end
d.categories = cell(1,size(d.Ymea,2));
for i = 1:size(d.Ymea,2)
    d.categories{i} = unique(d.Ymea(d.observed(:,i),i));
    if ismember(o.MeaSpecMatrix(i),[1 2])
        assert(numel(d.categories{i}) >= 2,'DCE:Hybrid:Categories',...
            'A categorical measurement needs at least two observed categories.');
    end
end
d.Ymea(~d.observed) = 0;
d.moments = hybridLatentMoments(d.Xstr,e,R);
d.E = reshape(e,size(e,1),R,N);
d.opt = o;
end
