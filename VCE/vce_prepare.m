function d = vce_prepare(model,data,o)
% Keep immutable inputs in respondent pages; MDCEV has one task per page.
d.opt = o; d.mixed = strcmp(model,'MMDCEV');
d.N = o.NP; d.T = o.NCT; d.R = 1;
if d.mixed, d.R = o.NRep; else, d.N = o.NP*o.NCT; d.T = 1; end
A = o.NAlt; K = o.NVarA;
assert(isa(data.Y,'double') && isequal(size(data.Y),[A,d.T*d.N]),'DCE:VCE:Shape','Incorrect quantity dimensions.');
assert(all(isfinite(data.Y),'all') && all(data.Y >= 0,'all') && all(any(data.Y > 0,1)),...
    'DCE:VCE:Quantities','Every task must have finite nonnegative quantities and at least one consumed good.');
assert(isequal(size(data.priceMat),size(data.Y)) && all(isfinite(data.priceMat),'all') && all(data.priceMat > 0,'all'),...
    'DCE:VCE:Prices','Prices must be finite and strictly positive.');
d.Y = data.Y; d.price = data.priceMat;
if d.mixed
    assert(size(data.Xa,1) == A*d.T && size(data.Xa,2) == K && size(data.Xa,3) == d.N,'DCE:VCE:Shape','Incorrect attribute pages.');
    assert(isequal(size(data.err),[K,d.R*d.N]),'DCE:VCE:Draws','Incorrect draw dimensions.');
    assert(all(ismember(o.Dist,[-1 0 1])),'DCE:VCE:Distribution','Supported distributions are fixed, normal and lognormal.');
    d.X = data.Xa; d.E = reshape(data.err,K,d.R,d.N);
    if o.NVarM > 0, d.M = reshape(data.Xm,o.NVarM,1,[]); end
else
    assert(isequal(size(data.Xa),[A*d.N,K]),'DCE:VCE:Shape','Incorrect attribute rows.');
    d.X = permute(reshape(data.Xa,A,d.N,K),[1 3 2]);
    d.E = zeros(K,1,d.N);
    if o.NVarM > 0, d.M = reshape(data.Xm',o.NVarM,1,[]); end
end
if o.NVarM == 0, d.M = zeros(0,1,d.N); end
assert(size(d.M,3) == d.N,'DCE:VCE:Shape','Incorrect mean covariate dimensions.');
if o.NVarU == 0, d.U = ones(d.T*d.N,1); else, d.U = data.Xu; end
assert(isequal(size(d.U),[d.T*d.N,1+o.NVarU]),'DCE:VCE:Shape','Incorrect profile covariate dimensions.');
for name = {'X','E','M','U'}
    assert((isnumeric(d.(name{1})) || islogical(d.(name{1}))) && isreal(d.(name{1})) && all(isfinite(d.(name{1})),'all'),...
        'DCE:VCE:Inputs','Attributes, covariates and draws must be finite real numbers.');
    d.(name{1}) = double(d.(name{1}));
end
assert(isequal(size(o.SpecProfile),[2,A]) && all(o.SpecProfile >= 0 & o.SpecProfile == fix(o.SpecProfile),'all'),...
    'DCE:VCE:Profile','Incorrect profile specification.');
for row = 1:2
    groups = unique(o.SpecProfile(row,o.SpecProfile(row,:) > 0));
    assert(isempty(groups) || isequal(groups,1:numel(groups)),'DCE:VCE:Profile','Profile group numbers must be consecutive.');
end
end
