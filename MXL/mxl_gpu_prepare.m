function data = mxl_gpu_prepare(C,needGradient)
% Upload immutable inputs once for the blocked GPU likelihood.
if nargin < 2, needGradient = true; end
opt = C.EstimOpt;
K = opt.NVarA;
NP = opt.NP;
R = opt.NRep;
T = opt.NAlt*opt.NCT;
MCount = opt.NVarM;
SCount = opt.NVarS;
assert(ismember(opt.FullCov,[0 1]) && all(ismember(opt.Dist,[-1 0 1])),...
    'The GPU path supports only normal/lognormal/fixed FullCov 0/1.');
assert(~isfield(opt,'NVarNLT') || opt.NVarNLT == 0,...
    'DCE:GPU:UnsupportedNLT',...
    'Nonlinear transformations are not supported by the GPU path.');
assert(~isfield(opt,'Johnson') || opt.Johnson == 0,...
    'Johnson distributions are not supported by the GPU path.');
assert(~isfield(opt,'ExpB') || isempty(opt.ExpB),...
    'DCE:GPU:UnsupportedExpB',...
    'ExpB is not supported by the GPU path.');
assert(numel(opt.Dist) == K && opt.WTP_space >= 0 && opt.WTP_space <= K,...
    'Unexpected distribution or WTP parameter layout.');
if opt.WTP_space > 0
    assert(numel(opt.WTP_matrix) == K-opt.WTP_space &&...
        all(ismember(opt.WTP_matrix,K-opt.WTP_space+1:K)),...
        'WTP mappings must target the final cost coefficient rows.');
end
if opt.FullCov == 1 && needGradient
    [i,j] = find(tril(ones(K)));
    assert(isfield(opt,'indx1') && isfield(opt,'indx2') &&...
        isequal(opt.indx1(:),i) && isequal(opt.indx2(:),j),...
        'The GPU path requires the standard Cholesky derivative order.');
end
if ~isfield(opt,'mCT'), opt.mCT = 0; end
if ~isfield(opt,'WTP_matrix'), opt.WTP_matrix = []; end
assert(isa(C.YY,'double') && isa(C.XXa,'double') && isa(C.err,'double') &&...
    (isempty(C.XXm) || isa(C.XXm,'double')) &&...
    (isempty(C.Xs) || isa(C.Xs,'double')),...
    'GPU inputs must be CPU double arrays.');
assert(numel(C.YY) == T*NP && numel(C.XXa) == T*K*NP &&...
    numel(C.err) == K*R*NP,'Unexpected choice data or draw dimensions.');
available = ~isnan(reshape(C.YY,[T,NP]));
chosen = reshape(C.YY,[T,NP]) == 1;
X = reshape(C.XXa,[T,K,NP]);
X(repmat(reshape(~available,[T,1,NP]),[1,K,1])) = 0;
chosenX = permute(sum(X.*reshape(chosen,[T,1,NP]),1),[2 1 3]);
missing = ~any(reshape(available,[opt.NAlt,opt.NCT,NP]),1);

if MCount == 0
    M = zeros(0,1,NP);
elseif opt.mCT ~= 0
    assert(numel(C.XXm) == MCount*T*NP,'Unexpected row-varying Xm dimensions.');
    M = reshape(C.XXm,[MCount,T*NP]);
    M(:,~available(:)) = 0;
    M = reshape(M,[MCount,T,NP]);
else
    assert(numel(C.XXm) == MCount*NP,'Unexpected respondent-specific Xm dimensions.');
    M = reshape(C.XXm,[MCount,1,NP]);
end
if SCount == 0
    Xs = zeros(T,0,NP);
else
    assert(numel(C.Xs) == T*NP*SCount,'Unexpected Xs dimensions.');
    Xs = reshape(C.Xs,[T*NP,SCount]);
    Xs(~available(:),:) = 0;
    Xs = permute(reshape(Xs,[T,NP,SCount]),[1 3 2]);
end
data = struct('opt',opt,'X',gpuArray(X),...
    'Y',gpuArray(reshape(chosen,[T,1,NP])),...
    'A',gpuArray(reshape(available,[T,1,NP])),...
    'missing',gpuArray(reshape(missing,[1,opt.NCT,1,NP])),...
    'E',gpuArray(reshape(C.err,[K,R,NP])),...
    'M',gpuArray(M),'Xs',gpuArray(Xs),'chosenX',gpuArray(chosenX));
end
