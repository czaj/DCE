function test_hlc_memory(useGPU)
if nargin < 1, useGPU = false; end
C = hybrid_fixture('HMXL',1,1,0,0,0,0,[2 2 4 7]); a = C.Args;
o = a{9}; o.NClass = 3; o.NVarC = 2; o.GPU = 'cpu'; o.HybridBlockSize = 1;
Xc = [ones(4,1),[-.4;.1;.3;.7]];
Y = repmat(reshape(a{1},3,[]),1,3);
X = reshape(permute(a{2},[1 3 2]),[],3);
B = [.1*sin((1:9)');.1*cos((1:8)');a{10}(end-7:end)];
args = {Y,X,Xc,a{5},a{6},a{7},a{8}(4:end,:),o,B};
f = LL_hlc(args{:}); args{8}.HybridBlockSize = 0; old = LL_hlc(args{:});
assert(max(abs(f-old)) < 1e-10);
args{8}.HybridBlockSize = 3; other = LL_hlc(args{:});
assert(max(abs(f-other)) < 1e-10);
if useGPU, checkGPU(args,f); end
% A category occurring in only one respondent must retain the global layout.
args{5} = [-2;1;4;1]; args{8}.MeaSpecMatrix = 2;
args{9} = [B(1:end-4);.15;-.2;-.3;log(.8)];
args{8}.HybridBlockSize = 1; f = LL_hlc(args{:});
args{8}.HybridBlockSize = 0; old = LL_hlc(args{:});
assert(max(abs(f-old)) < 1e-10);
if useGPU, checkGPU(args,f); end
fprintf('HLC pooled normalization, singleton blocks and ordered categories passed.\n');
end

function checkGPU(args,expected)
args{8}.GPU = 'gpu'; args{8}.HybridBlockSize = 1;
[actual,~,backend] = hybrid_likelihood('HLC',args,false);
assert(strcmp(backend,'gpu') && max(abs(actual-expected)) < 1e-8);
fprintf('HLC forced GPU passed.\n');
end
