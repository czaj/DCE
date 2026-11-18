function test_mxl_importance_sampling(outDir)
% Wiktor Budzinski's wide-scale WTP example, with the original SML start and draws.
if nargin < 1, outDir = tempname; end
if ~isfolder(outDir), mkdir(outDir); end
rng(20261117);
N = 300; T = 12; A = 3;
X = zeros(A,T,N,5);
X(3,:,:,1) = 1;
X(1:2,:,:,2) = randi([0 1],2,T,N);
X(1:2,:,:,3) = randi([0 1],2,T,N);
X(1:2,:,:,4) = randi([0 3],2,T,N);
cost = [1 2 4 6 8];
X(1:2,:,:,5) = -cost(randi(5,2,T,N));
wtp = [-2;2;-1.5;1] + [2;1.5;1;0.5].*randn(4,N);
scale = exp(log(0.5)+2.5*randn(1,N));
Y = zeros(A,T,N);
for n = 1:N
    v = reshape(X(:,:,n,:),A*T,5)*(scale(n)*[wtp(:,n);1]);
    [~,choice] = max(reshape(v,A,T)-log(-log(rand(A,T))),[],1);
    Y(sub2ind([A,T,N],choice,1:T,repmat(n,1,T))) = 1;
end
I = struct('Y',Y(:),'Xa',reshape(X,[],5));
E = struct('NP',N,'NCT',T,'NAlt',A,'Display',0,'OutputDir',outDir, ...
    'ProjectName','wide WTP','NamesA',{{'Opt-out';'x1';'x2';'x3';'-cost'}});
[I,old,E,O] = DataCleanDCE(I,E);
E.WTP_space = 1; E.FullCov = 0; E.Dist = [0 0 0 0 1];
E.NRep = 1000; E.Seed1 = 179424673; E.SobolLayout = 'global'; E.HaltonSkip = 1;
E.RobustStd = 0; E.HessEstFix = 5; E.EffectiveMoments = 0; E.Drawskeep = 0;
old.MXL_d.b0 = [-1.96797904;2.06329773;-1.79898506;.97510629;-.84921456; ...
    1.89221266;1.59785858;.98006346;.51424786;2.28389752];
base = MXL(I,old,E,O);
assert(~isfield(base,'IS') && abs(base.std(1)-.0404183)<1e-5)
E.ImportanceSampling = 1; E.Display = 1; E.OutputFilePrefix = 'wide';
both = MXL(I,old,E,O);
is = both.IS;
expected = [.201624;.169214;.147747;.067652;.151368;.181569;.161647;.161032;.089861;.209876];
assert(is.success && max(abs(is.std-expected))<2e-4)
assert(abs(is.LL+2982.210807)<2e-4 && abs(is.bhat(1)+1.678767)<2e-4)
assert(norm(is.hess*is.ihess-eye(10),'inf')<1e-9 && max(abs(is.std-sqrt(diag(is.ihess))))<1e-10)
assert(norm(sum(is.jacobian,1)'-is.g)<1e-7 && abs(is.LL+sum(is.LLdetailed))<1e-8)
assert(abs(both.LL-base.LL)<1e-7 && max(abs(both.std-base.std))<1e-6)
assert(~isfield(both,'err') && isfile(fullfile(outDir,'wide.xlsx')) && isfile(is.OutputFile))
assert(all(is.ESS>=1 & is.ESS<=is.NRep) && max(is.RefreshStepSE)<1)
sheet = readcell(is.OutputFile,'Sheet','Results');
assert(abs(sheet{4,5}-is.std(1))<1e-9 && any(strcmp(sheet(:,1),'Fresh IS draws')))
fprintf('WIKTOR_IS_PASS pilot SE %.8f, IS SE %.8f, IS LL %.8f, refresh %.4f SE\n',base.std(1),is.std(1),is.LL,max(is.RefreshStepSE));

E.ISSeed = uint32(4321); E.ISNRep = uint32(1001); E.Display = 0;
other = MXL(I,old,E,O);
assert(other.IS.success && other.IS.NRep==1001 && other.IS.Seed==4321)
assert(other.IS.LL~=is.LL && abs(other.IS.std(1)/.206369-1)<.15)
assert(~isfield(other.IS,'OutputFile'))
% Duplicating every respondent's likelihood contribution doubles information.
E.ISNRep = 1000; E.ISSeed = 1234;
Iw = I; Iw.W = 2*I.W;
weighted = MXL(Iw,old,E,O);
assert(weighted.IS.success && abs(weighted.IS.LL-2*is.LL)<1e-3)
assert(max(abs(weighted.IS.std*sqrt(2)-is.std))<1e-4)
E.FullCov = 1;
try
    MXL(I,old,E,O);
    error('test:ExpectedError','Unsupported covariance was accepted.');
catch ME
    assert(strcmp(ME.identifier,'MXL:ISSpecification'),ME.message)
end
% This log-panel likelihood is finite even though multiplying probabilities underflows.
t = 1200;
logP = mxl_log_panel(repmat([0;1],t,1),repmat([true;false],t,1),1,2);
assert(abs(logP+t*log1p(exp(1)))<1e-8)
fprintf('ALL_MXL_IS_TESTS_PASS\n');
end
