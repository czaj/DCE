function C = hybrid_fixture(model,wtp,fullCov,meanMode,scaleMode,missingMode,scaleLV,counts)
if nargin < 8, counts = [1 2 4 7]; end
K = 3;
A = 3;
T = 3;
LCount = counts(1);
CCount = counts(2);
NP = counts(3);
R = counts(4);
rows = A*T;
X = zeros(rows,K,NP);
Y = zeros(rows,NP);
for n = 1:NP
    j = (1:rows)';
    X(:,:,n) = [sin(j/2+n/3),-.3-.02*j-.01*n,-.7-.03*j-.02*n];
    for t = 1:T, Y((t-1)*A+mod(t+n-2,A)+1,n) = 1; end
end
if missingMode > 0
    Y(A+(1:A),1) = NaN;
    X(A+(1:A),:,1) = NaN;
    Y(2*A+(1:A),3) = NaN;
    X(2*A+(1:A),:,3) = NaN;
end
if missingMode > 1
    Y(1:A,4) = NaN;
    X(1:A,:,4) = NaN;
    for n = 1:2
        block = (n-1)*A+(1:A);
        unavailable = block(find(Y(block,n) == 0,1,'last'));
        Y(unavailable,n) = NaN;
        X(unavailable,:,n) = NaN;
    end
end
MCount = double(meanMode > 0);
meanLong = zeros(rows,NP);
scaleLong = zeros(rows,NP);
for n = 1:NP
    task = ceil((1:rows)'/A);
    alt = mod((1:rows)'-1,A)+1;
    meanLong(:,n) = .15+.03*n;
    scaleLong(:,n) = .10+.02*n;
    if meanMode >= 2, meanLong(:,n) = meanLong(:,n)+.07*task; end
    if meanMode == 3, meanLong(:,n) = meanLong(:,n)+.04*alt; end
    if scaleMode >= 2, scaleLong(:,n) = scaleLong(:,n)+.05*task; end
    if scaleMode == 3, scaleLong(:,n) = scaleLong(:,n)+.03*alt; end
end
if meanMode >= 2, meanLong(isnan(Y)) = NaN; end
scaleLong(isnan(Y)) = NaN;
if meanMode == 0
    Xm = zeros(0,NP);
elseif meanMode == 1
    Xm = meanLong(1,:);
elseif strcmp(model,'HMXL')
    Xm = reshape(meanLong,[rows,1,NP]);
else
    Xm = reshape(meanLong,[1,rows*NP]);
end
if scaleMode == 0, Xs = zeros(rows*NP,0); else, Xs = scaleLong(:); end
NVarS = size(Xs,2);
opt = struct('NAlt',A,'NCT',T,'NP',NP,'NRep',R,'NVarA',K,...
    'NVarM',MCount,'NVarS',NVarS,'FullCov',fullCov,'NumGrad',0,...
    'Dist',[0 1 0],'WTP_space',wtp,'WTP_matrix',[],...
    'Triang',[],'NVarNLT',0,'NLTVariables',[],'NLTType',[],...
    'Johnson',0,'RealMin',0,'ExpB',[],'mCT',double(meanMode >= 2));
if wtp == 1, opt.Dist = [0 0 1]; opt.WTP_matrix = [3 3]; end
if wtp == 2, opt.Dist = [0 1 1]; opt.WTP_matrix = 2; end
[ci,cj] = find(tril(ones(K)));
opt.indx1 = ci';
opt.indx2 = cj';
opt.DiagIndex = find(ci == cj);
available = reshape(~isnan(Y),[A,T,NP]);
count = reshape(sum(available,1),[T,NP]);
opt.MissingCT = count == 0;
opt.NCTMiss = sum(~opt.MissingCT,1)';
opt.NAltMiss = sum(count,1)'./opt.NCTMiss;
count(opt.MissingCT) = A;
opt.NAltMissInd = count;
opt.NAltMissIndExp = repelem(count,A,1);
opt.NLatent = LCount;
opt.NVarStr = 2;
opt.NVarMeaExp = 0;
opt.MeaMatrix = ones(LCount,1);
opt.MeaSpecMatrix = 0;
opt.MeaExpMatrix = 0;
opt.MissingIndMea = false(NP,1);
opt.indx3 = [];
opt.ScaleLV = scaleLV;
opt.GPU = 'cpu';
mu = [-.2;.1;-.3];
V = [.2 0 0;.03 .25 0;-.02 .04 .3];
if fullCov, covariance = V(sub2ind([K,K],ci,cj)); else, covariance = diag(V); end
if MCount, bm = [.08;-.04;.06]; else, bm = zeros(0,1); end
bs = .08*ones(NVarS,1);
if strcmp(model,'MXL')
    B = [mu;covariance;bm;bs];
    E = reshape(.9*sin((1:K*R*NP)/4),[K,R*NP]);
    args = {Y,X,Xm,Xs,E,opt,B};
    entry = 'LL_mxl';
elseif strcmp(model,'HMXL')
    opt.NVarS = NVarS+LCount*scaleLV;
    structuralX = [-.4;-.1;.3;.8];
    Xstr = [ones(NP,1),structuralX(1:NP)];
    measurementY = [-.2;.4;.7;-.1];
    Xmea = measurementY(1:NP);
    if missingMode == 2
        Xmea(2) = NaN;
        opt.MissingIndMea(2) = true;
    end
    E = reshape(.9*sin((1:(K+LCount)*R*NP)/4),[K+LCount,R*NP]);
    latentInteraction = [.1 -.07;-.08 .11;.12 .05];
    structuralB = [.15 -.1;.6 .45];
    measurementB = [.4;-.25];
    bl = latentInteraction(:,1:LCount);
    bstr = structuralB(:,1:LCount);
    latentScale = (.07-.02*(0:LCount-1)')*scaleLV;
    if ~scaleLV, latentScale = zeros(0,1); end
    B = [mu;covariance;bm;bl(:);bs;latentScale;bstr(:);...
        .1;measurementB(1:LCount);log(.8)];
    args = {Y,X,Xm,Xs,Xstr,Xmea,zeros(NP,0),E,opt,B};
    entry = 'LL_hmxl';
else
    opt.NClass = CCount;
    opt.NVarC = 2;
    opt.Dist = repmat(opt.Dist(:),[1,CCount]);
    opt.indx1 = reshape(ci+K*(0:CCount-1),1,[]);
    opt.indx2 = reshape(cj+K*(0:CCount-1),1,[]);
    classX = [-.3;.1;.4;.8];
    Xc = [ones(NP,1),classX(1:NP)];
    E = reshape(.9*cos((1:K*CCount*R*NP)/5),[K*CCount,R*NP]);
    classMeans = mu+.1*(0:CCount-1);
    classCovariance = covariance*(1+.1*(0:CCount-1));
    classScale = bs+.03*(0:CCount-1);
    membership = [.2;-.15]*(1+.25*(0:CCount-2));
    B = [classMeans(:);classCovariance(:);classScale(:);membership(:)];
    args = {Y,X,Xc,Xs,E,opt,B};
    entry = 'LL_lcmxl';
end
C = struct('Name',sprintf('%s_WTP%d_Cov%d_Xm%d_Xs%d_Missing%d_LVScale%d',...
    model,wtp,fullCov,meanMode,scaleMode,missingMode,scaleLV),...
    'Model',model,'Entry',entry,'Args',{args});
end
