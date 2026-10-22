function c = vce_fixture(model,profile,fullCov,u,m,NP,R,T,K)
if nargin < 6, NP = 4; end
if nargin < 7, R = 7; end
if nargin < 8, T = 3; end
if nargin < 9, K = 3; end
A = 3; N = NP*T;
o = struct('NP',NP,'NCT',T,'NAlt',A,'NRep',R,'NVarA',K,'NVarM',m,...
    'NVarU',u,'FullCov',fullCov,'Profile',profile,'RealMin',1,'GPU','cpu','VCEBlockSize',2,...
    'Dist',zeros(1,K),'NumGrad',0);
o.Dist(end) = 1;
if profile == 1, o.SpecProfile = [1 1 2;0 0 0];
elseif profile == 2, o.SpecProfile = [0 0 0;1 2 2];
else, o.SpecProfile = [1 0 1;0 1 0]; end
o.NVarP = 2;
d.Y = repmat([.7 0 .5;0 .8 .6;1 1.2 1],1,ceil(N/3)); d.Y = d.Y(:,1:N);
d.priceMat = [ones(1,N);1.2+.1*sin(1:N);2*ones(1,N)];
X = reshape(.2*sin(1:A*N*K),A*N,K); X(3:A:end,1) = 1;
d.Xu = zeros(N,0); if u, d.Xu = [ones(N,1),reshape(.2*cos(1:N*u),N,u)]; end
b = .1*cos((1:K)');
if strcmp(model,'MMDCEV')
    d.Xa = permute(reshape(X,A*T,NP,K),[1 3 2]);
    d.err = reshape(.6*sin((1:K*R*NP)/5),K,R*NP);
    d.Xm = reshape(.3*cos(1:NP*m),m,NP);
    V = tril(.03*ones(K)); V(1:K+1:end) = .2;
    [o.indx1,o.indx2] = find(tril(true(K))); o.indx1 = o.indx1'; o.indx2 = o.indx2';
    if fullCov, b = [b;V(tril(true(K)))]; else, b = [b;diag(V)]; end
else
    d.Xa = X; d.Xm = reshape(.3*cos(1:N*m),N,m);
end
b = [b;.1*sin((1:K*m)');.2*cos((1:2*(1+u))');-.1];
o.BActive = ones(1,numel(b));
d.W = repelem(1+.2*sin((1:NP)'),T);
c = struct('model',model,'data',d,'opt',o,'b',b,...
    'name',sprintf('%s_P%d_C%d_U%d_M%d_N%d_R%d_T%d_K%d',model,profile,fullCov,u,m,NP,R,T,K));
end
