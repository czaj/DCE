function test_mdcev_mmdcev
% Small likelihood/gradient regression check; no optimizer or fixtures.
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'MDCEV'),fullfile(repo,'MMDCEV'));

md.NAlt = 3; md.NCT = 2; md.NP = 2;
md.NVarA = 2; md.NVarM = 1; md.NVarU = 1; md.NVarP = 3;
md.Profile = 1; md.SpecProfile = [1 2 3;0 0 0]; md.RealMin = 1;
n = md.NCT*md.NP;
alt = repmat((1:md.NAlt)',n,1);
data.Xa = [alt == 3,linspace(-1,1,md.NAlt*n)'];
data.Y = [1 0 .5 0;0 1 .5 0;1 1 1 1];
data.priceMat = [1 1.1 .9 1.2;1.2 .9 1.1 1;1000 1000 1000 1000];
data.Xm = [0;1;0;1];
data.Xu = [ones(n,1),[-1;1;-1;1]];
b = [0.2;-0.1;0.3;-0.2;repmat([0.8;0.1],3,1);-0.2];
f0 = probs_mdcev(data,md,b);
b(3) = b(3)+0.2;
f1 = probs_mdcev(data,md,b);
assert(all(isfinite(f0)) && all(abs(f1(data.Xm == 0)-f0(data.Xm == 0)) < 1e-12));
assert(any(abs(f1(data.Xm == 1)-f0(data.Xm == 1)) > 1e-8));
rmd = struct('EstimOpt',md,'INPUT',data,'bhat',b);
rmd.EstimOpt.NSim = 2; rmd.INPUT.I = 2000*ones(n,1);
dmd = MDCEV_demand(rmd,[]);
assert(max(abs(sum(dmd.*data.priceMat,1)-rmd.INPUT.I')) < 1e-8);

mm.NAlt = 3; mm.NCT = 2; mm.NP = 3; mm.NRep = 8;
mm.NVarA = 2; mm.NVarM = 1; mm.NVarU = 0; mm.NVarP = 1;
mm.Profile = 1; mm.SpecProfile = [1 1 1;0 0 0]; mm.RealMin = 1;
mm.FullCov = 0; mm.Dist = [0 0]; mm.indx1 = []; mm.indx2 = [];
data2.Y = repmat([1 0;0 1;1 1],1,mm.NP);
data2.priceMat = repmat([1 1.1;1.2 .9;1000 1000],1,mm.NP);
data2.Xa = zeros(mm.NAlt*mm.NCT,mm.NVarA,mm.NP);
for i=1:mm.NP
    data2.Xa(:,:,i) = [zeros(mm.NAlt*mm.NCT,1),linspace(-1,1,mm.NAlt*mm.NCT)'];
    data2.Xa(3:mm.NAlt:end,1,i) = 1;
end
data2.Xm = [0 1 -1];
data2.Xu = zeros(mm.NCT*mm.NP,0);
data2.err = reshape(linspace(-1.5,1.5,mm.NVarA*mm.NRep*mm.NP),mm.NVarA,[]);
b2 = [0.2;-0.1;0.3;0.2;0.1;-0.2;0.8;-0.1];
[f,g] = probs_mmdcev(data2,mm,b2);
gn = zeros(size(g));
h = 1e-6;
for k=1:numel(b2)
    d = zeros(size(b2)); d(k) = h;
    gn(:,k) = (probs_mmdcev(data2,mm,b2+d)-probs_mmdcev(data2,mm,b2-d))/(2*h);
end
assert(all(isfinite(f)) && max(abs(g-gn),[],'all') < 2e-4);
rmm = struct('EstimOpt',mm,'INPUT',data2,'bhat',b2);
rmm.EstimOpt.NSim = 2; rmm.INPUT.I = 2000*ones(mm.NCT*mm.NP,1);
dm = MMDCEV_demand(rmm,[]);
assert(max(abs(sum(dm.*data2.priceMat,1)-rmm.INPUT.I')) < 1e-8);

mm.NP = 1; mm.NCT = 200;
data2.Y = repmat([1;0;1],1,mm.NCT);
data2.priceMat = repmat([1;1.2;1000],1,mm.NCT);
data2.Xa = repmat(data2.Xa(1:3,:,1),mm.NCT,1);
data2.Xa = reshape(data2.Xa,mm.NAlt*mm.NCT,mm.NVarA,1);
data2.Xm = 0; data2.Xu = zeros(mm.NCT,0);
data2.err = data2.err(:,1:mm.NRep);
assert(all(isfinite(probs_mmdcev(data2,mm,b2))));
disp('test_mdcev_mmdcev passed')
end
