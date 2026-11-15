function test_hessian_diagnostics(outDir)
% Hessian quality diagnostics and the Hessian settings they document.
%
%   test_hessian_diagnostics(outDir)
%
% Checks, in order:
%   1. hessianDiagnostics on matrices with known answers, including the
%      information matrix equality warning and its silence inside [0.8,1.25];
%   2. MNL on the demo data - the diagnostics are computed, stored in
%      Results.diagnostics and printed/saved next to the 'Hessian' row. MNL
%      simulates nothing, so SE(BHHH)/SE(H) has to be about 1;
%   3. MXL_d on the demo data - the same mechanics with a simulated model;
%   4. weights together with the analytical Hessian are an error, not a
%      silent reset of all weights to 1;
%   5. HessEstFix = 5 works in LC and LCMXL (it used to fail there).
%
% Standard errors and estimates are not affected by any of this; the
% diagnostics only read matrices the modules have already computed.

if nargin < 1 || isempty(outDir)
    outDir = fullfile(tempdir,'test_hessian_diagnostics');
end
if ~isfolder(outDir), mkdir(outDir); end

warnId = 'hessianDiagnostics:InformationMatrixEquality';


%% 1. the helper on matrices with known answers

n = 4;
H = diag([1 2 3 4]);
iH = inv(H);
J = sqrt(H);                       % J'*J == H, so SE(BHHH) == SE(H)
BActive = ones(n,1);

lastwarn('');
[D,TailRows] = hessianDiagnostics(H,iH,zeros(n,1),J,BActive,'unit');
[~,id] = lastwarn;
assert(abs(D.MedianSERatioBHHH - 1) < 1e-12,'ratio of equal standard errors')
assert(~D.SERatioOutsideBounds && isempty(id),'no warning at ratio 1')
assert(abs(D.MinEigHessian - 1) < 1e-12,'smallest eigenvalue')
assert(D.MaxAbsGradient == 0 && D.MaxNewtonStepSE == 0,'zero gradient')
assert(D.NNonFiniteSE == 0,'all standard errors finite')
assert(isequal(size(TailRows),[5 2]) && ischar(TailRows{1,2}),'tail rows')

lastwarn('');
D = hessianDiagnostics(H,iH,zeros(n,1),2*J,BActive,'unit');   % J'*J == 4*H
[~,id] = lastwarn;
assert(abs(D.MedianSERatioBHHH - 0.5) < 1e-12,'ratio of halved standard errors')
assert(D.SERatioOutsideBounds && strcmp(id,warnId),'warning below the lower bound')

lastwarn('');
D = hessianDiagnostics(H,iH,zeros(n,1),J/2,BActive,'unit');   % J'*J == H/4
[~,id] = lastwarn;
assert(abs(D.MedianSERatioBHHH - 2) < 1e-12 && D.SERatioOutsideBounds, ...
    'ratio of doubled standard errors')
assert(strcmp(id,warnId),'warning above the upper bound')

% a non-unit gradient: the Newton step is reported in standard errors
g = [1;0;0;0];
D = hessianDiagnostics(H,iH,g,J,BActive,'unit');
assert(abs(D.MaxNewtonStepSE - 1) < 1e-12,'Newton step in standard errors')

% a fixed parameter is excluded, a missing jacobian is reported, not faked
D = hessianDiagnostics(H(1:3,1:3),iH,g,[],[1 1 1 0],'unit');
assert(isnan(D.MedianSERatioBHHH) && D.NFree == 3,'restricted diagnostics')
D = hessianDiagnostics(H,iH,g,'reason text',BActive,'unit');
assert(isnan(D.MedianSERatioBHHH) && strcmp(D.Note,'reason text'),'skip reason')

% a singular Hessian is reported instead of throwing
D = hessianDiagnostics(zeros(n),zeros(n),g,J,BActive,'unit');
assert(D.MinEigHessian == 0 && D.NNonFiniteSE == n,'singular Hessian')
fprintf('HESS_DIAG_UNIT_PASS\n');


%% data for the module checks (the toolbox demo data set)

demo = fullfile(fileparts(fileparts(mfilename('fullpath'))),'DCE_demo','NEWFOREX_DCE_demo1.mat');
assert(exist(demo,'file') == 2,'demo data not found: %s',demo);
DATA = load(demo);

INPUT = struct();
INPUT.Y = DATA.Choice;
INPUT.Xa = [DATA.SQ,DATA.GOS,DATA.CEN,DATA.VIS == 2,DATA.VIS == 1,-DATA.FEE/4/10];

EstimOpt = struct();
EstimOpt.NamesA = {'Status quo';'Passive protection of commertial forests'; ...
    'Passive protection of natural regeneration forests';'Visitors 7500 per day'; ...
    'Visitors 5000 per day';'-Cost (10 EUR)'};
EstimOpt.NCT = 12;
EstimOpt.NAlt = 3;
EstimOpt.NP = numel(INPUT.Y)/EstimOpt.NCT/EstimOpt.NAlt;
EstimOpt.OutputDir = outDir;              % before DataCleanDCE: it starts the diary
EstimOpt.ProjectName = 'hess diagnostics test';

[INPUT,Results,EstimOpt,OptimOpt] = DataCleanDCE(INPUT,EstimOpt);
EstimOpt.NumGrad = 0;
EstimOpt.RobustStd = 0;
OptimOpt.Algorithm = 'quasi-newton';
OptimOpt.GradObj = 'on';
OptimOpt.Hessian = 'off';
OptimOpt.Display = 'off';


%% 2. MNL - nothing is simulated, so the two matrices have to agree

E = EstimOpt; E.HessEstFix = 5;
lastwarn('');
Results.MNL = MNL(INPUT,Results,E,OptimOpt);
[~,id] = lastwarn;
D = Results.MNL.diagnostics;
assert(isfinite(D.MedianSERatioBHHH) && isfinite(D.MinEigHessian) && ...
    isfinite(D.MaxAbsGradient) && isfinite(D.MaxNewtonStepSE),'MNL diagnostics finite')
assert(D.MedianSERatioBHHH > 0.8 && D.MedianSERatioBHHH < 1.25, ...
    'MNL median SE(BHHH)/SE(H) = %.4g outside [0.8, 1.25]',D.MedianSERatioBHHH)
assert(~D.SERatioOutsideBounds && ~strcmp(id,warnId),'MNL raised the simulator warning')
assert(D.MinEigHessian > 0 && D.NNonFiniteSE == 0,'MNL Hessian at a maximum')
assert(any(strcmp(Results.MNL.R_out(:,1),'Median SE(BHHH)/SE(H)')), ...
    'diagnostics missing from the MNL output')
fprintf('HESS_DIAG_MNL_PASS ratio %.4g, min eig %.4g, max|g| %.3g, Newton %.3g SE\n', ...
    D.MedianSERatioBHHH,D.MinEigHessian,D.MaxAbsGradient,D.MaxNewtonStepSE);


%% 3. MXL_d - the same mechanics with a simulated likelihood

E = EstimOpt;
E.HessEstFix = 5;
E.FullCov = 0;
E.NRep = 500;
E.Dist = zeros(1,EstimOpt.NVarA);
Results.MXL_d = MXL(INPUT,Results,E,OptimOpt);
D = Results.MXL_d.diagnostics;
assert(isfinite(D.MedianSERatioBHHH),'MXL_d median standard error ratio')
assert(isfinite(D.MinEigHessian) && isfinite(D.MaxAbsGradient),'MXL_d diagnostics finite')
assert(any(strcmp(Results.MXL_d.R_out(:,1),'Median SE(BHHH)/SE(H)')), ...
    'diagnostics missing from the MXL output')
fprintf('HESS_DIAG_MXL_PASS ratio %.4g, min eig %.4g, max|g| %.3g, Newton %.3g SE\n', ...
    D.MedianSERatioBHHH,D.MinEigHessian,D.MaxAbsGradient,D.MaxNewtonStepSE);


%% 4. weights with the analytical Hessian are an error, not a silent reset

Iw = INPUT;
Iw.W = 1 + mod((1:EstimOpt.NP)',2)/2;      % weights that are not all 1
Iw.W = EstimOpt.NP*Iw.W/sum(Iw.W);
E = EstimOpt; E.HessEstFix = 4; E.Display = 0;
try
    MXL(Iw,Results,E,OptimOpt);
    error('weights with HessEstFix = 4 did not raise an error');
catch ME
    assert(strcmp(ME.identifier,'MXL:WeightsWithAnalyticalHessian'), ...
        'unexpected error: %s (%s)',ME.message,ME.identifier)
end
fprintf('HESS_WEIGHTS_ERROR_PASS\n');


%% 5. HessEstFix = 5 in LC and LCMXL

O = OptimOpt;
E = EstimOpt;
E.HessEstFix = 5;
E.NClass = 2;
E.Display = 0;
Rlc = Results;
Rlc.LC.b0 = [1.1931370757907798;0.2815764715419146;0.5374717296012220; ...
    -0.3054383941964302;-0.1831477362870355;1.1001689372825088; ...
    -1.4131343817410376;0.5570029252510648;0.7351411877222224; ...
    0.3272949568749655;0.1471108267204569;0.6886643341302896;-0.4277287840231599];
R = LC(INPUT,Rlc,E,O);
assert(all(isfinite(R.hess(:))) && all(isfinite(R.std)),'LC HessEstFix = 5')
assert(min(real(eig((R.hess + R.hess')/2))) > 0,'LC Hessian not positive definite')
fprintf('HESS_LC_5_PASS min eig %.4g\n',min(real(eig((R.hess + R.hess')/2))));

E.NRep = 200;
E.FullCov = 0;
E.Dist = zeros(1,EstimOpt.NVarA);
Rlcmxl = Results;
Rlcmxl.LCMXL_d.b0 = [0.0380807426830999;0.2281848732223802;0.3174276128610416; ...
    -0.1256544028697331;-0.1623555534496395;-0.3226970208204383; ...
    -3.8239003331275563;1.4181961167717461;1.9364405524738382; ...
    0.9761073108398031;0.6015870106916039;0.3925777507002695; ...
    2.0568840041824745;0.4156842096340016;-0.3252185735872216; ...
    0.3949601108592588;0.3354961926182393;1.3919258768556675; ...
    2.1335591515587713;1.2754360293186153;1.2419153013633313; ...
    -0.8653219958227978;-0.1212910287944660;1.0419217222020118;0.4290906748324989];
R = LCMXL(INPUT,Rlcmxl,E,O);
assert(all(isfinite(R.hess(:))) && all(isfinite(R.std)),'LCMXL HessEstFix = 5')
assert(min(real(eig((R.hess + R.hess')/2))) > 0,'LCMXL Hessian not positive definite')
fprintf('HESS_LCMXL_5_PASS min eig %.4g\n',min(real(eig((R.hess + R.hess')/2))));


%% 6. with robust standard errors the diagnostics say which covariance they describe

E = EstimOpt; E.HessEstFix = 5; E.RobustStd = 1; E.Display = 0;
R = MNL(INPUT,Results,E,OptimOpt);
assert(isfinite(R.diagnostics.MedianSERatioBHHH),'robust run: diagnostics')
assert(contains(R.diagnostics.Note,'model-based'),'robust run: note')
fprintf('HESS_ROBUST_NOTE_PASS ratio %.4g\n',R.diagnostics.MedianSERatioBHHH);


%% 7. the compatibility overrides keep the estimator the user asked for

% WTP space has no analytical Hessian for the optimizer, so ApproxHess has to
% go; HessEstFix is the ex-post estimator of the covariance and must survive
O = OptimOpt; O.MaxIterations = 1;
E = EstimOpt;
E.WTP_space = 1;
E.FullCov = 0;
E.NRep = 100;
E.Dist = zeros(1,EstimOpt.NVarA);
E.ApproxHess = 0;
E.Display = 0;
for asked = [0 3 5]
    E.HessEstFix = asked;
    R = MXL(INPUT,Results,E,O);
    assert(R.EstimOpt.ApproxHess == 1,'ApproxHess not turned off')
    assert(R.EstimOpt.HessEstFix == asked, ...
        'HessEstFix changed from %d to %d',asked,R.EstimOpt.HessEstFix)
end
% 4 is the one value that cannot be kept - it is the missing analytical Hessian
E.HessEstFix = 4;
R = MXL(INPUT,Results,E,O);
assert(R.EstimOpt.ApproxHess == 1 && R.EstimOpt.HessEstFix == 5, ...
    'HessEstFix = 4 became %d, expected 5',R.EstimOpt.HessEstFix)
fprintf('HESS_OVERRIDE_KEEPS_SETTING_PASS\n');


%% 8. the demand models report the diagnostics too

settings = parallel.Settings; previous = settings.Pool.AutoCreate;
restore = onCleanup(@() setAuto(settings,previous)); %#ok<NASGU>
settings.Pool.AutoCreate = false;
for model = {'MDCEV','MMDCEV'}
    name = model{1};
    % 60 respondents: BHHH needs more observations than parameters, and the
    % MMDCEV jacobian has one row per respondent
    c = vce_fixture(name,3,1,1,1,60,16,3,3); e = c.opt;
    I = c.data; I.Y = I.Y(:); I.priceMat = I.priceMat(:);
    if strcmp(name,'MMDCEV')
        I.Xa = reshape(permute(I.Xa,[1 3 2]),e.NAlt*e.NCT*e.NP,e.NVarA);
        I.Xm = repelem(I.Xm',e.NAlt*e.NCT,1); I = rmfield(I,'err');
    else
        I.Xm = repelem(I.Xm,e.NAlt,1);
    end
    I.Xu = repelem(I.Xu(:,2:end),e.NAlt,1);
    e.Display = 0; e.ConstVarActive = 0; e.ApproxHess = 1; e.NSdSim = 100;
    e.Seed1 = 123; e.Seed2 = 456; e.Draws = 1; e.NObs = e.NCT*e.NP;
    e.ProjectName = ['hess diag ',name]; e.SaveTxtOutput = 0; e.OutputDir = outDir;
    e.HessEstFix = 5;
    old = struct(name,struct('b0',c.b));
    oo = optimset('Algorithm','quasi-newton','GradObj','on','Hessian','off', ...
        'MaxIter',1,'Display','off','FinDiffType','central');
    R = feval(name,I,old,e,oo);
    D = R.diagnostics;
    % one iteration only, so the point is not a maximum: what is checked here
    % is that the diagnostics are computed and that an unavailable ratio says
    % why instead of passing for a healthy 1
    assert(isfinite(D.MinEigHessian) && isfinite(D.MaxAbsGradient),'%s: diagnostics finite',name)
    assert(isfinite(D.MedianSERatioBHHH) || ~isempty(D.Note), ...
        '%s: an unavailable standard error ratio has to state a reason',name)
    fprintf('HESS_DIAG_%s_PASS ratio %.4g, min eig %.4g, max|g| %.4g %s\n', ...
        name,D.MedianSERatioBHHH,D.MinEigHessian,D.MaxAbsGradient,D.Note);
end


%% 9. every module is wired to the diagnostics

repo = fileparts(fileparts(mfilename('fullpath')));
modules = {'MNL/MNL.m','MXL/MXL.m','LC/LC.m','LCMXL/LCMXL.m','HMNL/HMNL.m','HMXL/HMXL.m', ...
    'GMXL/GMXL.m','HLC/HLC.m','MIMIC/MIMIC.m','LML/LML.m','MDCEV/MDCEV.m','MMDCEV/MMDCEV.m'};
for i = 1:numel(modules)
    src = fileread(fullfile(repo,modules{i}));
    assert(contains(src,'dceHessDiagnostics('),'%s does not compute the diagnostics',modules{i})
    assert(contains(src,'Tail = [Tail; DiagnosticsTail]'),'%s does not print them',modules{i})
end
fprintf('HESS_DIAG_WIRING_PASS %d modules\n',numel(modules));

fprintf('ALL_HESS_DIAGNOSTICS_TESTS_PASS\n');
end

function setAuto(settings,value)
settings.Pool.AutoCreate = value;
end
