function test_sim_diagnostics(outDir)
% Simulator diagnostics: effective draws per respondent and the draw-set check.
%
%   test_sim_diagnostics(outDir)
%
% Checks, in order:
%   1. simDiagnostics on inputs with known answers, including both warnings
%      and their silence at the thresholds;
%   2. mxl_drawESS integrates the same thing LL_mxl does - the log-likelihood
%      it recomputes has to equal the module's own, for a diagonal and a full
%      covariance, in preference space and in WTP space with a lognormal cost
%      coefficient;
%   3. the error paths of mxl_drawESS, and that an unsupported specification
%      leaves the model estimated and the rest of the diagnostics in place;
%   4. MXL on the demo data - both parts reach Results.diagnostics.Simulator
%      and the printed/saved output;
%   5. every module with draws computes the draw-set check.
%
% Estimates and standard errors are not affected by any of this.

if nargin < 1 || isempty(outDir)
    outDir = fullfile(tempdir,'test_sim_diagnostics');
end
if ~isfolder(outDir), mkdir(outDir); end

starvedId = 'simDiagnostics:StarvedDraws';
shiftId = 'simDiagnostics:DrawSetSpecific';


%% 1. the helper on inputs with known answers

rotOK = struct('MaxAbsGradient',1e-4,'MaxNewtonStepSE',0.01, ...
    'MedianNewtonStepSE',0.002,'Offset',50,'ReferenceMaxAbsGradient',1e-6);

lastwarn('');
[D,TailRows] = simDiagnostics(100*ones(40,1),1000,rotOK,'unit');
[~,id] = lastwarn;
assert(D.MedianESS == 100 && abs(D.MedianESSShare - 0.1) < 1e-12,'median effective draws')
assert(D.MinESS == 100 && D.NBelowMinDraws == 0 && D.ShareBelowMinDraws == 0,'nothing starved')
assert(~D.StarvedDraws && ~D.DrawSetSpecific && isempty(id),'no warning on a healthy simulator')
assert(isequal(size(TailRows),[4 2]) && ischar(TailRows{1,2}),'tail rows')
assert(D.RotationOffset == 50 && D.RotatedStepSE == 0.01,'rotation statistics stored')

% exactly 5% of respondents below 10 draws is not yet a warning, 10% is
ess = 100*ones(40,1);
ess(1:2) = 3;                                  % 5% of 40
lastwarn('');
D = simDiagnostics(ess,1000,rotOK,'unit');
[~,id] = lastwarn;
assert(D.ShareBelowMinDraws == 0.05 && ~D.StarvedDraws && isempty(id), ...
    'no warning at the threshold share')
ess(3:4) = 3;
lastwarn('');
D = simDiagnostics(ess,1000,rotOK,'unit');
[~,id] = lastwarn;
assert(D.StarvedDraws && strcmp(id,starvedId),'warning above the threshold share')
assert(D.NBelowMinDraws == 4 && D.MinESS == 3,'starved respondents counted')

% the draw-set check warns when rotating the blocks would move the estimate
rotBad = rotOK;
rotBad.MaxNewtonStepSE = 1.01;
rotBad.MedianNewtonStepSE = 0.3;
rotBad.MaxAbsGradient = 53;
lastwarn('');
D = simDiagnostics(100*ones(40,1),1000,rotBad,'unit');
[~,id] = lastwarn;
assert(D.DrawSetSpecific && strcmp(id,shiftId),'warning above the shift threshold')
lastwarn('');
rotEdge = rotOK; rotEdge.MaxNewtonStepSE = 1;
D = simDiagnostics(100*ones(40,1),1000,rotEdge,'unit');
[~,id] = lastwarn;
assert(~D.DrawSetSpecific && isempty(id),'no warning at the shift threshold')

% thresholds can be changed
lastwarn('');
D = simDiagnostics(ess,1000,rotOK,'unit',[2 0.05 1]);
[~,id] = lastwarn;
assert(~D.StarvedDraws && isempty(id),'the starvation threshold is settable')

% unavailable parts are reported, not invented
[D,TailRows] = simDiagnostics('no analytical jacobian',1000,'numerical gradient','unit');
assert(isnan(D.MedianESS) && isnan(D.RotatedStepSE),'nothing invented')
assert(strcmp(D.ESSNote,'no analytical jacobian') && ...
    strcmp(D.RotationNote,'numerical gradient'),'reasons kept')
assert(isequal(size(TailRows),[2 2]) && contains(TailRows{1,2},'not available'), ...
    'reasons reach the output')
fprintf('SIM_DIAG_UNIT_PASS\n');


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
EstimOpt.ProjectName = 'sim diagnostics test';

[INPUT,Results,EstimOpt,OptimOpt] = DataCleanDCE(INPUT,EstimOpt);
EstimOpt.NumGrad = 0;
EstimOpt.RobustStd = 0;
EstimOpt.HessEstFix = 5;
EstimOpt.NRep = 500;
EstimOpt.Drawskeep = 1;
OptimOpt.Algorithm = 'quasi-newton';
OptimOpt.GradObj = 'on';
OptimOpt.Hessian = 'off';
OptimOpt.Display = 'off';

NVarA = size(INPUT.Xa,2);
XXa = permute(reshape(INPUT.Xa,[EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP,NVarA]),[1 3 2]);
I2 = INPUT;
I2.XXa = XXa;
I2.YY = reshape(INPUT.Y,[EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP]);


%% 2. mxl_drawESS integrates the same thing LL_mxl does

cases = struct( ...
    'name',   {'diagonal, preference space','full covariance','WTP space, lognormal cost'}, ...
    'fullcov',{0,1,0}, ...
    'wtp',    {0,0,1}, ...
    'lognorm',{false,false,true});

for c = 1:numel(cases)
    E = EstimOpt;
    E.FullCov = cases(c).fullcov;
    E.WTP_space = cases(c).wtp;
    E.Dist = zeros(1,NVarA);
    if cases(c).lognorm
        E.Dist(end) = 1;      % lognormal cost coefficient (the attribute is -FEE,
    end                       % so its coefficient is positive)
    INPUTc = INPUT;
    I2c = I2;
    E.ProjectName = ['sim diagnostics ',cases(c).name];
    R = MXL(INPUTc,Results,E,OptimOpt);
    Eu = E;
    Eu.NVarA = NVarA;
    Eu.ExpB = [];
    [ess,ll] = mxl_drawESS(I2c,Eu,R.err,R.bhat);
    assert(abs(ll - R.LL) < 1e-6*max(1,abs(R.LL)), ...
        '%s: mxl_drawESS log-likelihood %.10f against the module''s %.10f', ...
        cases(c).name,ll,R.LL);
    assert(all(ess > 0 & ess <= E.NRep + 1e-9),'%s: effective draws out of range',cases(c).name)
    S = R.diagnostics.Simulator;
    assert(abs(median(ess) - S.MedianESS) < 1e-9,'%s: the module reports other numbers',cases(c).name)
    fprintf('SIM_DIAG_ESS_PASS %-28s LL %.6f, effective draws median %.1f of %d, min %.1f\n', ...
        cases(c).name,ll,S.MedianESS,E.NRep,S.MinESS);
end


%% 3. the error paths, and an unsupported specification

Eu = EstimOpt; Eu.NVarA = NVarA; Eu.FullCov = 0; Eu.WTP_space = 0;
Eu.Dist = zeros(1,NVarA); Eu.Dist(1) = 3;                 % triangular
try
    mxl_drawESS(I2,Eu,zeros(NVarA,EstimOpt.NRep*EstimOpt.NP),zeros(2*NVarA,1));
    error('a triangular distribution did not raise an error');
catch ME
    assert(strcmp(ME.identifier,'mxl_drawESS:Distribution'),'unexpected error: %s',ME.message)
end
Eu.Dist = zeros(1,NVarA);
try
    mxl_drawESS(I2,Eu,zeros(NVarA,7),zeros(2*NVarA,1));
    error('a wrong draw matrix did not raise an error');
catch ME
    assert(strcmp(ME.identifier,'mxl_drawESS:Draws'),'unexpected error: %s',ME.message)
end

% covariates of the means: the effective draws are not available, the model is
% estimated and the draw-set check still runs
E = EstimOpt;
E.FullCov = 0;
E.Dist = zeros(1,NVarA);
E.NVarM = 1;
E.ProjectName = 'sim diagnostics covariates';
Im = INPUT;
Im.Xm = DATA.AGE(:)/100;
R = MXL(Im,Results,E,OptimOpt);
S = R.diagnostics.Simulator;
assert(isnan(S.MedianESS) && ~isempty(S.ESSNote),'unsupported specification not reported')
assert(isfinite(S.RotatedStepSE),'the draw-set check did not run with covariates')
assert(isfinite(R.LL) && all(isfinite(R.std(R.std > 0))),'the model did not estimate')
fprintf('SIM_DIAG_FALLBACK_PASS "%s", draw-set shift %.3g SE\n',S.ESSNote,S.RotatedStepSE);


%% 4. the diagnostics reach the output

E = EstimOpt;
E.FullCov = 0;
E.Dist = zeros(1,NVarA);
E.ProjectName = 'sim diagnostics output';
R = MXL(INPUT,Results,E,OptimOpt);
S = R.diagnostics.Simulator;
assert(isfinite(S.MedianESS) && isfinite(S.MinESS) && isfinite(S.RotatedStepSE) && ...
    isfinite(S.RotatedMedianStepSE) && isfinite(S.RotatedMaxAbsGradient), ...
    'simulator diagnostics finite')
assert(S.RotatedMedianStepSE <= S.RotatedStepSE,'the median shift exceeds the maximum')
assert(S.RotationOffset == floor(EstimOpt.NP/2),'rotation offset')
assert(any(strcmp(R.R_out(:,1),'Effective draws')),'effective draws missing from the output')
assert(any(strcmp(R.R_out(:,1),'Shift, other draws')),'draw-set check missing from the output')
assert(isfield(R.diagnostics,'MedianSERatioBHHH'),'the Hessian diagnostics are still there')
fprintf('SIM_DIAG_OUTPUT_PASS effective draws median %.1f of %d, min %.1f, <10 draws %.1f%%, shift %.3g SE (max|g| %.3g against %.3g)\n', ...
    S.MedianESS,E.NRep,S.MinESS,100*S.ShareBelowMinDraws,S.RotatedStepSE, ...
    S.RotatedMaxAbsGradient,S.ReferenceMaxAbsGradient);


%% 5. every module with draws computes the check

root = fileparts(fileparts(mfilename('fullpath')));
mods = {'MXL','GMXL','LCMXL','HMNL','HMXL','HLC','MIMIC'};
for k = 1:numel(mods)
    f = fileread(fullfile(root,mods{k},[mods{k},'.m']));
    assert(contains(f,'dceSimDiagnostics('),'%s does not compute the simulator diagnostics',mods{k})
    assert(contains(f,'Tail = [Tail; SimDiagnosticsTail]'), ...
        '%s does not put them in the output',mods{k})
end
fprintf('SIM_DIAG_WIRING_PASS %d modules\n',numel(mods));

fprintf('\nALL_SIM_DIAGNOSTICS_TESTS_PASS\n');
end
