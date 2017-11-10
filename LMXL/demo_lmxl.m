% Different no. of alternatives DEMO
% DIFFERENT NO. FOR THE SAME INDIVIDUAL
% MXL simulation
% 2016-12-28

clear
clc

global B_backup
%% EstimOpt
n = 1.5e3; % no of respondents
rng(1000001);
EstimOpt.NP = n;
EstimOpt.NAlt = 4;
EstimOpt.NCT = 6;
EstimOpt.NRep =500; %1e3; % no of draws
EstimOpt.Draws = 6;
EstimOpt.HaltonSkip = 1;
EstimOpt.HaltonLeap = 0;
EstimOpt.HessEstFix = 1;
EstimOpt.Display = 1;
EstimOpt.ConstVarActive = 0;
EstimOpt.NObs = EstimOpt.NP*EstimOpt.NCT;
EstimOpt.NVarA = 4;
EstimOpt.NSDSIM = 10000;

OptimOpt = optimoptions('fminunc');
%OptimOpt.Algorithm = 'quasi-newton';
OptimOpt.GradObj = 'on';
% OptimOpt.Hessian = 'user-supplied'; % ('off'), only used by trust-region
OptimOpt.MaxIter = 1000;
OptimOpt.FunValCheck = 'on';
OptimOpt.Diagnostics = 'on';
OptimOpt.MaxFunEvals = 1e6; %Maximum number of function evaluations allowed (1000)
OptimOpt.OutputFcn = @outputf;
OptimOpt.TolFun = 1.e-7;
OptimOpt.TolX = 1.e-7;

EstimOpt.NamesA = {'SQ', 'Var1', 'Var2', 'Var3'};

EstimOpt.WTP_space = 0;

%% Generating data for MXL
Lognormal = 0;
EstimOpt.FullCov =0;

ba = [5, 4, 3, -2]';
bvar = [3,3,3,1];

if EstimOpt.FullCov == 0
    u = normrnd(0,3,3, EstimOpt.NP);
    u = [u; normrnd(0,1,1, EstimOpt.NP)];
else
    Sigma = [9, 3, 3, 2; 3, 9, 2, 2; 3, 2, 9, 2; 2, 2, 2, 1];
    u = mvnrnd(zeros(1,4), Sigma, EstimOpt.NP)';
end

bmxl = ba(:, ones(EstimOpt.NP,1)) +u;


if Lognormal > 0
    bmxl(4,:) =  exp(bmxl(4,:));
end

INPUT.Xa = 2*unifrnd(0,1,EstimOpt.NAlt*EstimOpt.NP*EstimOpt.NCT,4);
INPUT.Xa(:,4) = -INPUT.Xa(:,4);

epsmnl = unifrnd(0,1,EstimOpt.NP*EstimOpt.NCT*EstimOpt.NAlt,1);
epsmnl = -log(-log(epsmnl)); % Gumbel distribution

U = zeros(EstimOpt.NAlt*EstimOpt.NCT*EstimOpt.NP,1);
for i = 1:EstimOpt.NP
    U((i-1)*EstimOpt.NAlt*EstimOpt.NCT+1:i*EstimOpt.NAlt*EstimOpt.NCT) = INPUT.Xa((i-1)*EstimOpt.NAlt*EstimOpt.NCT+1:i*EstimOpt.NAlt*EstimOpt.NCT,:)*bmxl(:,i) ; %NALT*CT x 1
end
U = U + epsmnl;
U =  reshape(U, EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP);
INPUT.Y = zeros(EstimOpt.NAlt,EstimOpt.NP*EstimOpt.NCT);


for i = 1:EstimOpt.NP*EstimOpt.NCT
    maxU = max(U(:,i));
    k = find(U(:,i) == maxU);
    INPUT.Y(k,i) = 1;
end

INPUT.Y = reshape(INPUT.Y, EstimOpt.NP*EstimOpt.NCT*EstimOpt.NAlt,1);
INPUT.MissingInd = zeros(EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP);


%% Some transformations and estimation

EstimOpt.NumGrad = 1;
%EstimOpt.FullCov =0;
[INPUT, Results, EstimOpt, OptimOpt] = DataCleanDCE(INPUT,EstimOpt);
EstimOpt.Bounds = [round(min(bmxl,[],2)), round(max(bmxl,[],2))];
if Lognormal == 0
    EstimOpt.Dist = [0,0,0,0];
else
    EstimOpt.Dist = [0,0,0,1];
end

if EstimOpt.FullCov == 1
    EstimOpt.FullCov =0;
    Results.LMXL_d = LMXL(INPUT,Results,EstimOpt,OptimOpt);
    EstimOpt.FullCov =1;
    Results.LMXL = LMXL(INPUT,Results,EstimOpt,OptimOpt);
else
    Results.LMXL = LMXL(INPUT,Results,EstimOpt,OptimOpt);
end
B_backup = [ba; bvar'];
if EstimOpt.FullCov == 1
    EstimOpt.FullCov =0;
    Results.MXL_d = MXL(INPUT,Results,EstimOpt,OptimOpt);
    EstimOpt.FullCov =1;
    Results.MXL = MXL(INPUT,Results,EstimOpt,OptimOpt);
else
    Results.MXL = MXL(INPUT,Results,EstimOpt,OptimOpt);
end
Res = [mean(bmxl,2); std(bmxl,[],2)];
Res = [Res, [Results.LMXL.Means; Results.LMXL.Stds], [Results.MXL.DetailsA(:,1); Results.MXL.DetailsV(:,1)]]

figure
for i = 1:4
    subplot(4,2,(i-1)*2 + 1) 
    hist(bmxl(i,:)',100)
    subplot(4,2,(i-1)*2 + 2) 
    plot(Results.LMXL.B2_sort(i,:), Results.LMXL.P2_sort(i,:))
end