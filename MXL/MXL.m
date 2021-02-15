function Results = MXL(INPUT,Results_old,EstimOpt,OptimOpt)

% MXL creates Mixed Logit Model that requires from the user to specify the distribution of each parameter.
%
% Syntax:   MXL(INPUT,EstimOpt,OptimOpt)
%           MXL(INPUT,Results_old,EstimOpt,OptimOpt)
%
% Inputs:
%    INPUT - clean, updated INPUT data from DataCleanDCE
%    EstimOpt - Estimation Options (check below)
%    OptimOpt - Optimizer Options define how algorithms converges to the final result. They are set by default based on provided EstimOpt in DataCleanDCE, however, they are subject to change.
%    Results_old - here one can provide old results to use as starting
%    values
%
% EstimOpt Options:
% Set them by e.g. Estimopt.DataFile = 'Project'
%
% MXL parameter options:
% �	Dist = 0; distribution of random parameters, by default set to normal. Set in a vector of numbers, each corresponding to specific distribution:
% o	-1 - constant, 
% o	0 - normal, 
% o	1 - lognormal, 
% o	2 - spike, 
% o	3 - Triangular, 
% o	4 - Weibull, 
% o	5 - Sinh-Arcsinh, 
% o	6 - Johnson Sb, 
% o	7 - Johnson Su
% �	FullCov = 0; set to 1 for correlated random parameters, 0 if not
% �	EffectiveMoments = 0; set to 1 to calculate effective moments
% 
% 
% General basics:
% �	DataFile � path/name of the .mat data file
% �	Display � 1; shows output, set to 0 to hide it 
% �	ProjectName � Name of the project/model
% �	WTP_space � set to 1 for estimation in WTP space. If missing or set to 0, MNL uses Preference Space
% �	NCT - Number of choice tasks per person 
% �	NAlt - Number of alternatives
% �	NP � Number of respondents
% 
% 
% Variables options:
% �	NamesA � Names of variables in list e.g. {'-Opt out';�-Cost (EUR)'}
% �	NamesM � Names of variables of means of random parameters
% �	NamesS � Names of variables of Scale
% �	NLTVariables � vector specifying which attributes are to subject to non-linear transformations
% �	NLTType � Transformation for non-linear variables. By default it is set to Box-Cox transformation (1), set to 2 in order to use Yeo-Johnson transformation
% 
% Numbers of variables are set automatically; you can check them in the following fields:
% o	NVarA - Number of attributes
% o	NVarM - Number of covariates of means of random parameters
% o	NVarS - Number of covariates of scale
% o	NVarNLT - Number of non-linear variables to transform
% 
% 
% Parameters options:
% �	ExpB = vector of 0; for each parameter set it to 1 to use ExpB, otherwise 0
% �	BActive = vector of 0; for each parameter set it to 1 to constrain model parameters to their initial values
% �	ConstVarActive = 0; set to 1 to constrain model parameters to its initial values 
% 
% 
% Modelling options from DataCleanDCE:
% �	ApproxHess = 1; for user supplied hessians, 1 for BHHH, 0 for analytical
% �	RobustStd = 0; by default not using robust standard errors, set to 1 to use them
% �	NumGrad = 0; uses analytical gradient in calculations, set to 1 for numerical gradient
% �	HessEstFix = 0; Options: 
% o	0 - use optimization Hessian, 
% o	1 - use jacobian-based (BHHH) Hessian, 
% o	2 - use high-precision jacobian-based (BHHH) Hessian,
% o	3 - use numerical Hessian, 
% o	4 - use analytical Hessian
% 
% 
% For drawing and simulations:
% �	HaltonSkip = 1; specify no of rows in halton sequence to skip
% �	HaltonLeap = 0; specify no of rows in halton sequence to leap
% �	Draws = 6; specify draws type, by default Sobol with scrambling. Options: 
% o	1 - pseudo-random, 
% o	2 - Latin Hypercube, 
% o	3 - Halton, 
% o	4 - Halton RR scrambled, 
% o	5 - Sobol, 
% o	6 - Sobol MAO scrambled
% �	NRep = 1e3; specify no. of draws for numerical simulation
% �	RealMin = by default 0, can be set to 1
% �	NSdSim = 1e4; number of draws for simulating standard deviations
%  
% 
% Precision:
% �	eps = 1.e-6; overall precision level
% �	Otherwise:
% o	FunctionTolerance - df / gradient precision level
% o	TolX - step precision level
% o	OptimalityTolerance - dB precision level
% 
% 
% Seeds by default:
% �	Seed1 = 179424673
% �	Seed2 = 7521436817
%
%
% Example: 
%    Results.MXL = MXL(INPUT,Results,EstimOpt,OptimOpt);
%
% Author: Mikolaj Czajkowski, Professor
% University of Warsaw, Faculty of Economic Sciences
% email address: mik@czaj.org 
% Website: http://czaj.org/#

% save tmp_MXL
% return

global B_backup

tic

Results.bhat = [];
Results.R = [];
Results.R_out = {};
Results.stats = [];

NVarA = EstimOpt.NVarA;
NSdSim = EstimOpt.NSdSim;


%% Check data and inputs


if nargin < 3 % check no. of inputs
    error('Too few input arguments for MXL(INPUT,EstimOpt,OptimOpt)')
end

disp(' ');
disp('__________________________________________________________________________________________________________________');
disp(' ');

warning off MATLAB:mir_warning_maybe_uninitialized_temporary

format shortG;
format compact;

if any(INPUT.W ~= 1)
    cprintf('Black','Estimating '); cprintf('*Black','weighted '); cprintf('Black','MXL model...\n');
else
    disp('Estimating MXL model ...')
end

if isfield(EstimOpt,'FullCov') == 0
    EstimOpt.FullCov = 0;
end
if ~isfield(EstimOpt,'WTP_space')
    EstimOpt.WTP_space = 0;
    EstimOpt.WTP_matrix = [];
elseif EstimOpt.WTP_space == 0
    EstimOpt.WTP_matrix = [];
end

if EstimOpt.FullCov == 0
    disp('with non-correlated random parameters ...')
    if EstimOpt.WTP_space > 0
        disp('in WTP-space ...')
    else
        disp('in preference-space ...')
    end
else
    disp('with correlated random parameters ...')
    if EstimOpt.WTP_space > 0
        disp('in WTP-space ...')
    else
        disp('in preference-space ...')
    end
end

if isfield(EstimOpt,'NLTVariables') && ~isempty(EstimOpt.NLTVariables)
    EstimOpt.NLTVariables = EstimOpt.NLTVariables(:);
    EstimOpt.NVarNLT = length(unique(EstimOpt.NLTVariables));
    if ~ismember(unique(EstimOpt.NLTVariables),1:NVarA)
        error('Incorrect non-linear variable(s) specification')
    end
    if isfield(EstimOpt,'NLTType') == 0
        cprintf(rgb('DarkOrange'),'WARNING: Assuming Box-Cox transformation \n')
        EstimOpt.NLTType = 1;
    elseif EstimOpt.NLTType == 1
        disp('with Box-Cox transformed variable(s).')
    elseif EstimOpt.NLTType == 2
        disp('with Yeo-Johnson transformed variable(s)')
    else
        error('Incorrect transformation type')
    end
    if EstimOpt.NLTType == 1
        if any(INPUT.Xa(:,EstimOpt.NLTVariables) < 0)
            cprintf(rgb('DarkOrange'),'WARNING: Values of Box-Cox transformed variables < 0 \n')
        elseif any(INPUT.Xa(:,EstimOpt.NLTVariables) == 0) % not sure if this is stil necessary
            cprintf(rgb('DarkOrange'),'WARNING: Values of Box-Cox transformed variables including zeros shifted by 0.00001 \n') % tutaj raczej powinna byc miara relatywna, np abs(mean(Xa(:,n)))*0.00001
            for i = 1:EstimOpt.NVarNLT
                if any(INPUT.Xa(:,EstimOpt.NLTVariables(i)) == 0)
                    INPUT.Xa(:,EstimOpt.NLTVariables(i)) = INPUT.Xa(:,EstimOpt.NLTVariables(i)) + 0.00001;
                end
            end
        end
    end
else
    EstimOpt.NVarNLT = 0;
    EstimOpt.NLTVariables = [];
    EstimOpt.NLTType = [];
end

if isfield(EstimOpt,'Dist') == 0 || isempty(EstimOpt.Dist)
    EstimOpt.Dist = zeros(1,NVarA);
    if EstimOpt.WTP_space == 0
        cprintf(rgb('DarkOrange'),'WARNING: distributions for random parameters not specified - assuming normality \n')
    else
        cprintf(rgb('DarkOrange'),'WARNING: distributions for random parameters not specified - assuming normality (monetary parameters assumed log-normal) \n')
        EstimOpt.Dist(end-EstimOpt.WTP_space+1:end) = 1; % cost in WTP-space models log-normally distributed
    end
else
    if length(EstimOpt.Dist) == 1
        EstimOpt.Dist = EstimOpt.Dist.*ones(1,NVarA); % needed?
    elseif length(EstimOpt.Dist) == NVarA
        EstimOpt.Dist = EstimOpt.Dist(:)';
    else
        error('Incorrect no. of random parameters'' distributions provided')
    end
end
if isfield(EstimOpt,'Triang') == 0 || length(EstimOpt.Triang) ~= sum(EstimOpt.Dist == 3,2) % Needed only if any parameter has triangular distribution
    EstimOpt.Triang = zeros(1,sum(EstimOpt.Dist == 3,2));
elseif length(EstimOpt.Triang) == 1
    EstimOpt.Triang = EstimOpt.Triang*ones(1,sum(EstimOpt.Dist == 3,2));
else
    EstimOpt.Triang = EstimOpt.Triang(:)';
end

EstimOpt.Johnson = sum(EstimOpt.Dist >= 5);

if (sum(EstimOpt.Dist >= 3 & EstimOpt.Dist <= 5) > 0 && any(find(EstimOpt.Dist >= 3 & EstimOpt.Dist <= 5) > sum(EstimOpt.Dist >= 3 & EstimOpt.Dist <= 5))) && EstimOpt.FullCov == 1
    cprintf(rgb('DarkOrange'),'WARNING: It is recommended to put variables with random parameters with Triangular/Weibull/Sinh-Arcsinh distribution first \n')
end

disp(['Random parameters distributions: ',num2str(EstimOpt.Dist),' (-1 - constant, 0 - normal, 1 - lognormal, 2 - spike, 3 - Triangular, 4  - Weibull, 5 - Sinh-Arcsinh, 6 - Johnson Sb, 7 - Johnson Su)'])

if EstimOpt.WTP_space > 0 && sum(EstimOpt.Dist(end-EstimOpt.WTP_space+1:end)==1) > 0 && any(mean(INPUT.Xa(:,end-EstimOpt.WTP_space+1:end)) >= 0)
    cprintf(rgb('DarkOrange'),'WARNING: Cost attributes with log-normally distributed parameters should enter utility function with a ''-'' sign \n')
end

if isfield(INPUT,'Xs') == 0
    INPUT.Xs = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarS = size(INPUT.Xs,2); % Number of covariates of scale
if isfield(INPUT,'Xm') == 0
    INPUT.Xm = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarM = size(INPUT.Xm,2); % Number of covariates of means of random parameters

% This does not currently work:
% if isfield(INPUT, 'Xv') == 0
%     INPUT.Xv = zeros(size(INPUT.Y,1),0);
% end
% EstimOpt.NVarV = size(INPUT.Xv,2); % Number of covariates of variances of random parameters

if EstimOpt.WTP_space > 0
    if isfield(EstimOpt,'WTP_matrix') == 0
        WTP_att = (NVarA-EstimOpt.WTP_space)/EstimOpt.WTP_space;
        if rem(WTP_att,1) ~= 0
            error('EstimOpt.WTP_matrix associating attributes with cost parameters not provided')
        else
            if EstimOpt.WTP_space > 1
                disp(['EstimOpt.WTP_matrix associating attributes with cost parameters not provided - assuming equal shares for each of the ',num2str(EstimOpt.WTP_space),' monetary attributes'])
            end
            EstimOpt.WTP_matrix = NVarA - EstimOpt.WTP_space + kron(1:EstimOpt.WTP_space,ones(1,WTP_att));
        end
    elseif size(EstimOpt.WTP_matrix,2) ~= NVarA - EstimOpt.WTP_space
        error('Dimensions of EstimOpt.WTP_matrix not correct - for each non-monetary attribute provide no. of attribute to multiply it with')
    else
        EstimOpt.WTP_matrix = EstimOpt.WTP_matrix(:)';
    end
end

if isfield(EstimOpt,'Scores') == 0 || isempty(EstimOpt.Scores)
    EstimOpt.Scores = 0;
end

if isfield(EstimOpt,'NamesA') == 0 || isempty(EstimOpt.NamesA) || length(EstimOpt.NamesA) ~= NVarA
    EstimOpt.NamesA = (1:NVarA)';
    EstimOpt.NamesA = cellstr(num2str(EstimOpt.NamesA));
elseif size(EstimOpt.NamesA,1) ~= NVarA
    EstimOpt.NamesA = EstimOpt.NamesA';
end

if EstimOpt.NVarM > 0
    if isfield(EstimOpt,'NamesM') == 0 || isempty(EstimOpt.NamesM)|| length(EstimOpt.NamesM) ~= EstimOpt.NVarM
        EstimOpt.NamesM = (1:EstimOpt.NVarM)';
        EstimOpt.NamesM = cellstr(num2str(EstimOpt.NamesM));
    elseif size(EstimOpt.NamesM,1) ~= EstimOpt.NVarM
        EstimOpt.NamesM = EstimOpt.NamesM';
    end
end
if EstimOpt.NVarS > 0
    if isfield(EstimOpt,'NamesS') == 0 || isempty(EstimOpt.NamesS) || length(EstimOpt.NamesS) ~= EstimOpt.NVarS
        EstimOpt.NamesS = (1:EstimOpt.NVarS)';
        EstimOpt.NamesS = cellstr(num2str(EstimOpt.NamesS));
    elseif size(EstimOpt.NamesS,1) ~= EstimOpt.NVarS
        EstimOpt.NamesS = EstimOpt.NamesS';
    end
end

NVar = (~EstimOpt.FullCov)*(NVarA*2 + EstimOpt.NVarM*NVarA + EstimOpt.NVarS + EstimOpt.NVarNLT + 2*EstimOpt.Johnson)+EstimOpt.FullCov*(NVarA*(1+EstimOpt.NVarM) + sum(1:NVarA) + EstimOpt.NVarS + EstimOpt.NVarNLT + 2*EstimOpt.Johnson);
if ~isfield(EstimOpt,'ExpB')
    EstimOpt.ExpB = [];
elseif ~isempty(EstimOpt.ExpB)
    EstimOpt.ExpB = EstimOpt.ExpB(:);
    if size(EstimOpt.ExpB,1) == 1
        EstimOpt.ExpB = EstimOpt.ExpB.*ones(NVar,1);
    elseif size(EstimOpt.ExpB,1) ~= NVar
        error('Dimensions of ExpB not correct - provide ExpB indicator for each parameter in B')
    elseif any((EstimOpt.ExpB ~= 0) & (EstimOpt.ExpB ~= 1))
        error('ExpB must only include logical (0 or 1) values')
    end
    EstimOpt.ExpB = (1:NVar)'.* EstimOpt.ExpB;
    EstimOpt.ExpB(EstimOpt.ExpB == 0) = [];
end

%% Starting values
[Results_old, EstimOpt, b0] = setStartingValues(INPUT, Results_old, EstimOpt);

%% Optimization Options
[EstimOpt, OptimOpt] = setOptimizationOptions(EstimOpt, OptimOpt, b0);

%% Generate pseudo-random draws from standard normal distribution
err_mtx = generateRandomDraws(EstimOpt);

% set column value to 0 for constant parameters
err_mtx(:,EstimOpt.Dist == -1) = 0;


%% Display Options
[INPUT, EstimOpt, OptimOpt] = setDisplayOptions(INPUT, EstimOpt, OptimOpt);


%% Restructure data
INPUT = restructureInputData(INPUT,EstimOpt);


err_mtx = err_mtx';
% change err_mtx from NRep*NP x NVarA to NP*NRep x NVarA (increasing the no. of draws only adds new draws for each respondent, does not change all draws per individual)
% err_mtx = reshape(permute(reshape(err_mtx,EstimOpt.NP,EstimOpt.NRep,NVarA),[2,1,3]),EstimOpt.NP*EstimOpt.NRep,NVarA)';
% problem - look at the first NRep draws for NVarA=1... all are positive...
if isfield(EstimOpt,'Drawskeep') && ~isempty(EstimOpt.Drawskeep) && EstimOpt.Drawskeep == 1
    Results.err = err_mtx;
end

% Number of var-cov matrix parameters
VC = tril(ones(NVarA));
VC(VC == 1) = (1:(NVarA*(NVarA-1)/2+NVarA))';
EstimOpt.DiagIndex = diag(VC); % numbers on diagonal

% Creating indices for analytical gradient
EstimOpt.indx1 = [];
EstimOpt.indx2 = [];
if EstimOpt.NumGrad == 0 && EstimOpt.FullCov == 1
    for i = 1:NVarA
        EstimOpt.indx1 = [EstimOpt.indx1,i:NVarA];
        EstimOpt.indx2 = [EstimOpt.indx2,i*ones(1,NVarA+1-i)];
    end
end

% save tmp2
% return

% if EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4; %calculations needed for analitical Hessian
%     EstimOpt.XXX = permute(mmx('square',permute(INPUT.XXa,[2,4,1,3]),[]),[3,1,2,4])
%     EstimOpt.XXX = zeros(EstimOpt.NAlt*EstimOpt.NCT,NVarA,NVarA, EstimOpt.NP);
%     if EstimOpt.FullCov == 0
%         EstimOpt.VCx = zeros(NVarA,NVarA,EstimOpt.NRep,EstimOpt.NP);
%     else
%         EstimOpt.VCx = zeros(NVarA*(NVarA-1)/2+NVarA,NVarA*(NVarA-1)/2+NVarA,EstimOpt.NRep,EstimOpt.NP);
%     end
%     err_tmp = reshape(err_mtx,NVarA,EstimOpt.NRep,EstimOpt.NP);
%     for i = 1:EstimOpt.NP
%         for j = 1:EstimOpt.NAlt*EstimOpt.NCT
%             EstimOpt.XXX(j,:,:,i) = (INPUT.XXa(j,:,i)')*INPUT.XXa(j,:,i);
%         end
%         for j = 1:EstimOpt.NRep
%             if EstimOpt.FullCov == 0
%                 EstimOpt.VCx(:,:,j,i) = err_tmp(:,j,i)*err_tmp(:,j,i)';
%             else
%                 EstimOpt.VCx(:,:,j,i) = err_tmp(EstimOpt.indx2,j,i)*err_tmp(EstimOpt.indx2,j,i)';
%             end
%         end
%     end

% end


%% Estimation


LLfun = @(B) LL_mxl_MATlike(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,err_mtx,INPUT.W,EstimOpt,OptimOpt,B);

if EstimOpt.ConstVarActive == 0 % no equality constraints
    
    if EstimOpt.HessEstFix == 0
        [Results.bhat,LL,Results.exitf,Results.output,Results.g,Results.hess] = fminunc(LLfun,b0,OptimOpt);
    else
        [Results.bhat,LL,Results.exitf,Results.output,Results.g] = fminunc(LLfun,b0,OptimOpt);
    end

    %         options_tmp = optimset('MaxFunEvals',1e100,'MaxIter',1e3,'TolFun',1e-6,'TolX',1e-6,'OutputFcn',@outputf);
    %
    %         [Results.beta, LL,Results.exitf,Results.output] = fminsearch(LLfun, b0,options_tmp);
    %
    %         save tmp1
    
    %     [x,fval,exitflag,output,lambda,grad,hessian] = knitromatlab(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,extendedFeatures,options,KNITROOptions)
    %     [Results.bhat,LL,Results.exitf,Results.output,Results.lambda,Results.g,Results.hess] = knitromatlab(LLfun,b0,[],[],[],[],[],[],[],[],[],'knitro.opt'); %
    %     [Results.bhat,LL,Results.exitf,Results.output,Results.lambda,Results.g,Results.hess] = knitromatlab(LLfun,b0,[],[],[],[],[],[],[],[],OptimOpt,'knitro.opt'); %

elseif EstimOpt.ConstVarActive == 1 % equality constraints
    EstimOpt.CONS1 = diag(1 - EstimOpt.BActive);
    EstimOpt.CONS1(sum(EstimOpt.CONS1,1)==0,:) = [];
    EstimOpt.CONS2 = zeros(size(EstimOpt.CONS1,1),1);
    %     EstimOpt.CONS1 = sparse(EstimOpt.CONS1);
    %     EstimOpt.CONS2 = sparse(EstimOpt.CONS2);
    if EstimOpt.HessEstFix == 0
        [Results.bhat,LL,Results.exitf,Results.output,Results.lambda,Results.g,Results.hess] = fmincon(LLfun,b0,[],[],EstimOpt.CONS1,EstimOpt.CONS2,[],[],[],OptimOpt);
    else
        [Results.bhat,LL,Results.exitf,Results.output,Results.lambda,Results.g] = fmincon(LLfun,b0,[],[],EstimOpt.CONS1,EstimOpt.CONS2,[],[],[],OptimOpt);
    end
end


%% Output
Results.LL = -LL;
Results.b0_old = b0;

[Results, EstimOpt] = calculateResults(Results, Results_old, INPUT, err_mtx, EstimOpt, OptimOpt);

% save tmp1
% return


%% Tworzebnie templatek do printu


Template1 = {'DetailsA','DetailsV'};
Template2 = {'DetailsA','DetailsV'};
Names.DetailsA = EstimOpt.NamesA;
Heads.DetailsA = {'Means';'tc'};
Heads.DetailsV = {'Standard Deviations';'lb'};
ST = {};
if EstimOpt.NVarM > 0
    Template1 = [Template1,'DetailsM'];
    Temp = cell(1,size(Template2,2));
    Temp(1,1) = {'DetailsM'};
    Template2 = [Template2;Temp];
    Heads.DetailsM(:,2) = [EstimOpt.NamesM;{'lb'}];
    Heads.DetailsM(1:2,1) = {'Interactions of means';'lc'};
end

if EstimOpt.NVarNLT > 0
    Template1 = [Template1,'DetailsNLT0'];
    Temp = cell(1,size(Template2,2));
    Temp(1,1) = {'DetailsNLT0'};
    Template2 = [Template2;Temp];
    if EstimOpt.NLTType == 1
        Heads.DetailsNLT0 = {'Box-Cox transformation parameters';'tb'};
    elseif EstimOpt.NLTType == 2
        Heads.DetailsNLT0 = {'Yeo-Johnson transformation parameters';'tb'};
    end
end

if EstimOpt.Johnson > 0
    Heads.ResultsJ = {'Johnson location parameters';'Johnson scale parameters';'tb'}; %heads need to be written vertically
    Template1 = [Template1,'ResultsJ'];
    Temp = cell(1,size(Template2,2));
    Temp(1,1) = {'ResultsJ'};
    Template2 = [Template2;Temp];
end

if EstimOpt.NVarS > 0
    Temp = cell(1,size(Template1,2));
    Temp(1,1) = {'DetailsS'};
    Template1 = [Template1;Temp];
    Temp = cell(1,size(Template2,2));
    Temp(1,1) = {'DetailsS'};
    Template2 = [Template2;Temp];
    Names.DetailsS = EstimOpt.NamesS;
    Heads.DetailsS = {'Covariates of Scale';'tb'};
    ST = {'DetailsS'};
end


%% Tworzenie naglowka


Head = cell(1,2);
if EstimOpt.FullCov == 0
    Head(1,1) = {'MXL_d'};
else
    Head(1,1) = {'MXL'};
end

if EstimOpt.WTP_space > 0
    Head(1,2) = {'in WTP-space'};
else
    Head(1,2) = {'in preference-space'};
end


%% Tworzenie stopki


Tail = cell(17,2);
Tail(2,1) = {'Model diagnostics'};
Tail(3:17,1) = {'LL at convergence';'LL at constant(s) only';strcat('McFadden''s pseudo-R',char(178));strcat('Ben-Akiva-Lerman''s pseudo-R',char(178));'AIC/n';'BIC/n';'n (observations)';'r (respondents)';'k (parameters)';'';'Estimation method';'Simulation with';'Optimization method';'Gradient';'Hessian'};

if isfield(Results_old,'MNL0') && isfield(Results_old.MNL0,'LL')
    Tail(3:11,2) = num2cell(Results.stats);
end

if any(INPUT.W ~= 1)
    Tail(13,2) = {'weighted simulated maximum likelihood'};
else
    Tail(13,2) = {'simulated maximum likelihood'};
end

switch EstimOpt.Draws
    case 1
        Tail(14,2) = {[num2str(EstimOpt.NRep),' ','pseudo-random draws']};
    case 2
        Tail(14,2) = {[num2str(EstimOpt.NRep),' ','Latin Hypercube Sampling draws']};
    case  3
        Tail(14,2) = {[num2str(EstimOpt.NRep),' ','Halton draws (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),')']};
    case 4
        Tail(14,2) = {[num2str(EstimOpt.NRep),' ','Halton draws with reverse radix scrambling (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),')']};
    case 5
        Tail(14,2) = {[num2str(EstimOpt.NRep),' ','Sobol draws (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),')']};
    case 6
        Tail(14,2) = {[num2str(EstimOpt.NRep),' ','Sobol draws with random linear scramble and random digital shift (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),')']};
end

Tail(15,2) = {OptimOpt.Algorithm;};

if strcmp(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        Tail(16,2) = {'user-supplied, analytical'};
    else
        Tail(16,2) = {['user-supplied, numerical ',num2str(OptimOpt.FinDiffType)]};
    end
else
    Tail(16,2) = {['built-in, ',num2str(OptimOpt.FinDiffType)]};
end

if isequal(OptimOpt.Algorithm,'quasi-newton')
    outHessian='off, ';
    switch EstimOpt.HessEstFix
        case 0
            outHessian = [outHessian,'retained from optimization'];
        case 1
            outHessian = [outHessian,'ex-post calculated using BHHH'];
        case 2
            outHessian = [outHessian,'ex-post calculated using high-precision BHHH'];
        case 3
            outHessian = [outHessian,'ex-post calculated numerically'];
        case 4
            outHessian = [outHessian,'ex-post calculated analytically'];
    end
else
    if strcmp(OptimOpt.Hessian,'user-supplied')
        if EstimOpt.ApproxHess == 1
            outHessian = 'user-supplied, BHHH, ';
        else
            outHessian = 'user-supplied, analytical, ';
        end
    else
        outHessian = ['built-in, ',num2str(OptimOpt.HessUpdate),', '];
    end
    switch EstimOpt.HessEstFix
        case 0
            outHessian = [outHessian,'retained from optimization'];
        case 1
            outHessian = [outHessian,'ex-post calculated using BHHH'];
        case 2
            outHessian = [outHessian,'ex-post calculated using high-precision BHHH'];
        case 3
            outHessian = [outHessian,'ex-post calculated numerically'];
        case 4
            outHessian = [outHessian,'ex-post calculated analytically'];
    end
end
Tail(17,2) = {outHessian};


%% Tworzenie ResultsOut, drukowanie na ekran i do pliku .xls


if EstimOpt.Display~=0
    Results.R_out = genOutput(EstimOpt,Results,Head,Tail,Names,Template1,Template2,Heads,ST);
end


end
