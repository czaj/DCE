function Results = LC(INPUT,Results_old,EstimOpt,OptimOpt)
% LC creates Latent Class Model.
%
% Syntax:   LC(INPUT,EstimOpt,OptimOpt)
%           LC(INPUT,Results_old,EstimOpt,OptimOpt)
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
% LC model assumes that parameter vectors are distributed among individuals with discrete distribution. The analyst does not know from the data which observation is in which class, hence the name latent classes:
% •	NClass = 2; number of latent classes
% •	NamesC – names of classes
% •	BActiveClass vector of 0; for each class set it to 1 to constrain parameters of the attributes with zeros equal between classes
% 
% 
% General basics:
% •	DataFile – path/name of the .mat data file
% •	Display – 1; shows output, set to 0 to hide it 
% •	ProjectName – Name of the project/model
% •	WTP_space – set to 1 for estimation in WTP space. If missing or set to 0, MNL uses Preference Space
% •	NCT - Number of choice tasks per person 
% •	NAlt - Number of alternatives
% •	NP – Number of respondents
% 
% 
% Variables options:
% •	NamesA – Names of variables in list e.g. {'-Opt out';’-Cost (EUR)'}
% •	NamesS – Names of variables of Scale
% 
% Numbers of variables are set automatically; you can check them in the following fields:
% o	NVarA - Number of attributes
% o	NVarS - Number of covariates of scale
% 
% 
% Parameters options:
% •	BActive = vector of 0; for each parameter set it to 1 to constrain model parameters to their initial values
% •	ConstVarActive = 0; set to 1 to constrain model parameters to its initial values 
% 
% 
% Modelling options from DataCleanDCE:
% •	ApproxHess = 1; for user supplied hessians, 1 for BHHH, 0 for analytical
% •	RobustStd = 0; by default not using robust standard errors, set to 1 to use them
% •	NumGrad = 0; uses analytical gradient in calculations, set to 1 for numerical gradient
% •	HessEstFix = 0; Options: 
% o	0 - use optimization Hessian, 
% o	1 - use jacobian-based (BHHH) Hessian, 
% o	2 - use high-precision jacobian-based (BHHH) Hessian,
% o	3 - use numerical Hessian, 
% o	4 - use analytical Hessian
% 
% 
% For drawing and simulations:
% •	HaltonSkip = 1; specify no of rows in halton sequence to skip
% •	HaltonLeap = 0; specify no of rows in halton sequence to leap
% •	Draws = 6; specify draws type, by default Sobol with scrambling. Options: 
% o	1 - pseudo-random, 
% o	2 - Latin Hypercube, 
% o	3 - Halton, 
% o	4 - Halton RR scrambled, 
% o	5 - Sobol, 
% o	6 - Sobol MAO scrambled
% •	NRep = 1e3; specify no. of draws for numerical simulation
% •	RealMin = by default 0, can be set to 1
% •	NSdSim = 1e4; number of draws for simulating standard deviations
%  
% 
% Precision:
% •	eps = 1.e-6; overall precision level
% •	Otherwise:
% o	FunctionTolerance - df / gradient precision level
% o	TolX - step precision level
% o	OptimalityTolerance - dB precision level
% 
% 
% Seeds by default:
% •	Seed1 = 179424673
% •	Seed2 = 7521436817
% 
% Example: 
%    Results.LC = LC(INPUT,Results,EstimOpt,OptimOpt);
%
% Author: Mikolaj Czajkowski, Professor
% University of Warsaw, Faculty of Economic Sciences
% email address: mik@czaj.org 
% Website: http://czaj.org/#

% save tmp_LC
% return

tic

global B_backup

Results.bhat = [];
Results.R = [];
Results.R_out = {};
Results.stats = [];
NSdSim = EstimOpt.NSdSim;

%% Check data and inputs

if nargin < 3 % check no. of inputs
    error('Too few input arguments for LC(INPUT,EstimOpt,OptimOpt)')
end

warning off MATLAB:mir_warning_maybe_uninitialized_temporary

format shortG;
format compact;

if isfield(EstimOpt,'Display') == 0
    EstimOpt.Display = 1;
end

if isfield(EstimOpt,'NClass') == 0
    cprintf(rgb('DarkOrange'),'WARNING: Number of latent classes not specified - assuming 2 classes. \n')
    EstimOpt.NClass = 2;
end

if EstimOpt.Display == 1
    disp(' ');
    disp('__________________________________________________________________________________________________________________');
    disp(' ');
end
if EstimOpt.Display == 1
    if any(INPUT.W ~= 1)
        cprintf('Black','Estimating '); cprintf('*Black','weighted '); cprintf('Black','LC model with '); cprintf('Black',num2str(EstimOpt.NClass)); cprintf('Black',' classes ...\n');
    else
        disp(num2str(EstimOpt.NClass,'Estimating LC model with %1.0f classes ...'))
    end
end

if isfield(EstimOpt,'WTP_space') == 0
    EstimOpt.WTP_space = 0;
    EstimOpt.WTP_matrix = [];
elseif EstimOpt.WTP_space == 0
    EstimOpt.WTP_matrix = [];
end

if EstimOpt.Display == 1
    if EstimOpt.WTP_space > 0
        disp('in WTP-space.')
    else
        disp('in preference-space.')
    end
end

if isfield(INPUT,'Xc') == 0 || numel(INPUT.Xc) == 0
    INPUT.Xc = ones(size(INPUT.Y,1),1);
else
    INPUT.Xc = [ones(size(INPUT.Y,1),1),INPUT.Xc];
end
EstimOpt.NVarC = size(INPUT.Xc,2); % no. of variables explaining class probabilities

if isfield(INPUT,'Xs') == 0 || numel(INPUT.Xs) == 0
    INPUT.Xs = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarS = size(INPUT.Xs,2); % no. of variables explaining class probabilities

if isfield(EstimOpt,'ClassScores') == 0
    EstimOpt.ClassScores = 0;
end

if EstimOpt.WTP_space > 0
    if isfield(EstimOpt,'WTP_matrix') == 0
        WTP_att = (EstimOpt.NVarA - EstimOpt.WTP_space)/EstimOpt.WTP_space;
        if rem(WTP_att,1) ~= 0
            error('EstimOpt.WTP_matrix associating attributes with cost parameters not provided')
        else
            if EstimOpt.WTP_space > 1
                disp(['EstimOpt.WTP_matrix associating attributes with cost parameters not provided - assuming equal shares for each of the ',num2str(EstimOpt.WTP_space),' monetary attributes'])
            end
            EstimOpt.WTP_matrix = EstimOpt.NVarA - EstimOpt.WTP_space + kron(1:EstimOpt.WTP_space,ones(1,WTP_att));
        end
    elseif size(EstimOpt.WTP_matrix,2) ~= EstimOpt.NVarA - EstimOpt.WTP_space
        error('Dimensions of EstimOpt.WTP_matrix not correct - for each non-monetary attribute provide no. of attribute to multiply it with')
    else
        EstimOpt.WTP_matrix = EstimOpt.WTP_matrix(:)';
    end
end
if isfield(EstimOpt,'NamesA') == 0 || isempty(EstimOpt.NamesA) || length(EstimOpt.NamesA) ~= EstimOpt.NVarA
    EstimOpt.NamesA = (1:EstimOpt.NVarA)';
    EstimOpt.NamesA = cellstr(num2str(EstimOpt.NamesA));
elseif size(EstimOpt.NamesA,1) ~= EstimOpt.NVarA
    EstimOpt.NamesA = EstimOpt.NamesA';
end

if EstimOpt.NVarC > 1
    if isfield(EstimOpt,'NamesC') == 0 || isempty(EstimOpt.NamesC)|| length(EstimOpt.NamesC) ~= (EstimOpt.NVarC-1)
        EstimOpt.NamesC = (1:EstimOpt.NVarC-1)';
        EstimOpt.NamesC = cellstr(num2str(EstimOpt.NamesC));
    elseif size(EstimOpt.NamesC,1) ~= EstimOpt.NVarC - 1
        EstimOpt.NamesC = EstimOpt.NamesC';
    end
    EstimOpt.NamesC = [{'Cons'};EstimOpt.NamesC];
else
    EstimOpt.NamesC = {'Cons'};
end
if EstimOpt.NVarS > 0
    if isfield(EstimOpt,'NamesS') == 0 || isempty(EstimOpt.NamesS) || length(EstimOpt.NamesS) ~= EstimOpt.NVarS
        EstimOpt.NamesS = (1:EstimOpt.NVarS)';
        EstimOpt.NamesS = cellstr(num2str(EstimOpt.NamesS));
    elseif size(EstimOpt.NamesS,1) ~= EstimOpt.NVarS
        EstimOpt.NamesS = EstimOpt.NamesS';
    end
end

if ~isfield(EstimOpt,'BActiveClass') || isempty(EstimOpt.BActiveClass) || sum(EstimOpt.BActiveClass == 0) == 0
    EstimOpt.BActiveClass = ones(EstimOpt.NVarA,1);
end

% for later:
% - covariates of scale (class-specific?)
% - allow for non-constant scale parameters for classes

%% Starting values

EstimOpt.jitter1 = 1.1; % Jittering parameter (relative) for MNL starting values (attributes)
EstimOpt.jitter2 = 0.1; % Jittering parameter (absolute) for class probabilities starting values

if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == (EstimOpt.NVarA+EstimOpt.NVarS)*EstimOpt.NClass + (EstimOpt.NClass-1)*EstimOpt.NVarC
    b0 = B_backup(:);
    if EstimOpt.Display == 1
        disp('Using the starting values from Backup')
    end
elseif isfield(Results_old,'LC') && isfield(Results_old.LC,'b0') % starting values provided
    Results_old.LC.b0_old = Results_old.LC.b0(:);
    Results_old.LC = rmfield(Results_old.LC,'b0');
    if length(Results_old.LC.b0_old) ~= (EstimOpt.NVarA + EstimOpt.NVarS)*EstimOpt.NClass + (EstimOpt.NClass - 1)*EstimOpt.NVarC
        if EstimOpt.Display == 1
            cprintf(rgb('DarkOrange'),'WARNING: Incorrect no. of starting values or model specification \n')
        end
        Results_old.LC = rmfield(Results_old.LC,'b0_old');
    else
        b0 = Results_old.LC.b0_old(:);
    end
end
if ~exist('b0','var')
    if isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat')
        if EstimOpt.Display == 1
            disp('Using MNL results as starting values')
        end
        if length(Results_old.MNL.bhat((1:size(Results_old.MNL.bhat,1))'*ones(1,EstimOpt.NClass),:)) == length((EstimOpt.jitter1.*unifrnd(0,ones((EstimOpt.NVarA+EstimOpt.NVarS).*EstimOpt.NClass,1))))
            b0 = [Results_old.MNL.bhat((1:size(Results_old.MNL.bhat,1))'*ones(1,EstimOpt.NClass),:).*(EstimOpt.jitter1.*unifrnd(0,ones((EstimOpt.NVarA+EstimOpt.NVarS).*EstimOpt.NClass,1)));EstimOpt.jitter2+unifrnd(-1,ones(EstimOpt.NVarC.*(EstimOpt.NClass-1),1))];
        else
            error('Run correct MNL options');
        end
    else
        error('No starting values available - run MNL first')
    end
end

%% Optimization Options

if  isfield(EstimOpt,'BActive')
    EstimOpt.BActive = EstimOpt.BActive(:)';
end

if isfield(EstimOpt,'BActiveClass') && ~isempty(EstimOpt.BActiveClass) && sum(EstimOpt.BActiveClass == 0) > 0 && isfield(EstimOpt,'BActive') && ~isempty(EstimOpt.BActive) && sum(EstimOpt.BActive == 0) > 0
    disp('WARINING: The model does not currently support simultaneous BActive and BActiveClass constraints')
end

if isfield(EstimOpt,'BActiveClass')
    EstimOpt.BActiveClass = EstimOpt.BActiveClass(:);
end
if ~isfield(EstimOpt,'BActiveClass') || isempty(EstimOpt.BActiveClass) || sum(EstimOpt.BActiveClass == 0) == 0
    EstimOpt.BActiveClass = ones(EstimOpt.NVarA,1);
else
    if length(EstimOpt.BActiveClass) ~= EstimOpt.NVarA
        error('Check no. of constraints (BActiveClass)')
    else
        disp(['Parameters of the attributes with zeros constrained equal between classes: ' mat2str(EstimOpt.BActiveClass)])
    end
end

if EstimOpt.ConstVarActive == 1
    if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive) || sum(EstimOpt.BActive == 0) == 0
        error('Are there any constraints on model parameters (EstimOpt.ConstVarActive)? Constraints not provided (EstimOpt.BActive).')
    elseif length(b0) ~= length(EstimOpt.BActive)
        error('Check no. of constraints')
    end
    if EstimOpt.Display == 1
        disp(['Initial values: ' mat2str(b0',2)])
        disp(['Parameters with zeros are constrained to their initial values: ' mat2str(EstimOpt.BActive')])
    end
else
    if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive) || sum(EstimOpt.BActive == 0) == 0
        EstimOpt.BActive = ones(1,length(b0));
        if EstimOpt.Display == 1
            disp(['Initial values: ' mat2str(b0',2)])
        end
    else
        if length(b0) ~= length(EstimOpt.BActive)
            error('Check no. of constraints')
        else
            if EstimOpt.Display == 1
                disp(['Initial values: ' mat2str(b0',2)])
                disp(['Parameters with zeros are constrained to their initial values: ' mat2str(EstimOpt.BActive')])
            end
        end
    end
end

if any(EstimOpt.MissingCT(:) == 1) && EstimOpt.NumGrad == 0
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient to numerical - missing choice tasks not supported by analytical gradient \n')
end

if any(EstimOpt.MissingAlt(:) == 1) && EstimOpt.NumGrad == 0
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient to numerical - missing alternatives not supported by analytical gradient \n')
end

if ((isfield(EstimOpt,'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    OptimOpt.GradObj = 'on';
end

if any(EstimOpt.BActiveClass == 0) && EstimOpt.NumGrad == 0
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient to numerical - parameters constrained between classes not supported by analytical gradient \n')
end

if EstimOpt.RobustStd == 1 && (EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2)
    EstimOpt.RobustStd = 0;
    cprintf(rgb('DarkOrange'),'WARNING: Setting off robust standard errors, they do not matter for BHHH aproximation of hessian \n')
end

if (isfield(EstimOpt, 'ConstVarActive') == 0 || EstimOpt.ConstVarActive == 0) && isequal(OptimOpt.Algorithm,'quasi-newton') && isequal(OptimOpt.Hessian,'user-supplied')
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied Hessian off - quasi-newton algorithm does not use it anyway \n')
    OptimOpt.Hessian = 'off';
end

cprintf('Opmization algorithm: '); cprintf('*Black',[OptimOpt.Algorithm '\n'])

if strcmp(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        cprintf('Gradient: '); cprintf('*Black','user-supplied, analytical \n')
    else
        cprintf('Gradient: '); cprintf('*Black',['user-supplied, numerical, ' OptimOpt.FinDiffType '\n'])
    end
else
    cprintf('Gradient: '); cprintf('*Black',['built-in, ' OptimOpt.FinDiffType '\n'])
end

if isequal(OptimOpt.Algorithm,'quasi-newton')
    cprintf('Hessian: '); cprintf('*Black','off, ')
    switch EstimOpt.HessEstFix
        case 0
            cprintf('*Black','retained from optimization \n')
        case 1
            cprintf('*Black','ex-post calculated using BHHH \n')
        case 2
            cprintf('*Black','ex-post calculated using high-precision BHHH \n')
        case 3
            cprintf('*Black','ex-post calculated numerically \n')
        case 4
            cprintf('*Black','ex-post calculated analytically \n')
    end
else
    if strcmp(OptimOpt.Hessian,'user-supplied')
        if EstimOpt.ApproxHess == 1
            cprintf('Hessian: '); cprintf('*Black','user-supplied, BHHH, ')
        else
            cprintf('Hessian: '); cprintf('*Black','user-supplied, analytical, ')
        end
    else
        cprintf('Hessian: '); cprintf('*Black',['built-in, ' OptimOpt.HessUpdate ', '])
    end
    switch EstimOpt.HessEstFix
        case 0
            cprintf('*Black','retained from optimization \n')
        case 1
            cprintf('*Black','ex-post calculated using BHHH \n')
        case 2
            cprintf('*Black','ex-post calculated using high-precision BHHH \n')
        case 3
            cprintf('*Black','ex-post calculated numerically \n')
        case 4
            cprintf('*Black','ex-post calculated analytically \n')
    end
end

%% Restructuring data

INPUT.XXc = INPUT.Xc(1:EstimOpt.NCT*EstimOpt.NAlt:end,:); % NP x NVarC


INPUT.YY = reshape(INPUT.Y,[EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP]);
INPUT.YY = INPUT.YY(:,(1:size(INPUT.YY,2))'*ones(1,EstimOpt.NClass)); %NAlt x NCT*NP*NClass
INPUT.YY(isnan(INPUT.YY)) = 0;
INPUT.MissingInd = INPUT.MissingInd((1:size(INPUT.MissingInd,1))'*ones(1,EstimOpt.NClass),:);
INPUT.MissingInd = reshape(INPUT.MissingInd,[EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP*EstimOpt.NClass]); %NAlt x NCT*NP*NClass

if sum(EstimOpt.BActiveClass == 0,1) > 0
    bactive_1 = EstimOpt.BActive(1:EstimOpt.NClass*EstimOpt.NVarA);
    b0_1 = b0(1:EstimOpt.NClass*EstimOpt.NVarA);
    bactive_2 = EstimOpt.BActive(EstimOpt.NClass*EstimOpt.NVarA+1:end);
    b0_2 = b0(EstimOpt.NClass*EstimOpt.NVarA+1:end);
    b0 = b0_1(1:EstimOpt.NVarA);
    EstimOpt.BActive = bactive_1(1:EstimOpt.NVarA);
    
    for i = 2:EstimOpt.NClass
        b0x = b0_1((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA);
        bactivex = bactive_1((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA);
        b0 = [b0;b0x(EstimOpt.BActiveClass == 1)];
        EstimOpt.BActive = [EstimOpt.BActive,bactivex(EstimOpt.BActiveClass == 1)];
    end
    b0 = [b0;b0_2];
    EstimOpt.BActive = [EstimOpt.BActive,bactive_2];
    clear bactive_1 b0_1 bactive_2 b0_2 b0x bactivex
end

%% Estimation

LLfun = @(B) LL_lc_MATlike(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xs,INPUT.MissingInd,INPUT.W,EstimOpt,OptimOpt,B);
if EstimOpt.ConstVarActive == 0
    
    if EstimOpt.HessEstFix == 0
        [Results.bhat,LL,Results.exitf,Results.output,Results.g,Results.hess] = fminunc(LLfun,b0,OptimOpt);
    else
        [Results.bhat,LL,Results.exitf,Results.output,Results.g] = fminunc(LLfun,b0,OptimOpt);
    end
    
elseif EstimOpt.ConstVarActive == 1 % equality constraints
    EstimOpt.CONS1 = diag(1-EstimOpt.BActive);
    EstimOpt.CONS1(sum(EstimOpt.CONS1,1) == 0,:) = [];
    EstimOpt.CONS2 = zeros(size(EstimOpt.CONS1,1),1);
    if EstimOpt.HessEstFix == 0
        [Results.bhat,LL,Results.exitf,Results.output,Results.lambda,Results.g,Results.hess] = fmincon(LLfun,b0,[],[],EstimOpt.CONS1,EstimOpt.CONS2,[],[],[],OptimOpt);
    else
        [Results.bhat,LL,Results.exitf,Results.output,Results.lambda,Results.g] = fmincon(LLfun,b0,[],[],EstimOpt.CONS1,EstimOpt.CONS2,[],[],[],OptimOpt);
    end
    
end

%% Output

Results.LL = -LL;
Results.b0_old = b0;

Results.LLdetailed = LL_lc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xs,INPUT.MissingInd,EstimOpt,Results.bhat);
Results.LLdetailed = Results.LLdetailed.*INPUT.W;

if any(INPUT.MissingInd == 1) % In case of some missing data
    idx = sum(reshape(INPUT.MissingInd,[EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP])) == EstimOpt.NAlt;
    idx = sum(reshape(idx,[EstimOpt.NCT,EstimOpt.NP]),1)'; % no. of missing NCT for every respondent
    idx = EstimOpt.NCT - idx;
    R2 = mean(exp(-Results.LLdetailed./idx),1);
else
    R2 = mean(exp(-Results.LLdetailed/EstimOpt.NCT),1);
end

if EstimOpt.HessEstFix == 1
    f = LL_lc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xs,INPUT.MissingInd,EstimOpt,Results.bhat);
    Results.jacobian = numdiff(@(B) INPUT.W.*LL_lc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xs,INPUT.MissingInd,EstimOpt,B),INPUT.W.*f,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) INPUT.W.*LL_lc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xs,INPUT.MissingInd,EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.hess = hessian(@(B) sum(INPUT.W.*LL_lc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xs,INPUT.MissingInd,EstimOpt,B)),Results.bhat);
end

% if sum(EstimOpt.BActive == 0) > 0
%     if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
%         Results.jacobian = Results.jacobian(:, EstimOpt.BActive == 1);
%         Results.hess = Results.jacobian'*Results.jacobian;
%     elseif EstimOpt.HessEstFix == 0 || EstimOpt.HessEstFix == 3
%         Results.hess = Results.hess(EstimOpt.BActive == 1,EstimOpt.BActive == 1);
%     end
%     Results.ihess = inv(Results.hess);
%     Results.ihess = direcXpnd(Results.ihess,EstimOpt.BActive);
%     Results.ihess = direcXpnd(Results.ihess',EstimOpt.BActive);
%     Results.std = sqrt(diag(Results.ihess));
% 	Results.std(EstimOpt.BActive == 0) = NaN;
% else
% 	if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
%         Results.hess = Results.jacobian'*Results.jacobian;
%     end
%     Results.ihess = full(inv(sparse(Results.hess)));
%     Results.std = sqrt(diag(Results.ihess));
% end

if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
    Results.hess = Results.jacobian'*Results.jacobian;
end
EstimOpt.BLimit = (sum(Results.hess) == 0 & EstimOpt.BActive == 1);
EstimOpt.BActive(EstimOpt.BLimit == 1) = 0;
Results.hess = Results.hess(EstimOpt.BActive == 1,EstimOpt.BActive == 1);
Results.ihess = inv(Results.hess);
Results.ihess = direcXpnd(Results.ihess,EstimOpt.BActive);
Results.ihess = direcXpnd(Results.ihess',EstimOpt.BActive);
if EstimOpt.RobustStd == 1
    if EstimOpt.NumGrad == 0
        [~, Results.jacobian] = LL_lc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xs,INPUT.MissingInd,EstimOpt,Results.bhat);
        Results.jacobian = Results.jacobian.*INPUT.W;
    else
        Results.jacobian = numdiff(@(B) INPUT.W.*LL_lc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xs,INPUT.MissingInd,EstimOpt,B),Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
    end
    RobustHess = Results.jacobian'*Results.jacobian;
    Results.ihess = Results.ihess*RobustHess*Results.ihess;
end

Results.std = sqrt(diag(Results.ihess));
Results.std(EstimOpt.BActive == 0) = NaN;
Results.std(EstimOpt.BLimit == 1) = 0;
Results.std(imag(Results.std) ~= 0) = NaN;

Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;
Results.INPUT = INPUT;

EstimOpt.params = length(Results.bhat);
EstimOpt.params = EstimOpt.params - sum(EstimOpt.BActive == 0) + sum(EstimOpt.BLimit == 1);

if sum(EstimOpt.BActiveClass == 0,1) == 0
    for i = 1:EstimOpt.NClass
        Results.DetailsA(1:EstimOpt.NVarA,4*i-3) = Results.bhat((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA);
        Results.DetailsA(1:EstimOpt.NVarA,4*i-1:4*i) = [Results.std((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA),pv(Results.bhat((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA),Results.std((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA))];
    end
    l = EstimOpt.NClass*EstimOpt.NVarA;
else
    Bclass = Results.bhat(1:EstimOpt.NVarA)*ones(1,EstimOpt.NClass);
    stdclass = Results.std(1:EstimOpt.NVarA)*ones(1,EstimOpt.NClass);
    for i = 1:EstimOpt.NClass-1
        Bclass(EstimOpt.BActiveClass == 1,i+1) = Results.bhat(EstimOpt.NVarA+(i-1)*sum(EstimOpt.BActiveClass,1)+1:EstimOpt.NVarA+i*sum(EstimOpt.BActiveClass,1));
        stdclass(EstimOpt.BActiveClass == 1,i+1) = Results.std(EstimOpt.NVarA+(i-1)*sum(EstimOpt.BActiveClass,1)+1:EstimOpt.NVarA+i*sum(EstimOpt.BActiveClass,1));
    end
    stdclass(EstimOpt.BActiveClass == 0,2:end) = NaN;
    for i = 1:EstimOpt.NClass
        Results.DetailsA(1:EstimOpt.NVarA,4*i-3) = Bclass(:,i);
        Results.DetailsA(1:EstimOpt.NVarA,4*i-1:4*i) = [stdclass(:,i),pv(Bclass(:,i),stdclass(:,i))];
    end
    l = EstimOpt.NVarA + (EstimOpt.NClass - 1)*sum(EstimOpt.BActiveClass,1);
end

if EstimOpt.NVarS > 0
    for i = 1:EstimOpt.NClass
        Results.DetailsS(:,4*i-3) = Results.bhat(l+(i-1)*EstimOpt.NVarS+1:l+i*EstimOpt.NVarS);
        Results.DetailsS(:,4*i-1:4*i) = [Results.std(l+(i-1)*EstimOpt.NVarS+1:l+i*EstimOpt.NVarS),pv(Results.bhat(l+(i-1)*EstimOpt.NVarS+1:l+i*EstimOpt.NVarS),Results.std(l+(i-1)*EstimOpt.NVarS+1:l+i*EstimOpt.NVarS))];
    end
    l = l + EstimOpt.NClass*EstimOpt.NVarS;
end

for i = 1:EstimOpt.NClass-1
    for j = 1:EstimOpt.NVarC
        Results.DetailsV(j,4*i-3) = Results.bhat(l+j+EstimOpt.NVarC*(i-1));
        Results.DetailsV(j,4*i-1:4*i) = [Results.std(l+j+EstimOpt.NVarC*(i-1)),pv(Results.bhat(l+j+EstimOpt.NVarC*(i-1)),Results.std(l+j+EstimOpt.NVarC*(i-1)))];
    end
end

Results.DetailsV = [Results.DetailsV,zeros(EstimOpt.NVarC,1),NaN(EstimOpt.NVarC,3)];

if sum(EstimOpt.BActiveClass == 0,1) == 0
    bclass = reshape([Results.bhat((EstimOpt.NVarA+EstimOpt.NVarS)*EstimOpt.NClass+1:end);zeros(EstimOpt.NVarC,1)],[EstimOpt.NVarC,EstimOpt.NClass]);
    bclass_sim = reshape([mvnrnd(Results.bhat((EstimOpt.NVarA+EstimOpt.NVarS)*EstimOpt.NClass+1:end),Results.ihess((EstimOpt.NVarA+EstimOpt.NVarS)*EstimOpt.NClass+1:end,(EstimOpt.NVarA+EstimOpt.NVarS)*EstimOpt.NClass+1:end),NSdSim)';zeros(EstimOpt.NVarC,NSdSim)],[EstimOpt.NVarC,EstimOpt.NClass,NSdSim]);
else
    bclass = reshape([Results.bhat((EstimOpt.NClass-1)*sum(EstimOpt.BActiveClass,1)+EstimOpt.NVarA+EstimOpt.NVarS*EstimOpt.NClass+1:end);zeros(EstimOpt.NVarC,1)],[EstimOpt.NVarC,EstimOpt.NClass]);
    bclass_sim = reshape([mvnrnd(Results.bhat((EstimOpt.NClass-1)*sum(EstimOpt.BActiveClass,1)+EstimOpt.NVarA+EstimOpt.NVarS*EstimOpt.NClass+1:end),Results.ihess((EstimOpt.NClass-1)*sum(EstimOpt.BActiveClass,1)+EstimOpt.NVarA+EstimOpt.NVarS*EstimOpt.NClass+1:end,(EstimOpt.NClass-1)*sum(EstimOpt.BActiveClass,1)+EstimOpt.NVarA+EstimOpt.NVarS*EstimOpt.NClass+1:end),NSdSim)';zeros(EstimOpt.NVarC,NSdSim)],[EstimOpt.NVarC,EstimOpt.NClass,NSdSim]);
    Results.bhat = [Results.DetailsA(:,1);Results.DetailsV(:,1)];
end

V = exp(INPUT.XXc*bclass);
Results.PClass = zeros(1,4*EstimOpt.NClass);
Results.PClass(1,1:4:EstimOpt.NClass*4-3) = mean(V./sum(V,2),1);

PClass_mean = zeros(NSdSim,EstimOpt.NClass);
XXc = INPUT.XXc;
parfor i = 1:NSdSim
    bhat_sim_i = bclass_sim(:,:,i);
    V_i = exp(XXc*bhat_sim_i);
    PC_i = V_i./sum(V_i,2);
    PClass_mean(i,:) = mean(PC_i,1);
end
Results.PClass(1,3:4:EstimOpt.NClass*4-1) = std(PClass_mean);
Results.PClass(1,4:4:EstimOpt.NClass*4) = pv(Results.PClass(1,1:4:EstimOpt.NClass*4-3),Results.PClass(1,3:4:EstimOpt.NClass*4-1));
Results.PClass(1,1:2:end) = Results.PClass(1,1:2:end)*100;
Results.PClass95ci = [quantile(PClass_mean,0.025);quantile(PClass_mean,0.975)];

Results.stats = [Results.LL;Results_old.MNL0.LL;1-Results.LL/Results_old.MNL0.LL;R2;((2*EstimOpt.params-2*Results.LL))/EstimOpt.NObs;((log(EstimOpt.NObs)*EstimOpt.params-2*Results.LL))/EstimOpt.NObs;EstimOpt.NObs;EstimOpt.NP;EstimOpt.params];

%% Header

Head = cell(1,2);
Head(1,1) = {'LC'};
if EstimOpt.WTP_space > 0
    Head(1,2) = {'in WTP-space'};
else
    Head(1,2) = {'in preference-space'};
end

%% Results

Template1 = {'DetailsA'};
Template2 = {'DetailsA'};
Names.DetailsA = EstimOpt.NamesA;
ST = {'DetailsA'};

for i = 1:EstimOpt.NClass
    Heads.DetailsA{i,1} = num2str(i, 'Class %1.0f');
end
Heads.DetailsA(end+1) = {'tb'};

if EstimOpt.NVarS >0
    Template1 = [Template1;{'DetailsS'}];
    Template2 = [Template2;{'DetailsS'}];
    Names.DetailsS = EstimOpt.NamesS;
    Heads.DetailsS(:,2) = Heads.DetailsA;
    Heads.DetailsS(end,2) = {'tb'};
    Heads.DetailsS(1:2,1) = {'Explanatory variables of scale';'lb'};
    ST = [ST,{'DetailsS'}];
end

Template1 = [Template1;{'DetailsV'}];
Template2 = [Template2;{'DetailsV'}];
Names.DetailsV = EstimOpt.NamesC;
Heads.DetailsV(:,2) = Heads.DetailsA(1:end-1);
Heads.DetailsV(end+1,2) = {'tb'};
Heads.DetailsV(1:2,1) = {'Probability model';'lb'};
ST = [ST,{'DetailsV'}];

Template1 = [Template1;{'PClass'}];
Template2 = [Template2;{'PClass'}];
Names.PClass = {"(%)"};
Heads.PClass(:,2) = Heads.DetailsA;
Heads.PClass(1:2,1) = {'Average class probabilities';'lb'};
ST = [ST,{'PClass'}];

%% Footer

Tail = cell(17,2);
Tail(2,1) = {'Model diagnostics'};
Tail(3:17,1) = {'LL at convergence';'LL at constant(s) only';strcat('McFadden''s pseudo-R',char(178));strcat('Ben-Akiva-Lerman''s pseudo-R',char(178));'AIC/n';'BIC/n';'n (observations)';'r (respondents)';'k (parameters)';' ';'Estimation method';'Simulation with';'Optimization method';'Gradient';'Hessian'};

if isfield(Results_old,'MNL0') && isfield(Results_old.MNL0,'LL')
    Tail(3:11,2) = num2cell(Results.stats);
end

if any(INPUT.W ~= 1)
    Tail(13,2) = {'weighted maximum likelihood'};
else
    Tail(13,2) = {'maximum likelihood'};
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

%%  Print to screen and .xls

if EstimOpt.Display ~= 0
    Results.Dist = -ones(EstimOpt.NVarA,1);
    Results.R_out = genOutput(EstimOpt,Results,Head,Tail,Names,Template1,Template2,Heads,ST);
end
end
