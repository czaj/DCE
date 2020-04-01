function Results = MIMIC(INPUT,Results_old,EstimOpt,OptimOpt)
% MIMIC creates Multiple Indicators Multiple Causes Model.
%
% Syntax:   MIMIC(INPUT,EstimOpt,OptimOpt)
%           MIMIC(INPUT,Results_old,EstimOpt,OptimOpt)
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
% The user should define variables to structural and measurement equations by setting them in INPUT.Xstr and INPUT.Xmea accordingly.
% •	MeaMatrix – matrix of measurement equations, by default model is assuming that every latent variable is in every measurement equation
% •	MeaSpecMatrix - measurement specification matrix
% •	NLatent = 1; number of latent variables
% •	NClass = 2; number of latent classes
% •	NamesC – names of classes
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
%
% Example: 
%    Results.MIMIC = MIMIC(INPUT,Results,EstimOpt,OptimOpt);
%
% Author: Mikolaj Czajkowski, Professor
% University of Warsaw, Faculty of Economic Sciences
% email address: mik@czaj.org 
% Website: http://czaj.org/#


% save tmp_MIMIC
% return

global B_backup

tic

Results.bhat = [];
Results.R = [];
Results.R_out = {};
Results.stats = [];


%% Check data and inputs


if nargin < 3 % check no. of inputs
	error('Too few input arguments for MIMIC(INPUT,EstimOpt,OptimOpt)')
end

disp(' ');
disp('__________________________________________________________________________________________________________________');
disp(' ');

warning off MATLAB:mir_warning_maybe_uninitialized_temporary

format shortG;
format compact;

if isfield(EstimOpt,'NLatent') == 0
	EstimOpt.NLatent = 1;
	disp(['Assuming ',num2str(EstimOpt.NLatent)',' NLatent Variable(s)']);
end

if any(INPUT.W ~= 1)
    cprintf('Black','Estimating '); cprintf('*Black','weighted '); cprintf('Black','MIMIC with '); cprintf('Black',num2str(EstimOpt.NLatent)); cprintf('Black',' Latent Variable(s) ...\n');
else
    disp(num2str(EstimOpt.NLatent,'Estimating MIMIC with %1.0f Latent Variable(s)'))
end

if isfield(INPUT,'Xstr') == 0 || numel(INPUT.Xstr) == 0
% 	error('Define variables to structural equations')
	cprintf(rgb('DarkOrange'),'WARNING: Structural equations empty.\n')
    INPUT.Xstr = zeros(size(INPUT.Y,1),0);
end
if isfield(INPUT,'Xmea') == 0 || numel(INPUT.Xmea) == 0
% 	error('Define attitudes to measurment equations')
    cprintf(rgb('DarkOrange'),'WARNING: Measurement equations empty.\n')
    INPUT.Xmea = zeros(size(INPUT.Y,1),0);
end

if isfield(INPUT,'Xmea_exp') == 0 || numel(INPUT.Xmea_exp) == 0 % additional covariates for explaining Measurment equations
	INPUT.Xmea_exp = [];
    EstimOpt.MeaExpMatrix = zeros(1,size(INPUT.Xmea,2));
else
    if isfield(EstimOpt,'MeaExpMatrix') == 0 || length(EstimOpt.MeaExpMatrix) ~= size(INPUT.Xmea,2)
        EstimOpt.MeaExpMatrix = ones(1,size(INPUT.Xmea,2));
        cprintf(rgb('DarkOrange'),'WARNING: MeaExpMatrix not defined - assuming that every measurment equation is explained with additional covariates \n')
    else
        EstimOpt.MeaExpMatrix = EstimOpt.MeaExpMatrix(:)';
    end
end

if isfield(EstimOpt,'MeaMatrix')
    if size(EstimOpt.MeaMatrix,2) ~= size(INPUT.Xmea,2) || size(EstimOpt.MeaMatrix,1) ~= EstimOpt.NLatent
        error('Measurment equations erroneously defined (incorrect size of the MeaMatrix)')
    elseif ismember(0,(ismember(EstimOpt.MeaMatrix,[0,1])))
        error('Measurment equations erroneously defined (incorrect values in the MeaMatrix)')
    elseif any(any(EstimOpt.MeaMatrix == 1,2) == 0)
        error('Measurment equations erroneously defined (some Latent Variables without measurement equations in the MeaMatrix)')
    elseif any(any(EstimOpt.MeaMatrix == 1) == 0)
        error('Measurment equations erroneously defined (some measurement variables unused)')
    end
else        
    EstimOpt.MeaMatrix = ones(EstimOpt.NLatent,size(INPUT.Xmea,2));
    disp('Assuming every Latent Variable in every measurment equation')
end

if isfield(EstimOpt,'MeaSpecMatrix') == 0 || numel(EstimOpt.MeaSpecMatrix) == 0
	EstimOpt.MeaSpecMatrix = zeros(1,size(INPUT.Xmea,2));
	disp('Using OLS for every measurment equation')
end
if numel(EstimOpt.MeaSpecMatrix) == 1 % if only one number provided - replicate it for all Measurement variables
    EstimOpt.MeaSpecMatrix = EstimOpt.MeaSpecMatrix*ones(1,size(INPUT.Xmea,2));
end
if size(EstimOpt.MeaSpecMatrix,2) ~= size(INPUT.Xmea,2)
	error('Measurment specification (model) erroneously defined')
end

EstimOpt.NVarStr = size(INPUT.Xstr,2);  % no. of variables in structural equation
EstimOpt.NVarMea = sum(sum(EstimOpt.MeaMatrix)); % no of parameters for Measurments without couting cutoffs, constants etc
EstimOpt.NVarMeaExp = size(INPUT.Xmea_exp,2);

if isfield(EstimOpt,'MissingIndMea') == 0
    EstimOpt.MissingIndMea = zeros(size(INPUT.Xmea));
end

if any(size(EstimOpt.MissingIndMea) ~= size(INPUT.Xmea))
    error('Incorrect size of EstimOpt.MissingIndMea matrix (must be NALT*NCT*NP x NXmea)')
end

INPUT.Xmea(EstimOpt.MissingIndMea == 1) = NaN;


if isfield(INPUT,'Xm') == 0
    INPUT.Xm = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarM = size(INPUT.Xm,2); % Number of covariates of means of random parameters
if isfield(INPUT,'Xs') == 0
    INPUT.Xs = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarS = size(INPUT.Xs,2); % Number of covariates of scale

% for i=1:size(EstimOpt.MeaMatrix,2)
%     if numel(EstimOpt.MeaSpecMatrix(i) > 0) > 0
%         if EstimOpt.MeaSpecMatrix(i) > 0 && numel(unique(INPUT.Xmea(:,i))) > 10 && EstimOpt.MeaSpecMatrix(i) < 3
%             cprintf(rgb('DarkOrange'),'WARNING: There are over 10 levels for measurement variable %d \n',i)
%         end
%     end
%     if sum(isnan(INPUT.Xmea(:,i))) > 0
%         cprintf(rgb('DarkOrange'),'WARNING: Measurement variable %d contains NaN values \n',i)
%     end
% 	if sum(isinf(INPUT.Xmea(:,i))) > 0
%         cprintf(rgb('DarkOrange'),'WARNING: Measurement variable %d contains Inf values \n',i)
% 	end
% end
% if any(sum(EstimOpt.MeaMatrix,2) == 1) == 1
%     error('There must be more than one attitude for every Latent Variable for model identification')    
% end

for i = 1:size(EstimOpt.MeaMatrix,2)
    if sum(isnan(INPUT.Xmea(INPUT.MissingInd==0 & (EstimOpt.MissingIndMea(:,i) == 0),i))) > 0
        cprintf(rgb('DarkOrange'),'WARNING: Measurement variable %d contains NaN values - they will be treated as mising. \n', i)
        EstimOpt.MissingIndMea(isnan(INPUT.Xmea(:,i)),i) = 1; 
    end
    if sum(isinf(INPUT.Xmea(INPUT.MissingInd==0 & (EstimOpt.MissingIndMea(:,i) == 0),i))) > 0
        cprintf(rgb('DarkOrange'),'WARNING: Measurement variable %d contains Inf values - they will be treated as mising. \n', i)
        EstimOpt.MissingIndMea(isinf(INPUT.Xmea(:,i)),i) = 1; 
    end    
    if numel(EstimOpt.MeaSpecMatrix(i) > 0) > 0
        if EstimOpt.MeaSpecMatrix(i) > 0 && numel(unique(INPUT.Xmea(INPUT.MissingInd == 0 & (EstimOpt.MissingIndMea(:,i) == 0),i))) > 10
            cprintf(rgb('DarkOrange'),'WARNING: There are over 10 levels for measurement variable %d \n',i)
        end
    end
    if numel(EstimOpt.MeaSpecMatrix(i) > 0) > 0
        if EstimOpt.MeaSpecMatrix(i) == 1 % MNL
            INPUT.Xmea(INPUT.MissingInd == 0 & (EstimOpt.MissingIndMea(:,i) == 0),i) = INPUT.Xmea(INPUT.MissingInd == 0 & (EstimOpt.MissingIndMea(:,i) == 0),i) - min(unique(INPUT.Xmea(INPUT.MissingInd == 0 & (EstimOpt.MissingIndMea(:,i) == 0),i))) + 1; % make unique values positive (avoid problems with dummyvar)
        end
    end
end

% for i=1:size(EstimOpt.NVarStr,1)
%     if sum(isnan(INPUT.Xstr(:,i))) > 0
%         cprintf(rgb('DarkOrange'),'WARNING: Structural variable %d contains NaN values \n',i)
%     end
% 	if sum(isinf(INPUT.Xstr(:,i))) > 0
%         cprintf(rgb('DarkOrange'),'WARNING: Structural variable %d contains Inf values \n',i)
% 	end
% end

for i = 1:EstimOpt.NVarStr
    if sum(isnan(INPUT.Xstr(INPUT.MissingInd == 0,i))) > 0
        cprintf(rgb('DarkOrange'),'WARNING: Structural variable %d contains NaN values \n',i)
    end
    if sum(isinf(INPUT.Xstr(INPUT.MissingInd == 0,i))) > 0
        cprintf(rgb('DarkOrange'),'WARNING: Structural variable %d contains Inf values \n',i)
    end
end

if EstimOpt.NVarMeaExp > 0
    if isfield(EstimOpt,'NamesMeaExp') == 0 || isempty(EstimOpt.NamesMeaExp) || length(EstimOpt.NamesMeaExp) ~= EstimOpt.NVarMeaExp
        EstimOpt.NamesMeaExp = (1:EstimOpt.NVarMeaExp)';
        EstimOpt.NamesMeaExp = cellstr(num2str(EstimOpt.NamesMeaExp));
    elseif size(EstimOpt.NamesMeaExp,1) ~= EstimOpt.NVarMeaExp
        EstimOpt.NamesMeaExp = EstimOpt.NamesMeaExp';
    end
end

if EstimOpt.NLatent > 0
    if ~isfield(EstimOpt,'NamesLV') || isempty(EstimOpt.NamesLV) || length(EstimOpt.NamesLV) ~= EstimOpt.NLatent
        EstimOpt.NamesLV = {};
        for i = 1:EstimOpt.NLatent
            EstimOpt.NamesLV = [EstimOpt.NamesLV;cellfun(@(x)[x num2str(i)],{'LV '},'UniformOutput',0)];
        end
    else
        EstimOpt.NamesLV = EstimOpt.NamesLV(:);
    end
end

EstimOpt.Names = [];% Names of the models
EstimOpt.NVarcut = 0; % no of cutoffs for ordered probits + constants + variances for OLS 
EstimOpt.NVarcut0 = 0; % no of cutoffs for HMNL0
EstimOpt.CutMatrix = zeros(1,size(INPUT.Xmea,2)); 
EstimOpt.NamesMea_tmp = {};
for i = 1:size(INPUT.Xmea,2)
    if EstimOpt.MeaSpecMatrix(i) == 2 %Ordered probit: cutoffs
        EstimOpt.NVarcut = EstimOpt.NVarcut + length(unique(INPUT.Xmea(:,i))) - 1 + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0);
        EstimOpt.CutMatrix(i) = length(unique(INPUT.Xmea(:,i))) - 1 + sum(EstimOpt.MeaMatrix(:,i))+ EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0);
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end    
        if EstimOpt.MeaExpMatrix(i) ~= 0
           EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;EstimOpt.NamesMeaExp];
        end
        for n = 1:(length(unique(INPUT.Xmea(:,i))) - 1)
            EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;cellfun(@(x)[x num2str(n)],{'Cutoff '},'UniformOutput',0)];
        end
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + length(unique(INPUT.Xmea(:,i))) - 1;
        EstimOpt.Names = [EstimOpt.Names,'OP '];
    elseif EstimOpt.MeaSpecMatrix(i) == 0
        EstimOpt.NVarcut = EstimOpt.NVarcut + 2 + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0); %OLS: constant + sigma
        EstimOpt.CutMatrix(i) = 2 + sum(EstimOpt.MeaMatrix(:,i)) + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0);
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + 2;
        EstimOpt.Names = [EstimOpt.Names,'OLS ']; 
        EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;{'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end    
        if EstimOpt.MeaExpMatrix(i) ~= 0
           EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;EstimOpt.NamesMeaExp];
        end
        EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;{'Sigma'}];
    elseif EstimOpt.MeaSpecMatrix(i) == 1 % MNL 
        EstimOpt.NVarcut = EstimOpt.NVarcut + (length(unique(INPUT.Xmea(:,i))) - 2)*sum(EstimOpt.MeaMatrix(:,i)) + (length(unique(INPUT.Xmea(:,i))) - 1)*(1 + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0)); % constants + additional coefficients 
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + length(unique(INPUT.Xmea(:,i))) - 1;
        EstimOpt.Names = [EstimOpt.Names,'MNL '];
        EstimOpt.CutMatrix(i) = (1 + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0) + sum(EstimOpt.MeaMatrix(:,i)))*(length(unique(INPUT.Xmea(:,i))) - 1);
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for j = 1:(length(unique(INPUT.Xmea(:,i)))-1)
            EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;{'Cons.'}];
            for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
                EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
            end
            if EstimOpt.MeaExpMatrix(i) ~= 0
                EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;EstimOpt.NamesMeaExp];
            end
        end
    elseif EstimOpt.MeaSpecMatrix(i) == 3 % Poisson    
        EstimOpt.NVarcut = EstimOpt.NVarcut + 1 + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0); %Poiss: only constant
        EstimOpt.CutMatrix(i) = 1 + sum(EstimOpt.MeaMatrix(:,i)) + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0);
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + 1;
        EstimOpt.Names = [EstimOpt.Names,'POISS ']; 
        EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;{'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end    
        if EstimOpt.MeaExpMatrix(i) ~= 0
           EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;EstimOpt.NamesMeaExp];
        end
    elseif EstimOpt.MeaSpecMatrix(i) == 4 % Negative Binomial   
        EstimOpt.NVarcut = EstimOpt.NVarcut + 2 + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0); %Poiss: only constant
        EstimOpt.CutMatrix(i) = 2 + sum(EstimOpt.MeaMatrix(:,i)) + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0);
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + 2;
        EstimOpt.Names = [EstimOpt.Names,'NB ']; 
        EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;{'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end    
        if EstimOpt.MeaExpMatrix(i) ~= 0
           EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;EstimOpt.NamesMeaExp];
        end
        EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;{'Theta'}];
    elseif EstimOpt.MeaSpecMatrix(i) == 5 % ZIP   
        EstimOpt.NVarcut = EstimOpt.NVarcut + 2 + sum(EstimOpt.MeaMatrix(:,i)) + 2*EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0); %Poiss: only constant
        EstimOpt.CutMatrix(i) = 2 + 2*sum(EstimOpt.MeaMatrix(:,i)) + 2*EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0);
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + 2;
        EstimOpt.Names = [EstimOpt.Names,'ZIP ']; 
        EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;{'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end    
        if EstimOpt.MeaExpMatrix(i) ~= 0
           EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;EstimOpt.NamesMeaExp];
        end
        EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;{'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end    
        if EstimOpt.MeaExpMatrix(i) ~= 0
           EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;EstimOpt.NamesMeaExp];
        end
    elseif EstimOpt.MeaSpecMatrix(i) == 6 % ZINB  
        EstimOpt.NVarcut = EstimOpt.NVarcut + 3 + sum(EstimOpt.MeaMatrix(:,i)) + 2*EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0); %Poiss: only constant
        EstimOpt.CutMatrix(i) = 3 + 2*sum(EstimOpt.MeaMatrix(:,i)) + 2*EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0);
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + 3;
        EstimOpt.Names = [EstimOpt.Names,'ZINB ']; 
        EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;{'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end    
        if EstimOpt.MeaExpMatrix(i) ~= 0
           EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;EstimOpt.NamesMeaExp];
        end
        
        EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;{'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end    
        if EstimOpt.MeaExpMatrix(i) ~= 0
           EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;EstimOpt.NamesMeaExp];
        end
        EstimOpt.NamesMea_tmp = [EstimOpt.NamesMea_tmp;{'Theta'}];
    end
    
end

if isfield(EstimOpt,'NamesStr') == 0 || isempty(EstimOpt.NamesStr) || length(EstimOpt.NamesStr) ~= EstimOpt.NVarStr
    EstimOpt.NamesStr = (1:EstimOpt.NVarStr)';
    EstimOpt.NamesStr = cellstr(num2str(EstimOpt.NamesStr));
elseif size(EstimOpt.NamesStr,1) ~= EstimOpt.NVarStr
    EstimOpt.NamesStr = EstimOpt.NamesStr';
end

if isfield(EstimOpt,'NamesMea') == 0 || isempty(EstimOpt.NamesMea) || length(EstimOpt.NamesMea) ~= size(INPUT.Xmea,2)
    EstimOpt.NamesMea = (1:size(INPUT.Xmea,2))';
    EstimOpt.NamesMea = cellstr(num2str(EstimOpt.NamesMea));
elseif size(EstimOpt.NamesMea,1) ~= size(INPUT.Xmea,2)
    EstimOpt.NamesMea = EstimOpt.NamesMea';
end

disp(['Using following models for attitudes: ' char(EstimOpt.Names)]);


%% Rescructure data


% INPUT.Xstr = double(INPUT.Xstr(1:EstimOpt.NAlt*EstimOpt.NCT:end,:));
% if isfield(EstimOpt,'StrNorm')
%     if length(EstimOpt.StrNorm) ~= EstimOpt.NVarStr && length(EstimOpt.StrNorm) ~= 1
%         cprintf(rgb('DarkOrange'),'WARNING: Structural variables normalization options (EstimOpt.StrNorm) incorrect - variables not normalized \n')
%     else
%         if length(EstimOpt.StrNorm) == 1
%             EstimOpt.StrNorm = EstimOpt.StrNorm(ones(EstimOpt.NVarStr,1),1);
%         end       
%         INPUT.Xstr(:,EstimOpt.StrNorm == 1) = (INPUT.Xstr(:,EstimOpt.StrNorm == 1) - mean(INPUT.Xstr))./std(INPUT.Xstr);
%     end
% end

INPUT.Xstr = double(INPUT.Xstr(1:EstimOpt.NAlt*EstimOpt.NCT:end,:));
% normalize explanatory variables for structural equations:
if ~isfield(EstimOpt,'StrNorm')
    EstimOpt.StrNorm = ones(EstimOpt.NVarStr,1);
elseif length(EstimOpt.StrNorm) ~= EstimOpt.NVarStr && length(EstimOpt.StrNorm) ~= 1
    cprintf(rgb('DarkOrange'),'WARNING: Structural variables normalization options (EstimOpt.StrNorm) incorrect - variables not normalized \n')
elseif length(EstimOpt.StrNorm) == 1
    EstimOpt.StrNorm = EstimOpt.StrNorm(ones(EstimOpt.NVarStr,1),1);
end
if any(EstimOpt.StrNorm > 0)
    meanx = mean(INPUT.Xstr);
    stdx = std(INPUT.Xstr);
    INPUT.Xstr(:,EstimOpt.StrNorm == 1) = (INPUT.Xstr(:,EstimOpt.StrNorm == 1) - meanx(:,EstimOpt.StrNorm == 1))./stdx(:,EstimOpt.StrNorm == 1);
end

INPUT.Xmea = INPUT.Xmea(1:EstimOpt.NAlt*EstimOpt.NCT:end,:);
EstimOpt.MissingIndMea = EstimOpt.MissingIndMea(1:EstimOpt.NAlt*EstimOpt.NCT:end,:);

% if EstimOpt.NVarMeaExp > 0
%     INPUT.Xmea_exp = INPUT.Xmea_exp(1:EstimOpt.NAlt*EstimOpt.NCT:end,:);
% end
if EstimOpt.NVarMeaExp > 0
    INPUT.Xmea_exp = INPUT.Xmea_exp(1:EstimOpt.NAlt*EstimOpt.NCT:end,:);
    % normalize explanatory variables for measurement equations:
    if ~isfield(EstimOpt,'MeaExpNorm')
        EstimOpt.MeaExpNorm = ones(EstimOpt.NVarMeaExp,1);
    elseif length(EstimOpt.MeaExpNorm) ~= EstimOpt.NVarMeaExp && length(EstimOpt.MeaExpNorm) ~= 1
        cprintf(rgb('DarkOrange'),'WARNING: Explanatory measurement variables normalization options (EstimOpt.MeaExpNorm) incorrect - variables not normalized \n')
    elseif length(EstimOpt.MeaExpNorm) == 1
        EstimOpt.MeaExpNorm = EstimOpt.MeaExpNorm(ones(EstimOpt.NVarMeaExp,1),1);
    end
    if any(EstimOpt.MeaExpNorm > 0)
        meanx = mean(INPUT.Xmea_exp);
        stdx = std(INPUT.Xmea_exp);
        INPUT.Xmea_exp(:,EstimOpt.MeaExpNorm == 1) = (INPUT.Xmea_exp(:,EstimOpt.MeaExpNorm == 1) - meanx(:,EstimOpt.MeaExpNorm == 1))./stdx(:,EstimOpt.MeaExpNorm == 1);
    end
end


%% Estimating MIMIC0


OptimOpt_0 = optimoptions('fminunc');
OptimOpt_0.Algorithm = 'quasi-newton';
OptimOpt_0.GradObj = 'off';
OptimOpt_0.Hessian = 'off';
OptimOpt_0.Display = 'off';
OptimOpt_0.FunValCheck = 'off'; 
OptimOpt_0.Diagnostics = 'off';

LLfun0 = @(B) LL_hmnl0(INPUT.Xmea,EstimOpt,B);
[Results.MIMIC0.bhat,LL0] = fminunc(LLfun0,0.001*ones(EstimOpt.NVarcut0,1),OptimOpt_0);
Results.MIMIC0.LL = -LL0;


%% Starting values

if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut
    b0 = B_backup(:);
    disp('Using the starting values from Backup')
elseif isfield(Results_old,'MIMIC') && isfield(Results_old.MIMIC,'b0') % starting values provided
    Results_old.MIMIC.b0_old = Results_old.MIMIC.b0(:);
    Results_old.MIMIC = rmfield(Results_old.MIMIC,'b0');
    if length(Results_old.MIMIC.b0_old) ~= EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut
        cprintf(rgb('DarkOrange'),'WARNING: Incorrect no. of starting values or model specification \n')
        Results_old.MIMIC = rmfield(Results_old.MIMIC,'b0_old');
    else
        b0 = Results_old.MIMIC.b0_old(:);
    end
end
if ~exist('b0','var')
    b0 = zeros(EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut,1);
    l = EstimOpt.NVarStr*EstimOpt.NLatent;
    k = 0;
    for i = 1:size(INPUT.Xmea,2)
        if EstimOpt.MeaSpecMatrix(i) == 2 %Ordered probit: cutoffs
            b0(l+sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp+1:l+sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp+length(unique(INPUT.Xmea(:,i)))-1) = Results.MIMIC0.bhat(k+1:k+length(unique(INPUT.Xmea(:,i)))-1);            
            l = l + sum(EstimOpt.MeaMatrix(:,i),1) + EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp + length(unique(INPUT.Xmea(:,i))) - 1;
            k = k + length(unique(INPUT.Xmea(:,i))) - 1;
        elseif EstimOpt.MeaSpecMatrix(i) == 0
            b0(l+1) = Results.MIMIC0.bhat(k+1);
            b0(l+2+sum(EstimOpt.MeaMatrix(:,i)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp,1)) = Results.MIMIC0.bhat(k+2);
            k = k + 2;
            l = l + 2 + sum(EstimOpt.MeaMatrix(:,i),1) + EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp;        
        elseif EstimOpt.MeaSpecMatrix(i) == 1 % MNL 
            for n = 1:(length(unique(INPUT.Xmea(:,i)))-1)
                b0(l+1) = Results.MIMIC0.bhat(k+1);
                l = l + 1 + sum(EstimOpt.MeaMatrix(:,i)) + EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp;
                k = k + 1;
            end
        elseif EstimOpt.MeaSpecMatrix(i) == 3 % POISS 
            b0(l+1) = Results.MIMIC0.bhat(k+1);
            k = k + 1;
            l = l + 1 + sum(EstimOpt.MeaMatrix(:,i),1) + EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp;     
        elseif EstimOpt.MeaSpecMatrix(i) == 4 % NB
            b0(l+1) = Results.MIMIC0.bhat(k+1);
            b0(l+2+sum(EstimOpt.MeaMatrix(:,i)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp,1)) = Results.MIMIC0.bhat(k+2); % theta
            k = k + 2;
            l = l + 2 + sum(EstimOpt.MeaMatrix(:,i),1) + EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp;  
        elseif EstimOpt.MeaSpecMatrix(i) == 5 % ZIP
            b0(l+1) = Results.MIMIC0.bhat(k+1);
            b0(l+2+sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp) = Results.MIMIC0.bhat(k+2); % second constant
            k = k + 2;
            l = l + 2 + 2*sum(EstimOpt.MeaMatrix(:,i),1) + 2*EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp;     
        elseif EstimOpt.MeaSpecMatrix(i) == 6 % ZINB
            b0(l+1) = Results.MIMIC0.bhat(k+1);
            b0(l+2+sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp) = Results.MIMIC0.bhat(k+2); % Second constant
            b0(l+3+2*sum(EstimOpt.MeaMatrix(:,i),1)+2*EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp) = Results.MIMIC0.bhat(k+3); % Theta 
            k = k + 3;
            l = l + 3 + 2*sum(EstimOpt.MeaMatrix(:,i),1) + 2*EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp;         
        end
    end
end

%% Optimization Options

if  isfield(EstimOpt,'BActive')
	EstimOpt.BActive = EstimOpt.BActive(:)';
end

if isfield(EstimOpt,'StrMatrix') && size(EstimOpt.StrMatrix,1) == EstimOpt.NLatent && size(EstimOpt.StrMatrix,2) == EstimOpt.NVarStr && any(any(EstimOpt.StrMatrix ~= 1))
    if ~isfield(EstimOpt,'BActive')
        EstimOpt.BActive = ones(1,EstimOpt.NVarStr*EstimOpt.NLatent+EstimOpt.NVarMea+EstimOpt.NVarcut);
    end
    EstimOpt.BActive = StrSelect(EstimOpt,0);
end

if EstimOpt.ConstVarActive == 1
    if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive) || sum(EstimOpt.BActive == 0) == 0
        error ('Are there any constraints on model parameters (EstimOpt.ConstVarActive)? Constraints not provided (EstimOpt.BActive).')
    elseif length(b0) ~= length(EstimOpt.BActive)
        error('Check no. of constraints')
    end
    disp(['Initial values: ' mat2str(b0',2)])
    disp(['Parameters with zeros are constrained to their initial values: ' mat2str(EstimOpt.BActive')])
else    
    if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive) || sum(EstimOpt.BActive == 0) == 0
        EstimOpt.BActive = ones(1,length(b0));
        disp(['Initial values: ' mat2str(b0',2)])
    else        
        if length(b0) ~= length(EstimOpt.BActive)
            error('Check no. of constraints')
        else
            disp(['Initial values: ' mat2str(b0',2)])
            disp(['Parameters with zeros are constrained to their initial values: ' mat2str(EstimOpt.BActive')])
        end
    end
end

%% Generate pseudo-random draws

if isfield(EstimOpt,'Seed1') == 1
	rng(EstimOpt.Seed1);
end

cprintf('Simulation with ');
cprintf('*blue',[num2str(EstimOpt.NRep) ' ']); 

if EstimOpt.Draws == 1
    cprintf('*blue','Pseudo-random '); cprintf('draws \n');
    err_mtx = randn(EstimOpt.NP*EstimOpt.NRep,EstimOpt.NLatent+EstimOpt.NVarA+1); %to be cut down later   
elseif EstimOpt.Draws == 2 % LHS
    cprintf('*blue','Latin Hypercube Sampling '); cprintf('draws \n');
    err_mtx = lhsnorm(zeros((EstimOpt.NLatent+EstimOpt.NVarA+1)*EstimOpt.NP,1),diag(ones((EstimOpt.NLatent+EstimOpt.NVarA+1)*EstimOpt.NP,1)),EstimOpt.NRep); 
    err_mtx = reshape(err_mtx,[EstimOpt.NRep*EstimOpt.NP,EstimOpt.NLatent]);
elseif EstimOpt.Draws >= 3 % Quasi random draws 
    if EstimOpt.Draws == 3
        cprintf('*blue','Halton '); cprintf(['draws (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n']);
        hm1 = haltonset(EstimOpt.NLatent+EstimOpt.NVarA+1,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); % 
    elseif EstimOpt.Draws == 4 % apply reverse-radix scrambling
        cprintf('*blue','Halton '); cprintf(['draws with reverse radix scrambling (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n']);
        hm1 = haltonset(EstimOpt.NLatent+EstimOpt.NVarA+1,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); % 
        hm1 = scramble(hm1,'RR2');
    elseif EstimOpt.Draws == 5
        cprintf('*blue','Sobol '); cprintf(['draws (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n']);
        hm1 = sobolset(EstimOpt.NLatent+EstimOpt.NVarA+1,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); 
    elseif EstimOpt.Draws == 6
        cprintf('*blue','Sobol '); cprintf(['draws with random linear scramble and random digital shift (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n']);
        hm1 = sobolset(EstimOpt.NLatent+EstimOpt.NVarA+1,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); 
        hm1 = scramble(hm1,'MatousekAffineOwen');
    end
    
    err_mtx = net(hm1,EstimOpt.NP*EstimOpt.NRep); % this takes every point:
    clear hm1;
	err_mtx = err_mtx(:,2+EstimOpt.NVarA:end);
    if EstimOpt.NP*EstimOpt.NRep < 3e+7
        err_mtx= icdf('Normal',err_mtx,0,1); %to be cut down later
    else % this is for very large number of draws * variables
        for i=1:EstimOpt.NLatent
            err_mtx(:,i) = icdf('Normal',err_mtx(:,i),0,1); %to be cut down later
        end
    end    
end

err_sliced = err_mtx'; % NLatent x NRep * NP

if isfield(EstimOpt,'Drawskeep') && ~isempty(EstimOpt.Drawskeep) && EstimOpt.Drawskeep == 1
    Results.err = err_sliced;
end

%% Display Options

if ((isfield(EstimOpt,'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    OptimOpt.GradObj = 'on';
end

if any(EstimOpt.MissingAlt(:) == 1) && EstimOpt.NumGrad == 0
   EstimOpt.NumGrad = 1;
   cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient to numerical - missing alternatives not supported by analytical gradient \n')
end

if any(EstimOpt.MeaSpecMatrix >= 3) && EstimOpt.NumGrad == 0 && any(any(INPUT.Xmea(:,EstimOpt.MeaSpecMatrix >= 3) > 100))
   cprintf(rgb('DarkOrange'),'WARNING: it is recommended to switch to numerical gradient, as analitycal can be not precise when Xmea take large values for NB \n')
end

% if any(EstimOpt.MeaSpecMatrix == 5) && EstimOpt.NumGrad == 0 
%    EstimOpt.NumGrad = 1;
%    cprintf(rgb('DarkOrange'), 'WARNING: ZIP\n')
% end

if (isfield(EstimOpt,'ConstVarActive') == 0 || EstimOpt.ConstVarActive == 0) && isequal(OptimOpt.Algorithm,'quasi-newton') && isequal(OptimOpt.Hessian,'user-supplied')
	cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied Hessian off - quasi-newton algorithm does not use it anyway \n')
    OptimOpt.Hessian = 'off';
end
if EstimOpt.RobustStd == 1 && (EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2)
    EstimOpt.RobustStd = 0; 
    cprintf(rgb('DarkOrange'),'WARNING: Setting off robust standard errors, they do not matter for BHHH aproximation of hessian \n')
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

%% Estimation

LLfun = @(B) LL_mimic_MATlike(INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,INPUT.W,EstimOpt,OptimOpt,B);
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

if EstimOpt.HessEstFix == 1
	f = LL_mimic(INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,EstimOpt,Results.bhat);
	Results.jacobian = numdiff(@(B) INPUT.W.*LL_mimic(INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,EstimOpt,B),INPUT.W.*f,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
	Results.jacobian = jacobianest(@(B) INPUT.W.*LL_mimic(INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
	Results.jacobian = hessian(@(B) sum(INPUT.W.*LL_mimic(INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,EstimOpt,B),1),Results.bhat);
end
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
           [~,Results.jacobian] = LL_mimic(INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,EstimOpt,Results.bhat);
           Results.jacobian = Results.jacobian.*INPUT.W;
    else
        Results.LLdetailed = LL_mimic(INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,EstimOpt,Results.bhat);
        Results.jacobian = numdiff(@(B) INPUT.W.*LL_mimic(INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,EstimOpt,B),INPUT.W.*Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
    end
    RobustHess = Results.jacobian'*Results.jacobian;
    Results.ihess = Results.ihess*RobustHess*Results.ihess;
end

Results.std = sqrt(diag(Results.ihess));
Results.std(EstimOpt.BActive == 0) = NaN;
Results.std(EstimOpt.BLimit == 1) = 0;

for i = 1:EstimOpt.NLatent
    Results.DetailsS(:,4*i-3) = Results.bhat((EstimOpt.NVarStr)*(i-1)+1:(EstimOpt.NVarStr)*i);
    Results.DetailsS(:,4*i-1:4*i) = [Results.std((EstimOpt.NVarStr)*(i-1)+1:(EstimOpt.NVarStr)*i),pv(Results.bhat((EstimOpt.NVarStr)*(i-1)+1:(EstimOpt.NVarStr)*i),Results.std((EstimOpt.NVarStr)*(i-1)+1:(EstimOpt.NVarStr)*i))];
end

Results.DetailsM(:,1) = Results.bhat((EstimOpt.NVarStr)*EstimOpt.NLatent+1:end);
Results.DetailsM(:,3:4) = [Results.std((EstimOpt.NVarStr)*EstimOpt.NLatent+1:end),pv(Results.bhat((EstimOpt.NVarStr)*EstimOpt.NLatent+1:end),Results.std((EstimOpt.NVarStr)*EstimOpt.NLatent+1:end))];

EstimOpt.params = length(b0);
EstimOpt.NObs = EstimOpt.NP;
if isfield(EstimOpt,'BActive')
	EstimOpt.params = EstimOpt.params - sum(EstimOpt.BActive == 0);
end
Results.stats = [Results.LL;Results.MIMIC0.LL;1 - Results.LL/Results.MIMIC0.LL;NaN;((2*EstimOpt.params - 2*Results.LL))/EstimOpt.NObs;((log(EstimOpt.NObs)*EstimOpt.params-2*Results.LL))/EstimOpt.NObs;EstimOpt.NObs;EstimOpt.NP;EstimOpt.params];
Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;
Results.INPUT = INPUT;

%% Tworzenie naglowka

Head = cell(1,2);
Head(1,1) = {'MIMIC'};
if EstimOpt.NVarS > 0 
    Template1 = {'DetailsS'};
    Template2 = {'DetailsS'};
    Names.DetailsS = EstimOpt.NamesStr;
    Heads.DetailsS(:,2) = EstimOpt.NamesLV;
    Heads.DetailsS(1:2,1) = {'Structural equations';'lc'};
    Heads.DetailsS(EstimOpt.NLatent+1,2) = {'lb'};
    ST = {'DetailsS'};
else
    Template1 = {};
    Template2 = {};
    ST = {};
end

k = 0;

if any(EstimOpt.MeaSpecMatrix == 2)
    Results.DetailsOP = [];
end

for i = 1:size(INPUT.Xmea,2)
    model = model_name(EstimOpt.MeaSpecMatrix(i));
    Heads.(strcat('Xmea',num2str(i)))(1:2,1) = [{['Measurment equation for:',' ',char(EstimOpt.NamesMea(i)),' (',model,')']};'lb'];
    if EstimOpt.MeaSpecMatrix(i) == 2
        l = k + EstimOpt.CutMatrix(i) - length(unique(INPUT.Xmea(:,i))) + 1;
        Results.DetailsOP = vertcat(Results.DetailsOP,Results.DetailsM(l+1:l+length(unique(INPUT.Xmea(:,i)))-1,:));
        if length(unique(INPUT.Xmea(:,i))) > 2 % if attitude is not binary
            g = [Results.DetailsM(l+1,1);exp(Results.DetailsM(l+2:l+length(unique(INPUT.Xmea(:,i)))-1,1))];
            for n = 2:length(unique(INPUT.Xmea(:,i)))-1
                stdx = sqrt(g(1:n)'*Results.ihess((EstimOpt.NVarStr)*EstimOpt.NLatent+l+1:(EstimOpt.NVarStr)*EstimOpt.NLatent+l+n,(EstimOpt.NVarStr)*EstimOpt.NLatent+l+1:(EstimOpt.NVarStr)*EstimOpt.NLatent+l+n)*g(1:n));
                Results.DetailsM(l+n,1) = sum(g(1:n),1);
                Results.DetailsM(l+n,3:4) = [stdx,pv(sum(g(1:n),1),stdx)];
            end
        end
    end
    if EstimOpt.MeaSpecMatrix(i) == 5 || EstimOpt.MeaSpecMatrix(i) == 6
        Heads.(strcat('Xmea',num2str(i)))(1:2,2) = [{'Probability of Non-participation (logit)'};{'lb'}];
        Results.(strcat('Xmea',num2str(i)))(:,1:4) = Results.DetailsM(k+1:k+floor(EstimOpt.CutMatrix(i)/2),:);
        Results.(strcat('Xmea2',num2str(i)))(:,1:4) = Results.DetailsM(k+floor(EstimOpt.CutMatrix(i)/2)+1:k+EstimOpt.CutMatrix(i),:);
        Names.(strcat('Xmea',num2str(i))) = EstimOpt.NamesMea_tmp(k+1:k+floor(EstimOpt.CutMatrix(i)/2));
        Names.(strcat('Xmea2',num2str(i))) = EstimOpt.NamesMea_tmp(k+floor(EstimOpt.CutMatrix(i)/2)+1:k+EstimOpt.CutMatrix(i));
        if EstimOpt.MeaSpecMatrix(i) == 5
            Heads.(strcat('Xmea2',num2str(i))) = [{'Poisson model'};{'lb'}];
        else %if EstimOpt.MeaSpecMatrix(i) == 6
            Heads.(strcat('Xmea2',num2str(i))) = [{'Negative binomial model'};{'lb'}];
        end
        Temp1 = cell(2,size(Template1,2));
        Temp1(1,1) = {strcat('Xmea',num2str(i))};
        Temp1(2,1) = {strcat('Xmea2',num2str(i))};
        Template1 = [Template1;Temp1]; %#ok<AGROW>
        Temp2 = cell(2,size(Template2,2));
        Temp2(1,1) = {strcat('Xmea',num2str(i))};
        Temp2(2,1) = {strcat('Xmea2',num2str(i))};
        Template2 = [Template2;Temp2]; %#ok<AGROW>
        ST = [ST,{strcat('Xmea',num2str(i))},{strcat('Xmea2',num2str(i))}]; %#ok<AGROW>
    else
        Results.(strcat('Xmea',num2str(i)))(1:EstimOpt.CutMatrix(i),1:4) = Results.DetailsM(k+1:k+EstimOpt.CutMatrix(i),:);
        Names.(strcat('Xmea',num2str(i))) = EstimOpt.NamesMea_tmp(k+1:k+EstimOpt.CutMatrix(i));
        Temp1 = cell(1,size(Template1,2));
        Temp1(1,1) = {strcat('Xmea',num2str(i))};
        Template1 = [Template1;Temp1]; %#ok<AGROW>
        Temp2 = cell(1,size(Template2,2));
        Temp2(1,1) = {strcat('Xmea',num2str(i))};
        Template2 = [Template2;Temp2]; %#ok<AGROW>
        ST = [ST,{strcat('Xmea',num2str(i))}]; %#ok<AGROW>
    end
    k = k + EstimOpt.CutMatrix(i);
end

%% Tworzenie stopki

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

%% Tworzenie ResultsOut, drukowanie na ekran i do pliku .xls

if EstimOpt.Display ~= 0
    Results.Dist = zeros(size(INPUT.Xa,2),1);
    Results.R_out = genOutput(EstimOpt,Results,Head,Tail,Names,Template1,Template2,Heads,ST,'MIMIC');
end

end