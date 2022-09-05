function Results = MMDCEV(INPUT,Results_old, EstimOpt,OptimOpt)
% MDCEV creates Mixed Multiple Discrete-Continuous Extreme Value model.
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
% General basics:
% �	DataFile � path/name of the .mat data file
% �	Display � 1; shows output, set to 0 to hide it 
% �	ProjectName � Name of the project/model
% �	NCT - Number of choice tasks per person 
% �	NAlt - Number of alternatives
% �	NP � Number of respondents
% 
% 
% Variables options:
% �	NamesA � Names of variables in list e.g. {'-Opt out';�-Cost (EUR)'}
% 
% 
% Modelling options from DataCleanDCE:
% �	ApproxHess = 1; for user supplied hessians, 1 for BHHH, 0 for analytical
% �	NumGrad = 0; uses analytical gradient in calculations, set to 1 for numerical
%    gradient (only NumGrad=1 at the moment)
% �	HessEstFix = 0; Options: 
% o	0 - use optimization Hessian, 
% o	1 - use jacobian-based (BHHH) Hessian, 
% o	2 - use high-precision jacobian-based (BHHH) Hessian,
% o	3 - use numerical Hessian, 
%
% Example: 
%    Results.MDCEV = MDCEV(INPUT,Results,EstimOpt,OptimOpt);
%
% Author: Mikolaj Czajkowski, Professor
% University of Warsaw, Faculty of Economic Sciences
% email address: mik@czaj.org 
% Website: http://czaj.org/#

% save tmp_MNL
% return

global B_backup

tic

Results.bhat = [];
Results.R = [];
Results.R_out = {};
Results.stats = [];

EstimOpt.NVarA = size(INPUT.Xa,2);
NVarA = EstimOpt.NVarA;


%% Check data and inputs

if nargin < 3
    error('Too few input arguments for MMDCEV(INPUT,EstimOpt,OptimOpt)')
end

warning off MATLAB:mir_warning_maybe_uninitialized_temporary

format shortG;
format compact;

if isfield(EstimOpt,'Display') == 0
    EstimOpt.Display = 1;
end

% if isfield(EstimOpt,'WTP_space') == 0
%     EstimOpt.WTP_space = 0;
%     EstimOpt.WTP_matrix = [];
% elseif EstimOpt.WTP_space == 0
%     EstimOpt.WTP_matrix = [];
% end

if EstimOpt.Display ~= 0
    disp(' ');
    disp('__________________________________________________________________________________________________________________');
    disp(' ');
    
    disp('Estimating Mixed MDCEV model ...')
    
    disp('in preference-space ...')
end
if isfield(EstimOpt,'FullCov') == 0
    EstimOpt.FullCov = 0;
end

if isfield(EstimOpt,'NRep') == 0
    EstimOpt.NRep = 1000;
    cprintf(rgb('DarkOrange'), 'WARNING: Setting no. of draws (EstimOpt.NRep) to 1000. \n')
end
if isfield(EstimOpt,'Draws') == 0
    EstimOpt.Draws = 6;
end
if isfield(EstimOpt,'HaltonSkip') == 0
    EstimOpt.HaltonSkip = 1; % specify no of rows in halton sequence to skip (default=1)
end
if isfield(EstimOpt,'HaltonLeap') == 0
    EstimOpt.HaltonLeap = 0; % specify no of rows in halton sequence to leap (default=0)
end
if isfield(EstimOpt,'Seed1') == 0
    EstimOpt.Seed1 = 179424673;
end
if isfield(EstimOpt,'RealMin') == 0
    EstimOpt.RealMin = 0;
end
if isfield(EstimOpt,'NSdSim') == 0
    EstimOpt.NSdSim = 10000;
end
if isfield(EstimOpt,'SpecProfile') == 0 || isempty(EstimOpt.SpecProfile) || size(EstimOpt.SpecProfile,1) ~= 2 || size(EstimOpt.SpecProfile,2) ~= EstimOpt.NAlt
    if isfield(EstimOpt,'Profile') == 0 || isempty(EstimOpt.Profile) || ~ismember(EstimOpt.Profile, [1, 2])
        if EstimOpt.Display ~= 0
            cprintf(rgb('DarkOrange'), 'WARNING: Assuming alpha-profile. For gamma-profile specify EstimOpt.Profile = 2 or use EstimOpt.SpecProfile.\n')
        end
        EstimOpt.Profile = 1;
        EstimOpt.SpecProfile = [1:EstimOpt.NAlt; zeros(1,EstimOpt.NAlt)]; % Alpha profile
    else
        if EstimOpt.Profile == 1
            EstimOpt.SpecProfile = [1:EstimOpt.NAlt; zeros(1,EstimOpt.NAlt)];
        else
            EstimOpt.SpecProfile = [zeros(1,EstimOpt.NAlt); 1:EstimOpt.NAlt];
        end
    end
else
    disp('Using profile specification from EstimOpt.SpecProfile')
    if sum(EstimOpt.SpecProfile(1,:)) ~= 0 && sum(EstimOpt.SpecProfile(2,:)) ~= 0 % Both
        EstimOpt.Profile = 3; 
    elseif sum(EstimOpt.SpecProfile(1,:)) ~= 0 % Alpha
        EstimOpt.Profile = 1; 
    elseif sum(EstimOpt.SpecProfile(2,:)) ~= 0 % Gamma
        EstimOpt.Profile = 2; 
    end
end

if isfield(EstimOpt,'Dist') == 0 || isempty(EstimOpt.Dist)
    EstimOpt.Dist = zeros(1,EstimOpt.NVarA);
    cprintf(rgb('DarkOrange'),'WARNING: distributions for random parameters not specified - assuming normality \n')
else
    if length(EstimOpt.Dist) == 1
        EstimOpt.Dist = EstimOpt.Dist.*ones(1,EstimOpt.NVarA); % needed?
    elseif length(EstimOpt.Dist) == EstimOpt.NVarA
        EstimOpt.Dist = EstimOpt.Dist(:)';
    else
        error('Incorrect no. of random parameters'' distributions provided')
    end
end
disp(['Random parameters distributions: ',num2str(EstimOpt.Dist),' (-1 - constant, 0 - normal, 1 - lognormal)'])

if isfield(INPUT,'Xm') == 0
    INPUT.Xm = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarM = size(INPUT.Xm,2); % Number of covariates of means of random parameters

if isfield(EstimOpt,'NamesA') == 0 || isempty(EstimOpt.NamesA) || length(EstimOpt.NamesA) ~= EstimOpt.NVarA
    EstimOpt.NamesA = (1:EstimOpt.NVarA)';
    EstimOpt.NamesA = cellstr(num2str(EstimOpt.NamesA));
elseif size(EstimOpt.NamesA,1) ~= EstimOpt.NVarA
    EstimOpt.NamesA = EstimOpt.NamesA';
end

if isfield(EstimOpt,'NamesM') == 0 || isempty(EstimOpt.NamesM) || length(EstimOpt.NamesM) ~= EstimOpt.NVarM
    EstimOpt.NamesM = (1:EstimOpt.NVarM)';
    EstimOpt.NamesM = cellstr(num2str(EstimOpt.NamesM));
elseif size(EstimOpt.NamesM,1) ~= EstimOpt.NVarM
    EstimOpt.NamesM = EstimOpt.NamesM';
end

% Alpha and gamma parameters
EstimOpt.NVarP = length(unique(EstimOpt.SpecProfile(1,EstimOpt.SpecProfile(1,:) ~= 0))) + length(unique(EstimOpt.SpecProfile(2,EstimOpt.SpecProfile(2,:) ~= 0)));

if isfield(EstimOpt,'NamesAlt') == 0 || isempty(EstimOpt.NamesAlt) || length(EstimOpt.NamesAlt) ~= EstimOpt.NVarP
    EstimOpt.NamesAlt = (1:EstimOpt.NVarP)';
    EstimOpt.NamesAlt = cellstr(num2str(EstimOpt.NamesAlt));
elseif size(EstimOpt.NamesAlt,1) ~= EstimOpt.NVarP
    EstimOpt.NamesAlt = EstimOpt.NamesAlt';
end

%% Starting values
if EstimOpt.FullCov == 0
    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == 2*EstimOpt.NVarA + EstimOpt.NVarP*(1+EstimOpt.NVarM) + 1
        b0 = B_backup(:);
        if EstimOpt.Display ~= 0
            disp('Using the starting values from Backup')
        end
    end
    if ~exist('b0','var')
        if (isfield(Results_old,'MDCEV') && isfield(Results_old.MDCEV,'bhat') && length(Results_old.MDCEV.bhat) == (EstimOpt.NVarA + EstimOpt.NVarP*(1+EstimOpt.NVarM) + 1))
            disp('Using MDCEV for starting values')
            b0 = [Results_old.MDCEV.bhat(1:EstimOpt.NVarA); 0.1*ones(EstimOpt.NVarA,1); Results_old.MDCEV.bhat(EstimOpt.NVarA+1:end)];
        else
            error('No starting values - Estimate MDCEV first')
        end   
    end
else
    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVarA + sum(1:EstimOpt.NVarA) + EstimOpt.NVarP*(1+EstimOpt.NVarM) + 1
        b0 = B_backup(:);
        if EstimOpt.Display ~= 0
            disp('Using the starting values from Backup')
        end
    end
    if ~exist('b0','var')
        if (isfield(Results_old,'MMDCEV_d') && isfield(Results_old.MMDCEV_d,'bhat') && length(Results_old.MMDCEV_d.bhat) == (2*EstimOpt.NVarA + EstimOpt.NVarP*(1+EstimOpt.NVarM) + 1))
            disp('Using MMDCEV_d for starting values')
            vc_tmp = abs(diag(Results_old.MMDCEV_d.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2)));
            b0 = [Results_old.MMDCEV_d.bhat(1:EstimOpt.NVarA); vc_tmp(tril(ones(size(vc_tmp))) == 1); Results_old.MMDCEV_d.bhat(2*EstimOpt.NVarA+1:end)];
        else
            error('No starting values - Estimate MMDCEV_d first')
        end   
    end
end


%% Optimization Options

if isfield(EstimOpt,'BActive')
    EstimOpt.BActive = EstimOpt.BActive(:)';
else
    EstimOpt.BActive = ones(length(b0),1);
end

if sum(EstimOpt.Dist == -1) > 0
%     if isfield(EstimOpt,'BActive') == 0 || isempty(EstimOpt.BActive)
%         EstimOpt.BActive = ones(1,length(b0));
%     end
    if EstimOpt.FullCov == 0
        EstimOpt.BActive(NVarA+find(EstimOpt.Dist == -1)) = 0;
    elseif EstimOpt.FullCov == 1
        Vt = tril(ones(NVarA));
        Vt(EstimOpt.Dist == -1,:) = 0;
        EstimOpt.BActive(NVarA+1:NVarA+sum(1:NVarA)) = EstimOpt.BActive(NVarA+1:NVarA+sum(1:NVarA)).*(Vt(tril(ones(size(Vt)))~=0)');
    end
end

if isfield(EstimOpt,'ConstVarActive') == 0 
    EstimOpt.ConstVarActive = 0;
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


if isfield(EstimOpt,'HessEstFix') == 0 
    EstimOpt.HessEstFix = 0;
end 
INPUT.W = ones(EstimOpt.NP,1); % Weights not supported for now

if isfield(EstimOpt,'NumGrad') == 0 || EstimOpt.NumGrad ~= 1
    EstimOpt.NumGrad = 1;
    cprintf(rgb('DarkOrange'), 'WARNING: Setting gradient to numerical. Analytical gradient is not supported.\n')
end
if EstimOpt.Display ~= 0
    fprintf('\n')
    cprintf('Opmization algorithm: '); cprintf('*Black',[OptimOpt.Algorithm '\n'])
    if strcmp(OptimOpt.GradObj,'on')
        if EstimOpt.NumGrad == 0
            cprintf('*Red', 'Analytical gradient not supported at the moment.')
            return
            % cprintf('Gradient: '); cprintf('*Black','user-supplied, analytical \n')
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
    fprintf('\n')
end

%% Generate pseudo-random draws from standard normal distribution
err_mtx = generateRandomDraws(EstimOpt);

% set column value to 0 for constant parameters
err_mtx(:,EstimOpt.Dist == -1) = 0;
INPUT.err = err_mtx';
%% Restucturing Data 

INPUT.Y = reshape(INPUT.Y, [EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP]);
INPUT.priceMat = reshape(INPUT.priceMat, [EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP]);

INPUT.Xa = reshape(INPUT.Xa,[EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP, EstimOpt.NVarA]);
INPUT.Xa = permute(INPUT.Xa,[1 3 2]);
if EstimOpt.NVarM > 0
    INPUT.Xm = INPUT.Xm(1:EstimOpt.NAlt:end,:);
    INPUT.Xm = [ones(size(INPUT.Xm,1),1), INPUT.Xm];
end


%% Estimation

% =========================================================================
% INPUT must contain Xa, priceMat, y (dependent variable; demands) ========
% =========================================================================

LLfun = @(B) LL_mmdcev_MATlike(INPUT, EstimOpt,OptimOpt,B);

if EstimOpt.ConstVarActive == 0
    if EstimOpt.HessEstFix == 0
        [Results.bhat,LL,Results.exitf,Results.output,Results.g,Results.hess] = fminunc(LLfun,b0,OptimOpt);
    else
        [Results.bhat,LL,Results.exitf,Results.output,Results.g] = fminunc(LLfun,b0,OptimOpt);
    end
end


%% Output


% save tmp_MNL_output

Results.LL = -LL;
Results.b0_old = b0;

if isfield(EstimOpt,'R2type') == 0
    EstimOpt.R2type = 0;
end

Results.LLdetailed = probs_mmdcev(INPUT, EstimOpt, Results.bhat);

if EstimOpt.HessEstFix == 1
    f = probs_mmdcev(INPUT, EstimOpt,Results.bhat);
    Results.jacobian = numdiff(@(B) -probs_mmdcev(INPUT, EstimOpt,B),-f,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) -probs_mmdcev(INPUT, EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.hess = hessian(@(B) -sum(probs_mmdcev(INPUT,EstimOpt,B),1),Results.bhat);
% elseif EstimOpt.HessEstFix == 4 % analytical - missing
%     Results.hess = hessian(@(B) -sum(probs_mdcev(INPUT,EstimOpt,B),1),Results.bhat);
end


if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
    Results.hess = Results.jacobian'*Results.jacobian;
end
% Results.ihess = inv(Results.hess);

EstimOpt.BLimit = (sum(Results.hess) == 0 & EstimOpt.BActive == 1);
EstimOpt.BActive(EstimOpt.BLimit == 1) = 0;
Results.hess = Results.hess(EstimOpt.BActive == 1,EstimOpt.BActive == 1);
Results.ihess = inv(Results.hess);
Results.ihess = direcXpnd(Results.ihess,EstimOpt.BActive);
Results.ihess = direcXpnd(Results.ihess',EstimOpt.BActive);


Results.std = sqrt(diag(Results.ihess));

if sum(EstimOpt.BActive == 0) > 0
    Results.std(EstimOpt.BActive == 0) = NaN;
end

Results.std(imag(Results.std) ~= 0) = NaN;
Results.R = [Results.bhat,Results.std,pv(Results.bhat,Results.std)];

EstimOpt.Params = length(b0);

Results.INPUT = INPUT;
Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;

Results.DetailsA(1:EstimOpt.NVarA,1) = Results.bhat(1:EstimOpt.NVarA);
Results.DetailsA(1:EstimOpt.NVarA,3:4) = [Results.std(1:EstimOpt.NVarA),pv(Results.bhat(1:EstimOpt.NVarA),Results.std(1:EstimOpt.NVarA))];
if EstimOpt.FullCov == 0
    Results.DetailsV(1:EstimOpt.NVarA,1) = abs(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2));
    Results.DetailsV(1:EstimOpt.NVarA,3:4) = ...
        [Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*2), ...
         pv(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2), ...
         Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*2))];
    l = 2*EstimOpt.NVarA;
else
    Results.DetailsV = ...
        sdtri(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2), ...
              Results.ihess(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2, ...
                            EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2), ...
              EstimOpt);
    Results.DetailsV = [Results.DetailsV(:,1),zeros(EstimOpt.NVarA,1),Results.DetailsV(:,2:3)];
    l = EstimOpt.NVarA*(EstimOpt.NVarA+3)/2; 
end
% alphas or gammas
if EstimOpt.NVarM == 0
    if EstimOpt.Profile == 1 % alpha profile (alphas = 1 - exp(-alphas))
        Results.DetailsProfile(1:EstimOpt.NAlt, 1) = ...
            1 - exp(-Results.bhat(l+1:l+EstimOpt.NAlt));
        Tmp = exp(-Results.bhat(l+1:l+EstimOpt.NAlt));
        Results.DetailsProfile(1:EstimOpt.NAlt,3:4) = ...
            [Results.std(l+1:l+EstimOpt.NAlt).*Tmp, ...
            pv(1 - exp(-Results.bhat(l+1:l+EstimOpt.NAlt)),...
            Results.std(l+1:l+EstimOpt.NAlt).*Tmp)];
    elseif EstimOpt.Profile == 2 % gamma profile (gammas = exp(gammas))
        Results.DetailsProfile(1:EstimOpt.NAlt, 1) = ...
        exp(Results.bhat(l+1:l+EstimOpt.NAlt));
        Tmp = exp(Results.bhat(l+1:l+EstimOpt.NAlt));
    Results.DetailsProfile(1:EstimOpt.NAlt,3:4) = ...
        [Results.std(l+1:l+EstimOpt.NAlt).*Tmp, ...
         pv(exp(Results.bhat(l+1:l+EstimOpt.NAlt)),...
         Results.std(l+1:l+EstimOpt.NAlt).*Tmp)];
    elseif EstimOpt.Profile == 3
        indx = length(unique(EstimOpt.SpecProfile(1,EstimOpt.SpecProfile(1,:) ~= 0)));
        b1 = Results.bhat(l+1:l+EstimOpt.NVarP);
        b2 = b1(indx+1:end);
        b1 = b1(1:indx);
        std1 = Results.std(l+1:l+EstimOpt.NVarP);
        std2 = std1(indx+1:end);
        std1 = std1(1:indx);
    
        Results.DetailsProfile(1:indx, 1) = 1 - exp(-b1);
        Tmp = exp(-b1);
        Results.DetailsProfile(1:indx,3:4) = [std1.*Tmp, pv(1 - exp(-b1), std1.*Tmp)];
        Results.DetailsProfile(indx+1:EstimOpt.NVarP, 1) = exp(b2);
        Tmp = exp(b2);
        Results.DetailsProfile(indx+1:EstimOpt.NVarP, 3:4) = [std2.*Tmp, pv(exp(b2),std2.*Tmp)];
    end
else

    Results.DetailsProfile = [];
    for i=1:EstimOpt.NVarP
        Results.DetailsProfile(1:(EstimOpt.NVarM+1),4*i-3) = Results.bhat(l+(EstimOpt.NVarM+1)*(i-1)+1:l+(EstimOpt.NVarM+1)*(i));
        Results.DetailsProfile(1:(EstimOpt.NVarM+1),4*i-1:4*i) = ...
            [Results.std(l+(EstimOpt.NVarM+1)*(i-1)+1:l+(EstimOpt.NVarM+1)*(i)), ...
             pv(Results.bhat(l+(EstimOpt.NVarM+1)*(i-1)+1:l+(EstimOpt.NVarM+1)*(i)), ...
             Results.std(l+(EstimOpt.NVarM+1)*(i-1)+1:l+(EstimOpt.NVarM+1)*(i)))];
    end
end
Results.DetailsScale(1, 1) = exp(Results.bhat(end));
Results.DetailsScale(1, 3:4) = ...
    [Results.std(end)*exp(Results.bhat(end)), ...
     pv(exp(Results.bhat(end)),...
     Results.std(end)*exp(Results.bhat(end)))];

%% Template filling

% Template1 = {'DetailsA'};
% Template2 = {'DetailsA'};
% Names.DetailsA = EstimOpt.NamesA;
% Heads.DetailsA = {'';'tb'};

Template1 = {'DetailsA','DetailsV'};
Template2 = {'DetailsA','DetailsV'};
Names.DetailsA = EstimOpt.NamesA;
Heads.DetailsA = {'Means';'tc'};
Heads.DetailsV = {'Standard Deviations';'lb'};
ST = {};

% Template1 = [Template1,{'DetailsProfile'}];
Temp = cell(1,size(Template2,2));
Temp(1,1) = {'DetailsProfile'};
Template1 = [Template1;Temp];
Template2 = [Template2;Temp];
if EstimOpt.NVarM == 0
    Names.DetailsProfile = EstimOpt.NamesAlt;
    if EstimOpt.Profile == 1
        Heads.DetailsProfile = {'Coefficients of Alpha-profile'};
    else
        Heads.DetailsProfile = {'Coefficients of Gamma-profile'};
    end
    Heads.DetailsProfile(2,:) = {'tb'};
%      ST = [ST,'DetailsProfile'];
else
    Heads.DetailsProfile(:,2) = [EstimOpt.NamesAlt;{'lb'}];
    Heads.DetailsProfile(1:2,1) = {'Alpha/Gamma profile coef.';'lc'};
    Names.DetailsProfile = [{'Const.'}; EstimOpt.NamesM];
%      ST = [ST,'DetailsProfile'];
end
ST = {'DetailsProfile'};

Temp = cell(1,size(Template2,2));
Temp(1,1) = {'DetailsScale'};
Template1 = [Template1;Temp];
Template2 = [Template2;Temp];
Names.DetailsScale = {'Scale'};
Heads.DetailsScale = {'Coefficients of scale'};
Heads.DetailsScale(2,:) = {'tb'};
ST = [ST;{'DetailsScale'}];

%% Header

Head = cell(1,2);
if EstimOpt.FullCov == 0
    Head(1,1) = {'Mixed MDCEV_d'};
else
    Head(1,1) = {'Mixed MDCEV'};
end
Head(1,2) = {'in preference-space'}; % Probably not needed

%% Footer
EstimOpt.NObs = EstimOpt.NP*EstimOpt.NCT; % May need to be changed if there are missing obs.

Results.stats = [Results.LL; NaN; NaN; NaN; ((2*EstimOpt.Params - 2*Results.LL))/EstimOpt.NObs;((log(EstimOpt.NObs)*EstimOpt.Params - 2*Results.LL))/EstimOpt.NObs;EstimOpt.NObs;EstimOpt.NP;EstimOpt.Params];
Tail = cell(16,2);
Tail(2,1) = {'Model diagnostics'};
Tail(3:16,1) = {'LL at convergence';'LL at constant(s) only';strcat('McFadden''s pseudo-R',char(178));strcat('Ben-Akiva-Lerman''s pseudo-R',char(178));'AIC/n';'BIC/n';'n (observations)';'r (respondents)';'k (parameters)';' ';'Estimation method';'Optimization method';'Gradient';'Hessian'};
Tail(3:11,2) = num2cell(Results.stats);

if any(INPUT.W ~= 1)
    Tail(13,2) = {'weighted maximum likelihood'};
else
    Tail(13,2) = {'maximum likelihood'};
end

Tail(14,2) = {OptimOpt.Algorithm;};

if strcmp(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        Tail(15,2) = {'user-supplied, analytical'};
    else
        Tail(15,2) = {['user-supplied, numerical ',num2str(OptimOpt.FinDiffType)]};
    end
else
    Tail(15,2) = {['built-in, ',num2str(OptimOpt.FinDiffType)]};
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

Tail(16,2) = {outHessian};


%%  Print to screen and .xls

if EstimOpt.Display ~= 0
     Results.Dist = EstimOpt.Dist;
    Results.R_out = genOutput(EstimOpt,Results,Head,Tail,Names,Template1,Template2,Heads,ST);
end
end
