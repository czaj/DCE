function Results = MDCEV(INPUT,EstimOpt,OptimOpt)
% MDCEV creates Multiple Discrete-Continuous Extreme Value model.
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

save tmp_MDCEV
% return

global B_backup

tic

Results.bhat = [];
Results.R = [];
Results.R_out = {};
Results.stats = [];


%% Check data and inputs

if nargin < 3
    error('Too few input arguments for MDCEV(INPUT,EstimOpt,OptimOpt)')
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
    
    disp('Estimating MDCEV model ...')
    
    disp('in preference-space ...')
end
EstimOpt.NVarA = size(INPUT.Xa,2);

if isfield(EstimOpt,'Profile') == 0 || isempty(EstimOpt.Profile) || ~ismember(EstimOpt.Profile, [1, 2])
    if EstimOpt.Display ~= 0
        cprintf(rgb('DarkOrange'), 'WARNING: Assuming alpha-profile. For gamma-profile specify EstimOpt.Profile = 2.\n')
    end
    EstimOpt.Profile = 1;

end
if isfield(EstimOpt,'NamesA') == 0 || isempty(EstimOpt.NamesA) || length(EstimOpt.NamesA) ~= EstimOpt.NVarA
    EstimOpt.NamesA = (1:EstimOpt.NVarA)';
    EstimOpt.NamesA = cellstr(num2str(EstimOpt.NamesA));
elseif size(EstimOpt.NamesA,1) ~= EstimOpt.NVarA
    EstimOpt.NamesA = EstimOpt.NamesA';
end
if isfield(EstimOpt,'NamesAlt') == 0 || isempty(EstimOpt.NamesAlt) || length(EstimOpt.NamesAlt) ~= EstimOpt.NAlt
    EstimOpt.NamesAlt = (1:EstimOpt.NAlt)';
    EstimOpt.NamesAlt = cellstr(num2str(EstimOpt.NamesAlt));
elseif size(EstimOpt.NamesAlt,1) ~= EstimOpt.NAlt
    EstimOpt.NamesAlt = EstimOpt.NamesAlt';
end


%% Starting values
if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVarA + EstimOpt.NAlt + 1
    b0 = B_backup(:);
    if EstimOpt.Display ~= 0
        disp('Using the starting values from Backup')
    end
end

if ~exist('b0','var')
    if EstimOpt.Display ~= 0
        disp('Using linear regression estimates as starting values')
    end
    % [covariates parameters; alphas or gammas; scale]
    b0 = [regress(INPUT.Y,INPUT.Xa); 0.5*ones(EstimOpt.NAlt,1); 1];
    
%     % Maybe add a transformation for alpha or gamma profile
%     if EstimOpt.Profile == 1 % alpha profile (alpha = 1 - exp(-alpha))
%         b0(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NAlt) = ...
%             -log(1 - b0(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NAlt));
%     else % gamma profile (gamma = exp(gamma))
%         b0(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NAlt) = ...
%             log(b0(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NAlt));
%     end
end


%% Optimization Options

if isfield(EstimOpt,'BActive')
    EstimOpt.BActive = EstimOpt.BActive(:)';
else
    EstimOpt.BActive = ones(length(b0),1);
end
if isfield(EstimOpt,'ConstVarActive') == 0 
    EstimOpt.ConstVarActive = 0;
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


%% Restucturing Data 

INPUT.Y = reshape(INPUT.Y, [EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP]);
INPUT.priceMat = reshape(INPUT.priceMat, [EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP]);

% idx = sum(reshape(INPUT.MissingInd,[EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP])) == EstimOpt.NAlt;
% idx = reshape(idx(ones(EstimOpt.NAlt,1),:),[EstimOpt.NAlt*EstimOpt.NCT*EstimOpt.NP,1]);
% INPUT.Y = INPUT.Y(idx == 0);
% INPUT.Xa(INPUT.MissingInd == 1,:) = NaN;
% 
% %Xa has NaNs where there is missing alternative
% 
% INPUT.Xa = INPUT.Xa(idx == 0,:);

%% Estimation

% =========================================================================
% INPUT must contain Xa, priceMat, y (dependent variable; demands) ========
% =========================================================================

LLfun = @(B) LL_mdcev_MATlike(INPUT, EstimOpt,OptimOpt,B);

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

Results.LLdetailed = probs_mdcev(INPUT, EstimOpt, Results.bhat);

if EstimOpt.HessEstFix == 1
    f = probs_mdcev(INPUT, EstimOpt,Results.bhat);
    Results.jacobian = numdiff(@(B) -probs_mdcev(INPUT, EstimOpt,B),-f,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) -probs_mdcev(INPUT, EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.hess = hessian(@(B) -sum(probs_mdcev(INPUT,EstimOpt,B),1),Results.bhat);
% elseif EstimOpt.HessEstFix == 4 % analytical - missing
%     Results.hess = hessian(@(B) -sum(probs_mdcev(INPUT,EstimOpt,B),1),Results.bhat);
end


if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
    Results.hess = Results.jacobian'*Results.jacobian;
end
Results.ihess = inv(Results.hess);
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

% alphas or gammas
if EstimOpt.Profile == 1 % alpha profile (alphas = 1 - exp(-alphas))
    Results.DetailsProfile(1:EstimOpt.NAlt, 1) = ...
        1 - exp(-Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NAlt));
    Tmp = exp(-Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NAlt));
    Results.DetailsProfile(1:EstimOpt.NAlt,3:4) = ...
        [Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NAlt).*Tmp, ...
        pv(1 - exp(-Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NAlt)),...
        Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NAlt).*Tmp)];
elseif EstimOpt.Profile == 2 % gamma profile (gammas = exp(gammas))
    Results.DetailsProfile(1:EstimOpt.NAlt, 1) = ...
    exp(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NAlt));
    Tmp = exp(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NAlt));
Results.DetailsProfile(1:EstimOpt.NAlt,3:4) = ...
    [Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NAlt).*Tmp, ...
     pv(exp(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NAlt)),...
     Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NAlt).*Tmp)];
end
 
Results.DetailsScale(1, 1) = exp(Results.bhat(end));
Results.DetailsScale(1, 3:4) = ...
    [Results.std(end)*exp(Results.bhat(end)), ...
     pv(exp(Results.bhat(end)),...
     Results.std(end)*exp(Results.bhat(end)))];

%% Template filling

Template1 = {'DetailsA'};
Template2 = {'DetailsA'};
Names.DetailsA = EstimOpt.NamesA;
Heads.DetailsA = {'';'tb'};
ST = {'DetailsA'};

Template1 = [Template1;'DetailsProfile'];
Template2 = [Template2;'DetailsProfile'];
Names.DetailsProfile = EstimOpt.NamesAlt;
if EstimOpt.Profile == 1
    Heads.DetailsProfile = {'Coefficients of Alpha-profile'};
else
    Heads.DetailsProfile = {'Coefficients of Gamma-profile'};
end
Heads.DetailsProfile(2,:) = {'tb'};
ST = [ST,'DetailsProfile'];

Template1 = [Template1;'DetailsScale'];
Template2 = [Template2;'DetailsScale'];
Names.DetailsScale = {'Scale'};
Heads.DetailsScale = {'Coefficients of scale'};
Heads.DetailsScale(2,:) = {'tb'};
ST = [ST,'DetailsScale'];

%% Header

Head = cell(1,2);
Head(1,1) = {'MDCEV'};
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
    Results.Dist = -ones(EstimOpt.NVarA,1);
    Results.R_out = genOutput(EstimOpt,Results,Head,Tail,Names,Template1,Template2,Heads,ST);
end
end
