function Results = MIMIC(INPUT,Results_old,EstimOpt,OptimOpt)


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
	disp(['Assuming ', num2str(EstimOpt.NLatent)',' NLatent Variable(s)']);
end
disp(num2str(EstimOpt.NLatent,'Estimating MIMIC with %1.0f NLatent Variable(s)'))


if isfield(INPUT, 'Xstr') == 0 || numel(INPUT.Xstr) == 0
	error('Define variables to structural equations')
end
if isfield(INPUT, 'Xmea') == 0 || numel(INPUT.Xmea) == 0
	error('Define attitudes to measurment equations')
end

if isfield(INPUT, 'Xmea_exp') == 0 || numel(INPUT.Xmea_exp) == 0 % additional covariates for explaining Measurment equations
	INPUT.Xmea_exp = [];
    EstimOpt.MeaExpMatrix = zeros(1,size(INPUT.Xmea,2));
else
    if isfield(EstimOpt,'MeaExpMatrix') == 0 || length(EstimOpt.MeaExpMatrix) ~= size(INPUT.Xmea,2)
        EstimOpt.MeaExpMatrix = ones(1, size(INPUT.Xmea,2));
        cprintf(rgb('DarkOrange'), 'WARNING: MeaExpMatrix not defined - assuming that every measurment equation is explained with additional covariates \n')
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
	EstimOpt.MeaSpecMatrix = zeros(1, size(INPUT.Xmea,2));
	disp('Using OLS for every measurment equation')
end
if numel(EstimOpt.MeaSpecMatrix) == 1 % if only one number provided - replicate it for all Measurement variables
    EstimOpt.MeaSpecMatrix = EstimOpt.MeaSpecMatrix * ones(1,size(INPUT.Xmea,2));
end
if size(EstimOpt.MeaSpecMatrix,2) ~= size(INPUT.Xmea,2)
	error('Measurment specification (model) erroneously defined')
end

EstimOpt.NVarstr = size(INPUT.Xstr,2);  % no. of variables in structural equation
EstimOpt.NVarmea = sum(sum(EstimOpt.MeaMatrix)); % no of parameters for Measurments without couting cutoffs, constants etc
EstimOpt.NVarmea_exp = size(INPUT.Xmea_exp,2);

for i=1:size(EstimOpt.MeaMatrix,2)
    if numel(EstimOpt.MeaSpecMatrix(i) > 0) > 0
        if EstimOpt.MeaSpecMatrix(i) > 0 && numel(unique(INPUT.Xmea(:,i))) > 10 && EstimOpt.MeaSpecMatrix(i) < 3
            cprintf(rgb('DarkOrange'), 'WARNING: There are over 10 levels for measurement variable %d \n', i)
        end
    end
    if sum(isnan(INPUT.Xmea(:,i))) > 0
        cprintf(rgb('DarkOrange'), 'WARNING: Measurement variable %d contains NaN values \n', i)
    end
	if sum(isinf(INPUT.Xmea(:,i))) > 0
        cprintf(rgb('DarkOrange'), 'WARNING: Measurement variable %d contains Inf values \n', i)
	end
end
if any(sum(EstimOpt.MeaMatrix,2) == 1) == 1
    error('There must be more than one attitude for every Latent Variable for model identification')    
end

for i=1:size(EstimOpt.NVarstr,1)
    if sum(isnan(INPUT.Xstr(:,i))) > 0
        cprintf(rgb('DarkOrange'), 'WARNING: Structural variable %d contains NaN values \n', i)
    end
	if sum(isinf(INPUT.Xstr(:,i))) > 0
        cprintf(rgb('DarkOrange'), 'WARNING: Structural variable %d contains Inf values \n', i)
	end
end

if EstimOpt.NVarmea_exp > 0
    if isfield(EstimOpt,'NamesMeaExp') == 0 || isempty(EstimOpt.NamesMeaExp) || length(EstimOpt.NamesMeaExp) ~= EstimOpt.NVarmea_exp
        EstimOpt.NamesMeaExp = (1:EstimOpt.NVarmea_exp)';
        EstimOpt.NamesMeaExp = cellstr(num2str(EstimOpt.NamesMeaExp));
    elseif size(EstimOpt.NamesMeaExp,1) ~= EstimOpt.NVarmea_exp
        EstimOpt.NamesMeaExp = EstimOpt.NamesMeaExp';
    end
end

EstimOpt.Names = [];% Names of the models
EstimOpt.NVarcut = 0; % no of cutoffs for ordered probits + constants + variances for OLS 
EstimOpt.NVarcut0 = 0; % no of cutoffs for HMNL0
EstimOpt.CutMatrix = zeros(1, size(INPUT.Xmea,2)); 
EstimOpt.NamesLV = {};
for i = 1:size(INPUT.Xmea,2)
    if EstimOpt.MeaSpecMatrix(i) == 2 %Ordered probit: cutoffs
        EstimOpt.NVarcut = EstimOpt.NVarcut + length(unique(INPUT.Xmea(:,i))) - 1 + EstimOpt.NVarmea_exp*(EstimOpt.MeaExpMatrix(i) ~=0);
        EstimOpt.CutMatrix(i) = length(unique(INPUT.Xmea(:,i))) - 1 + sum(EstimOpt.MeaMatrix(:,i))+ EstimOpt.NVarmea_exp*(EstimOpt.MeaExpMatrix(i) ~=0);
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV; cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end    
        if EstimOpt.MeaExpMatrix(i) ~=0
           EstimOpt.NamesLV = [EstimOpt.NamesLV; EstimOpt.NamesMeaExp];
        end
        for n = 1:(length(unique(INPUT.Xmea(:,i))) - 1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV; cellfun(@(x)[x num2str(n)],{'Cutoff '},'UniformOutput',0)];
        end
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + length(unique(INPUT.Xmea(:,i))) - 1;
        EstimOpt.Names = [EstimOpt.Names, 'OP '];
    elseif EstimOpt.MeaSpecMatrix(i) == 0
        EstimOpt.NVarcut = EstimOpt.NVarcut + 2+EstimOpt.NVarmea_exp*(EstimOpt.MeaExpMatrix(i) ~=0); %OLS: constant + sigma
        EstimOpt.CutMatrix(i) = 2+sum(EstimOpt.MeaMatrix(:,i))+ EstimOpt.NVarmea_exp*(EstimOpt.MeaExpMatrix(i) ~=0);
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + 2;
        EstimOpt.Names = [EstimOpt.Names, 'OLS ']; 
        EstimOpt.NamesLV = [EstimOpt.NamesLV; {'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV; cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end    
        if EstimOpt.MeaExpMatrix(i) ~=0
           EstimOpt.NamesLV = [EstimOpt.NamesLV; EstimOpt.NamesMeaExp];
        end
        EstimOpt.NamesLV = [EstimOpt.NamesLV; {'Sigma'}];
    elseif EstimOpt.MeaSpecMatrix(i) == 1 % MNL 
        EstimOpt.NVarcut = EstimOpt.NVarcut + (length(unique(INPUT.Xmea(:,i))) - 2)*sum(EstimOpt.MeaMatrix(:,i)) + (length(unique(INPUT.Xmea(:,i))) - 1)*(1+ EstimOpt.NVarmea_exp*(EstimOpt.MeaExpMatrix(i) ~=0)); % constants + additional coefficients 
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + length(unique(INPUT.Xmea(:,i)))-1;
        EstimOpt.Names = [EstimOpt.Names, 'MNL '];
        EstimOpt.CutMatrix(i) = (1+ EstimOpt.NVarmea_exp*(EstimOpt.MeaExpMatrix(i) ~=0)+sum(EstimOpt.MeaMatrix(:,i)))*(length(unique(INPUT.Xmea(:,i)))-1);
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for j = 1:(length(unique(INPUT.Xmea(:,i)))-1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV; {'Cons.'}];
            for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
                EstimOpt.NamesLV = [EstimOpt.NamesLV; cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
            end
            if EstimOpt.MeaExpMatrix(i) ~=0
                EstimOpt.NamesLV = [EstimOpt.NamesLV; EstimOpt.NamesMeaExp];
            end
        end
    elseif EstimOpt.MeaSpecMatrix(i) == 3 % Poisson    
        EstimOpt.NVarcut = EstimOpt.NVarcut +1+EstimOpt.NVarmea_exp*(EstimOpt.MeaExpMatrix(i) ~=0); %Poiss: only constant
        EstimOpt.CutMatrix(i) = 1+sum(EstimOpt.MeaMatrix(:,i))+ EstimOpt.NVarmea_exp*(EstimOpt.MeaExpMatrix(i) ~=0);
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + 1;
        EstimOpt.Names = [EstimOpt.Names, 'POISS ']; 
        EstimOpt.NamesLV = [EstimOpt.NamesLV; {'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV; cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end    
        if EstimOpt.MeaExpMatrix(i) ~=0
           EstimOpt.NamesLV = [EstimOpt.NamesLV; EstimOpt.NamesMeaExp];
        end
    elseif EstimOpt.MeaSpecMatrix(i) == 4 % Negative Binomial   
        EstimOpt.NVarcut = EstimOpt.NVarcut +2+EstimOpt.NVarmea_exp*(EstimOpt.MeaExpMatrix(i) ~=0); %Poiss: only constant
        EstimOpt.CutMatrix(i) = 2+sum(EstimOpt.MeaMatrix(:,i))+ EstimOpt.NVarmea_exp*(EstimOpt.MeaExpMatrix(i) ~=0);
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + 2;
        EstimOpt.Names = [EstimOpt.Names, 'NB ']; 
        EstimOpt.NamesLV = [EstimOpt.NamesLV; {'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV; cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end    
        if EstimOpt.MeaExpMatrix(i) ~=0
           EstimOpt.NamesLV = [EstimOpt.NamesLV; EstimOpt.NamesMeaExp];
        end
        EstimOpt.NamesLV = [EstimOpt.NamesLV; {'Theta'}];
    elseif EstimOpt.MeaSpecMatrix(i) == 5 % ZIP   
        EstimOpt.NVarcut = EstimOpt.NVarcut +2 +sum(EstimOpt.MeaMatrix(:,i)) +2*EstimOpt.NVarmea_exp*(EstimOpt.MeaExpMatrix(i) ~=0); %Poiss: only constant
        EstimOpt.CutMatrix(i) = 2+2*sum(EstimOpt.MeaMatrix(:,i))+ 2*EstimOpt.NVarmea_exp*(EstimOpt.MeaExpMatrix(i) ~=0);
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + 2;
        EstimOpt.Names = [EstimOpt.Names, 'ZIP ']; 
        EstimOpt.NamesLV = [EstimOpt.NamesLV; {'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV; cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end    
        if EstimOpt.MeaExpMatrix(i) ~=0
           EstimOpt.NamesLV = [EstimOpt.NamesLV; EstimOpt.NamesMeaExp];
        end
        EstimOpt.NamesLV = [EstimOpt.NamesLV; {'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV; cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end    
        if EstimOpt.MeaExpMatrix(i) ~=0
           EstimOpt.NamesLV = [EstimOpt.NamesLV; EstimOpt.NamesMeaExp];
        end
    end
    
end

if isfield(EstimOpt,'NamesStr') == 0 || isempty(EstimOpt.NamesStr) || length(EstimOpt.NamesStr) ~= EstimOpt.NVarstr
    EstimOpt.NamesStr = (1:EstimOpt.NVarstr)';
    EstimOpt.NamesStr = cellstr(num2str(EstimOpt.NamesStr));
elseif size(EstimOpt.NamesStr,1) ~= EstimOpt.NVarstr
    EstimOpt.NamesStr = EstimOpt.NamesStr';
end

if isfield(EstimOpt,'NamesMea') == 0 || isempty(EstimOpt.NamesMea) || length(EstimOpt.NamesMea) ~= size(INPUT.Xmea,2)
    EstimOpt.NamesMea = (1:size(INPUT.Xmea,2))';
    EstimOpt.NamesMea = cellstr(num2str(EstimOpt.NamesMea));
elseif size(EstimOpt.NamesMea,1) ~= size(INPUT.Xmea,2)
    EstimOpt.NamesMea = EstimOpt.NamesMea';
end

disp(['Using following models for attitudes: '  char(EstimOpt.Names)]);


%% Rescructure data


INPUT.Xstr = double(INPUT.Xstr(1:EstimOpt.NAlt*EstimOpt.NCT:end,:));
if isfield(EstimOpt, 'StrNorm')
    if length(EstimOpt.StrNorm) ~= EstimOpt.NVarstr && length(EstimOpt.StrNorm) ~= 1
        cprintf(rgb('DarkOrange'), 'WARNING: Structural variables normalization options (EstimOpt.StrNorm) incorrect - variables not normalized \n')
    else
        if length(EstimOpt.StrNorm) == 1
            EstimOpt.StrNorm = EstimOpt.StrNorm(ones(EstimOpt.NVarstr,1),1);
        end       
        meanx = mean(INPUT.Xstr);
        stdx = std(INPUT.Xstr);
        INPUT.Xstr(:,EstimOpt.StrNorm == 1) = (INPUT.Xstr(:,EstimOpt.StrNorm == 1) - meanx(ones(size(INPUT.Xstr,1),1),EstimOpt.StrNorm == 1))./stdx(ones(size(INPUT.Xstr,1),1),EstimOpt.StrNorm == 1);
    end
end

INPUT.Xmea = INPUT.Xmea(1:EstimOpt.NAlt*EstimOpt.NCT:end,:);
if EstimOpt.NVarmea_exp > 0
    INPUT.Xmea_exp = INPUT.Xmea_exp(1:EstimOpt.NAlt*EstimOpt.NCT:end,:);
end


%% Estimating MIMIC0


OptimOpt_0 = optimoptions('fminunc');
OptimOpt_0.Algorithm = 'quasi-newton';
OptimOpt_0.GradObj = 'off';
OptimOpt_0.Hessian = 'off';
OptimOpt_0.Display = 'off';
OptimOpt_0.FunValCheck = 'off'; 
OptimOpt_0.Diagnostics = 'off';

LLfun0 = @(B) LL_hmnl0(INPUT.Xmea, EstimOpt, B);
[Results.MIMIC0.bhat, LL0] = fminunc(LLfun0, 0.001*ones(EstimOpt.NVarcut0,1),OptimOpt_0);
Results.MIMIC0.LL = -LL0;


%% Starting values

if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVarstr*EstimOpt.NLatent + EstimOpt.NVarmea + EstimOpt.NVarcut
    b0 = B_backup(:);
    disp('Using the starting values from Backup')
elseif isfield(Results_old,'MIMIC') && isfield(Results_old.MIMIC,'b0') % starting values provided
    Results_old.MIMIC.b0_old = Results_old.MIMIC.b0(:);
    Results_old.MIMIC = rmfield(Results_old.MIMIC,'b0');
    if length(Results_old.MIMIC.b0_old) ~= EstimOpt.NVarstr*EstimOpt.NLatent + EstimOpt.NVarmea + EstimOpt.NVarcut
        cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
        Results_old.MIMIC = rmfield(Results_old.MIMIC,'b0_old');
    else
        b0 = Results_old.MIMIC.b0_old(:);
    end
end
if  ~exist('b0','var')
    b0 = zeros(EstimOpt.NVarstr*EstimOpt.NLatent + EstimOpt.NVarmea + EstimOpt.NVarcut,1);
    l = EstimOpt.NVarstr*EstimOpt.NLatent;
    k = 0;
    for i = 1:size(INPUT.Xmea,2)
        if EstimOpt.MeaSpecMatrix(i) == 2 %Ordered probit: cutoffs
            b0(l+ sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp + 1:l+ sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp + length(unique(INPUT.Xmea(:,i))) -1) = Results.MIMIC0.bhat(k+1:k+length(unique(INPUT.Xmea(:,i)))-1);
            l = l + sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp + length(unique(INPUT.Xmea(:,i)))-1;
            k = k+ length(unique(INPUT.Xmea(:,i)))-1;
        elseif EstimOpt.MeaSpecMatrix(i) == 0
            b0(l + 1) = Results.MIMIC0.bhat(k+1);
            b0(l+2+sum(EstimOpt.MeaMatrix(:,i)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp,1)) = Results.MIMIC0.bhat(k+2);
            k = k+2;
            l = l+2 + sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp;        
        elseif EstimOpt.MeaSpecMatrix(i) == 1 % MNL 
            for n = 1:(length(unique(INPUT.Xmea(:,i)))-1)
                b0(l + 1) = Results.MIMIC0.bhat(k+1);
                l = l+ 1+ sum(EstimOpt.MeaMatrix(:,i))+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp;
                k = k+ 1;
            end
        elseif EstimOpt.MeaSpecMatrix(i) == 3 % POISS 
            b0(l + 1) = Results.MIMIC0.bhat(k+1);
            k = k+1;
            l = l+1 + sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp;     
        elseif EstimOpt.MeaSpecMatrix(i) == 4 % NB
            b0(l + 1) = Results.MIMIC0.bhat(k+1);
            b0(l + 2+ sum(EstimOpt.MeaMatrix(:,i)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp,1)) = Results.MIMIC0.bhat(k+2); % theta
            k = k+2;
            l = l+2 + sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp;  
        elseif EstimOpt.MeaSpecMatrix(i) == 5 % ZIP
            b0(l + 1) = Results.MIMIC0.bhat(k+1);
            b0(l + 2+ sum(EstimOpt.MeaMatrix(:,i)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp,1)) = Results.MIMIC0.bhat(k+2); % second constant
            k = k+2;
            l = l+2 + 2*sum(EstimOpt.MeaMatrix(:,i),1)+2*EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarmea_exp;     
        end
    end
end


%% Optimization Options


if  isfield(EstimOpt,'BActive')
	EstimOpt.BActive = EstimOpt.BActive(:)';
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
    err_mtx = randn(EstimOpt.NP*EstimOpt.NRep, EstimOpt.NLatent+EstimOpt.NVarA+1); %to be cut down later   
elseif EstimOpt.Draws == 2 % LHS
    cprintf('*blue','Latin Hypercube Sampling '); cprintf('draws \n');
    err_mtx=lhsnorm(zeros((EstimOpt.NLatent+EstimOpt.NVarA+1)*EstimOpt.NP,1),diag(ones((EstimOpt.NLatent+EstimOpt.NVarA+1)*EstimOpt.NP,1)),EstimOpt.NRep); 
    err_mtx = reshape(err_mtx, EstimOpt.NRep*EstimOpt.NP, EstimOpt.NLatent);
elseif EstimOpt.Draws >= 3 % Quasi random draws 
    if EstimOpt.Draws == 3
        cprintf('*blue','Halton '); cprintf('draws (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = haltonset(EstimOpt.NLatent+EstimOpt.NVarA+1,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); % 
    elseif EstimOpt.Draws == 4 % apply reverse-radix scrambling
        cprintf('*blue','Halton '); cprintf('draws with reverse radix scrambling (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = haltonset(EstimOpt.NLatent+EstimOpt.NVarA+1,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); % 
        hm1 = scramble(hm1,'RR2');
    elseif EstimOpt.Draws == 5
        cprintf('*blue','Sobol '); cprintf('draws (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = sobolset(EstimOpt.NLatent+EstimOpt.NVarA+1,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); 
    elseif EstimOpt.Draws == 6
        cprintf('*blue','Sobol '); cprintf('draws with random linear scramble and random digital shift (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
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

if isfield(EstimOpt, 'Drawskeep') && ~isempty(EstimOpt.Drawskeep) && EstimOpt.Drawskeep == 1
    Results.err = err_sliced;
end


%% Display Options


if ((isfield(EstimOpt, 'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    OptimOpt.GradObj = 'on';
end

if any(EstimOpt.MissingAlt(:) == 1) && EstimOpt.NumGrad == 0
   EstimOpt.NumGrad = 1;
   cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - missing alternatives not supported by analytical gradient \n')
end

if any(EstimOpt.MeaSpecMatrix >= 3) && EstimOpt.NumGrad == 0 && any(any(INPUT.Xmea(:, EstimOpt.MeaSpecMatrix >=3) > 100))
   cprintf(rgb('DarkOrange'), 'WARNING: it is recommended to switch to numerical gradient, as analitycal can be not precise when Xmea take large values for NB \n')
end

% if any(EstimOpt.MeaSpecMatrix == 5) && EstimOpt.NumGrad == 0 
%    EstimOpt.NumGrad = 1;
%    cprintf(rgb('DarkOrange'), 'WARNING: ZIP\n')
% end

if (isfield(EstimOpt, 'ConstVarActive') == 0 || EstimOpt.ConstVarActive == 0) && isequal(OptimOpt.Algorithm,'quasi-newton') && isequal(OptimOpt.Hessian,'user-supplied')
	cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied Hessian off - quasi-newton algorithm does not use it anyway \n')
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


%% Estimation


LLfun = @(B) LL_mimic_MATlike(INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp, err_sliced,EstimOpt,OptimOpt,B);
if EstimOpt.ConstVarActive == 0  
    
    if EstimOpt.HessEstFix == 0
        [Results.bhat, LL, Results.exitf, Results.output, Results.g, Results.hess] = fminunc(LLfun, b0, OptimOpt);
    else
        [Results.bhat, LL, Results.exitf, Results.output, Results.g] = fminunc(LLfun, b0, OptimOpt);
    end      
    
elseif EstimOpt.ConstVarActive == 1 % equality constraints
        
    EstimOpt.CONS1 = diag(1 - EstimOpt.BActive);
    EstimOpt.CONS1(sum(EstimOpt.CONS1,1)==0,:)=[];
    EstimOpt.CONS2 = zeros(size(EstimOpt.CONS1,1),1);
%     EstimOpt.CONS1 = sparse(EstimOpt.CONS1);
%     EstimOpt.CONS2 = sparse(EstimOpt.CONS2);
    if EstimOpt.HessEstFix == 0
        [Results.bhat, LL, Results.exitf, Results.output, Results.lambda, Results.g, Results.hess] = fmincon(LLfun,b0,[],[],EstimOpt.CONS1,EstimOpt.CONS2,[],[],[],OptimOpt);
    else
        [Results.bhat, LL, Results.exitf, Results.output, Results.lambda, Results.g] = fmincon(LLfun,b0,[],[],EstimOpt.CONS1,EstimOpt.CONS2,[],[],[],OptimOpt);
    end

end


%% Output


Results.LL = -LL;
Results.b0_old = b0;

if EstimOpt.HessEstFix == 1
	f = LL_mimic(INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,Results.bhat);...
	Results.jacobian = numdiff(@(B) LL_mimic(INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,B), f, Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
	Results.jacobian = jacobianest(@(B) LL_mimic(INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
	Results.jacobian = hessian(@(B) sum(LL_mimic(INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,B),1), Results.bhat);
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
Results.std = sqrt(diag(Results.ihess));
Results.std(EstimOpt.BActive == 0) = NaN;
Results.std(EstimOpt.BLimit == 1) = 0;

Results.DetailsS = [Results.bhat(1:(EstimOpt.NVarstr)*EstimOpt.NLatent), Results.std(1:(EstimOpt.NVarstr)*EstimOpt.NLatent), pv(Results.bhat(1:(EstimOpt.NVarstr)*EstimOpt.NLatent), Results.std(1:(EstimOpt.NVarstr)*EstimOpt.NLatent))];
Results.DetailsM = [Results.bhat((EstimOpt.NVarstr)*EstimOpt.NLatent+1:end), Results.std((EstimOpt.NVarstr)*EstimOpt.NLatent+1:end), pv(Results.bhat((EstimOpt.NVarstr)*EstimOpt.NLatent+1:end), Results.std((EstimOpt.NVarstr)*EstimOpt.NLatent+1:end))];

EstimOpt.params = length(b0);
EstimOpt.NObs = EstimOpt.NP;
if isfield(EstimOpt,'BActive')
	EstimOpt.params = EstimOpt.params - sum(EstimOpt.BActive == 0);
end
Results.stats = [Results.MIMIC0.LL; Results.LL; 1-Results.LL/Results.MIMIC0.LL; ((2*EstimOpt.params-2*Results.LL) + 2*EstimOpt.params*(EstimOpt.params+1)/(EstimOpt.NObs-EstimOpt.params-1))/EstimOpt.NObs; EstimOpt.NObs; EstimOpt.params];

Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;
Results.INPUT = INPUT;

l = 0;
for i = 1:EstimOpt.NLatent
    disp(' ')
    disp(num2str(i,'Structural equation of Latent Variable %1.0f'));
    disp(['var.', blanks(size(char(EstimOpt.NamesStr),2)-2) ,'coef.      st.err.  p-value'])
    disp([char(EstimOpt.NamesStr) ,blanks(EstimOpt.NVarstr)',num2str(Results.DetailsS(l+1:l+EstimOpt.NVarstr,1),'%8.4f'), star_sig(Results.DetailsS(l+1:l+EstimOpt.NVarstr,3)), num2str(Results.DetailsS(l+1:l+EstimOpt.NVarstr,2:3),'%8.4f %8.4f')])
    l = l + EstimOpt.NVarstr;
end

l = 0;

for i = 1:size(INPUT.Xmea,2)
    tmp = EstimOpt.NVarmea_exp*(EstimOpt.MeaExpMatrix(i) ~=0);
    disp(' ')
    disp([num2str('Measurment equation for '), char(EstimOpt.NamesMea(i))]);
    if EstimOpt.MeaSpecMatrix(i) == 0
        disp('Estimated using OLS')
        disp('var.   coef.     st.err.  p-value')
        Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1:3) = [exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)), Results.DetailsM(l+tmp+sum(EstimOpt.MeaMatrix(:,i))+2,2)*exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)),pv(exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)), Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,2)*exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)))];
        
        disp([char(EstimOpt.NamesLV(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%7.4f %8.4f')])
        l = l+sum(EstimOpt.MeaMatrix(:,i))+2+tmp;
    elseif EstimOpt.MeaSpecMatrix(i) == 1
        disp('Estimated using MNL')
        for n = 1:(length(unique(INPUT.Xmea(:,i)))-1)
            disp(num2str(n+1, 'Parameters for %1.0f alternative'))
            disp('var.   coef.     st.err.  p-value')
            disp([char(EstimOpt.NamesLV(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1), '%11.4f'), star_sig(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%8.4f %8.4f')])
            l = l + 1 + sum(EstimOpt.MeaMatrix(:,i),1)+tmp;
        end
    elseif EstimOpt.MeaSpecMatrix(i) == 2
        disp('Estimated using Ordered Probit')
        disp('var.       coef.     st.err.  p-value')
        disp([char(EstimOpt.NamesLV(l+1:l+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,repmat(blanks(sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',1,6),num2str(Results.DetailsM(l+1:l+sum(EstimOpt.MeaMatrix(:,i))+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+sum(EstimOpt.MeaMatrix(:,i))+tmp,3)), num2str(Results.DetailsM(l+1:l+sum(EstimOpt.MeaMatrix(:,i))+tmp,2:3),'%8.4f %8.4f')])
        l = l+sum(EstimOpt.MeaMatrix(:,i))+tmp;
        disp([num2str([1, Results.DetailsM(l+1,1)],'Cutoff %1.0f %7.4f'), star_sig(Results.DetailsM(l+1,3)), num2str(Results.DetailsM(l+1,2:3),'%8.4f %8.4f')])
        if length(unique(INPUT.Xmea(:,i))) > 2 % if attitude is not binary 
            g = [Results.DetailsM(l+1,1) ; exp(Results.DetailsM(l+2:l+length(unique(INPUT.Xmea(:,i)))-1,1))];
            for n = 2:length(unique(INPUT.Xmea(:,i)))-1
                stdx = sqrt(g(1:n)'*Results.ihess((EstimOpt.NVarstr)*EstimOpt.NLatent+l+1:(EstimOpt.NVarstr)*EstimOpt.NLatent+l+n,(EstimOpt.NVarstr)*EstimOpt.NLatent+ l+1:(EstimOpt.NVarstr)*EstimOpt.NLatent+l+n)*g(1:n));
                Results.DetailsM(l+n,1:3) = [sum(g(1:n),1), stdx, pv(sum(g(1:n),1), stdx)];
                disp([num2str([n, Results.DetailsM(l+n,1)],'Cutoff %1.0f %7.4f'), star_sig(Results.DetailsM(l+n,3)), num2str(Results.DetailsM(l+n,2:3),'%8.4f %8.4f')])
            end
        end       
        l = l+length(unique(INPUT.Xmea(:,i)))-1;
    elseif EstimOpt.MeaSpecMatrix(i) == 3
        disp('Estimated using Poisson regression')
        disp('var.   coef.     st.err.  p-value')
        
        disp([char(EstimOpt.NamesLV(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%7.4f %8.4f')])
        l = l+sum(EstimOpt.MeaMatrix(:,i))+1+tmp;
    elseif EstimOpt.MeaSpecMatrix(i) == 4
        disp('Estimated using Negative Binomial regression')
        disp('var.   coef.     st.err.  p-value')
        Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1:3) = [exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)), Results.DetailsM(l+tmp+sum(EstimOpt.MeaMatrix(:,i))+2,2)*exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)),pv(exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)), Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,2)*exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)))];

        disp([char(EstimOpt.NamesLV(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%7.4f %8.4f')])
        l = l+sum(EstimOpt.MeaMatrix(:,i))+2+tmp;
    elseif EstimOpt.MeaSpecMatrix(i) == 5
        disp('Estimated using Zero Inflated Poisson regression')
        disp('Probability of Non-participation (logit)')
        disp('var.   coef.     st.err.  p-value')

        disp([char(EstimOpt.NamesLV(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%7.4f %8.4f')])
        l = l+sum(EstimOpt.MeaMatrix(:,i))+1+tmp;
        
        disp('Poisson model')
        disp('var.   coef.     st.err.  p-value')
        disp([char(EstimOpt.NamesLV(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%7.4f %8.4f')])        
        l = l+sum(EstimOpt.MeaMatrix(:,i))+1+tmp;
    end
    
end

Results.R = [Results.DetailsS;Results.DetailsM]; 
Results.R_out  = cell(3+EstimOpt.NVarstr + 3*size(INPUT.Xmea,2)+EstimOpt.NVarcut+EstimOpt.NVarmea+ 2 + 6, 1+3*EstimOpt.NLatent);  
Results.R_out(1,1) = {'MIMIC'};
head = {'var.' , 'coef.', 'st.err.' , 'p-value'};
headx = [head, repmat(head(1,2:4),1,EstimOpt.NLatent-1)];
LVlist = {'LV 1' , 'LV 2' , 'LV 3' , 'LV 4' , 'LV 5' , 'LV 6' , 'LV 7' , 'LV 8' ,'LV 9' ,'LV 10'};
if EstimOpt.NLatent<=10
    Results.R_out(2,2:3:(2+(EstimOpt.NLatent-1)*3)) = LVlist(1,1:EstimOpt.NLatent);
end
Results.R_out(3,:) = headx;
Results.R_out(4:(EstimOpt.NVarstr+3),1) = EstimOpt.NamesStr;
for i = 1: EstimOpt.NLatent
    Results.R_out(4:(EstimOpt.NVarstr+3),2 + 3*(i-1):1+3*i) = num2cell(Results.DetailsS((i-1)*EstimOpt.NVarstr+1:i*EstimOpt.NVarstr,:));
end
l = 0; % this is for indexing in R_out
k = 0;
for i = 1:size(INPUT.Xmea,2)
    Results.R_out(EstimOpt.NVarstr+3+l+1,1) =  cellfun(@(x)[x char(EstimOpt.NamesMea(i))],{'Measurment equation for '},'UniformOutput',0);
    if EstimOpt.MeaSpecMatrix(i) == 0 
        model = 'OLS';
    elseif EstimOpt.MeaSpecMatrix(i) == 1
        model = 'MNL';
    elseif EstimOpt.MeaSpecMatrix(i) == 2
        model = 'OP';
    elseif EstimOpt.MeaSpecMatrix(i) == 3
        model = 'POISS';
    elseif EstimOpt.MeaSpecMatrix(i) == 4
        model = 'NB';
    elseif EstimOpt.MeaSpecMatrix(i) == 5
        model = 'ZIP';
    end
    Results.R_out(EstimOpt.NVarstr+3+l+2,1) = cellfun(@(x)[x model],{'Estimated using '},'UniformOutput',0);
    Results.R_out(EstimOpt.NVarstr+3+l+3,1:4) = head;
    Results.R_out(EstimOpt.NVarstr+3+l+4:EstimOpt.NVarstr+3+l+3 + EstimOpt.CutMatrix(i) ,1:4) = [EstimOpt.NamesLV(k+1:k+EstimOpt.CutMatrix(i)), num2cell(Results.DetailsM(k+1:k+EstimOpt.CutMatrix(i),:))];
    l = l+3+EstimOpt.CutMatrix(i);
    k = k + EstimOpt.CutMatrix(i);
end
Results.R_out(EstimOpt.NVarstr+3+l+2,1) = {'Model characteristics'};
Results.R_out(EstimOpt.NVarstr+3+l+3:end,1) = {'LL0'; 'LL' ; 'Pseudo R2'; 'AIC/n' ; 'n'; 'k'};
Results.R_out(EstimOpt.NVarstr+3+l+3:end,2) = num2cell(Results.stats);

disp(' ')
disp(['LL at convergence: ',num2str(Results.LL,'%8.4f')])
disp(' ');
Results.clocknote = clock;
Results.tocnote = toc;
[~,DayName] = weekday(now,'long');
disp(['Estimation completed on ' DayName ', ' num2str(Results.clocknote(1)) '-' sprintf('%02.0f',Results.clocknote(2)) '-' sprintf('%02.0f',Results.clocknote(3)) ' at ' sprintf('%02.0f',Results.clocknote(4)) ':' sprintf('%02.0f',Results.clocknote(5)) ':' sprintf('%02.0f',Results.clocknote(6))])
disp(['Estimation took ' num2str(Results.tocnote) ' seconds ('  num2str(floor(Results.tocnote/(60*60))) ' hours ' num2str(floor(rem(Results.tocnote,60*60)/60)) ' minutes ' num2str(rem(Results.tocnote,60)) ' seconds).']);
disp(' ');