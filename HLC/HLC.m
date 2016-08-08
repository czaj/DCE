function Results = HLC(INPUT,Results_old,EstimOpt,OptimOpt)


global B_backup;

% save tmp_HLC
% return

tic

Results.bhat = [];
Results.R = [];
Results.R_out = {};
Results.stats = [];


%% Check data and inputs


if nargin < 3 % check no. of inputs
    error('Too few input arguments for HLC(INPUT,EstimOpt,OptimOpt)')
end

disp(' ');
disp('__________________________________________________________________________________________________________________');
disp(' ');

warning off MATLAB:mir_warning_maybe_uninitialized_temporary

format shortG;
format compact;

if isfield(EstimOpt,'NClass') == 0;
    EstimOpt.NClass = 2;
end

if isfield(EstimOpt,'NLatent') == 0
    EstimOpt.NLatent = 1;
    disp(['Assuming ', num2str(EstimOpt.NLatent)',' Latent Variable(s)']);
end


if any(INPUT.W ~= 1)
    cprintf('Black','Estimating '); cprintf('*Black','weighted '); cprintf('Black','HLC model with '); cprintf('Black',num2str(EstimOpt.NLatent)); cprintf('Black',' Latent Variable(s) and '); cprintf('Black',num2str(EstimOpt.NClass)); cprintf('Black',' classes...\n');
else
    disp(num2str([EstimOpt.NLatent, EstimOpt.NClass],'Estimating HLC model with %1.0f Latent Variable(s) and %1.0f classes ...'))
end

if isfield(EstimOpt, 'WTP_space') == 0
    EstimOpt.WTP_space = 0;
    EstimOpt.WTP_matrix = [];
elseif EstimOpt.WTP_space == 0;
    EstimOpt.WTP_matrix = [];
end

if EstimOpt.WTP_space > 0
    disp('in WTP-space.')
else
    disp('in preference-space.')
end

if EstimOpt.WTP_space > 0
    if isfield(EstimOpt, 'WTP_matrix') == 0
        WTP_att = (EstimOpt.NVarA-EstimOpt.WTP_space)/EstimOpt.WTP_space;
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

if isfield(INPUT, 'Xc') == 0 || numel(INPUT.Xc) == 0
    INPUT.Xc = ones(size(INPUT.Y,1),1);
else
    INPUT.Xc = [ones(size(INPUT.Y,1),1), INPUT.Xc];
end

EstimOpt.NVarC = size(INPUT.Xc,2); % no. of variables explaining class probabilities (without counting LV)

if isfield(INPUT, 'Xstr') == 0 || numel(INPUT.Xstr) == 0
    % 	error('Define variables to structural equations')
    cprintf(rgb('DarkOrange'), 'WARNING: Structural equations empty. \n')
    INPUT.Xstr = zeros(size(INPUT.Y,1),0);
end
if isfield(INPUT, 'Xmea') == 0 || numel(INPUT.Xmea) == 0
    % 	error('Define attitudes to measurment equations')
    cprintf(rgb('DarkOrange'), 'WARNING: Measurement equations empty.\n')
    INPUT.Xmea = zeros(size(INPUT.Y,1),0);
end

if isfield(INPUT, 'Xmea_exp') == 0 || numel(INPUT.Xmea_exp) == 0 % additional covariates for explaining Measurment equations
    INPUT.Xmea_exp = [];
    EstimOpt.MeaExpMatrix = zeros(1,size(INPUT.Xmea,2));
else
    if isfield(EstimOpt,'MeaExpMatrix') == 0 || length(EstimOpt.MeaExpMatrix) ~= size(INPUT.Xmea,2)
        EstimOpt.MeaExpMatrix = ones(1, size(INPUT.Xmea,2));
        cprintf(rgb('DarkOrange'), 'WARNING: MeaExpMatrix not defined - assuming that every measurment equation is explained with additional covariates /n')
    else
        EstimOpt.MeaExpMatrix = EstimOpt.MeaExpMatrix(:)';
    end
end

if isfield(EstimOpt,'MeaMatrix')
    if size(EstimOpt.MeaMatrix,2) ~= size(INPUT.Xmea,2) || size(EstimOpt.MeaMatrix,1) ~= EstimOpt.NLatent
        error('Measurment equations erroneously defined (incorrect size of the MeaMatrix)')
        %         EstimOpt = rmfield(EstimOpt,'MeaMatrix');
    elseif ismember(0,(ismember(EstimOpt.MeaMatrix,[0,1])))
        error('Measurment equations erroneously defined (incorrect values in the MeaMatrix)')
        %         EstimOpt = rmfield(EstimOpt,'MeaMatrix');
    elseif any(any(EstimOpt.MeaMatrix == 1,2) == 0)
        %         error('Measurment equations erroneously defined (some Latent Variables without measurement equations in the MeaMatrix)')
        cprintf(rgb('DarkOrange'), 'WARNING: Some of the LV not associated with any measurement equations.\n')
        %         EstimOpt = rmfield(EstimOpt,'MeaMatrix');
    elseif any(any(EstimOpt.MeaMatrix == 1) == 0)
        error('Measurment equations erroneously defined (some measurement variables unused)')
        %         EstimOpt = rmfield(EstimOpt,'MeaMatrix');
    end
elseif size(INPUT.Xmea,2) == 0
    EstimOpt.MeaMatrix = ones(EstimOpt.NLatent,size(INPUT.Xmea,2));
else
    EstimOpt.MeaMatrix = ones(EstimOpt.NLatent,size(INPUT.Xmea,2));
    disp('Assuming every Latent Variable in every measurment equation')
end

if isfield(EstimOpt,'MeaSpecMatrix') == 0 || numel(EstimOpt.MeaSpecMatrix) == 0
    EstimOpt.MeaSpecMatrix = zeros(1, size(INPUT.Xmea,2));
    if size(INPUT.Xmea,2) > 0
        disp('Using OLS for every measurment equation')
    end
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
        if EstimOpt.MeaSpecMatrix(i) > 0 && numel(unique(INPUT.Xmea(INPUT.MissingInd==0,i))) > 10
            cprintf(rgb('DarkOrange'), 'WARNING: There are over 10 levels for measurement variable %d \n', i)
        end
    end
    if sum(isnan(INPUT.Xmea(INPUT.MissingInd==0,i))) > 0
        cprintf(rgb('DarkOrange'), 'WARNING: Measurement variable %d contains NaN values \n', i)
    end
    if sum(isinf(INPUT.Xmea(INPUT.MissingInd==0,i))) > 0
        cprintf(rgb('DarkOrange'), 'WARNING: Measurement variable %d contains Inf values \n', i)
    end
end

for i=1:EstimOpt.NVarstr
    if sum(isnan(INPUT.Xstr(INPUT.MissingInd==0,i))) > 0
        cprintf(rgb('DarkOrange'), 'WARNING: Structural variable %d contains NaN values \n', i)
    end
    if sum(isinf(INPUT.Xstr(INPUT.MissingInd==0,i))) > 0
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
    elseif size(EstimOpt.NamesC,1) ~= EstimOpt.NVarC
        EstimOpt.NamesC = EstimOpt.NamesC';
    end
    EstimOpt.NamesC = [{'Cons'};EstimOpt.NamesC];
else
    EstimOpt.NamesC = {'Cons'};
end
for i = 1:EstimOpt.NLatent
    EstimOpt.NamesC =  [EstimOpt.NamesC; cellfun(@(x)[x num2str(i)],{'LV '},'UniformOutput',0)];
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

if size(INPUT.Xmea,2) > 0
    disp(['Using following models for attitudes: '  char(EstimOpt.Names)]);
end

%% Rescructure data


INPUT.YYY = reshape(INPUT.Y,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP);
idx = sum(reshape(INPUT.MissingInd,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP)) == EstimOpt.NAlt; ...
    INPUT.YYY(idx(ones(EstimOpt.NAlt,1),:) == 1) = NaN; % replace YYY in missing choice-tasks with NaN
INPUT.YY = reshape(INPUT.YYY(:,:,ones(EstimOpt.NClass,1)),EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP*EstimOpt.NClass,1);

INPUT.Xa(INPUT.MissingInd == 1,:) = NaN;

INPUT.Xstr = double(INPUT.Xstr(1:EstimOpt.NAlt*EstimOpt.NCT:end,:));
% normalize explanatory variables for structural equations:
if ~isfield(EstimOpt,'StrNorm')
    EstimOpt.StrNorm = ones(EstimOpt.NVarStr,1);
elseif length(EstimOpt.StrNorm) ~= EstimOpt.NVarStr && length(EstimOpt.StrNorm) ~= 1
    cprintf(rgb('DarkOrange'), 'WARNING:  Structural variables normalization options (EstimOpt.StrNorm) incorrect - variables not normalized \n')
elseif length(EstimOpt.StrNorm) == 1
    EstimOpt.StrNorm = EstimOpt.StrNorm(ones(EstimOpt.NVarStr,1),1);
end
if any(EstimOpt.StrNorm > 0)
    meanx = mean(INPUT.Xstr);
    stdx = std(INPUT.Xstr);
    INPUT.Xstr(:,EstimOpt.StrNorm == 1) = (INPUT.Xstr(:,EstimOpt.StrNorm == 1) - meanx(ones(size(INPUT.Xstr,1),1),EstimOpt.StrNorm == 1))./stdx(ones(size(INPUT.Xstr,1),1),EstimOpt.StrNorm == 1);
end

INPUT.Xmea = INPUT.Xmea(1:EstimOpt.NAlt*EstimOpt.NCT:end,:);
INPUT.XXc = INPUT.Xc(1:EstimOpt.NCT*EstimOpt.NAlt:end,:); % NP x NVarC
if EstimOpt.NVarMeaExp > 0
    INPUT.Xmea_exp = INPUT.Xmea_exp(1:EstimOpt.NAlt*EstimOpt.NCT:end,:);  
    % normalize explanatory variables for measurement equations:
    if ~isfield(EstimOpt, 'MeaExpNorm')
        EstimOpt.MeaExpNorm = ones(EstimOpt.NVarMeaExp,1);
    elseif length(EstimOpt.MeaExpNorm) ~= EstimOpt.NVarMeaExp && length(EstimOpt.MeaExpNorm) ~= 1
        cprintf(rgb('DarkOrange'), 'WARNING: Explanatory measurement variables normalization options (EstimOpt.MeaExpNorm) incorrect - variables not normalized \n')
    elseif length(EstimOpt.MeaExpNorm) == 1
        EstimOpt.MeaExpNorm = EstimOpt.MeaExpNorm(ones(EstimOpt.NVarMeaExp,1),1);
    end
    if any(EstimOpt.MeaExpNorm > 0)
        meanx = mean(INPUT.Xmea_exp);
        stdx = std(INPUT.Xmea_exp);
        INPUT.Xmea_exp(:,EstimOpt.MeaExpNorm == 1) = (INPUT.Xmea_exp(:,EstimOpt.MeaExpNorm == 1) - meanx(ones(size(INPUT.Xmea_exp,1),1),EstimOpt.MeaExpNorm == 1))./stdx(ones(size(INPUT.Xmea_exp,1),1),EstimOpt.MeaExpNorm == 1);
    end
end


%% Estimating MIMIC0


if size(INPUT.Xmea,2) > 0
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
else
    Results.MIMIC0.bhat = [];
    Results.MIMIC0.LL = 0;
end


%% Starting values


EstimOpt.jitter1 = 0.8; % Jittering parameter (relative) for HMNL starting values (attributes)
EstimOpt.jitter2 = 0.1; % Jittering parameter (absolute) for class probabilities starting values

if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == (EstimOpt.NVarA*EstimOpt.NClass + (EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass - 1) + EstimOpt.NVarstr*EstimOpt.NLatent + EstimOpt.NVarmea + EstimOpt.NVarcut)
    b0 = B_backup(:);
    disp('Using the starting values from Backup')
elseif isfield(Results_old,'HLC') && isfield(Results_old.HLC,'b0') % starting values provided
    Results_old.HLC.b0_old = Results_old.HLC.b0(:);
    Results_old.HLC = rmfield(Results_old.HLC,'b0');
    if length(Results_old.HLC.b0_old) ~= (EstimOpt.NVarA*EstimOpt.NClass + (EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass - 1) + EstimOpt.NVarstr*EstimOpt.NLatent + EstimOpt.NVarmea + EstimOpt.NVarcut)
        cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
        Results_old.HLC = rmfield(Results_old.HLC,'b0_old');
    else
        b0 = Results_old.HLC.b0_old(:);
    end
end

if  ~exist('b0','var')
    if isfield(Results_old,'MIMIC') && isfield(Results_old.MIMIC,'bhat')
        disp('Using MIMIC results as starting values')
        Results_old.MIMIC.bhat = Results_old.MIMIC.bhat(:);
        bstr = reshape(Results_old.MIMIC.bhat(1:EstimOpt.NVarstr*EstimOpt.NLatent), EstimOpt.NVarstr, EstimOpt.NLatent);
        LV = INPUT.Xstr*bstr;
        LV = LV(:,:, ones(EstimOpt.NCT*EstimOpt.NAlt,1)); % NP x NLatent x NCT*NAlt
        LV = permute(LV, [3 1 2]);
        LV = reshape(LV, EstimOpt.NCT*EstimOpt.NAlt*EstimOpt.NP,EstimOpt.NLatent);
        disp('Estimating additional MNL for starting values')
        EstimOpt0 = EstimOpt;
        EstimOpt0.DISPLAY = 0;
        EstimOpt0.BActive =  ones(1,EstimOpt.NVarA);
        EstimOpt0.NClass = EstimOpt.NClass;
        EstimOpt0.WTP_space = EstimOpt.WTP_space;
        EstimOpt0.WTP_matrix = EstimOpt.WTP_matrix;
        EstimOpt0.NVarA = EstimOpt.NVarA;
        Results.MNL = MNL(INPUT,[],EstimOpt0,OptimOpt);
        disp('Estimating additional LC for starting values')
        if size(INPUT.Xc,2) > 1
            INPUT.Xc = [INPUT.Xc(:,2:end), LV];
        else
            INPUT.Xc = LV;
        end
        if isfield(EstimOpt, 'BActive') == 1 && numel(EstimOpt.BActive) ~= 0
            EstimOpt0.BActive =  EstimOpt.BActive(1:EstimOpt.NVarA*EstimOpt.NClass+(EstimOpt.NClass-1)*(EstimOpt.NVarC+EstimOpt.NLatent));
        else
            EstimOpt0.BActive =  ones(1,EstimOpt.NVarA*EstimOpt.NClass+(EstimOpt.NClass-1)*(EstimOpt.NVarC+EstimOpt.NLatent));
        end
        Results.LC_LV = LC(INPUT,Results_old,EstimOpt0,OptimOpt);
        b0 = [Results.LC_LV.bhat; Results_old.MIMIC.bhat];
    elseif isfield(Results_old,'HMNL') && isfield(Results_old.HMNL,'bhat')
        disp('Using HMNL results as starting values')
        Results_old.HMNL.bhat = Results_old.HMNL.bhat(:);
        if EstimOpt.WTP_space > 0
            Results_old.HMNL.bhat(1:(EstimOpt.NVarA-EstimOpt.WTP_space)) = Results_old.HMNL.bhat(1:(EstimOpt.NVarA-EstimOpt.WTP_space))./Results_old.HMNL.bhat(EstimOpt.WTP_matrix);
        end
        b0 = [Results_old.HMNL.bhat((1:EstimOpt.NVarA)'*ones(1,EstimOpt.NClass),:) .* ...
            (EstimOpt.jitter1 .* unifrnd(0,ones(EstimOpt.NVarA.*EstimOpt.NClass,1))) ; ...
            EstimOpt.jitter2 + unifrnd(-1,ones((EstimOpt.NVarC+EstimOpt.NLatent).*(EstimOpt.NClass-1),1)); ...
            Results_old.HMNL.bhat(EstimOpt.NVarA*(EstimOpt.NLatent+1)+1:end)];
    else
        error('No starting values available - run MIMIC or HMNL first')
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
    err_mtx = randn(EstimOpt.NP*EstimOpt.NRep, EstimOpt.NLatent); %to be cut down later
elseif EstimOpt.Draws == 2 % LHS
    cprintf('*blue','Latin Hypercube Sampling '); cprintf('draws \n');
    err_mtx=lhsnorm(zeros((EstimOpt.NLatent)*EstimOpt.NP,1),diag(ones((EstimOpt.NLatent)*EstimOpt.NP,1)),EstimOpt.NRep);
    err_mtx = reshape(err_mtx, EstimOpt.NRep*EstimOpt.NP, EstimOpt.NLatent);
elseif EstimOpt.Draws >= 3 % Quasi random draws
    if EstimOpt.Draws == 3
        cprintf('*blue','Halton '); cprintf('draws (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = haltonset(EstimOpt.NLatent,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); %
    elseif EstimOpt.Draws == 4 % apply reverse-radix scrambling
        cprintf('*blue','Halton '); cprintf('draws with reverse radix scrambling (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = haltonset(EstimOpt.NLatent,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); %
        hm1 = scramble(hm1,'RR2');
    elseif EstimOpt.Draws == 5
        cprintf('*blue','Sobol '); cprintf('draws (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = sobolset(EstimOpt.NLatent,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
    elseif EstimOpt.Draws == 6
        cprintf('*blue','Sobol '); cprintf('draws with random linear scramble and random digital shift (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = sobolset(EstimOpt.NLatent,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
        hm1 = scramble(hm1,'MatousekAffineOwen');
    end
    
    err_mtx = net(hm1,EstimOpt.NP*EstimOpt.NRep); % this takes every point:
    clear hm1;
    if EstimOpt.NP*EstimOpt.NRep < 3e+7
        err_mtx= icdf('Normal',err_mtx,0,1); %to be cut down later
    else % this is for very large number of draws * variables
        for i=1:EstimOpt.NLatent
            err_mtx(:,i) = icdf('Normal',err_mtx(:,i),0,1); %to be cut down later
        end
    end
    
end

err_sliced = err_mtx'; % NLatent x NRep * NP


%% Display Options


if EstimOpt.NumGrad == 0
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical -  analytical gradient not supported in the HLC model \n')
end

if ((isfield(EstimOpt, 'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    OptimOpt.GradObj = 'on';
end

if (isfield(EstimOpt, 'ConstVarActive') == 0 || EstimOpt.ConstVarActive == 0) && isequal(OptimOpt.Algorithm,'quasi-newton') && isequal(OptimOpt.Hessian,'user-supplied')
    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied Hessian off - quasi-newton algorithm does not use it anyway \n')
    OptimOpt.Hessian = 'off';
end
if EstimOpt.RobustStd == 1 && (EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2)
    EstimOpt.RobustStd = 0;
    cprintf(rgb('DarkOrange'), 'WARNING: Setting off robust standard errors, they do not matter for BHHH aproximation of hessian \n')
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


LLfun = @(B) LL_hlc_MATlike(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp, err_sliced,INPUT.W,EstimOpt,OptimOpt,B);
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


% save tmp_HLC

Results.LL = -LL;
R2 = R2_hybrid(INPUT.YY,INPUT.Xa,[],INPUT.Xstr,INPUT.XXc,INPUT.MissingInd,err_sliced,EstimOpt,Results.bhat,1);

Results.b0_old = b0;

if EstimOpt.HessEstFix == 1
    f = LL_hlc(INPUT.YY,INPUT.Xa, INPUT.XXc, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,Results.bhat);...
        Results.jacobian = numdiff(@(B) INPUT.W.*LL_hlc(INPUT.YY,INPUT.Xa, INPUT.XXc, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,B), INPUT.W.*f, Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) INPUT.W.*LL_hlc(INPUT.YY,INPUT.Xa, INPUT.XXc, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.jacobian = hessian(@(B) sum(INPUT.W.*LL_hlc(INPUT.YY,INPUT.Xa, INPUT.XXc, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,B),1), Results.bhat);
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
        [~, Results.jacobian] = LL_hlc(INPUT.YY,INPUT.Xa, INPUT.XXc, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,Results.bhat);
        Results.jacobian = Results.jacobian.*INPUT.W(:, ones(1,size(Results.jacobian,2)));
    else
        Results.LLdetailed = LL_hlc(INPUT.YY,INPUT.Xa, INPUT.XXc, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,Results.bhat);
        Results.jacobian = numdiff(@(B) INPUT.W.*LL_hlc(INPUT.YY,INPUT.Xa, INPUT.XXc, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,B) ,INPUT.W.*Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
    end
    RobustHess = Results.jacobian'*Results.jacobian;
    Results.ihess = Results.ihess*RobustHess*Results.ihess;
end

Results.std = sqrt(diag(Results.ihess));
Results.std(EstimOpt.BActive == 0) = NaN;
Results.std(EstimOpt.BLimit == 1) = 0;

Results.DetailsA = [Results.bhat(1:EstimOpt.NClass*EstimOpt.NVarA), Results.std(1:EstimOpt.NClass*EstimOpt.NVarA), pv(Results.bhat(1:EstimOpt.NClass*EstimOpt.NVarA), Results.std(1:EstimOpt.NClass*EstimOpt.NVarA))];
l = EstimOpt.NClass*EstimOpt.NVarA;
Results.DetailsV = [Results.bhat(l+1:l+(EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1)), Results.std(l+1:l+(EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1)), pv(Results.bhat(l+1:l+(EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1)), Results.std(l+1:l+(EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1)))];
l = l+(EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1);
Results.DetailsS = [Results.bhat(l+1:l+EstimOpt.NLatent*EstimOpt.NVarstr), Results.std(l+1:l+EstimOpt.NLatent*EstimOpt.NVarstr), pv(Results.bhat(l+1:l+EstimOpt.NLatent*EstimOpt.NVarstr), Results.std(l+1:l+EstimOpt.NLatent*EstimOpt.NVarstr))];
l = l + EstimOpt.NLatent*EstimOpt.NVarstr;
Results.DetailsM = [Results.bhat(l+1:end), Results.std(l+1:end), pv(Results.bhat(l+1:end), Results.std(l+1:end))];

Results.LL0 = Results.MIMIC0.LL + Results_old.MNL0.LL;
EstimOpt.params = length(Results.bhat);
if isfield(EstimOpt,'BActive')
    EstimOpt.params = EstimOpt.params - sum(EstimOpt.BActive == 0);
end
Results.stats = [Results.LL0; Results.LL; 1-Results.LL/Results.LL0;R2; ((2*EstimOpt.params-2*Results.LL) + 2*EstimOpt.params*(EstimOpt.params+1)/(EstimOpt.NObs-EstimOpt.params-1))/EstimOpt.NObs; EstimOpt.NObs; EstimOpt.params];

Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;
Results.INPUT = INPUT;

for j = 1:EstimOpt.NClass
    disp(' ')
    disp(num2str(j,'Latent class parameters %1.0f'));
    disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
    disp([char(EstimOpt.NamesA) ,blanks(EstimOpt.NVarA)',num2str(Results.DetailsA(((j-1)*EstimOpt.NVarA+1):j*EstimOpt.NVarA,1),'%8.4f'), star_sig(Results.DetailsA(((j-1)*EstimOpt.NVarA+1):j*EstimOpt.NVarA,3)), num2str(Results.DetailsA(((j-1)*EstimOpt.NVarA+1):j*EstimOpt.NVarA,2:3),'%8.4f %8.4f')])
end
l = 0;
for j = 1:EstimOpt.NClass -1
    disp(' ')
    disp(num2str(j,'Latent class probability model %1.0f'))
    disp('var.  coef.     st.err.  p-value')
    disp([char(EstimOpt.NamesC) ,blanks(EstimOpt.NVarC + EstimOpt.NLatent)', num2str(Results.DetailsV(l+1:l+EstimOpt.NVarC+EstimOpt.NLatent,1),'%8.4f'), star_sig(Results.DetailsV(l+1:l+EstimOpt.NVarC+EstimOpt.NLatent,3)), num2str(Results.DetailsV(l+1:l+EstimOpt.NVarC+EstimOpt.NLatent,2:3),'%8.4f %8.4f')])
    l = l+EstimOpt.NVarC+EstimOpt.NLatent;
end

bclass = reshape([Results.bhat(EstimOpt.NVarA*EstimOpt.NClass +1 :EstimOpt.NVarA*EstimOpt.NClass + (EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass - 1)); zeros(EstimOpt.NVarC+EstimOpt.NLatent,1)], EstimOpt.NVarC+EstimOpt.NLatent, EstimOpt.NClass);
bstr = reshape(Results.bhat(EstimOpt.NVarA*EstimOpt.NClass + (EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass - 1)+1:EstimOpt.NVarA*EstimOpt.NClass + (EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass - 1) + EstimOpt.NLatent*EstimOpt.NVarstr), EstimOpt.NVarstr, EstimOpt.NLatent);
LV = INPUT.Xstr*bstr;
V = exp([INPUT.XXc, LV]*bclass);% NP x NClass
Vsum = sum(V,2);
Results.PNClass = mean(V./Vsum(:,ones(EstimOpt.NClass,1)),1); %1 x NClass
disp(' ')
disp(['Avarage class probabilities: ', num2str(Results.PNClass)])

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


%
Results.R = [Results.DetailsA; Results.DetailsV ; Results.DetailsS;Results.DetailsM];

Results.R_out  = cell(3+EstimOpt.NVarA + 2+ EstimOpt.NVarC + EstimOpt.NLatent+ 3 +4+EstimOpt.NVarstr +  3*size(INPUT.Xmea,2)+EstimOpt.NVarcut+EstimOpt.NVarmea+ 2 + 7, 4+3*max(EstimOpt.NClass-1,EstimOpt.NLatent -1));
Results.R_out(1,1) = {'HLC'};
head = {'var.' , 'coef.', 'st.err.' , 'p-value'};
NClasses = {'NClass 1','NClass 2', 'NClass 3', 'NClass 4', 'NClass 5', 'NClass 6', 'NClass 7', 'NClass 8', 'NClass 9','NClass 10'};
LVlist = {'LV 1' , 'LV 2' , 'LV 3' , 'LV 4' , 'LV 5' , 'LV 6' , 'LV 7' , 'LV 8' ,'LV 9' ,'LV 10'};
if EstimOpt.NClass<=10
    Results.R_out(2,2:3:(2+(EstimOpt.NClass-1)*3)) = NClasses(1,1:EstimOpt.NClass);
end
Results.R_out(3,1:1+3*EstimOpt.NClass) = [head, repmat(head(1,2:4),1,EstimOpt.NClass-1)];
Results.R_out(4:(EstimOpt.NVarA+3),1) = EstimOpt.NamesA;
for i = 1: EstimOpt.NClass
    Results.R_out(4:(EstimOpt.NVarA+3),2 + 3*(i-1):1+3*i) = num2cell(Results.DetailsA((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA,:));
end
Results.R_out(EstimOpt.NVarA+4,1) = {'Latent class probability models'};
Results.R_out(EstimOpt.NVarA+5,1:4+3*(EstimOpt.NClass-2)) = [head, repmat(head(1,2:4),1,EstimOpt.NClass-2)];
Results.R_out(EstimOpt.NVarA+6:EstimOpt.NVarA+5+EstimOpt.NVarC+EstimOpt.NLatent,1) = EstimOpt.NamesC;
for i = 1: EstimOpt.NClass-1
    Results.R_out(EstimOpt.NVarA+6:EstimOpt.NVarA+5+EstimOpt.NVarC+EstimOpt.NLatent,2+ 3*(i-1):1+3*i) = num2cell(Results.DetailsV((i-1)*(EstimOpt.NVarC+EstimOpt.NLatent)+1:i*(EstimOpt.NVarC+EstimOpt.NLatent),:));
end
Results.R_out(EstimOpt.NVarA+5+EstimOpt.NVarC+EstimOpt.NLatent+2,1) = {'Average class probabilities'};
Results.R_out(EstimOpt.NVarA+5+EstimOpt.NVarC+EstimOpt.NLatent+3,1:EstimOpt.NClass) = num2cell(Results.PNClass);

Results.R_out(EstimOpt.NVarA+5+EstimOpt.NVarC+EstimOpt.NLatent+5,1) = {'Structural equations'};
for i = 1:EstimOpt.NLatent
    Results.R_out(EstimOpt.NVarA+5+EstimOpt.NVarC+EstimOpt.NLatent+6,2:3:(2+(EstimOpt.NLatent-1)*3)) = LVlist(1,1:EstimOpt.NLatent);
end
Results.R_out(EstimOpt.NVarA+5+EstimOpt.NVarC+EstimOpt.NLatent+7,1:(1+3*EstimOpt.NLatent)) = [head, repmat(head(1,2:4),1,EstimOpt.NLatent-1)];
Results.R_out(EstimOpt.NVarA+5+EstimOpt.NVarC+EstimOpt.NLatent+8:EstimOpt.NVarA+5+EstimOpt.NVarC+EstimOpt.NLatent+7+EstimOpt.NVarstr,1) = EstimOpt.NamesStr;
for i = 1: EstimOpt.NLatent
    Results.R_out(EstimOpt.NVarA+5+EstimOpt.NVarC+EstimOpt.NLatent+8:EstimOpt.NVarA+5+EstimOpt.NVarC+EstimOpt.NLatent+7+EstimOpt.NVarstr,2 + 3*(i-1):1+3*i) = num2cell(Results.DetailsS((i-1)*EstimOpt.NVarstr+1:i*EstimOpt.NVarstr,:));
end
l = EstimOpt.NVarA+5+EstimOpt.NVarC+EstimOpt.NLatent+4; % this is for indexing in R_out
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
Results.R_out(EstimOpt.NVarstr+3+l+3:end,1) = {'LL0'; 'LL' ; 'McFadden R2';'Ben-Akiva R2' ;'AIC/n' ; 'n'; 'k'};
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