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

if isfield(EstimOpt,'NClass') == 0
    EstimOpt.NClass = 2;
end

if isfield(EstimOpt,'NLatent') == 0
    EstimOpt.NLatent = 1;
    disp(['Assuming ',num2str(EstimOpt.NLatent)',' Latent Variable(s)']);
end


if any(INPUT.W ~= 1)
    cprintf('Black','Estimating '); cprintf('*Black','weighted '); cprintf('Black','HLC model with '); cprintf('Black',num2str(EstimOpt.NLatent)); cprintf('Black',' Latent Variable(s) and '); cprintf('Black',num2str(EstimOpt.NClass)); cprintf('Black',' classes...\n');
else
    disp(num2str([EstimOpt.NLatent, EstimOpt.NClass],'Estimating HLC model with %1.0f Latent Variable(s) and %1.0f classes ...'))
end

if isfield(EstimOpt,'WTP_space') == 0
    EstimOpt.WTP_space = 0;
    EstimOpt.WTP_matrix = [];
elseif EstimOpt.WTP_space == 0
    EstimOpt.WTP_matrix = [];
end

if EstimOpt.WTP_space > 0
    disp('in WTP-space.')
else
    disp('in preference-space.')
end

if EstimOpt.WTP_space > 0
    if isfield(EstimOpt,'WTP_matrix') == 0
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

if isfield(INPUT,'Xc') == 0 || numel(INPUT.Xc) == 0
    INPUT.Xc = ones(size(INPUT.Y,1),1);
else
    INPUT.Xc = [ones(size(INPUT.Y,1),1),INPUT.Xc];
end

EstimOpt.NVarC = size(INPUT.Xc,2); % no. of variables explaining class probabilities (without counting LV)

if isfield(INPUT,'Xstr') == 0 || numel(INPUT.Xstr) == 0
    cprintf(rgb('DarkOrange'),'WARNING: Structural equations empty. \n')
    INPUT.Xstr = zeros(size(INPUT.Y,1),0);
end
if isfield(INPUT,'Xmea') == 0 || numel(INPUT.Xmea) == 0
    cprintf(rgb('DarkOrange'),'WARNING: Measurement equations empty.\n')
    INPUT.Xmea = zeros(size(INPUT.Y,1),0);
elseif size(INPUT.Xmea,1) ~= size(INPUT.Xa)
    error('Incorrect number of observations for measurement equations.')
end

if isfield(INPUT,'Xmea_exp') == 0 || numel(INPUT.Xmea_exp) == 0 % additional covariates for explaining Measurment equations
    INPUT.Xmea_exp = [];
    EstimOpt.MeaExpMatrix = zeros(1,size(INPUT.Xmea,2));
else
    if isfield(EstimOpt,'MeaExpMatrix') == 0 || length(EstimOpt.MeaExpMatrix) ~= size(INPUT.Xmea,2)
        EstimOpt.MeaExpMatrix = ones(1,size(INPUT.Xmea,2));
        cprintf(rgb('DarkOrange'),'WARNING: MeaExpMatrix not defined - assuming that every measurment equation is explained with additional covariates /n')
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
        cprintf(rgb('DarkOrange'),'WARNING: Some of the LV not associated with any measurement equations.\n')
    elseif any(any(EstimOpt.MeaMatrix == 1) == 0)
        error('Measurment equations erroneously defined (some measurement variables unused)')
    end
elseif size(INPUT.Xmea,2) == 0
    EstimOpt.MeaMatrix = ones(EstimOpt.NLatent,size(INPUT.Xmea,2));
else
    EstimOpt.MeaMatrix = ones(EstimOpt.NLatent,size(INPUT.Xmea,2));
    disp('Assuming every Latent Variable in every measurment equation')
end

if isfield(EstimOpt,'MeaSpecMatrix') == 0 || numel(EstimOpt.MeaSpecMatrix) == 0
    EstimOpt.MeaSpecMatrix = zeros(1,size(INPUT.Xmea,2));
    if size(INPUT.Xmea,2) > 0
        disp('Using OLS for every measurment equation')
    end
end
if numel(EstimOpt.MeaSpecMatrix) == 1 % if only one number provided - replicate it for all Measurement variables
    EstimOpt.MeaSpecMatrix = EstimOpt.MeaSpecMatrix*ones(1,size(INPUT.Xmea,2));
end
if size(EstimOpt.MeaSpecMatrix,2) ~= size(INPUT.Xmea,2)
    error('Measurment specification (model) erroneously defined')
end

EstimOpt.NVarStr = size(INPUT.Xstr,2);  % no. of variables in structural equation
EstimOpt.NVarMea = sum(sum(EstimOpt.MeaMatrix)); % no. of parameters for Measurments without couting cutoffs, constants etc
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

for i = 1:EstimOpt.NVarStr
    if sum(isnan(INPUT.Xstr(INPUT.MissingInd == 0,i))) > 0
        cprintf(rgb('DarkOrange'),'WARNING: Structural variable %d contains NaN values \n', i)
    end
    if sum(isinf(INPUT.Xstr(INPUT.MissingInd == 0,i))) > 0
        cprintf(rgb('DarkOrange'),'WARNING: Structural variable %d contains Inf values \n', i)
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

EstimOpt.Names = [];% Names of the models
EstimOpt.NVarcut = 0; % no of cutoffs for ordered probits + constants + variances for OLS
EstimOpt.NVarcut0 = 0; % no of cutoffs for HMNL0
EstimOpt.CutMatrix = zeros(1,size(INPUT.Xmea,2));
EstimOpt.NamesLV = {};
for i = 1:size(INPUT.Xmea,2)
    if EstimOpt.MeaSpecMatrix(i) == 2 %Ordered probit: cutoffs
        EstimOpt.NVarcut = EstimOpt.NVarcut + length(unique(INPUT.Xmea(:,i))) - 1 + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0);
        EstimOpt.CutMatrix(i) = length(unique(INPUT.Xmea(:,i))) - 1 + sum(EstimOpt.MeaMatrix(:,i))+ EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0);
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end
        if EstimOpt.MeaExpMatrix(i) ~=0
            EstimOpt.NamesLV = [EstimOpt.NamesLV;EstimOpt.NamesMeaExp];
        end
        for n = 1:(length(unique(INPUT.Xmea(:,i))) - 1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV; cellfun(@(x)[x num2str(n)],{'Cutoff '},'UniformOutput',0)];
        end
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + length(unique(INPUT.Xmea(:,i))) - 1;
        EstimOpt.Names = [EstimOpt.Names,'OP '];
    elseif EstimOpt.MeaSpecMatrix(i) == 0
        EstimOpt.NVarcut = EstimOpt.NVarcut + 2 + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0); %OLS: constant + sigma
        EstimOpt.CutMatrix(i) = 2 + sum(EstimOpt.MeaMatrix(:,i))+ EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0);
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + 2;
        EstimOpt.Names = [EstimOpt.Names,'OLS '];
        EstimOpt.NamesLV = [EstimOpt.NamesLV;{'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end
        if EstimOpt.MeaExpMatrix(i) ~= 0
            EstimOpt.NamesLV = [EstimOpt.NamesLV;EstimOpt.NamesMeaExp];
        end
        EstimOpt.NamesLV = [EstimOpt.NamesLV;{'Sigma'}];
    elseif EstimOpt.MeaSpecMatrix(i) == 1 % MNL
        EstimOpt.NVarcut = EstimOpt.NVarcut + (length(unique(INPUT.Xmea(:,i))) - 2)*sum(EstimOpt.MeaMatrix(:,i)) + (length(unique(INPUT.Xmea(:,i))) - 1)*(1 + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0)); % constants + additional coefficients
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + length(unique(INPUT.Xmea(:,i))) - 1;
        EstimOpt.Names = [EstimOpt.Names,'MNL '];
        EstimOpt.CutMatrix(i) = (1 + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0) + sum(EstimOpt.MeaMatrix(:,i)))*(length(unique(INPUT.Xmea(:,i))) - 1);
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for j = 1:(length(unique(INPUT.Xmea(:,i))) - 1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV;{'Cons.'}];
            for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
                EstimOpt.NamesLV = [EstimOpt.NamesLV;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
            end
            if EstimOpt.MeaExpMatrix(i) ~= 0
                EstimOpt.NamesLV = [EstimOpt.NamesLV;EstimOpt.NamesMeaExp];
            end
        end
    elseif EstimOpt.MeaSpecMatrix(i) == 3 % Poisson
        EstimOpt.NVarcut = EstimOpt.NVarcut + 1 + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0); %Poiss: only constant
        EstimOpt.CutMatrix(i) = 1 + sum(EstimOpt.MeaMatrix(:,i)) + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0);
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + 1;
        EstimOpt.Names = [EstimOpt.Names,'POISS '];
        EstimOpt.NamesLV = [EstimOpt.NamesLV;{'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end
        if EstimOpt.MeaExpMatrix(i) ~= 0
            EstimOpt.NamesLV = [EstimOpt.NamesLV;EstimOpt.NamesMeaExp];
        end
    elseif EstimOpt.MeaSpecMatrix(i) == 4 % Negative Binomial
        EstimOpt.NVarcut = EstimOpt.NVarcut + 2 + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0); %Poiss: only constant
        EstimOpt.CutMatrix(i) = 2 + sum(EstimOpt.MeaMatrix(:,i)) + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0);
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + 2;
        EstimOpt.Names = [EstimOpt.Names,'NB '];
        EstimOpt.NamesLV = [EstimOpt.NamesLV;{'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end
        if EstimOpt.MeaExpMatrix(i) ~= 0
            EstimOpt.NamesLV = [EstimOpt.NamesLV;EstimOpt.NamesMeaExp];
        end
        EstimOpt.NamesLV = [EstimOpt.NamesLV;{'Theta'}];
    elseif EstimOpt.MeaSpecMatrix(i) == 5 % ZIP
        EstimOpt.NVarcut = EstimOpt.NVarcut + 2 + sum(EstimOpt.MeaMatrix(:,i)) + 2*EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0); %Poiss: only constant
        EstimOpt.CutMatrix(i) = 2 + 2*sum(EstimOpt.MeaMatrix(:,i)) + 2*EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0);
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + 2;
        EstimOpt.Names = [EstimOpt.Names,'ZIP '];
        EstimOpt.NamesLV = [EstimOpt.NamesLV;{'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end
        if EstimOpt.MeaExpMatrix(i) ~= 0
            EstimOpt.NamesLV = [EstimOpt.NamesLV;EstimOpt.NamesMeaExp];
        end
        EstimOpt.NamesLV = [EstimOpt.NamesLV;{'Cons.'}];
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV;cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end
        if EstimOpt.MeaExpMatrix(i) ~= 0
            EstimOpt.NamesLV = [EstimOpt.NamesLV;EstimOpt.NamesMeaExp];
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
    EstimOpt.NamesC =  [EstimOpt.NamesC;cellfun(@(x)[x num2str(i)],{'LV '},'UniformOutput',0)];
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

if size(INPUT.Xmea,2) > 0
    disp(['Using following models for attitudes: ' char(EstimOpt.Names)]);
end

%% Rescructure data

INPUT.YYY = reshape(INPUT.Y,[EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP]);
idx = sum(reshape(INPUT.MissingInd,[EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP])) == EstimOpt.NAlt;
INPUT.YYY(idx(ones(EstimOpt.NAlt,1),:) == 1) = NaN; % replace YYY in missing choice-tasks with NaN
INPUT.YY = reshape(INPUT.YYY(:,:,ones(EstimOpt.NClass,1)),[EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP*EstimOpt.NClass,1]);

INPUT.Xa(INPUT.MissingInd == 1,:) = NaN;

INPUT.Xstr = double(INPUT.Xstr(1:EstimOpt.NAlt*EstimOpt.NCT:end,:));
% normalize explanatory variables for structural equations:
if ~isfield(EstimOpt,'StrNorm')
    EstimOpt.StrNorm = ones(EstimOpt.NVarStr,1);
elseif length(EstimOpt.StrNorm) ~= EstimOpt.NVarStr && length(EstimOpt.StrNorm) ~= 1
    cprintf(rgb('DarkOrange'),'WARNING:  Structural variables normalization options (EstimOpt.StrNorm) incorrect - variables not normalized \n')
elseif length(EstimOpt.StrNorm) == 1
    EstimOpt.StrNorm = EstimOpt.StrNorm(ones(EstimOpt.NVarStr,1),1);
end
if any(EstimOpt.StrNorm > 0)
    INPUT.Xstr(:,EstimOpt.StrNorm == 1) = (INPUT.Xstr(:,EstimOpt.StrNorm == 1) - mean(INPUT.Xstr))./std(INPUT.Xstr);
end

INPUT.Xmea = INPUT.Xmea(1:EstimOpt.NAlt*EstimOpt.NCT:end,:);
EstimOpt.MissingIndMea = EstimOpt.MissingIndMea(1:EstimOpt.NAlt*EstimOpt.NCT:end,:);

INPUT.XXc = INPUT.Xc(1:EstimOpt.NCT*EstimOpt.NAlt:end,:); % NP x NVarC
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
        INPUT.Xmea_exp(:,EstimOpt.MeaExpNorm == 1) = (INPUT.Xmea_exp(:,EstimOpt.MeaExpNorm == 1) - mean(INPUT.Xmea_exp))./std(INPUT.Xmea_exp);
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
    LLfun0 = @(B) LL_hmnl0(INPUT.Xmea,EstimOpt,B);
    [Results.MIMIC0.bhat,LL0] = fminunc(LLfun0,0.001*ones(EstimOpt.NVarcut0,1),OptimOpt_0);
    Results.MIMIC0.LL = -LL0;
else
    Results.MIMIC0.bhat = [];
    Results.MIMIC0.LL = 0;
end

%% Starting values

EstimOpt.jitter1 = 0.8; % Jittering parameter (relative) for HMNL starting values (attributes)
EstimOpt.jitter2 = 0.1; % Jittering parameter (absolute) for class probabilities starting values

if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == (EstimOpt.NVarA*EstimOpt.NClass + (EstimOpt.NVarC + EstimOpt.NLatent)*(EstimOpt.NClass - 1) + EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut)
    b0 = B_backup(:);
    disp('Using the starting values from Backup')
elseif isfield(Results_old,'HLC') && isfield(Results_old.HLC,'b0') % starting values provided
    Results_old.HLC.b0_old = Results_old.HLC.b0(:);
    Results_old.HLC = rmfield(Results_old.HLC,'b0');
    if length(Results_old.HLC.b0_old) ~= (EstimOpt.NVarA*EstimOpt.NClass + (EstimOpt.NVarC + EstimOpt.NLatent)*(EstimOpt.NClass - 1) + EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut)
        cprintf(rgb('DarkOrange'),'WARNING: Incorrect no. of starting values or model specification \n')
        Results_old.HLC = rmfield(Results_old.HLC,'b0_old');
    else
        b0 = Results_old.HLC.b0_old(:);
    end
end

if ~exist('b0','var')
    if isfield(Results_old,'MIMIC') && isfield(Results_old.MIMIC,'bhat')
        disp('Using MIMIC results as starting values')
        Results_old.MIMIC.bhat = Results_old.MIMIC.bhat(:);
        bstr = reshape(Results_old.MIMIC.bhat(1:EstimOpt.NVarStr*EstimOpt.NLatent),[EstimOpt.NVarStr,EstimOpt.NLatent]);
        LV = INPUT.Xstr*bstr;
        LV = permute(LV(:,:,ones(EstimOpt.NCT*EstimOpt.NAlt,1)),[3 1 2]);
        LV = reshape(LV,[EstimOpt.NCT*EstimOpt.NAlt*EstimOpt.NP,EstimOpt.NLatent]);
        disp('Estimating additional MNL for starting values')
        EstimOpt0 = EstimOpt;
        EstimOpt0.DISPLAY = 0;
        EstimOpt0.BActive = ones(1,EstimOpt.NVarA);
        EstimOpt0.NClass = EstimOpt.NClass;
        EstimOpt0.WTP_space = EstimOpt.WTP_space;
        EstimOpt0.WTP_matrix = EstimOpt.WTP_matrix;
        EstimOpt0.NVarA = EstimOpt.NVarA;
        Results.MNL = MNL(INPUT,[],EstimOpt0,OptimOpt);
        disp('Estimating additional LC for starting values')
        if size(INPUT.Xc,2) > 1
            INPUT.Xc = [INPUT.Xc(:,2:end),LV];
        else
            INPUT.Xc = LV;
        end
        if isfield(EstimOpt,'BActive') == 1 && numel(EstimOpt.BActive) ~= 0
            EstimOpt0.BActive = EstimOpt.BActive(1:EstimOpt.NVarA*EstimOpt.NClass+(EstimOpt.NClass-1)*(EstimOpt.NVarC+EstimOpt.NLatent));
        else
            EstimOpt0.BActive = ones(1,EstimOpt.NVarA*EstimOpt.NClass+(EstimOpt.NClass-1)*(EstimOpt.NVarC+EstimOpt.NLatent));
        end
        Results.LC_LV = LC(INPUT,Results_old,EstimOpt0,OptimOpt);
        b0 = [Results.LC_LV.bhat;Results_old.MIMIC.bhat];
    elseif isfield(Results_old,'HMNL') && isfield(Results_old.HMNL,'bhat')
        disp('Using HMNL results as starting values')
        Results_old.HMNL.bhat = Results_old.HMNL.bhat(:);
        if EstimOpt.WTP_space > 0
            Results_old.HMNL.bhat(1:(EstimOpt.NVarA-EstimOpt.WTP_space)) = Results_old.HMNL.bhat(1:(EstimOpt.NVarA-EstimOpt.WTP_space))./Results_old.HMNL.bhat(EstimOpt.WTP_matrix);
        end
        b0 = [Results_old.HMNL.bhat((1:EstimOpt.NVarA)'*ones(1,EstimOpt.NClass),:).* ...
             (EstimOpt.jitter1.*unifrnd(0,ones(EstimOpt.NVarA.*EstimOpt.NClass,1))); ...
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
    err_mtx = randn(EstimOpt.NP*EstimOpt.NRep,EstimOpt.NLatent); %to be cut down later
elseif EstimOpt.Draws == 2 % LHS
    cprintf('*blue','Latin Hypercube Sampling '); cprintf('draws \n');
    err_mtx = lhsnorm(zeros((EstimOpt.NLatent)*EstimOpt.NP,1),diag(ones((EstimOpt.NLatent)*EstimOpt.NP,1)),EstimOpt.NRep);
    err_mtx = reshape(err_mtx,[EstimOpt.NRep*EstimOpt.NP,EstimOpt.NLatent]);
elseif EstimOpt.Draws >= 3 % Quasi random draws
    if EstimOpt.Draws == 3
        cprintf('*blue','Halton '); cprintf(['draws (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n'])
        hm1 = haltonset(EstimOpt.NLatent,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); %
    elseif EstimOpt.Draws == 4 % apply reverse-radix scrambling
        cprintf('*blue','Halton '); cprintf(['draws with reverse radix scrambling (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n'])
        hm1 = haltonset(EstimOpt.NLatent,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); %
        hm1 = scramble(hm1,'RR2');
    elseif EstimOpt.Draws == 5
        cprintf('*blue','Sobol '); cprintf(['draws (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n'])
        hm1 = sobolset(EstimOpt.NLatent,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
    elseif EstimOpt.Draws == 6
        cprintf('*blue','Sobol '); cprintf(['draws with random linear scramble and random digital shift (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n'])
        hm1 = sobolset(EstimOpt.NLatent,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
        hm1 = scramble(hm1,'MatousekAffineOwen');
    end

    err_mtx = net(hm1,EstimOpt.NP*EstimOpt.NRep); % this takes every point:
    clear hm1;
    if EstimOpt.NP*EstimOpt.NRep < 3e+7
        err_mtx = icdf('Normal',err_mtx,0,1); %to be cut down later
    else % this is for very large number of draws * variables
        for i = 1:EstimOpt.NLatent
            err_mtx(:,i) = icdf('Normal',err_mtx(:,i),0,1); %to be cut down later
        end
    end

end

err_sliced = err_mtx'; % NLatent x NRep * NP


%% Display Options


if EstimOpt.NumGrad == 0
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient to numerical -  analytical gradient not supported in the HLC model \n')
end

if ((isfield(EstimOpt,'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    OptimOpt.GradObj = 'on';
end

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

LLfun = @(B) LL_hlc_MATlike(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,INPUT.W,EstimOpt,OptimOpt,B);
if EstimOpt.ConstVarActive == 0
    if EstimOpt.HessEstFix == 0
        [Results.bhat,LL,Results.exitf,Results.output,Results.g,Results.hess] = fminunc(LLfun,b0,OptimOpt);
    else
        [Results.bhat,LL,Results.exitf,Results.output,Results.g] = fminunc(LLfun,b0,OptimOpt);
    end
elseif EstimOpt.ConstVarActive == 1 % equality constraints
    EstimOpt.CONS1 = diag(1 - EstimOpt.BActive);
    EstimOpt.CONS1(sum(EstimOpt.CONS1,1) == 0,:) = [];
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
R2 = R2_hybrid(INPUT.YY,INPUT.Xa,INPUT.Xstr,INPUT.XXc,[],[],INPUT.MissingInd,err_sliced,EstimOpt,Results.bhat,1);

Results.b0_old = b0;

if EstimOpt.HessEstFix == 1
    f = LL_hlc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,EstimOpt,Results.bhat);
    Results.jacobian = numdiff(@(B) INPUT.W.*LL_hlc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,EstimOpt,B),INPUT.W.*f,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) INPUT.W.*LL_hlc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.jacobian = hessian(@(B) sum(INPUT.W.*LL_hlc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,EstimOpt,B),1),Results.bhat);
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
        [~,Results.jacobian] = LL_hlc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,EstimOpt,Results.bhat);
        Results.jacobian = Results.jacobian.*INPUT.W;
    else
        Results.LLdetailed = LL_hlc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,EstimOpt,Results.bhat);
        Results.jacobian = numdiff(@(B) INPUT.W.*LL_hlc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,EstimOpt,B),INPUT.W.*Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
    end
    RobustHess = Results.jacobian'*Results.jacobian;
    Results.ihess = Results.ihess*RobustHess*Results.ihess;
end

Results.std = sqrt(diag(Results.ihess));
Results.std(EstimOpt.BActive == 0) = NaN;
Results.std(EstimOpt.BLimit == 1) = 0;
Results.std(imag(Results.std) ~= 0) = NaN;
Results.DetailsMea = [];

for i = 1:EstimOpt.NClass
    Results.DetailsA(1:EstimOpt.NVarA,4*i-3) = Results.bhat((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA);
    Results.DetailsA(1:EstimOpt.NVarA,4*i-1:4*i) = [Results.std((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA),pv(Results.bhat((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA),Results.std((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA))];
    l = EstimOpt.NClass*EstimOpt.NVarA;
    if i ~= EstimOpt.NClass
        for j = 1:EstimOpt.NVarC+EstimOpt.NLatent
            Results.DetailsV(j,4*i-3) = Results.bhat(l+j+(EstimOpt.NVarC+EstimOpt.NLatent)*(i-1));
            Results.DetailsV(j,4*i-1:4*i) = [Results.std(l+j+(EstimOpt.NVarC+EstimOpt.NLatent)*(i-1)),pv(Results.bhat(l+j+(EstimOpt.NVarC+EstimOpt.NLatent)*(i-1)),Results.std(l+j+(EstimOpt.NVarC+EstimOpt.NLatent)*(i-1)))];
        end
    end
    l = l+(EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1);
    if EstimOpt.NVarStr > 0
        Results.DetailsS(:,4*i-3) = Results.bhat(l+(i-1)*EstimOpt.NVarStr+1:l+i*EstimOpt.NVarStr);
        Results.DetailsS(:,4*i-1:4*i)  = [Results.std(l+(i-1)*EstimOpt.NVarStr+1:l+i*EstimOpt.NVarStr),pv(Results.bhat(l+(i-1)*EstimOpt.NVarStr+1:l+i*EstimOpt.NVarStr),Results.std(l+(i-1)*EstimOpt.NVarStr+1:l+i*EstimOpt.NVarStr))];
        l = l + EstimOpt.NClass*EstimOpt.NVarStr;
    end
%     Results.DetailsM((i-1)*EstimOpt.NVarMea+1:i*EstimOpt.NVarMea,1) = Results.bhat(l+(i-1)*EstimOpt.NVarMea+1:l+i*EstimOpt.NVarMea);
%     Results.DetailsM((i-1)*EstimOpt.NVarMea+1:i*EstimOpt.NVarMea,3:4) = [Results.std(l+(i-1)*EstimOpt.NVarMea+1:l+i*EstimOpt.NVarMea),pv(Results.bhat(l+(i-1)*EstimOpt.NVarMea+1:l+i*EstimOpt.NVarMea),Results.std(l+(i-1)*EstimOpt.NVarMea+1:l+i*EstimOpt.NVarMea))];
end

if sum(EstimOpt.CutMatrix) > 0
    Results.DetailsM(1:sum(EstimOpt.CutMatrix),1) = Results.bhat(l+1:end);
    Results.DetailsM(1:sum(EstimOpt.CutMatrix),3:4) = [Results.std(l+1:end),pv(Results.bhat(l+1:end), Results.std(l+1:end))];
end

Results.DetailsV = [Results.DetailsV,zeros(EstimOpt.NVarC+EstimOpt.NLatent,1),NaN(EstimOpt.NVarC+EstimOpt.NLatent,3)];
Results.LL0 = Results.MIMIC0.LL + Results_old.MNL0.LL;
EstimOpt.params = length(Results.bhat);
if isfield(EstimOpt,'BActive')
    EstimOpt.params = EstimOpt.params - sum(EstimOpt.BActive == 0);
end
Results.stats = [Results.LL;Results_old.MNL0.LL;1 - Results.LL/Results_old.MNL0.LL;R2;((2*EstimOpt.params - 2*Results.LL))/EstimOpt.NObs;((log(EstimOpt.NObs)*EstimOpt.params - 2*Results.LL))/EstimOpt.NObs;EstimOpt.NObs;EstimOpt.NP;EstimOpt.params];

Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;
Results.INPUT = INPUT;

NSdSim = 10000;
bclass = reshape([Results.bhat(EstimOpt.NVarA*EstimOpt.NClass+1:EstimOpt.NVarA*EstimOpt.NClass+(EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1));zeros(EstimOpt.NVarC+EstimOpt.NLatent,1)],[EstimOpt.NVarC+EstimOpt.NLatent,EstimOpt.NClass]);
bstr = reshape(Results.bhat(EstimOpt.NVarA*EstimOpt.NClass+(EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1)+1:EstimOpt.NVarA*EstimOpt.NClass+(EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1)+EstimOpt.NLatent*EstimOpt.NVarStr),[EstimOpt.NVarStr,EstimOpt.NLatent]);
bclass_sim = reshape([mvnrnd(Results.bhat(EstimOpt.NVarA*EstimOpt.NClass+1:EstimOpt.NVarA*EstimOpt.NClass+(EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1)),Results.ihess(EstimOpt.NVarA*EstimOpt.NClass+1:EstimOpt.NVarA*EstimOpt.NClass+(EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1),EstimOpt.NVarA*EstimOpt.NClass+1:EstimOpt.NVarA*EstimOpt.NClass+(EstimOpt.NVarC+EstimOpt.NLatent)*(EstimOpt.NClass-1)),NSdSim)';zeros(EstimOpt.NVarC+EstimOpt.NLatent,NSdSim)],[EstimOpt.NVarC+EstimOpt.NLatent,EstimOpt.NClass,NSdSim]);

LV = INPUT.Xstr*bstr;
V = exp([INPUT.XXc,LV]*bclass);% NP x NClass
Vsum = sum(V,2);
Results.PClass = zeros(1,4*EstimOpt.NClass);
Results.PClass(1,1:4:EstimOpt.NClass*4-3) = mean(V./Vsum,1);

PClass_mean = zeros(NSdSim,EstimOpt.NClass);

parfor i = 1:NSdSim
    bhat_sim_i = bclass_sim(:,:,i);
    V_i = exp([INPUT.XXc,LV]*bhat_sim_i);
    Vsum_i = sum(V_i,2);
    PC_i = V_i./Vsum_i;
    PClass_mean(i,:) = mean(PC_i,1);
end
Results.PClass(1,3:4:EstimOpt.NClass*4-1) = std(PClass_mean);
Results.PClass(1,4:4:EstimOpt.NClass*4) = pv(Results.PClass(1,1:4:EstimOpt.NClass*4-3),Results.PClass(1,3:4:EstimOpt.NClass*4-1));
Results.PClass(1,1:2:end) = Results.PClass(1,1:2:end)*100;
Results.PClass95ci = [quantile(PClass_mean,0.025);quantile(PClass_mean,0.975)];

%% Header

Head = cell(1,2);
Head(1,1) = {'HLC'};
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
    Heads.DetailsA{i,1} = num2str(i,'Class %1.0f');
end
Heads.DetailsA(end+1) = {'tb'};

if EstimOpt.NVarStr > 0
    Template1 = [Template1;{'DetailsS'}];
    Template2 = [Template2;{'DetailsS'}];
    Names.DetailsS = EstimOpt.NamesStr;
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
Heads.DetailsV(1:2,1) = {'Latent class probability model';'lb'};
ST = [ST,{'DetailsV'}];

Template1 = [Template1;{'PClass'}];
Template2 = [Template2;{'PClass'}];
Names.PClass = {"(%)"};
Heads.PClass(:,2) = Heads.DetailsA;
Heads.PClass(1:2,1) = {'Average class probabilities';'lb'};
ST = [ST,{'PClass'}];

if any(EstimOpt.MeaSpecMatrix == 2)
    Results.DetailsOP = [];
end

k = 0;

if any(EstimOpt.MeaSpecMatrix == 2)
    Results.DetailsOP = [];
end

for i = 1:size(INPUT.Xmea,2)
    tmp = EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~= 0);
    model = model_name(EstimOpt.MeaSpecMatrix(i));
    Heads.(strcat('Xmea',num2str(i)))(1:2,1) = [{['Measurment equation for:',' ',char(EstimOpt.NamesMea(i)),' (',model,')']};'lb'];

    if EstimOpt.MeaSpecMatrix(i) == 2
        l = l + sum(EstimOpt.MeaMatrix(:,i)) + tmp;
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
        Names.(strcat('Xmea',num2str(i))) = EstimOpt.NamesLV(k+1:k+floor(EstimOpt.CutMatrix(i)/2));
        Names.(strcat('Xmea2',num2str(i))) = EstimOpt.NamesLV(k+floor(EstimOpt.CutMatrix(i)/2)+1:k+EstimOpt.CutMatrix(i));
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
        Names.(strcat('Xmea',num2str(i))) = EstimOpt.NamesLV(k+1:k+EstimOpt.CutMatrix(i));
        Temp1 = cell(1,size(Template1,2));
        Temp1(1,1) = {strcat('Xmea',num2str(i))};
        Template1 = [Template1;Temp1]; %#ok<AGROW>
        Temp2 = cell(1,size(Template2,2));
        Temp2(1,1) = {strcat('Xmea',num2str(i))};
        Template2 = [Template2;Temp2]; %#ok<AGROW>
        ST = [ST, {strcat('Xmea',num2str(i))}]; %#ok<AGROW>
    end
    k = k + EstimOpt.CutMatrix(i);
end
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
        Tail(14,2) = {[num2str(EstimOpt.NRep),' ','Halton draws with reverse radix scrambling (skip = ',num2str(EstimOpt.HaltonSkip), '; leap = ', num2str(EstimOpt.HaltonLeap),')']};
    case 5
        Tail(14,2) = {[num2str(EstimOpt.NRep),' ','Sobol draws (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),')']};
    case 6
        Tail(14,2) = {[num2str(EstimOpt.NRep),' ','Sobol draws with random linear scramble and random digital shift (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),')']};
end

Tail(15,2) = {OptimOpt.Algorithm};

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
