function Results = HMXL(INPUT,Results_old,EstimOpt,OptimOpt)


global B_backup;

% save tmp_HMXL
% return

tic

Results.bhat = [];
Results.R = [];
Results.R_out = {};
Results.stats = [];


%% Check data and inputs


if nargin < 3 % check no. of inputs
    error('Too few input arguments for HMXL(INPUT,EstimOpt,OptimOpt)')
end

disp(' ');
disp('__________________________________________________________________________________________________________________');
disp(' ');

warning off MATLAB:mir_warning_maybe_uninitialized_temporary

format shortG;
format compact;

if any(INPUT.W ~= 1)
    cprintf('Black','Estimating '); cprintf('*Black','weighted '); cprintf('Black','HMXL model...\n');
else
    disp('Estimating HMXL model ...')
end

if isfield(EstimOpt,'FullCov') == 0;
    EstimOpt.FullCov = 0;
end
if isfield(EstimOpt, 'WTP_space') == 0
    EstimOpt.WTP_space = 0;
    EstimOpt.WTP_matrix = [];
elseif EstimOpt.WTP_space == 0;
    EstimOpt.WTP_matrix = [];
end

if EstimOpt.FullCov == 0
    disp('with non-correlated random parameters ...')
    if EstimOpt.WTP_space > 0
        disp('in WTP-space.')
    else
        disp('in preference-space.')
    end
elseif EstimOpt.FullCov == 1
    disp('with correlated random parameters ...')
    if EstimOpt.WTP_space > 0
        disp('in WTP-space.')
    else
        disp('in preference-space.')
    end
elseif EstimOpt.FullCov == 2
    disp('with correlated random parameters and latent variables ...')
    if EstimOpt.WTP_space > 0
        disp('in WTP-space.')
    else
        disp('in preference-space.')
    end
end

if isfield(EstimOpt,'Dist') == 0 || isempty(EstimOpt.Dist)
    EstimOpt.Dist = zeros(1,EstimOpt.NVarA);
    if EstimOpt.WTP_space == 0
        cprintf(rgb('DarkOrange'), 'WARNING: distributions for random parameters not specified - assuming normality \n')
    else
        cprintf(rgb('DarkOrange'), 'WARNING: distributions for random parameters not specified - assuming normality (monetary parameters assumed log-normal) \n')
        EstimOpt.Dist(end-EstimOpt.WTP_space+1:end) = 1; % cost in WTP-space models log-normally distributed
    end
else
    if length(EstimOpt.Dist) == 1
        EstimOpt.Dist = EstimOpt.Dist.*ones(1,EstimOpt.NVarA);
    elseif length(EstimOpt.Dist) == EstimOpt.NVarA
        EstimOpt.Dist = EstimOpt.Dist(:)';
    else
        error('Incorrect no. of random parameters'' distributions provided')
    end
end
disp(['Random parameters distributions: ', num2str(EstimOpt.Dist),' (-1 - constant, 0 - normal, 1 - lognormal, 2 - spike, 5 - Weibull)'])

if EstimOpt.WTP_space > 0 && sum(EstimOpt.Dist(end-EstimOpt.WTP_space+1:end)==1) > 0 && any(mean(INPUT.Xa(:,end-EstimOpt.WTP_space+1:end)) >= 0)
    cprintf(rgb('DarkOrange'), 'WARNING: Cost attributes with log-normally distributed parameters should enter utility function with a ''-'' sign \n')
end

if isfield(EstimOpt,'NLatent') == 0
    EstimOpt.NLatent = 1;
    disp(['Assuming ', num2str(EstimOpt.NLatent)',' Latent Variable(s)']);
else
    disp(num2str(EstimOpt.NLatent,'The model has %1.0f Latent Variable(s)'))
end

if isfield(INPUT, 'Xm') == 0
    INPUT.Xm = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarM = size(INPUT.Xm,2); % Number of covariates of means of random parameters

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

if isfield(INPUT, 'Xstr') == 0 || numel(INPUT.Xstr) == 0
    % 	error('Define variables to structural equations')
    cprintf(rgb('DarkOrange'), 'WARNING: Structural equations empty.\n')
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
    if isfield(EstimOpt,'MeaExpMatrix') == 0
        cprintf(rgb('DarkOrange'), 'WARNING: MeaExpMatrix not defined - assuming that every measurment equation is explained with additional covariates \n')
        EstimOpt.MeaExpMatrix = ones(1, size(INPUT.Xmea,2));
    elseif length(EstimOpt.MeaExpMatrix) ~= size(INPUT.Xmea,2)
        error('Additional covariates of measurment equations erroneously defined (incorrect size of the MeaExpMatrix)')
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
        %         error('Measurment equations erroneously defined (some Latent Variables without measurement equations in the MeaMatrix)')
        cprintf(rgb('DarkOrange'), 'WARNING: Some of the LV not associated with any measurement equations.\n')
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

EstimOpt.NVarStr = size(INPUT.Xstr,2);  % no. of variables in structural equation
EstimOpt.NVarMea = sum(sum(EstimOpt.MeaMatrix)); % no of parameters for Measurments without couting cutoffs, constants etc
EstimOpt.NVarMeaExp = size(INPUT.Xmea_exp,2);

for i=1:size(EstimOpt.MeaMatrix,2)
    if sum(isnan(INPUT.Xmea(INPUT.MissingInd==0,i))) > 0
        cprintf(rgb('DarkOrange'), 'WARNING:  Measurement variable %d contains NaN values \n' ,i)
    end
    if sum(isnan(INPUT.Xmea(INPUT.MissingInd==0,i))) > 0
        cprintf(rgb('DarkOrange'), 'WARNING:  Measurement variable %d contains Inf values \n',i)
    end
    if numel(EstimOpt.MeaSpecMatrix(i) > 0) > 0
        if EstimOpt.MeaSpecMatrix(i) > 0 && numel(unique(INPUT.Xmea(INPUT.MissingInd==0,i))) > 10
            cprintf(rgb('DarkOrange'), 'WARNING: There are over 10 levels for measurement variable %d \n', i)
        end
    end
    if numel(EstimOpt.MeaSpecMatrix(i) > 0) > 0
        if EstimOpt.MeaSpecMatrix(i) == 1 % MNL
            INPUT.Xmea(INPUT.MissingInd==0,i) = INPUT.Xmea(INPUT.MissingInd==0,i) - min(unique(INPUT.Xmea(INPUT.MissingInd==0,i))) + 1; % make unique values positive (avoid problems with dummyvar)
            %             cprintf(rgb('DarkOrange'), 'WARNING: There are over 10 levels for measurement variable %d \n', i)
        end
    end
end

for i = 1:EstimOpt.NVarStr
    if sum(isnan(INPUT.Xstr(INPUT.MissingInd==0,i))) > 0
        cprintf(rgb('DarkOrange'), 'WARNING:  Structural variable %d contains NaN values \n',i)
    end
    if sum(isinf(INPUT.Xstr(INPUT.MissingInd==0,i))) > 0
        cprintf(rgb('DarkOrange'), 'WARNING:  Structural variable %d contains Inf values \n',i)
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
EstimOpt.CutMatrix = zeros(1, size(INPUT.Xmea,2));
EstimOpt.NamesLV = {};
for i = 1:size(INPUT.Xmea,2)
    if EstimOpt.MeaSpecMatrix(i) == 2 %Ordered probit: cutoffs
        EstimOpt.NVarcut = EstimOpt.NVarcut + length(unique(INPUT.Xmea(INPUT.MissingInd==0,i))) - 1 + EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0);
        EstimOpt.CutMatrix(i) = length(unique(INPUT.Xmea(INPUT.MissingInd==0,i))) - 1 + sum(EstimOpt.MeaMatrix(:,i))+ EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0);
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV; cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
        end
        if EstimOpt.MeaExpMatrix(i) ~=0
            EstimOpt.NamesLV = [EstimOpt.NamesLV; EstimOpt.NamesMeaExp];
        end
        for n = 1:(length(unique(INPUT.Xmea(INPUT.MissingInd==0,i))) - 1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV; cellfun(@(x)[x num2str(n)],{'Cutoff '},'UniformOutput',0)];
        end
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + length(unique(INPUT.Xmea(INPUT.MissingInd==0,i))) - 1;
        EstimOpt.Names = [EstimOpt.Names, 'OP '];
    elseif EstimOpt.MeaSpecMatrix(i) == 0
        EstimOpt.NVarcut = EstimOpt.NVarcut + 2+EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0); %OLS: constant + sigma
        EstimOpt.CutMatrix(i) = 2+sum(EstimOpt.MeaMatrix(:,i))+ EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0);
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
        EstimOpt.NVarcut = EstimOpt.NVarcut + (length(unique(INPUT.Xmea(INPUT.MissingInd==0,i))) - 2)*sum(EstimOpt.MeaMatrix(:,i)) + (length(unique(INPUT.Xmea(INPUT.MissingInd==0,i))) - 1)*(1+ EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0)); % constants + additional coefficients
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + length(unique(INPUT.Xmea(INPUT.MissingInd==0,i)))-1;
        EstimOpt.Names = [EstimOpt.Names, 'MNL '];
        EstimOpt.CutMatrix(i) = (1+ EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0)+sum(EstimOpt.MeaMatrix(:,i)))*(length(unique(INPUT.Xmea(INPUT.MissingInd==0,i)))-1);
        k = find(EstimOpt.MeaMatrix(:,i) == 1);
        for j = 1:(length(unique(INPUT.Xmea(INPUT.MissingInd==0,i)))-1)
            EstimOpt.NamesLV = [EstimOpt.NamesLV; {'Cons.'}];
            for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
                EstimOpt.NamesLV = [EstimOpt.NamesLV; cellfun(@(x)[x num2str(k(n))],{'LV '},'UniformOutput',0)];
            end
            if EstimOpt.MeaExpMatrix(i) ~=0
                EstimOpt.NamesLV = [EstimOpt.NamesLV; EstimOpt.NamesMeaExp];
            end
        end
    elseif EstimOpt.MeaSpecMatrix(i) == 3 % Poisson
        EstimOpt.NVarcut = EstimOpt.NVarcut +1+EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0); %Poiss: only constant
        EstimOpt.CutMatrix(i) = 1+sum(EstimOpt.MeaMatrix(:,i))+ EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0);
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
        EstimOpt.NVarcut = EstimOpt.NVarcut +2+EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0); %Poiss: only constant
        EstimOpt.CutMatrix(i) = 2+sum(EstimOpt.MeaMatrix(:,i))+ EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0);
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
        EstimOpt.NVarcut = EstimOpt.NVarcut +2 +sum(EstimOpt.MeaMatrix(:,i)) +2*EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0); %Poiss: only constant
        EstimOpt.CutMatrix(i) = 2+2*sum(EstimOpt.MeaMatrix(:,i))+ 2*EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0);
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
    elseif EstimOpt.MeaSpecMatrix(i) == 6 % ZINB
        EstimOpt.NVarcut = EstimOpt.NVarcut +3 +sum(EstimOpt.MeaMatrix(:,i)) +2*EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0); %Poiss: only constant
        EstimOpt.CutMatrix(i) = 3+2*sum(EstimOpt.MeaMatrix(:,i))+ 2*EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0);
        EstimOpt.NVarcut0 = EstimOpt.NVarcut0 + 3;
        EstimOpt.Names = [EstimOpt.Names, 'ZINB '];
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
        EstimOpt.NamesLV = [EstimOpt.NamesLV; {'Theta'}];
    end
end

if isfield(EstimOpt,'NamesA') == 0 || isempty(EstimOpt.NamesA) || length(EstimOpt.NamesA) ~= EstimOpt.NVarA
    EstimOpt.NamesA = (1:EstimOpt.NVarA)';
    EstimOpt.NamesA = cellstr(num2str(EstimOpt.NamesA));
elseif size(EstimOpt.NamesA,1) ~= EstimOpt.NVarA
    EstimOpt.NamesA = EstimOpt.NamesA';
end
if EstimOpt.NVarM > 0
    if isfield(EstimOpt,'NamesM') == 0 || isempty(EstimOpt.NamesM)|| length(EstimOpt.NamesM) ~= EstimOpt.NVarM
        EstimOpt.NamesM = (1:EstimOpt.NVarM)';
        EstimOpt.NamesM = cellstr(num2str(EstimOpt.NamesM));
    elseif size(EstimOpt.NamesM,1) ~= EstimOpt.NVarM
        EstimOpt.NamesM = EstimOpt.NamesM';
    end
else
    EstimOpt.NamesM = {};
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
    disp(['Using following models for attitudes: '  char(EstimOpt.Names)]);
end


%% Rescructure data


%INPUT.YY = reshape(INPUT.Y,EstimOpt.NAlt,EstimOpt.NCT,1,EstimOpt.NP);
%INPUT.YY = INPUT.YY(:,:,ones(1,1,EstimOpt.NRep,1),:); % NAlt x NCT x NRep x NP

INPUT.YYY = reshape(INPUT.Y,EstimOpt.NAlt,EstimOpt.NCT,EstimOpt.NP);
idx = sum(reshape(INPUT.MissingInd,EstimOpt.NAlt,EstimOpt.NCT,EstimOpt.NP)) == EstimOpt.NAlt; ...
    INPUT.YYY(idx(ones(EstimOpt.NAlt,1),:,:) == 1) = NaN; % replace YYY in missing choice-tasks with NaN
INPUT.YY = reshape(INPUT.YYY,EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP);


INPUT.Xa(INPUT.MissingInd == 1,:) = NaN;
INPUT.XXa = reshape(INPUT.Xa',EstimOpt.NVarA, EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP);
INPUT.XXa = permute(INPUT.XXa, [2 1 3]); % NAlt*NCT x NVarA x NP

INPUT.Xstr = double(INPUT.Xstr(1:EstimOpt.NAlt*EstimOpt.NCT:end,:));
if isfield(EstimOpt, 'StrNorm')
    if length(EstimOpt.StrNorm) ~= EstimOpt.NVarStr && length(EstimOpt.StrNorm) ~= 1
        cprintf(rgb('DarkOrange'), 'WARNING:  Structural variables normalization options (EstimOpt.StrNorm) incorrect - variables not normalized \n')
    else
        if length(EstimOpt.StrNorm) == 1
            EstimOpt.StrNorm = EstimOpt.StrNorm(ones(EstimOpt.NVarStr,1),1);
        end
        meanx = mean(INPUT.Xstr);
        stdx = std(INPUT.Xstr);
        INPUT.Xstr(:,EstimOpt.StrNorm == 1) = (INPUT.Xstr(:,EstimOpt.StrNorm == 1) - meanx(ones(size(INPUT.Xstr,1),1),EstimOpt.StrNorm == 1))./stdx(ones(size(INPUT.Xstr,1),1),EstimOpt.StrNorm == 1);
    end
end

INPUT.Xmea = INPUT.Xmea(1:EstimOpt.NAlt*EstimOpt.NCT:end,:);
INPUT.Xm = INPUT.Xm(1:EstimOpt.NAlt*EstimOpt.NCT:end,:)'; % NVarM x NP
if EstimOpt.NVarMeaExp > 0
    INPUT.Xmea_exp = INPUT.Xmea_exp(1:EstimOpt.NAlt*EstimOpt.NCT:end,:);
end
EstimOpt.indx1 = [];
EstimOpt.indx2 = [];
EstimOpt.indx3 = [];
if EstimOpt.NumGrad == 0 && EstimOpt.FullCov > 0
    for i = 1:EstimOpt.NVarA
        EstimOpt.indx1 = [EstimOpt.indx1, i:EstimOpt.NVarA];
        EstimOpt.indx2 = [EstimOpt.indx2, i*ones(1,EstimOpt.NVarA+1-i)];
    end
    if EstimOpt.FullCov == 2
        aa = 1;
        ab = sum(1:EstimOpt.NVarA);
        for i = 1:EstimOpt.NVarA
            EstimOpt.indx3 = [EstimOpt.indx3, aa:(aa+EstimOpt.NVarA-i), ab+i];
            aa = aa+EstimOpt.NVarA-i+1;
        end
    end
end

%% Estimating MIMIC0


if size(INPUT.Xmea,2) > 0
    if isfield(Results_old,'HMNL') && isfield(Results_old.HMNL,'MIMIC0') && isfield(Results_old.HMNL.MIMIC0,'LL') && ~isempty(Results_old.HMNL.MIMIC0.LL)
        Results.MIMIC0.LL = Results_old.HMNL.MIMIC0.LL;
    elseif isfield(Results_old,'HMXL_d') && isfield(Results_old.HMXL_d,'MIMIC0') && isfield(Results_old.HMXL_d.MIMIC0,'LL') && ~isempty(Results_old.HMXL_d.MIMIC0.LL)
        Results.MIMIC0.LL = Results_old.HMXL_d.MIMIC0.LL;
    else
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
    end
else
    Results.MIMIC0.bhat = [];
    Results.MIMIC0.LL = 0;
end


%% Starting values

if EstimOpt.FullCov == 0
    
    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == (EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM) + EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut)
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'HMXL_d') && isfield(Results_old.HMXL_d,'b0') % starting values provided
        Results_old.HMXL_d.b0_old = Results_old.HMXL_d.b0(:);
        Results_old.HMXL_d = rmfield(Results_old.HMXL_d,'b0');
        if length(Results_old.HMXL_d.b0_old) ~= (EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM) + EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut)
            cprintf(rgb('DarkOrange'), 'WARNING:  Incorrect no. of starting values or model specification \n')
            Results_old.HMXL_d = rmfield(Results_old.HMXL_d,'b0_old');
        else
            b0 = Results_old.HMXL_d.b0_old(:);
        end
    end
    if  ~exist('b0','var')
        if isfield(Results_old,'HMNL') && isfield(Results_old.HMNL,'bhat')
            if isfield(Results_old,'MXL_d') && isfield(Results_old.MXL_d,'bhat')
                disp('Using HMNL and MXL_d results as starting values')
                Results_old.HMNL.bhat = Results_old.HMNL.bhat(:);
                Results_old.MXL_d.bhat = Results_old.MXL_d.bhat(:);
                % b0 = [Results_old.MXL_d.bhat(1:(2+EstimOpt.NVarM)*EstimOpt.NVarA); Results_old.HMNL.bhat(EstimOpt.NVarA+EstimOpt.NVarM+1:end)];
                b0 = [Results_old.MXL_d.bhat(1:2*EstimOpt.NVarA); Results_old.HMNL.bhat(EstimOpt.NVarA+1:end)];
            else
                disp('Using HMNL results as starting values')
                Results_old.HMNL.bhat = Results_old.HMNL.bhat(:);
                b0 = [Results_old.HMNL.bhat(1:EstimOpt.NVarA); max(1,sqrt(abs(Results_old.HMNL.bhat(1:EstimOpt.NVarA)))); Results_old.HMNL.bhat(EstimOpt.NVarA+1:end)];
                if sum(EstimOpt.Dist == 1) > 0 && ~any(Results_old.HMNL.EstimOpt.MNLDist==1)
                    b0(EstimOpt.Dist == 1) = log(b0(EstimOpt.Dist == 1));
                end
            end
        else
            error('No starting values available - run HMNL first')
        end
    end
elseif EstimOpt.FullCov == 1
    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == (EstimOpt.NVarA*(1+EstimOpt.NLatent+EstimOpt.NVarM) + sum(1:EstimOpt.NVarA) +EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut)
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'HMXL') && isfield(Results_old.HMXL,'b0') % starting values provided
        Results_old.HMXL.b0_old = Results_old.HMXL.b0(:);
        Results_old.HMXL = rmfield(Results_old.HMXL,'b0');
        if length(Results_old.HMXL.b0_old) ~= (EstimOpt.NVarA*(1+EstimOpt.NLatent+EstimOpt.NVarM) + sum(1:EstimOpt.NVarA) +EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut)
            cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.HMXL = rmfield(Results_old.HMXL,'b0_old');
        else
            b0 = Results_old.HMXL.b0_old(:);
        end
    end
    if  ~exist('b0','var')
        if isfield(Results_old,'HMXL_d') && isfield(Results_old.HMXL_d,'bhat')
            disp('Using HMXL_d results as starting values')
            Results_old.HMXL_d.bhat = Results_old.HMXL_d.bhat(:);
            vc_tmp = diag(Results_old.HMXL_d.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2)).^2;
            b0 = [Results_old.HMXL_d.bhat(1:EstimOpt.NVarA); vc_tmp(tril(ones(size(vc_tmp)))==1);Results_old.HMXL_d.bhat(EstimOpt.NVarA*2+1:end)];
        elseif isfield(Results_old,'HMNL') && isfield(Results_old.HMNL,'bhat')
            disp('Using HMNL results as starting values')
            Results_old.HMNL.bhat = Results_old.HMNL.bhat(:);
            b0 = [Results_old.HMNL.bhat(1:EstimOpt.NVarA); zeros(sum(1:EstimOpt.NVarA,2),1); Results_old.HMNL.bhat(EstimOpt.NVarA+1:end)];
            if sum(EstimOpt.Dist == 1) > 0 && ~isfield(Results_old.HMNL.EstimOpt,'XDist')
                b0(EstimOpt.Dist == 1) = log(b0(EstimOpt.Dist == 1));
            end
        else
            error('No starting values available - run HMNL or HMXL_d first')
        end
    end
    
elseif EstimOpt.FullCov == 2 % allowing for correlation between random terms and LV
    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == (EstimOpt.NVarA*(1+EstimOpt.NLatent+EstimOpt.NVarM) + sum(1:(EstimOpt.NVarA+EstimOpt.NLatent))-EstimOpt.NLatent +EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut)
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'HMXL_e') && isfield(Results_old.HMXL_e,'b0') % starting values provided
        Results_old.HMXL_e.b0_old = Results_old.HMXL_e.b0(:);
        Results_old.HMXL_e = rmfield(Results_old.HMXL_e,'b0');
        if length(Results_old.HMXL_e.b0_old) ~= (EstimOpt.NVarA*(1+EstimOpt.NLatent+EstimOpt.NVarM) + sum(1:(EstimOpt.NVarA+EstimOpt.NLatent))-EstimOpt.NLatent +EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut)
            cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.HMXL_e = rmfield(Results_old.HMXL_e,'b0_old');
        else
            b0 = Results_old.HMXL.b0_old(:);
        end
    end
    if  ~exist('b0','var')
        if isfield(Results_old,'HMXL') && isfield(Results_old.HMXL,'bhat')
            disp('Using HMXL results as starting values')
            Results_old.HMXL.bhat = Results_old.HMXL.bhat(:);
            VC = tril(ones(EstimOpt.NVarA));
            VC(VC==1)= Results_old.HMXL.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA));
            VC = [VC, zeros(EstimOpt.NVarA,EstimOpt.NLatent)];
            VC = [VC; zeros(EstimOpt.NLatent, EstimOpt.NVarA+EstimOpt.NLatent)];
            VCtmp = tril(ones(size(VC)));
            VCtmp2 = diag(ones(EstimOpt.NLatent+EstimOpt.NVarA,1));
            VCtmp2(1:EstimOpt.NVarA, 1:EstimOpt.NVarA) = 0;
            VCtmp(VCtmp2 == 1) = 0;
            b0 = [Results_old.HMXL.bhat(1:EstimOpt.NVarA); VC(VCtmp==1);Results_old.HMXL.bhat(EstimOpt.NVarA+sum(1:EstimOpt.NVarA)+1:end)];
        else
            error('No starting values available - run HMXL')
        end
    end
end


%% Optimization Options


if  isfield(EstimOpt,'BActive')
    EstimOpt.BActive = EstimOpt.BActive(:)';
end

if isfield(EstimOpt, 'StrMatrix')
    if size(EstimOpt.StrMatrix,1) == EstimOpt.NLatent && (size(EstimOpt.StrMatrix,2) == 1 || size(EstimOpt.StrMatrix,2) == EstimOpt.NVarStr) && any(any(EstimOpt.StrMatrix ~= 1))
        if size(EstimOpt.StrMatrix,2) ~= EstimOpt.NVarStr
            EstimOpt.StrMatrix = repmat(EstimOpt.StrMatrix,[1,EstimOpt.NVarStr]);
        end                 
        if ~isfield(EstimOpt,'BActive')
            EstimOpt.BActive = ones(1,size(b0,1));
        end
        EstimOpt.BActive = StrSelect(EstimOpt,EstimOpt.FullCov+2);
    else
        error('Structural variables matrix (EstimOpt.StrMatrix) erroneously defined')
    end
end

if sum(EstimOpt.Dist == -1) > 0
    if isfield(EstimOpt,'BActive') == 0 || isempty(EstimOpt.BActive)
        EstimOpt.BActive = ones(1,length(b0));
    end
    if EstimOpt.FullCov == 0
        EstimOpt.BActive(EstimOpt.NVarA+find(EstimOpt.Dist == -1)) = 0;
    elseif EstimOpt.FullCov == 1
        Vt = tril(ones(EstimOpt.NVarA));
        Vt(EstimOpt.Dist ==-1,:) = 0;
        Vt(:,EstimOpt.Dist ==-1) = 0;
        EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA)) = EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA)) .* (Vt(tril(ones(size(Vt)))~=0)');
    elseif EstimOpt.FullCov == 2
        % TO BE DONE
        Vt = tril(ones(EstimOpt.NVarA+EstimOpt.NLatent));
        Vt(find(EstimOpt.Dist ==-1),:) = 0;
        Vt(:,find(EstimOpt.Dist ==-1)) = 0;
        VCtmp = tril(ones(size(Vt)));
        VCtmp2 = diag(ones(EstimOpt.NLatent+EstimOpt.NVarA,1));
        VCtmp2(1:EstimOpt.NVarA, 1:EstimOpt.NVarA) = 0;
        VCtmp(VCtmp2 == 1) = 0;
        EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent) = EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent) .* (Vt(VCtmp~=0)');
    end
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
    err_mtx = err_mtx(:,2:end);
    if EstimOpt.NP*EstimOpt.NRep < 3e+7
        err_mtx = icdf('Normal',err_mtx,0,1); %to be cut down later
    else % this is for very large number of draws * variables
        for i = 1:EstimOpt.NLatent+EstimOpt.NVarA
            err_mtx(:,i) = icdf('Normal',err_mtx(:,i),0,1); %to be cut down later
        end
    end
end

% err_mtx(:,EstimOpt.NLatent + find(EstimOpt.Dist == -1)) = 0;
% err_mtx(:,find(EstimOpt.Dist == -1)) = 0;
err_mtx(:,EstimOpt.Dist == -1) = 0;
err_sliced = err_mtx'; % NVarA + NLatent x NRep * NP

if isfield(EstimOpt, 'Drawskeep') && ~isempty(EstimOpt.Drawskeep) && EstimOpt.Drawskeep == 1
    Results.err = err_sliced;
end


%% Display Options


if ((isfield(EstimOpt, 'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    OptimOpt.GradObj = 'on';
end

% if any(EstimOpt.MissingAlt(:) == 1) && EstimOpt.NumGrad == 0
%    EstimOpt.NumGrad = 1;
%    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - missing alternatives not supported by analytical gradient \n')
% end

if EstimOpt.FullCov == 2 && EstimOpt.NumGrad == 0 && (any(EstimOpt.MissingAlt(:) == 1) || EstimOpt.NLatent > 1 || any(EstimOpt.MeaSpecMatrix > 0))
    EstimOpt.NumGrad = 1;
    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - correlation of random parameters and LV not supported by analytical gradient with missing, more than 1 LV and not OLS measurment equations\n')
end
if any(EstimOpt.MeaSpecMatrix >= 3) && EstimOpt.NumGrad == 0 && any(any(INPUT.Xmea(:, EstimOpt.MeaSpecMatrix >=3) > 100))
    cprintf(rgb('DarkOrange'), 'WARNING: it is recommended to switch to numerical gradient, as analitycal can be not precise when Xmea take large values for NB \n')
end

if any(EstimOpt.Dist > 1) && EstimOpt.NumGrad == 0
    EstimOpt.NumGrad = 1;
    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - analytical gradient available for normally or lognormally distributed parameters only \n')
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


LLfun = @(B) LL_hmxl_MATlike(INPUT.YY,INPUT.XXa,INPUT.Xm,INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp, err_sliced,INPUT.W,EstimOpt,OptimOpt,B);
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


% save tmp_HMXL


Results.LL = -LL;
R2 = R2_hybrid(INPUT.YY,INPUT.XXa,INPUT.Xm,INPUT.Xstr,[],INPUT.MissingInd,err_sliced,EstimOpt,Results.bhat,2);

Results.b0_old = b0;
if EstimOpt.HessEstFix == 1
    if isequal(OptimOpt.GradObj,'on') && EstimOpt.NumGrad == 0
        [~, Results.jacobian] = LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm, INPUT.Xstr, INPUT.Xmea, INPUT.Xmea_exp, err_sliced, EstimOpt, Results.bhat);
    else
        f = LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm, INPUT.Xstr, INPUT.Xmea, INPUT.Xmea_exp, err_sliced, EstimOpt, Results.bhat);
        Results.jacobian = numdiff(@(B) LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp,  err_sliced, EstimOpt,B), f, Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
    end
    Results.jacobian = INPUT.W(:,ones(size(Results.jacobian,2),1)).*Results.jacobian;
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) INPUT.W.*LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp,  err_sliced, EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.hess = hessian(@(B) sum(INPUT.W.*LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp,  err_sliced, EstimOpt,B),1), Results.bhat);
end
if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
    Results.hess = Results.jacobian(:,EstimOpt.BActive == 1)'*Results.jacobian(:,EstimOpt.BActive == 1);
    EstimOpt.BLimit = (sum(Results.hess) == 0);
    EstimOpt.BLimit = direcXpnd(EstimOpt.BLimit,EstimOpt.BActive);
else
    EstimOpt.BLimit = (sum(Results.hess) == 0 & EstimOpt.BActive == 1);
    EstimOpt.BActive(EstimOpt.BLimit == 1) = 0;
    Results.hess = Results.hess(EstimOpt.BActive == 1,EstimOpt.BActive == 1);
end
Results.ihess = inv(Results.hess);
Results.ihess = direcXpnd(Results.ihess,EstimOpt.BActive);
Results.ihess = direcXpnd(Results.ihess',EstimOpt.BActive);
if EstimOpt.RobustStd == 1
    if EstimOpt.NumGrad == 0
        [~, Results.jacobian] = LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm, INPUT.Xstr, INPUT.Xmea, INPUT.Xmea_exp, err_sliced, EstimOpt, Results.bhat);
        Results.jacobian = Results.jacobian.*INPUT.W(:, ones(1,size(Results.jacobian,2)));
    else
        Results.LLdetailed = LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm, INPUT.Xstr, INPUT.Xmea, INPUT.Xmea_exp, err_sliced, EstimOpt, Results.bhat);
        Results.jacobian = numdiff(@(B) INPUT.W.*LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm, INPUT.Xstr, INPUT.Xmea, INPUT.Xmea_exp, err_sliced, EstimOpt, B) ,INPUT.W.*Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
    end
    RobustHess = Results.jacobian'*Results.jacobian;
    Results.ihess = Results.ihess*RobustHess*Results.ihess;
end

Results.std = sqrt(diag(Results.ihess));
Results.std(EstimOpt.BActive == 0) = NaN;
Results.std(EstimOpt.BLimit == 1) = 0;
Results.std(imag(Results.std) ~= 0) = NaN;


if EstimOpt.FullCov == 0
    Results.DetailsA = [Results.bhat(1:EstimOpt.NVarA), Results.std(1:EstimOpt.NVarA), pv(Results.bhat(1:EstimOpt.NVarA), Results.std(1:EstimOpt.NVarA))];
    %Results.DetailsV = [Results.bhat(EstimOpt.NVarA+1:2*EstimOpt.NVarA).^2, 2*abs(Results.bhat(EstimOpt.NVarA+1:2*EstimOpt.NVarA)).*Results.std(EstimOpt.NVarA+1:2*EstimOpt.NVarA), pv(Results.bhat(EstimOpt.NVarA+1:2*EstimOpt.NVarA).^2, 2*abs(Results.bhat(EstimOpt.NVarA+1:2*EstimOpt.NVarA)).*Results.std(EstimOpt.NVarA+1:2*EstimOpt.NVarA))];
    Results.DetailsV = [abs(Results.bhat(EstimOpt.NVarA+1:2*EstimOpt.NVarA)), Results.std(EstimOpt.NVarA+1:2*EstimOpt.NVarA), pv(abs(Results.bhat(EstimOpt.NVarA+1:2*EstimOpt.NVarA)), Results.std(EstimOpt.NVarA+1:2*EstimOpt.NVarA))];
    
    if EstimOpt.NVarM > 0
        Results.DetailsCM = [Results.bhat(2*EstimOpt.NVarA+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)), Results.std(2*EstimOpt.NVarA+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)), pv(Results.bhat(2*EstimOpt.NVarA+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)), Results.std(2*EstimOpt.NVarA+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)))];
    else
        Results.DetailsCM = [];
    end
    Results.DetailsL = [Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NVarM)+1:EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)), Results.std(EstimOpt.NVarA*(2+EstimOpt.NVarM)+1:EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)), pv(Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NVarM)+1:EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)), Results.std(EstimOpt.NVarA*(2+EstimOpt.NVarM)+1:EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)))];
    Results.DetailsS = [Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA), Results.std(EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA), pv(Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA), Results.std(EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA))];
    Results.DetailsM = [Results.bhat((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA+1:end), Results.std((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA+1:end), pv(Results.bhat((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA+1:end), Results.std((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA+1:end))];
elseif EstimOpt.FullCov == 1
    Results.DetailsA = [Results.bhat(1:EstimOpt.NVarA), Results.std(1:EstimOpt.NVarA), pv(Results.bhat(1:EstimOpt.NVarA), Results.std(1:EstimOpt.NVarA))];
    Results.DetailsV = sdtri(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2), Results.ihess(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2,EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2),EstimOpt);
    %   This is to calculate DetailsV with delta method - used for endogeneity study
    %     H = jacobianest(@(b) sdtriHe2(b, EstimOpt, -1), Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2)); % standard deviations
    %     covtmp = Results.ihess(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2,EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2);
    %     covtmp = covtmp(EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2), EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2));
    %     VarTmp = H*covtmp*H';
    %     Results.DetailsV = zeros(EstimOpt.NVarA,3);
    %     Results.DetailsV(EstimOpt.Dist ~= -1,:) = [sdtriHe2(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2),EstimOpt,-1) , sqrt(diag(VarTmp)), pv(sdtriHe2(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2), EstimOpt, -1), sqrt(diag(VarTmp)))];
    %
    l = EstimOpt.NVarA*(EstimOpt.NVarA+3)/2;
    if EstimOpt.NVarM > 0
        Results.DetailsCM = [Results.bhat(l+1:l+EstimOpt.NVarA*EstimOpt.NVarM), Results.std(l+1:l+EstimOpt.NVarA*EstimOpt.NVarM), pv(Results.bhat(l+1:l+EstimOpt.NVarA*EstimOpt.NVarM), Results.std(l+1:l+EstimOpt.NVarA*EstimOpt.NVarM))];
    else
        Results.DetailsCM = [];
    end
    l = l + EstimOpt.NVarM;
    Results.DetailsL = [Results.bhat(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent), Results.std(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent), pv(Results.bhat(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent), Results.std(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent))];
    Results.DetailsS = [Results.bhat(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent), Results.std(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent), pv(Results.bhat(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent), Results.std(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent))];
    Results.DetailsM = [Results.bhat(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end), Results.std(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end), pv(Results.bhat(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end), Results.std(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end))];
    
elseif EstimOpt.FullCov == 2
    Results.DetailsA = [Results.bhat(1:EstimOpt.NVarA), Results.std(1:EstimOpt.NVarA), pv(Results.bhat(1:EstimOpt.NVarA), Results.std(1:EstimOpt.NVarA))];
    %     VC = tril(ones(EstimOpt.NLatent+EstimOpt.NVarA));
    %     VC(EstimOpt.NVarA+1:end,:) = 0;
    %     VCtmp = tril(ones(size(VC)));
    %     VCtmp2 = diag(ones(EstimOpt.NLatent+EstimOpt.NVarA,1));
    %     VCtmp2(1:EstimOpt.NVarA, 1:EstimOpt.NVarA) = 0;
    %     VCtmp(VCtmp2 == 1) = 0;
    %     Indx = VC(VCtmp ==1);
    %     bhattmp = Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent);
    %     bhattmp = bhattmp(Indx == 1);
    %     covtmp =  Results.ihess(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent,EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent);
    %     covtmp = covtmp(Indx == 1,Indx==1);
    %     BActivetmp = EstimOpt.BActive;
    %     EstimOpt.BActive = EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent);
    %     EstimOpt.BActive = EstimOpt.BActive(Indx == 1);
    %     EstimOpt.BActive = [ones(1,EstimOpt.NVarA), EstimOpt.BActive];
    %     Results.DetailsV = sdtri(bhattmp, covtmp,EstimOpt);
    %     EstimOpt.BActive = BActivetmp;
    l = EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent;
    if EstimOpt.NVarM > 0
        Results.DetailsCM = [Results.bhat(l+1:l+EstimOpt.NVarA*EstimOpt.NVarM), Results.std(l+1:l+EstimOpt.NVarA*EstimOpt.NVarM), pv(Results.bhat(l+1:l+EstimOpt.NVarA*EstimOpt.NVarM), Results.std(l+1:l+EstimOpt.NVarA*EstimOpt.NVarM))];
    else
        Results.DetailsCM = [];
    end
    l = l + EstimOpt.NVarM;
    Results.DetailsL = [Results.bhat(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent), Results.std(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent), pv(Results.bhat(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent), Results.std(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent))];
    Results.DetailsS = [Results.bhat(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent), Results.std(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent), pv(Results.bhat(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent), Results.std(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent))];
    Results.DetailsM = [Results.bhat(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end), Results.std(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end), pv(Results.bhat(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end), Results.std(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end))];
    if any(EstimOpt.Dist == -1)
        VC = tril(ones(EstimOpt.NLatent+EstimOpt.NVarA));
        VC(EstimOpt.Dist == -1,:) = 0;
        VC(:,EstimOpt.Dist == -1) = 0;
        VCtmp = tril(ones(size(VC)));
        VCtmp2 = diag(ones(EstimOpt.NLatent+EstimOpt.NVarA,1));
        VCtmp2(1:EstimOpt.NVarA, 1:EstimOpt.NVarA) = 0;
        VCtmp(VCtmp2 == 1) = 0;
        Indx = VC(VCtmp ==1);
        bhattmp = Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent);
        bhattmp = bhattmp(Indx == 1);
        covtmp =  Results.ihess(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent,EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent);
        covtmp = covtmp(Indx == 1,Indx==1);
        %         [Results.DetailsCOV, Results.DetailsCORR] = sdtriHe(bhattmp, covtmp,EstimOpt);
        
        % Using delta method instead simulation
        H = jacobianest(@(b) sdtriHe2(b, EstimOpt, 0), bhattmp); % Covariance
        VarTmp = H*covtmp*H';
        Results.DetailsCOV = [sdtriHe2(bhattmp, EstimOpt, 0), sqrt(diag(VarTmp)), pv(sdtriHe2(bhattmp, EstimOpt, 0), sqrt(diag(VarTmp)))];
        H = jacobianest(@(b) sdtriHe2(b, EstimOpt, 1), bhattmp); % Correlation
        VarTmp = H*covtmp*H';
        Results.DetailsCORR = [sdtriHe2(bhattmp, EstimOpt, 1), sqrt(diag(VarTmp)), pv(sdtriHe2(bhattmp, EstimOpt, 1), sqrt(diag(VarTmp)))];
        
        H = jacobianest(@(b) sdtriHe2(b, EstimOpt, 2), bhattmp); % standard deviations
        VarTmp = H*covtmp*H';
        Results.DetailsV = zeros(EstimOpt.NVarA,3);
        Results.DetailsV(EstimOpt.Dist ~= -1,:) = [sdtriHe2(bhattmp,EstimOpt,2) , sqrt(diag(VarTmp)), pv(sdtriHe2(bhattmp, EstimOpt, 2), sqrt(diag(VarTmp)))];
    else
        
        bhattmp = Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent);
        covtmp =  Results.ihess(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent,EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent);
        % Using delta method instead simulation
        H = jacobianest(@(b) sdtriHe2(b, EstimOpt, 0), bhattmp); % Covariance
        VarTmp = H*covtmp*H';
        Results.DetailsCOV = [sdtriHe2(bhattmp, EstimOpt, 0), sqrt(diag(VarTmp)), pv(sdtriHe2(bhattmp, EstimOpt, 0), sqrt(diag(VarTmp)))];
        H = jacobianest(@(b) sdtriHe2(b, EstimOpt, 1), bhattmp); % Correlation
        VarTmp = H*covtmp*H';
        Results.DetailsCORR = [sdtriHe2(bhattmp, EstimOpt, 1), sqrt(diag(VarTmp)), pv(sdtriHe2(bhattmp, EstimOpt, 1), sqrt(diag(VarTmp)))];
        %[Results.DetailsCOV, Results.DetailsCORR] = sdtriHe(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent), Results.ihess(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent,EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent),EstimOpt);
        
        H = jacobianest(@(b) sdtriHe2(b, EstimOpt, 2), bhattmp); % standard deviations
        VarTmp = H*covtmp*H';
        Results.DetailsV = [sdtriHe2(bhattmp,EstimOpt,2) , sqrt(diag(VarTmp)), pv(sdtriHe2(bhattmp, EstimOpt, 2), sqrt(diag(VarTmp)))];
        
        
    end
end

disp(' ')
disp('Means');
disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
disp([char(EstimOpt.NamesA) ,blanks(EstimOpt.NVarA)', num2str(Results.DetailsA(:,1),'%8.4f') star_sig(Results.DetailsA(:,3)) num2str(Results.DetailsA(:,2:3),'%8.4f %8.4f')])

disp(' ')
disp('Standard deviations');
disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
disp([char(EstimOpt.NamesA) ,blanks(EstimOpt.NVarA)', num2str(Results.DetailsV(:,1),'%8.4f') star_sig(Results.DetailsV(:,3)) num2str(Results.DetailsV(:,2:3),'%8.4f %8.4f')])

if EstimOpt.NVarM > 0
    for i = 1:EstimOpt.NVarM
        disp(' ');
        disp(['Explanatory variable of random parameters'' means - ', char(EstimOpt.NamesM(i))]);
        disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
        disp([char(EstimOpt.NamesA),blanks(EstimOpt.NVarA)',num2str(Results.DetailsCM((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA,1),'%8.4f'), star_sig(Results.DetailsCM((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA,3)), num2str(Results.DetailsCM((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA,2:3),'%8.4f %8.4f')])
    end
end

for i = 1:EstimOpt.NLatent
    disp(' ')
    disp(num2str(i,'Interactions with Latent Variable %1.0f'));
    disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
    disp([char(EstimOpt.NamesA),blanks(EstimOpt.NVarA)',num2str(Results.DetailsL((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA,1),'%8.4f'), star_sig(Results.DetailsL((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA,3)), num2str(Results.DetailsL((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA,2:3),'%8.4f %8.4f')])
end

for i = 1:EstimOpt.NLatent
    disp(' ')
    disp(num2str(i,'Structural equation of Latent Variable %1.0f'));
    disp(['var.', blanks(size(char(EstimOpt.NamesStr),2)-2) ,'coef.      st.err.  p-value'])
    disp([char(EstimOpt.NamesStr),blanks(EstimOpt.NVarStr)',num2str(Results.DetailsS((i-1)*EstimOpt.NVarStr+1:i*EstimOpt.NVarStr,1),'%8.4f'), star_sig(Results.DetailsS((i-1)*EstimOpt.NVarStr+1:i*EstimOpt.NVarStr,3)), num2str(Results.DetailsS((i-1)*EstimOpt.NVarStr+1:i*EstimOpt.NVarStr,2:3),'%8.4f %8.4f')])
end

l = 0;
for i = 1:size(INPUT.Xmea,2)
    tmp = EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0);
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
                stdx = sqrt(g(1:n)'*Results.ihess((EstimOpt.NVarStr)*EstimOpt.NLatent+l+1:(EstimOpt.NVarStr)*EstimOpt.NLatent+l+n,(EstimOpt.NVarStr)*EstimOpt.NLatent+l+1:(EstimOpt.NVarStr)*EstimOpt.NLatent+l+n)*g(1:n));
                stdx(~isreal(stdx)) = NaN;
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
    elseif EstimOpt.MeaSpecMatrix(i) == 6
        disp('Estimated using Zero Inflated Negative Binomial regression')
        disp('Probability of Non-participation (logit)')
        disp('var.   coef.     st.err.  p-value')
        disp([char(EstimOpt.NamesLV(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%7.4f %8.4f')])
        l = l+sum(EstimOpt.MeaMatrix(:,i))+1+tmp;
        
        disp('Negative binomial model')
        disp('var.   coef.     st.err.  p-value')
        Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1:3) = [exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)), Results.DetailsM(l+tmp+sum(EstimOpt.MeaMatrix(:,i))+2,2)*exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)),pv(exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)), Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,2)*exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)))];
        disp([char(EstimOpt.NamesLV(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%7.4f %8.4f')])
        l = l+sum(EstimOpt.MeaMatrix(:,i))+2+tmp;
    end
end

Results.R = [Results.DetailsA; Results.DetailsV ; Results.DetailsCM ; Results.DetailsL ; Results.DetailsS;Results.DetailsM];
Results.LL0 = Results.MIMIC0.LL + Results_old.MNL0.LL;
EstimOpt.params = length(Results.bhat);
if isfield(EstimOpt,'BActive')
    EstimOpt.params = EstimOpt.params - sum(EstimOpt.BActive == 0);
end

Results.stats = [Results.LL0; Results.LL; 1-Results.LL/Results.LL0;R2; ((2*EstimOpt.params-2*Results.LL) + 2*EstimOpt.params*(EstimOpt.params+1)/(EstimOpt.NObs-EstimOpt.params-1))/EstimOpt.NObs; EstimOpt.NObs; EstimOpt.params];

Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;
Results.INPUT = INPUT;

Results.R_out  = cell(3+EstimOpt.NVarStr +EstimOpt.NVarA+3+ 3*size(INPUT.Xmea,2)+EstimOpt.NVarcut+EstimOpt.NVarMea+ 2 + 7, 7+3*(EstimOpt.NLatent+EstimOpt.NVarM));
if EstimOpt.FullCov == 0
    Results.R_out(1,1) = {'HMXL_d'};
else
    Results.R_out(1,1) = {'HMXL'};
end
head = {'var.' , 'coef.', 'st.err.' , 'p-value'};
headx = [head, repmat(head(1,2:4),1,1+EstimOpt.NLatent+EstimOpt.NVarM)];
LVlist = {'LV 1' , 'LV 2' , 'LV 3' , 'LV 4' , 'LV 5' , 'LV 6' , 'LV 7' , 'LV 8' ,'LV 9' ,'LV 10'};
if EstimOpt.NLatent<=10
    Results.R_out(2,2:3:(2+3+EstimOpt.NLatent*3+EstimOpt.NVarM*3)) = [{'Means', 'Standard Deviations'}, LVlist(1,1:EstimOpt.NLatent), EstimOpt.NamesM'];
end
Results.R_out(3,:) = headx;
Results.R_out(4:(EstimOpt.NVarA+3),1:7) = [EstimOpt.NamesA,num2cell(Results.DetailsA),num2cell(Results.DetailsV)];
for i = 1:EstimOpt.NLatent
    Results.R_out(4:(EstimOpt.NVarA+3),8+(i-1)*3:7+i*3) = num2cell(Results.DetailsL(1+(i-1)*EstimOpt.NVarA:i*EstimOpt.NVarA,:));
end
if EstimOpt.NVarM > 0
    for i = 1:EstimOpt.NVarM
        Results.R_out(4:(EstimOpt.NVarA+3),7+3*EstimOpt.NLatent+1+(i-1)*3:7+3*EstimOpt.NLatent+i*3) = num2cell(Results.DetailsCM((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA,:));
    end
end
Results.R_out(EstimOpt.NVarA+4,1) = {'Structural equations'};
Results.R_out(EstimOpt.NVarA+5,2:3:(2+(EstimOpt.NLatent-1)*3)) = LVlist(1,1:EstimOpt.NLatent);
Results.R_out(EstimOpt.NVarA+6,1:1+3*EstimOpt.NLatent) = headx(1:1+3*EstimOpt.NLatent);
Results.R_out(EstimOpt.NVarA+7:EstimOpt.NVarA+6+EstimOpt.NVarStr,1) = EstimOpt.NamesStr;
for i = 1: EstimOpt.NLatent
    Results.R_out(EstimOpt.NVarA+7:EstimOpt.NVarA+6+EstimOpt.NVarStr,2 + 3*(i-1):1+3*i) = num2cell(Results.DetailsS((i-1)*EstimOpt.NVarStr+1:i*EstimOpt.NVarStr,:));
end
l = EstimOpt.NVarA+3; % this is for indexing in R_out
k = 0;
for i = 1:size(INPUT.Xmea,2)
    Results.R_out(EstimOpt.NVarStr+3+l+1,1) =  cellfun(@(x)[x char(EstimOpt.NamesMea(i))],{'Measurment equation for '},'UniformOutput',0);
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
    elseif EstimOpt.MeaSpecMatrix(i) == 6
        model = 'ZINB';
    end
    Results.R_out(EstimOpt.NVarStr+3+l+2,1) = cellfun(@(x)[x model],{'Estimated using '},'UniformOutput',0);
    Results.R_out(EstimOpt.NVarStr+3+l+3,1:4) = head;
    Results.R_out(EstimOpt.NVarStr+3+l+4:EstimOpt.NVarStr+3+l+3 + EstimOpt.CutMatrix(i) ,1:4) = [EstimOpt.NamesLV(k+1:k+EstimOpt.CutMatrix(i)), num2cell(Results.DetailsM(k+1:k+EstimOpt.CutMatrix(i),:))];
    l = l+3+EstimOpt.CutMatrix(i);
    k = k + EstimOpt.CutMatrix(i);
end
Results.R_out(EstimOpt.NVarStr+3+l+2,1) = {'Model characteristics'};
Results.R_out(EstimOpt.NVarStr+3+l+3:end,1) = {'LL0'; 'LL' ; 'McFadden R2';'Ben-Akiva R2' ;'AIC/n' ; 'n'; 'k'};
Results.R_out(EstimOpt.NVarStr+3+l+3:end,2) = num2cell(Results.stats);


disp(' ')
disp(['LL at convergence: ',num2str(Results.LL,'%8.4f')])
disp(' ');
Results.clocknote = clock;
Results.tocnote = toc;
[~,DayName] = weekday(now,'long');
disp(['Estimation completed on ' DayName ', ' num2str(Results.clocknote(1)) '-' sprintf('%02.0f',Results.clocknote(2)) '-' sprintf('%02.0f',Results.clocknote(3)) ' at ' sprintf('%02.0f',Results.clocknote(4)) ':' sprintf('%02.0f',Results.clocknote(5)) ':' sprintf('%02.0f',Results.clocknote(6))])
disp(['Estimation took ' num2str(Results.tocnote) ' seconds ('  num2str(floor(Results.tocnote/(60*60))) ' hours ' num2str(floor(rem(Results.tocnote,60*60)/60)) ' minutes ' num2str(rem(Results.tocnote,60)) ' seconds).']);
disp(' ');