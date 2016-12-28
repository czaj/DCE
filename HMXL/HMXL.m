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
if isfield(INPUT, 'Xs') == 0
    INPUT.Xs = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarS = size(INPUT.Xs,2); % Number of covariates of scale

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

if isfield(EstimOpt,'CountCut') == 0
    EstimOpt.CountCut = 80; 
end
if any(EstimOpt.MeaSpecMatrix == 3 | EstimOpt.MeaSpecMatrix == 4 | EstimOpt.MeaSpecMatrix == 5 | EstimOpt.MeaSpecMatrix == 6)
   IndxTmp =  EstimOpt.MeaSpecMatrix == 3 | EstimOpt.MeaSpecMatrix == 4 | EstimOpt.MeaSpecMatrix == 5 | EstimOpt.MeaSpecMatrix == 6;
   XmeaTmp = INPUT.Xmea(:, IndxTmp);
   if any(XmeaTmp(:) > EstimOpt.CountCut)
        XmeaTmp(XmeaTmp > EstimOpt.CountCut) = EstimOpt.CountCut;
        INPUT.Xmea(:, IndxTmp) = XmeaTmp;
        cprintf(rgb('DarkOrange'), ['WARNING: Values of count data measurment equations are too high, they were censored to ', num2str(EstimOpt.CountCut), '\n'])
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
        % rescale Xmea MNL so that the minimum value is 1
        if min(unique(INPUT.Xmea(INPUT.MissingInd==0,i))) <= 0
            INPUT.Xmea(INPUT.MissingInd==0,i) = INPUT.Xmea(INPUT.MissingInd==0,i) - min(unique(INPUT.Xmea(INPUT.MissingInd==0,i))) + 1;
        end
        % test if only integer valus
        if ~all(ceil(INPUT.Xmea(INPUT.MissingInd==0,i)) == floor(INPUT.Xmea(INPUT.MissingInd==0,i)))
            error(num2str(i,'Measurement variable %1.0f is modelled using MNL - all values must be integer'));
        end   
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
if EstimOpt.NVarS > 0
    if isfield(EstimOpt,'NamesS') == 0 || isempty(EstimOpt.NamesS) || length(EstimOpt.NamesS) ~= EstimOpt.NVarS
        EstimOpt.NamesS = (1:EstimOpt.NVarS)';
        EstimOpt.NamesS = cellstr(num2str(EstimOpt.NamesS));
    elseif size(EstimOpt.NamesS,1) ~= EstimOpt.NVarS
        EstimOpt.NamesS = EstimOpt.NamesS';
    end
end

if isfield(EstimOpt, 'ScaleLV') == 0
    EstimOpt.ScaleLV = 0;
elseif EstimOpt.ScaleLV == 1
    EstimOpt.NVarS = EstimOpt.NVarS+EstimOpt.NLatent;
end    
if EstimOpt.ScaleLV == 1
    if EstimOpt.NVarS == EstimOpt.NLatent
        EstimOpt.NamesS = [];
    end
   for i = 1:EstimOpt.NLatent
        EstimOpt.NamesS = [EstimOpt.NamesS;cellfun(@(x)[x num2str(i)],{'LV '},'UniformOutput',0)];
   end
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
INPUT.Xm = INPUT.Xm(1:EstimOpt.NAlt*EstimOpt.NCT:end,:)'; % NVarM x NP
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
   if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == (EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM) + EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut+EstimOpt.NVarS)
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'HMXL_d') && isfield(Results_old.HMXL_d,'b0') % starting values provided
        Results_old.HMXL_d.b0_old = Results_old.HMXL_d.b0(:);
        Results_old.HMXL_d = rmfield(Results_old.HMXL_d,'b0');
        if length(Results_old.HMXL_d.b0_old) ~= (EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM) + EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut+EstimOpt.NVarS)
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
                    b0(EstimOpt.Dist == 1) = log(abs(b0(EstimOpt.Dist == 1)));
                end
            end
        else
            error('No starting values available - run HMNL first')
        end
    end
elseif EstimOpt.FullCov == 1
    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == (EstimOpt.NVarA*(1+EstimOpt.NLatent+EstimOpt.NVarM) + sum(1:EstimOpt.NVarA) +EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut+EstimOpt.NVarS)
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'HMXL') && isfield(Results_old.HMXL,'b0') % starting values provided
        Results_old.HMXL.b0_old = Results_old.HMXL.b0(:);
        Results_old.HMXL = rmfield(Results_old.HMXL,'b0');
        if length(Results_old.HMXL.b0_old) ~= (EstimOpt.NVarA*(1+EstimOpt.NLatent+EstimOpt.NVarM) + sum(1:EstimOpt.NVarA) +EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut+EstimOpt.NVarS)
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
                b0(EstimOpt.Dist == 1) = log(abs(b0(EstimOpt.Dist == 1)));
            end
        else
            error('No starting values available - run HMNL or HMXL_d first')
        end
    end
    
elseif EstimOpt.FullCov == 2 % allowing for correlation between random terms and LV
    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == (EstimOpt.NVarA*(1+EstimOpt.NLatent+EstimOpt.NVarM) + sum(1:(EstimOpt.NVarA+EstimOpt.NLatent))-EstimOpt.NLatent +EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut+EstimOpt.NVarS)
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'HMXL_e') && isfield(Results_old.HMXL_e,'b0') % starting values provided
        Results_old.HMXL_e.b0_old = Results_old.HMXL_e.b0(:);
        Results_old.HMXL_e = rmfield(Results_old.HMXL_e,'b0');
        if length(Results_old.HMXL_e.b0_old) ~= (EstimOpt.NVarA*(1+EstimOpt.NLatent+EstimOpt.NVarM) + sum(1:(EstimOpt.NVarA+EstimOpt.NLatent))-EstimOpt.NLatent +EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut+EstimOpt.NVarS)
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
if EstimOpt.ScaleLV == 1 && EstimOpt.WTP_space > 0
    cprintf(rgb('DarkOrange'), 'WARNING: LV in scale cannot be identified in WTP-space with log-normal cost parameter unless there are some further constraints in the model\n')
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


LLfun = @(B) LL_hmxl_MATlike(INPUT.YY,INPUT.XXa,INPUT.Xm,INPUT.Xs,INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp,err_sliced,INPUT.W,EstimOpt,OptimOpt,B);
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
R2 = R2_hybrid(INPUT.YY,INPUT.XXa,INPUT.Xstr,[],INPUT.Xm,INPUT.Xs,INPUT.MissingInd,err_sliced,EstimOpt,Results.bhat,2);

Results.b0_old = b0;
if EstimOpt.HessEstFix == 1
    if isequal(OptimOpt.GradObj,'on') && EstimOpt.NumGrad == 0
        [~, Results.jacobian] = LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm,INPUT.Xs, INPUT.Xstr, INPUT.Xmea, INPUT.Xmea_exp, err_sliced, EstimOpt, Results.bhat);
    else
        f = LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm,INPUT.Xs, INPUT.Xstr, INPUT.Xmea, INPUT.Xmea_exp, err_sliced, EstimOpt, Results.bhat);
        Results.jacobian = numdiff(@(B) LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm,INPUT.Xs, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp,  err_sliced, EstimOpt,B), f, Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
    end
    Results.jacobian = INPUT.W(:,ones(size(Results.jacobian,2),1)).*Results.jacobian;
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) INPUT.W.*LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm,INPUT.Xs, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp,  err_sliced, EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.hess = hessian(@(B) sum(INPUT.W.*LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm,INPUT.Xs, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp,  err_sliced, EstimOpt,B),1), Results.bhat);
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
        [~, Results.jacobian] = LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm,INPUT.Xs, INPUT.Xstr, INPUT.Xmea, INPUT.Xmea_exp, err_sliced, EstimOpt, Results.bhat);
        Results.jacobian = Results.jacobian.*INPUT.W(:, ones(1,size(Results.jacobian,2)));
    else
        Results.LLdetailed = LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm,INPUT.Xs, INPUT.Xstr, INPUT.Xmea, INPUT.Xmea_exp, err_sliced, EstimOpt, Results.bhat);
        Results.jacobian = numdiff(@(B) INPUT.W.*LL_hmxl(INPUT.YY,INPUT.XXa, INPUT.Xm,INPUT.Xs, INPUT.Xstr, INPUT.Xmea, INPUT.Xmea_exp, err_sliced, EstimOpt, B) ,INPUT.W.*Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
    end
    RobustHess = Results.jacobian'*Results.jacobian;
    Results.ihess = Results.ihess*RobustHess*Results.ihess;
end

Results.std = sqrt(diag(Results.ihess));
Results.std(EstimOpt.BActive == 0) = NaN;
Results.std(EstimOpt.BLimit == 1) = 0;
Results.std(imag(Results.std) ~= 0) = NaN;


if EstimOpt.FullCov == 0
    Results.DetailsA(:,1) = Results.bhat(1:EstimOpt.NVarA);
	Results.DetailsA(:,3:4) = [Results.std(1:EstimOpt.NVarA), pv(Results.bhat(1:EstimOpt.NVarA), Results.std(1:EstimOpt.NVarA))];
    %Results.DetailsV = [Results.bhat(EstimOpt.NVarA+1:2*EstimOpt.NVarA).^2, 2*abs(Results.bhat(EstimOpt.NVarA+1:2*EstimOpt.NVarA)).*Results.std(EstimOpt.NVarA+1:2*EstimOpt.NVarA), pv(Results.bhat(EstimOpt.NVarA+1:2*EstimOpt.NVarA).^2, 2*abs(Results.bhat(EstimOpt.NVarA+1:2*EstimOpt.NVarA)).*Results.std(EstimOpt.NVarA+1:2*EstimOpt.NVarA))];
    Results.DetailsV(:,1) = abs(Results.bhat(EstimOpt.NVarA+1:2*EstimOpt.NVarA));
    Results.DetailsV(:,3:4) = [Results.std(EstimOpt.NVarA+1:2*EstimOpt.NVarA), pv(abs(Results.bhat(EstimOpt.NVarA+1:2*EstimOpt.NVarA)), Results.std(EstimOpt.NVarA+1:2*EstimOpt.NVarA))]; 
	
    if EstimOpt.NVarM > 0
			Results.DetailsCM(:,1) = Results.bhat(2*EstimOpt.NVarA+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)); 
			Results.DetailsCM(:,3:4) = [Results.std(2*EstimOpt.NVarA+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)), pv(Results.bhat(2*EstimOpt.NVarA+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)), Results.std(2*EstimOpt.NVarA+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)))];
	else
        Results.DetailsCM = [];
    end
    Results.DetailsL(:,1) = Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NVarM)+1:EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM));
    Results.DetailsL(:,3:4) = [Results.std(EstimOpt.NVarA*(2+EstimOpt.NVarM)+1:EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)), pv(Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NVarM)+1:EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)), Results.std(EstimOpt.NVarA*(2+EstimOpt.NVarM)+1:EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)))];
	if EstimOpt.NVarS > 0
        Results.DetailsScale(:,1) = Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NVarM+EstimOpt.NLatent)+1:EstimOpt.NVarA*(2+EstimOpt.NVarM+EstimOpt.NLatent)+EstimOpt.NVarS);
		Results.DetailsScale(:,3:4) = [Results.std(EstimOpt.NVarA*(2+EstimOpt.NVarM+EstimOpt.NLatent)+1:EstimOpt.NVarA*(2+EstimOpt.NVarM+EstimOpt.NLatent)+EstimOpt.NVarS), pv(Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NVarM+EstimOpt.NLatent)+1:EstimOpt.NVarA*(2+EstimOpt.NVarM+EstimOpt.NLatent)+EstimOpt.NVarS), Results.std(EstimOpt.NVarA*(2+EstimOpt.NVarM+EstimOpt.NLatent)+1:EstimOpt.NVarA*(2+EstimOpt.NVarM+EstimOpt.NLatent)+EstimOpt.NVarS))];
    
	else
        Results.DetailsScale = [];
    end
    Results.DetailsS(:,1) = Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)+EstimOpt.NVarS+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA+EstimOpt.NVarS);
	Results.DetailsS(:,3:4) = [Results.std(EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)+EstimOpt.NVarS+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA+EstimOpt.NVarS), pv(Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)+EstimOpt.NVarS+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA+EstimOpt.NVarS), Results.std(EstimOpt.NVarA*(2+EstimOpt.NLatent+EstimOpt.NVarM)+EstimOpt.NVarS+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA+EstimOpt.NVarS))];
	Results.DetailsM(:,1) = Results.bhat((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA+EstimOpt.NVarS+1:end);
	Results.DetailsM(:,3:4) = [Results.std((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA+EstimOpt.NVarS+1:end), pv(Results.bhat((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA+EstimOpt.NVarS+1:end), Results.std((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+(2+EstimOpt.NVarM)*EstimOpt.NVarA+EstimOpt.NVarS+1:end))];

elseif EstimOpt.FullCov == 1
    Results.DetailsA(:,1) = Results.bhat(1:EstimOpt.NVarA);
	Results.DetailsA(:,3:4) = [Results.std(1:EstimOpt.NVarA), pv(Results.bhat(1:EstimOpt.NVarA), Results.std(1:EstimOpt.NVarA))];
    Results.DetailsV = sdtri(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2), Results.ihess(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2,EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2),EstimOpt);
    Results.DetailsV = [Results.DetailsV(:,1),zeros(EstimOpt.NVarA,1),Results.DetailsV(:,2:3)];
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
		for i=1:EstimOpt.NVarM
	%         Results.DetailsCM(1:EstimOpt.NVarA,4*i-3) = Results.bhat(EstimOpt.NVarA*(EstimOpt.NVarA/2+0.5+i)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+i)); 
	%         Results.DetailsCM(1:EstimOpt.NVarA,4*i-1:4*i) = [Results.std(EstimOpt.NVarA*(EstimOpt.NVarA/2+0.5+i)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+i)),pv(Results.bhat(EstimOpt.NVarA*(EstimOpt.NVarA/2+0.5+i)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+i)),Results.std(EstimOpt.NVarA*(EstimOpt.NVarA/2+0.5+i)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+i)))]; 
			Results.DetailsCM(1:EstimOpt.NVarA,4*i-3) = Results.bhat(EstimOpt.NVarA*i+1:EstimOpt.NVarA*(i+1)); 
			Results.DetailsCM(1:EstimOpt.NVarA,4*i-1:4*i) = [Results.std(EstimOpt.NVarA*i+1:EstimOpt.NVarA*(i+1)),pv(Results.bhat(EstimOpt.NVarA*i+1:EstimOpt.NVarA*(i+1)),Results.std(EstimOpt.NVarA*i+1:EstimOpt.NVarA*(i+1)))]; 

		
		end
	else
		Results.DetailsCM = [];
	end

	
	if EstimOpt.NVarM > 0
        for i=1:EstimOpt.NVarM
			Results.DetailsCM(:,4*i-3) = Results.bhat(EstimOpt.NVarA*(EstimOpt.NVarA+3+2*(i-1))/2+1:EstimOpt.NVarA*(EstimOpt.NVarA+3+2*i)/2);
			Results.DetailsCM(:,4*i-1:4*i) = [Results.std(EstimOpt.NVarA*(EstimOpt.NVarA+3+2*(i-1))/2+1:EstimOpt.NVarA*(EstimOpt.NVarA+3+2*i)/2), pv(Results.bhat(EstimOpt.NVarA*(EstimOpt.NVarA+3+2*(i-1))/2+1:EstimOpt.NVarA*(EstimOpt.NVarA+3+2*i)/2), Results.std(EstimOpt.NVarA*(EstimOpt.NVarA+3+2*(i-1))/2+1:EstimOpt.NVarA*(EstimOpt.NVarA+3+2*i)/2))];
		end
	else
        Results.DetailsCM = [];
    end
    l = l + EstimOpt.NVarM;
    if EstimOpt.NVarS > 0
        Results.DetailsScale = Results.bhat(l+1:l+EstimOpt.NVarS);
		Results.DetailsScale = [Results.std(l+1:l+EstimOpt.NVarS), pv(Results.bhat(l+1:l+EstimOpt.NVarS), Results.std(l+1:l+EstimOpt.NVarS))];
	else
        Results.DetailsScale = [];
    end
    l = l + EstimOpt.NVarS;
    Results.DetailsL(:,1) = Results.bhat(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent);
    Results.DetailsL(:,3:4) = [Results.std(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent), pv(Results.bhat(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent), Results.std(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent))];
	Results.DetailsS(:,1) = Results.bhat(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent);
	Results.DetailsS(:,3:4) = [Results.std(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent), pv(Results.bhat(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent), Results.std(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent))];
    Results.DetailsM(:,1) = Results.bhat(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end);
	Results.DetailsM(:,3:4) = [Results.std(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end), pv(Results.bhat(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end), Results.std(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end))];
    
elseif EstimOpt.FullCov == 2
    Results.DetailsA(:,1) = Results.bhat(1:EstimOpt.NVarA);
	Results.DetailsA(:,3:4) = [Results.std(1:EstimOpt.NVarA), pv(Results.bhat(1:EstimOpt.NVarA), Results.std(1:EstimOpt.NVarA))];
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
        Results.DetailsCM(:,1) = Results.bhat(l+1:l+EstimOpt.NVarA*EstimOpt.NVarM); 
		Results.DetailsCM(:,3:4) = [Results.std(l+1:l+EstimOpt.NVarA*EstimOpt.NVarM), pv(Results.bhat(l+1:l+EstimOpt.NVarA*EstimOpt.NVarM), Results.std(l+1:l+EstimOpt.NVarA*EstimOpt.NVarM))];
	else
        Results.DetailsCM = [];
    end
    l = l + EstimOpt.NVarM;
    if EstimOpt.NVarS > 0
        Results.DetailsScale(:,1) = Results.bhat(l+1:l+EstimOpt.NVarS);
		Results.DetailsScale(:,3:4) = [Results.std(l+1:l+EstimOpt.NVarS), pv(Results.bhat(l+1:l+EstimOpt.NVarS), Results.std(l+1:l+EstimOpt.NVarS))];
	else
        Results.DetailsScale = [];
    end
    l = l + EstimOpt.NVarS;
    Results.DetailsL(:,1) = Results.bhat(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent);
    Results.DetailsL(:,3:4) = [Results.std(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent), pv(Results.bhat(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent), Results.std(l+1:l+EstimOpt.NVarA*EstimOpt.NLatent))];
	Results.DetailsS(:,1) = Results.bhat(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent);
    Results.DetailsS(:,3:4) = [Results.std(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent), pv(Results.bhat(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent), Results.std(l+EstimOpt.NVarA*EstimOpt.NLatent+1:l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent))];
	Results.DetailsM(:,1) = Results.bhat(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end);
    Results.DetailsM(:,3:4) = [Results.std(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end), pv(Results.bhat(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end), Results.std(l+(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+1:end))];
	if any(EstimOpt.Dist == -1)
        VC = tril(ones(EstimOpt.NLatent+EstimOpt.NVarA));
        VC(EstimOpt.Dist == -1,:) = 0;
        VC(:,EstimOpt.Dist == -1) = 0;
        VCtmp = tril(ones(size(VC)));
        VCtmp2 = diag(ones(EstimOpt.NLatent+EstimOpt.NVarA,1));
        VCtmp2(1:EstimOpt.NVarA, 1:EstimOpt.NVarA) = 0;
        VCtmp(VCtmp2 == 1) = 0;
        Indx = VC(VCtmp == 1);
        bhattmp = Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent);
        bhattmp = bhattmp(Indx == 1);
        covtmp =  Results.ihess(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent,EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent);
        covtmp = covtmp(Indx == 1,Indx==1);
        %         [Results.DetailsCOV, Results.DetailsCORR] = sdtriHe(bhattmp, covtmp,EstimOpt);
        
        % Using delta method instead simulation
        H = jacobianest(@(b) sdtriHe2(b, EstimOpt, 0), bhattmp); % Covariance
        VarTmp = H*covtmp*H';
        Results.DetailsCOV = [sdtriHe2(bhattmp, EstimOpt, 0), zeros(EstimOpt.NVarA,1), sqrt(diag(VarTmp)), pv(sdtriHe2(bhattmp, EstimOpt, 0), sqrt(diag(VarTmp)))];
        H = jacobianest(@(b) sdtriHe2(b, EstimOpt, 1), bhattmp); % Correlation
        VarTmp = H*covtmp*H';
        Results.DetailsCORR = [sdtriHe2(bhattmp, EstimOpt, 1), zeros(EstimOpt.NVarA,1), sqrt(diag(VarTmp)), pv(sdtriHe2(bhattmp, EstimOpt, 1), sqrt(diag(VarTmp)))];
        
        H = jacobianest(@(b) sdtriHe2(b, EstimOpt, 2), bhattmp); % standard deviations
        VarTmp = H*covtmp*H';
        Results.DetailsV = zeros(EstimOpt.NVarA,4);
        Results.DetailsV(EstimOpt.Dist ~= -1,:) = [sdtriHe2(bhattmp,EstimOpt,2), zeros(EstimOpt.NVarA,1), sqrt(diag(VarTmp)), pv(sdtriHe2(bhattmp, EstimOpt, 2), sqrt(diag(VarTmp)))];
    else
        
        bhattmp = Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent);
        covtmp =  Results.ihess(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent,EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent);
        % Using delta method instead simulation
        H = jacobianest(@(b) sdtriHe2(b, EstimOpt, 0), bhattmp); % Covariance
        VarTmp = H*covtmp*H';
        Results.DetailsCOV = [sdtriHe2(bhattmp, EstimOpt, 0), zeros(EstimOpt.NVarA,1), sqrt(diag(VarTmp)), pv(sdtriHe2(bhattmp, EstimOpt, 0), sqrt(diag(VarTmp)))];
        H = jacobianest(@(b) sdtriHe2(b, EstimOpt, 1), bhattmp); % Correlation
        VarTmp = H*covtmp*H';
        Results.DetailsCORR = [sdtriHe2(bhattmp, EstimOpt, 1), zeros(EstimOpt.NVarA,1), sqrt(diag(VarTmp)), pv(sdtriHe2(bhattmp, EstimOpt, 1), sqrt(diag(VarTmp)))];
        %[Results.DetailsCOV, Results.DetailsCORR] = sdtriHe(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent), Results.ihess(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent,EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA+EstimOpt.NLatent)-EstimOpt.NLatent),EstimOpt);
        
        H = jacobianest(@(b) sdtriHe2(b, EstimOpt, 2), bhattmp); % standard deviations
        VarTmp = H*covtmp*H';
        Results.DetailsV = [sdtriHe2(bhattmp,EstimOpt,2), zeros(EstimOpt.NVarA,1), sqrt(diag(VarTmp)), pv(sdtriHe2(bhattmp, EstimOpt, 2), sqrt(diag(VarTmp)))];
        
        
    end
end



Results.R = [Results.DetailsA; Results.DetailsV ; Results.DetailsCM ; Results.DetailsL ; Results.DetailsS;Results.DetailsM];
Results.LL0 = Results.MIMIC0.LL + Results_old.MNL0.LL;
EstimOpt.params = length(Results.bhat);
if isfield(EstimOpt,'BActive')
    EstimOpt.params = EstimOpt.params - sum(EstimOpt.BActive == 0);
end

Results.stats = [Results.LL; Results.LL0;  1-Results.LL/Results.LL0;R2; ((2*EstimOpt.params-2*Results.LL))/EstimOpt.NObs; ((log(EstimOpt.NObs)*EstimOpt.params-2*Results.LL))/EstimOpt.NObs ;EstimOpt.NObs; EstimOpt.NP; EstimOpt.params];

Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;
Results.INPUT = INPUT;

Results.R_out  = cell(3+EstimOpt.NVarStr +EstimOpt.NVarA+3+ 3*size(INPUT.Xmea,2)+EstimOpt.NVarcut+EstimOpt.NVarMea+ 2 + 7+(2 + EstimOpt.NVarS)*(EstimOpt.NVarS>0), 7+3*(EstimOpt.NLatent+EstimOpt.NVarM));


Head = cell(1,2);

if EstimOpt.FullCov == 0
    Head(1,1) = {'HMXL_d'};
else
    Head(1,1) = {'HMXL'};
end

if EstimOpt.WTP_space > 0
    Head(1,2) = {'in WTP-space'};
else
    Head(1,2) = {'in preference-space'};
end

Template1 = {'DetailsA', 'DetailsV'};
Template2 = {'DetailsA', 'DetailsV'};
Names.DetailsA = EstimOpt.NamesA;
Heads.DetailsA = {'Means'};
Heads.DetailsV = {'Standard Deviations'};
ST = [];
LVlist = {'LV 1' , 'LV 2' , 'LV 3' , 'LV 4' , 'LV 5' , 'LV 6' , 'LV 7' , 'LV 8' ,'LV 9' ,'LV 10'};
LVlist2 = {'LV1' , 'LV2' , 'LV3' , 'LV4' , 'LV5' , 'LV6' , 'LV7' , 'LV8' ,'LV9' ,'LV10'};
% head = {'var.' , 'coef.', 'st.err.' , 'p-value'};
% headx = [head, repmat(head(1,2:4),1,1+EstimOpt.NLatent+EstimOpt.NVarM)];
Temp = [];
for i = 1:EstimOpt.NLatent
    Results.(LVlist2{i}) = Results.DetailsL(1+(i-1)*EstimOpt.NVarA:i*EstimOpt.NVarA,:);
    Heads.(LVlist2{i}) = LVlist(i);
    Template1 = [Template1, LVlist2(i)];
    Temp = [Temp, LVlist2(i)];
    %Results.R_out(4:(EstimOpt.NVarA+3),5+(i-1)*3:4+i*3) = num2cell(Results.DetailsL(1+(i-1)*EstimOpt.NVarA:i*EstimOpt.NVarA,:));
end
Template2 = cell(2,max(size(Temp,2),2));
Template2(1,1) = {'DetailsA'};
Template2(1,2) = {'DetailsV'};
Template2(end,:) = Temp;

% for i = 1:EstimOpt.NLatent
    % Results.R_out(4:(EstimOpt.NVarA+3),8+(i-1)*3:7+i*3) = num2cell(Results.DetailsL(1+(i-1)*EstimOpt.NVarA:i*EstimOpt.NVarA,:));
% end
if EstimOpt.NVarM > 0
    Template1 = [Template1, 'DetailsCM'];
	Temp2 = cell(1, size(Template2,2));
    Temp2(1,1) = {'DetailsCM'};
    Template2 = [Template2; Temp2];
	Heads.DetailsCM = EstimOpt.NamesM;
end


if EstimOpt.NVarS >0
    Temp1 = cell(1, size(Template1,2));
    Temp1(1,1) = {'DetailsScale'};
    Template1 = [Template1; Temp1];
    Temp2 = cell(1, size(Template2,2));
    Temp2(1,1) = {'DetailsScale'};
    Template2 = [Template2; Temp2];
    Names.DetailsScale = EstimOpt.NamesS;
    Heads.DetailsScale = {'Covariates of Scale'};
    ST = {'DetailsScale'};
end


for i = 1: EstimOpt.NLatent 
    if EstimOpt.NVarStr > 0
        Results.SE(1:EstimOpt.NVarStr,4*i-3:4*i) = Results.DetailsS((i-1)*EstimOpt.NVarStr+1:i*EstimOpt.NVarStr,:);
        Heads.SE{i,1} = 'Structural equations';
        Heads.SE(i,2) = LVlist(i);
    end

%Results.R_out(EstimOpt.NVarA+7+(2 + EstimOpt.NVarS)*(EstimOpt.NVarS>0):EstimOpt.NVarA+6+EstimOpt.NVarStr+(2 + EstimOpt.NVarS)*(EstimOpt.NVarS>0),2 + 3*(i-1):1+3*i) = num2cell(Results.DetailsS((i-1)*EstimOpt.NVarStr+1:i*EstimOpt.NVarStr,:));
end

if isfield(Results, 'SE')
        Temp1 = cell(1, size(Template1,2));
        Temp1(1,1) = {'SE'};
        Template1 = [Template1; Temp1];
        Temp2 = cell(1, size(Template2,2));
        Temp2(1,1) = {'SE'};
        Template2 = [Template2; Temp2];
        Names.SE = EstimOpt.NamesStr;
        ST = [ST, {'SE'}];
    end
%l = EstimOpt.NVarA+3+(2 + EstimOpt.NVarS)*(EstimOpt.NVarS>0); % this is for indexing in R_out
k = 0;
for i = 1:size(INPUT.Xmea,2)
    Heads.(strcat('Xmea',num2str(i))){1,1} = ['Measurment equation for',' ',char(EstimOpt.NamesMea(i))];
    model = model_name(EstimOpt.MeaSpecMatrix(i));
    Heads.(strcat('Xmea',num2str(i))){1,2} = ['Estimated using',' ', model];
   
    Results.(strcat('Xmea',num2str(i)))(1:EstimOpt.CutMatrix(i),1:4) = Results.DetailsM(k+1:k+EstimOpt.CutMatrix(i),:); 
    Names.(strcat('Xmea',num2str(i))) = EstimOpt.NamesLV(k+1:k+EstimOpt.CutMatrix(i));
	
	%l = l+3+EstimOpt.CutMatrix(i);
    k = k + EstimOpt.CutMatrix(i);
	
	Temp1 = cell(1, size(Template1,2));
    Temp1(1,1) = {strcat('Xmea',num2str(i))};
    Template1 = [Template1; Temp1];
    Temp2 = cell(1, size(Template2,2));
    Temp2(1,1) = {strcat('Xmea',num2str(i))};
    Template2 = [Template2; Temp2];
    ST = [ST, {strcat('Xmea',num2str(i))}];
end
%% Tworzenie stopki
Tail = cell(17,2);
Tail(2,1) = {'Model diagnostics'};
Tail(3:17,1) = {'LL at convergence' ;'LL at constant(s) only';  strcat('McFadden''s pseudo-R',char(178));strcat('Ben-Akiva-Lerman''s pseudo-R',char(178))  ;'AIC/n' ;'BIC/n'; 'n (observations)'; 'r (respondents)';'k (parameters)';' ';'Estimation method';'Simulation with';'Optimization method';'Gradient';'Hessian'};

if isfield(Results_old,'MNL0') && isfield(Results_old.MNL0,'LL')
    Tail(3:11,2) = num2cell(Results.stats);
end

if any(INPUT.W ~= 1)
    Tail(13,2) = {'weighted'};
else
    Tail(13,2) = {'maximum likelihood'};
end

switch EstimOpt.Draws
     case 1
     Tail(14,2) = {[num2str(EstimOpt.NRep),' ','pseudo-random draws']};
     case 2
     Tail(14,2) = {[num2str(EstimOpt.NRep),' ','Latin Hypercube Sampling draws']};
     case  3
     Tail(14,2) = {[num2str(EstimOpt.NRep),' ','Halton draws (skip = ', num2str(EstimOpt.HaltonSkip), '; leap = ', num2str(EstimOpt.HaltonLeap),')']};
     case 4 
     Tail(14,2) = {[num2str(EstimOpt.NRep),' ','Halton draws with reverse radix scrambling (skip = ', num2str(EstimOpt.HaltonSkip), '; leap = ', num2str(EstimOpt.HaltonLeap),')']};
     case 5
     Tail(14,2) = {[num2str(EstimOpt.NRep),' ','Sobol draws (skip = ', num2str(EstimOpt.HaltonSkip), '; leap = ', num2str(EstimOpt.HaltonLeap),')']};
     case 6
     Tail(14,2) = {[num2str(EstimOpt.NRep),' ','Sobol draws with random linear scramble and random digital shift (skip = ', num2str(EstimOpt.HaltonSkip), '; leap = ', num2str(EstimOpt.HaltonLeap),')']};    
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

outHessian = [];
if isequal(OptimOpt.Algorithm,'quasi-newton')
    outHessian='off, ';
    switch EstimOpt.HessEstFix
        case 0
            outHessian = [outHessian, 'retained from optimization'];
        case 1
            outHessian = [outHessian, 'ex-post calculated using BHHH'];
        case 2
            outHessian = [outHessian, 'ex-post calculated using high-precision BHHH'];
        case 3
            outHessian = [outHessian, 'ex-post calculated numerically'];
        case 4
            outHessian = [outHessian, 'ex-post calculated analytically'];
    end
else
    if strcmp(OptimOpt.Hessian,'user-supplied')
        if EstimOpt.ApproxHess == 1
            outHessian = 'user-supplied, BHHH, ';
        else
            outHessian = 'user-supplied, analytical, ';
        end
    else
        outHessian = ['built-in, ', num2str(OptimOpt.HessUpdate), ', '];
    end
    switch EstimOpt.HessEstFix
        case 0
            outHessian = [outHessian, 'retained from optimization'];
        case 1
            outHessian = [outHessian, 'ex-post calculated using BHHH'];
        case 2
            outHessian = [outHessian, 'ex-post calculated using high-precision BHHH'];
        case 3
            outHessian = [outHessian, 'ex-post calculated numerically'];
        case 4
            outHessian = [outHessian, 'ex-post calculated analytically'];
    end
end

Tail(17,2) = {outHessian};
%% Tworzenie ResultsOut, drukowanie na ekran i do pliku .xls
if EstimOpt.Display~=0

    Results.Dist = EstimOpt.MNLDist;
    Results.R_out = genOutput(EstimOpt, Results, Head, Tail, Names, Template1, Template2, Heads, ST);
    fullOrgTemplate = which('template.xls');   
    currFld = pwd;
    
	if EstimOpt.FullCov == 0
		if isfield(EstimOpt,'ProjectName')
			fullSaveName = strcat(currFld,'\HMXL_d_results_',EstimOpt.ProjectName,'.xls');
		else
			fullSaveName = strcat(currFld,'\HMXL_d_results.xls');
		end
    else
		if isfield(EstimOpt,'ProjectName')
			fullSaveName = strcat(currFld,'\HMXL_results_',EstimOpt.ProjectName,'.xls');
		else
			fullSaveName = strcat(currFld,'\HMXL_results.xls');
		end
	end
    copyfile(fullOrgTemplate, 'templateTMP.xls')
    fullTMPTemplate = which('templateTMP.xls');

    excel = actxserver('Excel.Application');
    excelWorkbook = excel.Workbooks.Open(fullTMPTemplate);
    excel.Visible = 1;
    excel.DisplayAlerts = 0;
    excelSheets = excel.ActiveWorkbook.Sheets;
    excelSheet1 = excelSheets.get('Item',1);
    excelSheet1.Activate;
    column = size(Results.R_out,2);
    columnName = [];    
        while column > 0
            modulo = mod(column - 1,26);
            columnName = [char(65 + modulo) , columnName];
            column = floor(((column - modulo) / 26));
        end

    rangeE = strcat('A1:',columnName,num2str(size(Results.R_out,1)));
    excelActivesheetRange = get(excel.Activesheet,'Range',rangeE);
    excelActivesheetRange.Value = Results.R_out;
    if isfield(EstimOpt,'xlsOverwrite') && EstimOpt.xlsOverwrite == 0
        i = 1;
        while exist(fullSaveName, 'file') == 2
            if isempty(strfind(fullSaveName, '('))
                pos = strfind(fullSaveName, '.xls');
                fullSaveName = strcat(fullSaveName(1:pos-1),'(',num2str(i),').xls');
            else
                pos = strfind(fullSaveName, '(');
                fullSaveName = strcat(fullSaveName(1:pos),num2str(i),').xls');
            end
            i = i+1;
        end
    end
    excelWorkbook.ConflictResolution = 2;
    SaveAs(excelWorkbook,fullSaveName);
    excel.DisplayAlerts = 0;
    excelWorkbook.Saved = 1;
    Close(excelWorkbook)
    Quit(excel)
    delete(excel)
    delete(fullTMPTemplate)
end