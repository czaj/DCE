function Results = HMNL(INPUT,Results_old,EstimOpt,OptimOpt)


global B_backup;

% save tmp_HMNL
% return

tic

Results.bhat = [];
Results.R = [];
Results.R_out = {};
Results.stats = [];


%% Check data and inputs


if nargin < 3 % check no. of inputs
    error('Too few input arguments for HMNL(INPUT,EstimOpt,OptimOpt)')
end

disp(' ');
disp('__________________________________________________________________________________________________________________');
disp(' ');

warning off MATLAB:mir_warning_maybe_uninitialized_temporary

format shortG;
format compact;
if isfield(EstimOpt,'NLatent') == 0
    EstimOpt.NLatent = 1;
    disp(['Assuming ', num2str(EstimOpt.NLatent)',' Latent Variable(s)']);
end
if any(INPUT.W ~= 1)
    cprintf('Black','Estimating '); cprintf('*Black','weighted '); cprintf('Black','HMNL with '); cprintf('Black',num2str(EstimOpt.NLatent)); cprintf('Black','  Latent Variable(s) ...\n');
else
    disp(num2str(EstimOpt.NLatent,'Estimating HMNL with %1.0f Latent Variable(s)'))
end

if isfield(EstimOpt, 'WTP_space') == 0
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

if (isfield(EstimOpt,'MeaSpecMatrix') == 0 || numel(EstimOpt.MeaSpecMatrix) == 0)
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

if isfield(INPUT, 'Xm') == 0
    INPUT.Xm = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarM = size(INPUT.Xm,2); % Number of covariates of means of random parameters
if isfield(INPUT, 'Xs') == 0
    INPUT.Xs = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarS = size(INPUT.Xs,2); % Number of covariates of scale


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

for i=1:EstimOpt.NVarStr
    if sum(isnan(INPUT.Xstr(INPUT.MissingInd==0,i))) > 0
        cprintf(rgb('DarkOrange'), 'WARNING: Structural variable %d contains NaN values \n', i)
    end
    if sum(isinf(INPUT.Xstr(INPUT.MissingInd==0,i))) > 0
        cprintf(rgb('DarkOrange'), 'WARNING: Structural variable %d contains Inf values \n', i)
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

if isfield(EstimOpt,'MNLDist') && any(EstimOpt.MNLDist ~= 0)
    EstimOpt.MNLDist = EstimOpt.MNLDist(:);
    if size(EstimOpt.MNLDist,1) ~= EstimOpt.NVarA
        error('Incorrect no. of random parameters'' distributions provided')
    end
else
    EstimOpt.MNLDist = zeros(EstimOpt.NVarA,1); % no transformations
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
    elseif EstimOpt.MeaSpecMatrix(i) == 0 % OLS
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

if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == (EstimOpt.NVarA*(1+EstimOpt.NLatent+ EstimOpt.NVarM) + EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut+EstimOpt.NVarS)
    b0 = B_backup(:);
    disp('Using the starting values from Backup')
elseif isfield(Results_old,'HMNL') && isfield(Results_old.HMNL,'b0') % starting values provided
    Results_old.HMNL.b0_old = Results_old.HMNL.b0(:);
    Results_old.HMNL = rmfield(Results_old.HMNL,'b0');
    if length(Results_old.HMNL.b0_old) ~= (EstimOpt.NVarA*(1+EstimOpt.NLatent+ EstimOpt.NVarM) + EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut+EstimOpt.NVarS)
        cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
        Results_old.HMNL = rmfield(Results_old.HMNL,'b0_old');
    else
        b0 = Results_old.HMNL.b0_old(:);
    end
end
if  ~exist('b0','var')
    if isfield(Results_old,'MIMIC') && isfield(Results_old.MIMIC,'bhat')
        disp('Using MIMIC results as starting values')
        Results_old.MIMIC.bhat = Results_old.MIMIC.bhat(:);
        bstr = reshape(Results_old.MIMIC.bhat(1:EstimOpt.NVarStr*EstimOpt.NLatent), EstimOpt.NVarStr, EstimOpt.NLatent);
        LV = INPUT.Xstr*bstr;
        mLV = mean(LV,1);
        sLV = std(LV); ...
            LV = (LV - mLV(ones(1,size(LV,1)),:))./sLV(ones(1,size(LV,1)),:); ... % normalilzing for 0 mean and std
            LV = LV(:,:, ones(EstimOpt.NCT*EstimOpt.NAlt,1)); % NP x NLatent x NCT*NAlt
        LV = permute(LV, [3 1 2]);
        LV = reshape(LV, EstimOpt.NCT*EstimOpt.NAlt*EstimOpt.NP,EstimOpt.NLatent);
        disp('Estimating additional MNL for starting values')
        INPUT0.Xa = INPUT.Xa;
        %             for i = 1:EstimOpt.NLatent
        %                 LVi = LV(:,i);
        %                INPUT0.Xa = [INPUT0.Xa, INPUT.Xa.*LVi(:,ones(1,EstimOpt.NVarA))];
        %             end
        INPUT0.Y = INPUT.Y;
        INPUT0.Xm = [INPUT.Xm, LV];
        INPUT0.W = ones(EstimOpt.NP,1);
        INPUT0.MissingInd = INPUT.MissingInd;
        INPUT0.Xs = INPUT.Xs;
        if EstimOpt.ScaleLV == 1
            INPUT0.Xs = [INPUT0.Xs, LV];
        end
        EstimOpt0 = EstimOpt;
        EstimOpt0.DISPLAY = 0;
        EstimOpt0.WTP_space = 0;
        %            EstimOpt0.NVarA = size(INPUT0.Xa,2);
        if isfield(EstimOpt, 'BActive') == 1
            EstimOpt0.BActive =  EstimOpt.BActive(1:EstimOpt.NVarA*(1+EstimOpt.NLatent));
        else
            EstimOpt0.BActive =  ones(1,EstimOpt.NVarA*(1+EstimOpt.NLatent));
        end
        Results.MNL_LV = MNL(INPUT0,[],EstimOpt0,OptimOpt);
        %             if EstimOpt.WTP_space > 0
        %                 Results.MNL_LV.bhat(1:(EstimOpt.NVarA-EstimOpt.WTP_space)) = Results.MNL_LV.bhat(1:(EstimOpt.NVarA-EstimOpt.WTP_space)).*Results.MNL_LV.bhat(EstimOpt.WTP_matrix);
        %             end
        b0 = [Results.MNL_LV.bhat; Results_old.MIMIC.bhat];
        
    elseif isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat')
        disp('Using MNL and MIMIC0 results as starting values')
        Results_old.MNL.bhat = Results_old.MNL.bhat(:);
        b0 = zeros(EstimOpt.NVarA*(1+EstimOpt.NLatent+EstimOpt.NVarM) + EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut+EstimOpt.NVarS,1);
        b0(1:EstimOpt.NVarA*(1+EstimOpt.NVarM)) = Results_old.MNL.bhat(1:EstimOpt.NVarA*(1+EstimOpt.NVarM));
        if EstimOpt.ScaleLV == 0
            b0(EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)+1:EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)+EstimOpt.NVarS) = Results_old.MNL.bhat(EstimOpt.NVarA*(1+EstimOpt.NVarM)+1:EstimOpt.NVarA*(1+EstimOpt.NVarM)+EstimOpt.NVarS);
        else
            b0(EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)+1:EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)+EstimOpt.NVarS-EstimOpt.NLatent) = Results_old.MNL.bhat(EstimOpt.NVarA*(1+EstimOpt.NVarM)+1:EstimOpt.NVarA*(1+EstimOpt.NVarM)+EstimOpt.NVarS-EstimOpt.NLatent);
        end
        l = EstimOpt.NVarStr*EstimOpt.NLatent + (EstimOpt.NLatent + 1+EstimOpt.NVarM)*EstimOpt.NVarA+EstimOpt.NVarS;
        k = 0;
        for i = 1:size(INPUT.Xmea,2)
            if EstimOpt.MeaSpecMatrix(i) == 2 %Ordered probit: cutoffs
                b0(l+ sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp + 1:l+ sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp + length(unique(INPUT.Xmea(:,i))) -1) = Results.MIMIC0.bhat(k+1:k+length(unique(INPUT.Xmea(:,i)))-1);
                l = l + sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp + length(unique(INPUT.Xmea(:,i)))-1;
                k = k+ length(unique(INPUT.Xmea(:,i)))-1;
            elseif EstimOpt.MeaSpecMatrix(i) == 0
                b0(l + 1) = Results.MIMIC0.bhat(k+1);
                b0(l+2+sum(EstimOpt.MeaMatrix(:,i)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp,1)) = Results.MIMIC0.bhat(k+2);
                k = k+2;
                l = l+2 + sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp;
            elseif EstimOpt.MeaSpecMatrix(i) == 1 % MNL
                for n = 1:(length(unique(INPUT.Xmea(:,i)))-1)
                    b0(l + 1) = Results.MIMIC0.bhat(k+1);
                    l = l+ 1+ sum(EstimOpt.MeaMatrix(:,i))+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp;
                    k = k+ 1;
                end
            elseif EstimOpt.MeaSpecMatrix(i) == 3 % POISS
                b0(l + 1) = Results.MIMIC0.bhat(k+1);
                k = k+1;
                l = l+1 + sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp;
            elseif EstimOpt.MeaSpecMatrix(i) == 4 % NB
                b0(l + 1) = Results.MIMIC0.bhat(k+1);
                b0(l + 2+ sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp) = Results.MIMIC0.bhat(k+2); % theta
                k = k+2;
                l = l+2 + sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp;
            elseif EstimOpt.MeaSpecMatrix(i) == 5 % ZIP
                b0(l + 1) = Results.MIMIC0.bhat(k+1);
                b0(l + 2+ sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp) = Results.MIMIC0.bhat(k+2); % second constant
                k = k+2;
                l = l+2 + 2*sum(EstimOpt.MeaMatrix(:,i),1)+2*EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp;
            elseif EstimOpt.MeaSpecMatrix(i) == 6 % ZINB
                b0(l + 1) = Results.MIMIC0.bhat(k+1);
                b0(l + 2+ sum(EstimOpt.MeaMatrix(:,i),1)+EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp) = Results.MIMIC0.bhat(k+2); % Second constant
                b0(l + 3+ 2*sum(EstimOpt.MeaMatrix(:,i),1)+2*EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp) = Results.MIMIC0.bhat(k+3); % Theta
                k = k+3;
                l = l+3 + 2*sum(EstimOpt.MeaMatrix(:,i),1)+2*EstimOpt.MeaExpMatrix(i)*EstimOpt.NVarMeaExp;
            end
        end
    else
        error('No starting values available - run MNL or MIMIC first')
    end
end
% it must be after starting values
INPUT.Xm = INPUT.Xm(1:EstimOpt.NCT*EstimOpt.NAlt:end,:);


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
            EstimOpt.BActive = ones(1,EstimOpt.NVarA*(1+EstimOpt.NLatent+ EstimOpt.NVarM) + EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut);
        end
        EstimOpt.BActive = StrSelect(EstimOpt,1);
    else
        error('Structural variables matrix (EstimOpt.StrMatrix) erroneously defined')
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
    %rng(EstimOpt.Seed1);
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

if any(EstimOpt.MNLDist == 2) && EstimOpt.NumGrad == 0
    EstimOpt.NumGrad = 1;
    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - spike distribution not supported by analytical gradient \n')
end

if any(EstimOpt.MeaSpecMatrix >= 3) && EstimOpt.NumGrad == 0 && any(any(INPUT.Xmea(:, EstimOpt.MeaSpecMatrix >=3) > 100))
    cprintf(rgb('DarkOrange'), 'WARNING: it is recommended to switch to numerical gradient, as analitycal can be not precise when Xmea take large values for NB \n')
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

LLfun = @(B) LL_hmnl_MATlike(INPUT.YY,INPUT.XXa,INPUT.Xm,INPUT.Xs,INPUT.Xstr,INPUT.Xmea,INPUT.Xmea_exp, err_sliced,INPUT.W,EstimOpt,OptimOpt,B);
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

% save tmp_HMNL

Results.LL = -LL;

%R2 = R2_hybrid(INPUT.YY,INPUT.XXa,INPUT.Xm,INPUT.Xstr,[],INPUT.MissingInd, err_sliced,EstimOpt,Results.bhat, 0);
R2 = R2_hybrid(INPUT.YY,INPUT.XXa,INPUT.Xstr,[],INPUT.Xm,INPUT.Xs, INPUT.MissingInd, err_sliced,EstimOpt,Results.bhat, 0);
Results.b0_old = b0;

if EstimOpt.HessEstFix == 1
    f = LL_hmnl(INPUT.YY,INPUT.XXa,INPUT.Xm,INPUT.Xs, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,Results.bhat);...
        Results.jacobian = numdiff(@(B) INPUT.W.*LL_hmnl(INPUT.YY,INPUT.XXa,INPUT.Xm ,INPUT.Xs, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,B), INPUT.W.*f, Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) INPUT.W.*LL_hmnl(INPUT.YY,INPUT.XXa,INPUT.Xm,INPUT.Xs, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.hess = hessian(@(B) sum(INPUT.W.*LL_hmnl(INPUT.YY,INPUT.XXa,INPUT.Xm,INPUT.Xs, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,B),1), Results.bhat);
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
        [~, Results.jacobian] = LL_hmnl(INPUT.YY,INPUT.XXa,INPUT.Xm,INPUT.Xs, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,Results.bhat);
        Results.jacobian = Results.jacobian.*INPUT.W(:, ones(1,size(Results.jacobian,2)));
    else
        Results.LLdetailed = LL_hmnl(INPUT.YY,INPUT.XXa,INPUT.Xm,INPUT.Xs, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, EstimOpt,Results.bhat);
        Results.jacobian = numdiff(@(B) INPUT.W.*LL_hmnl(INPUT.YY,INPUT.XXa,INPUT.Xm,INPUT.Xs, INPUT.Xstr, INPUT.Xmea,INPUT.Xmea_exp, err_sliced, B) ,INPUT.W.*Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
    end
    RobustHess = Results.jacobian'*Results.jacobian;
    Results.ihess = Results.ihess*RobustHess*Results.ihess;
end

Results.std = sqrt(diag(Results.ihess));
Results.std(EstimOpt.BActive == 0) = NaN;
Results.std(EstimOpt.BLimit == 1) = 0;
Results.std(imag(Results.std) ~= 0) = NaN;

%%1
Results.DetailsA(1:EstimOpt.NVarA,1) = Results.bhat(1:EstimOpt.NVarA);
Results.DetailsA(1:EstimOpt.NVarA,3:4) = [Results.std(1:EstimOpt.NVarA), pv(Results.bhat(1:EstimOpt.NVarA), Results.std(1:EstimOpt.NVarA))];
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


Results.DetailsL(1:EstimOpt.NVarA*EstimOpt.NLatent,1) = Results.bhat(EstimOpt.NVarA*(1+EstimOpt.NVarM)+1:EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent));
Results.DetailsL(1:EstimOpt.NVarA*EstimOpt.NLatent,3:4) = [Results.std(EstimOpt.NVarA*(1+EstimOpt.NVarM)+1:EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)), pv(Results.bhat(EstimOpt.NVarA*(1+EstimOpt.NVarM)+1:EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)), Results.std(EstimOpt.NVarA*(1+EstimOpt.NVarM)+1:EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)))];


if EstimOpt.NVarS > 0
    Results.DetailsScale(1:EstimOpt.NVarS,1) = Results.bhat(EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)+1:EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)+EstimOpt.NVarS);
    Results.DetailsScale(1:EstimOpt.NVarS,3:4) = [Results.std(EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)+1:EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)+EstimOpt.NVarS), pv(Results.bhat(EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)+1:EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)+EstimOpt.NVarS), Results.std(EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)+1:EstimOpt.NVarA*(1+EstimOpt.NVarM+EstimOpt.NLatent)+EstimOpt.NVarS))];
    
end

if EstimOpt.NVarStr > 0
    Results.DetailsS(1:(EstimOpt.NVarStr)*EstimOpt.NLatent,1) = Results.bhat(EstimOpt.NVarA*(1+EstimOpt.NLatent+EstimOpt.NVarM)+EstimOpt.NVarS+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+EstimOpt.NVarA*(1+EstimOpt.NVarM)+EstimOpt.NVarS);
    Results.DetailsS(1:(EstimOpt.NVarStr)*EstimOpt.NLatent,3:4) = [Results.std(EstimOpt.NVarA*(1+EstimOpt.NLatent+EstimOpt.NVarM)+EstimOpt.NVarS+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+EstimOpt.NVarA*(1+EstimOpt.NVarM)+EstimOpt.NVarS), pv(Results.bhat(EstimOpt.NVarA*(1+EstimOpt.NLatent+EstimOpt.NVarM)+EstimOpt.NVarS+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+EstimOpt.NVarA*(1+EstimOpt.NVarM)+EstimOpt.NVarS), Results.std(EstimOpt.NVarA*(1+EstimOpt.NLatent+EstimOpt.NVarM)+EstimOpt.NVarS+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+EstimOpt.NVarA*(1+EstimOpt.NVarM)+EstimOpt.NVarS))];
end

if sum(EstimOpt.CutMatrix) > 0
    Results.DetailsM(1:sum(EstimOpt.CutMatrix),1) = Results.bhat((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+EstimOpt.NVarA*(1+EstimOpt.NVarM)+EstimOpt.NVarS+1:end);
    Results.DetailsM(1:sum(EstimOpt.CutMatrix),3:4) = [Results.std((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+EstimOpt.NVarA*(1+EstimOpt.NVarM)+EstimOpt.NVarS+1:end), pv(Results.bhat((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+EstimOpt.NVarA*(1+EstimOpt.NVarM)+EstimOpt.NVarS+1:end), Results.std((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.NLatent+EstimOpt.NVarA*(1+EstimOpt.NVarM)+EstimOpt.NVarS+1:end))];
end






%Results.R = [Results.DetailsA; Results.DetailsL ; Results.DetailsS; Results.DetailsM];
Results.LL0 = Results.MIMIC0.LL + Results_old.MNL0.LL;
EstimOpt.params = length(Results.bhat);
if isfield(EstimOpt,'BActive')
    EstimOpt.params = EstimOpt.params - sum(EstimOpt.BActive == 0);
end

Results.stats = [Results.LL; Results.LL0;  1-Results.LL/Results.LL0;R2; ((2*EstimOpt.params-2*Results.LL))/EstimOpt.NObs; ((log(EstimOpt.NObs)*EstimOpt.params-2*Results.LL))/EstimOpt.NObs ;EstimOpt.NObs; EstimOpt.NP; EstimOpt.params];

Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;
Results.INPUT = INPUT;

%Results.R_out = cell(3+EstimOpt.NVarStr +EstimOpt.NVarA+3+ 3*size(INPUT.Xmea,2)+EstimOpt.NVarcut+EstimOpt.NVarMea+ 2 + 7 + (2 + EstimOpt.NVarS)*(EstimOpt.NVarS>0), 4+3*EstimOpt.NLatent);
%Results.R_out(1,1) = {'HMNL'};
%head = {'var.' , 'coef.', 'st.err.' , 'p-value'};
%headx = [head, repmat(head(1,2:4),1,EstimOpt.NLatent)];



Template1 = {'DetailsA'};

Names.DetailsA = EstimOpt.NamesA;
Heads.DetailsA = {'MNL'};


ST = [];
LVlist = {'LV 1' , 'LV 2' , 'LV 3' , 'LV 4' , 'LV 5' , 'LV 6' , 'LV 7' , 'LV 8' ,'LV 9' ,'LV 10'};
LVlist2 = {'LV1' , 'LV2' , 'LV3' , 'LV4' , 'LV5' , 'LV6' , 'LV7' , 'LV8' ,'LV9' ,'LV10'};
%if EstimOpt.NLatent<=10
%   Results.R_out(2,2:3:(2+EstimOpt.NLatent*3)) = [{'MNL'}, LVlist(1,1:EstimOpt.NLatent)];
%end
%Results.R_out(3,:) = headx;
%Results.R_out(4:(EstimOpt.NVarA+3),1:4) = [EstimOpt.NamesA,num2cell(Results.DetailsA)];
Temp = [];
for i = 1:EstimOpt.NLatent
    Results.(LVlist2{i}) = Results.DetailsL(1+(i-1)*EstimOpt.NVarA:i*EstimOpt.NVarA,:);
    Heads.(LVlist2{i}) = LVlist(i);
    Template1 = [Template1, LVlist2(i)];
    Temp = [Temp, LVlist2(i)];
    %Results.R_out(4:(EstimOpt.NVarA+3),5+(i-1)*3:4+i*3) = num2cell(Results.DetailsL(1+(i-1)*EstimOpt.NVarA:i*EstimOpt.NVarA,:));
end
Template2 = cell(2,size(Temp,2));
Template2(1,1) = {'DetailsA'};
Template2(end,:) = Temp;

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
    
    %     Results.R_out(EstimOpt.NVarA+4,1) = {'Covariates of scale'};
    %     Results.R_out(EstimOpt.NVarA+5,1:4) = head;
    %     Results.R_out(EstimOpt.NVarA+6:EstimOpt.NVarA+5+EstimOpt.NVarS,1) = EstimOpt.NamesS;
    %     Results.R_out(EstimOpt.NVarA+6:EstimOpt.NVarA+5+EstimOpt.NVarS,2:4) = num2cell(Results.DetailsScale);
end

%Results.R_out(EstimOpt.NVarA+5+(2 + EstimOpt.NVarS)*(EstimOpt.NVarS>0),2:3:(2+(EstimOpt.NLatent-1)*3)) = LVlist(1,1:EstimOpt.NLatent);
% Results.R_out(EstimOpt.NVarA+6+(2 + EstimOpt.NVarS)*(EstimOpt.NVarS>0),1:end-3) = headx(1:end-3);
% Results.R_out(EstimOpt.NVarA+7+(2 + EstimOpt.NVarS)*(EstimOpt.NVarS>0):EstimOpt.NVarA+6+EstimOpt.NVarStr+(2 + EstimOpt.NVarS)*(EstimOpt.NVarS>0),1) = EstimOpt.NamesStr;

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
    
    
    %Results.R_out(EstimOpt.NVarStr+3+l+1,1) =  cellfun(@(x)[x char(EstimOpt.NamesMea(i))],{'Measurment equation for '},'UniformOutput',0);
    Results.(strcat('Xmea',num2str(i)))(1:EstimOpt.CutMatrix(i),1:4) = Results.DetailsM(k+1:k+EstimOpt.CutMatrix(i),:);
    %Results.R_out(EstimOpt.NVarStr+3+l+3,1:4) = head;
    %Results.R_out(EstimOpt.NVarStr+3+l+4:EstimOpt.NVarStr+3+l+3 + EstimOpt.CutMatrix(i) ,1:4) = [EstimOpt.NamesLV(k+1:k+EstimOpt.CutMatrix(i)), num2cell(Results.DetailsM(k+1:k+EstimOpt.CutMatrix(i),:))];
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


%% Tworzenie naglowka
Head = cell(1,2);
Head(1,1) = {'HMNL'};
if EstimOpt.WTP_space > 0
    Head(1,2) = {'in WTP-space'};
else
    Head(1,2) = {'in preference-space'};
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
    EstimOpt.Dist = -ones(1,EstimOpt.NVarA+1);
    Results.Dist = transpose(EstimOpt.Dist(:,2:end));
    Results.R_out = genOutput(EstimOpt, Results, Head, Tail, Names, Template1, Template2, Heads, ST);
    fullOrgTemplate = which('template.xls');
    currFld = pwd;
    
    
    if isfield(EstimOpt,'ProjectName')
        fullSaveName = strcat(currFld,'\HMNL_results_',EstimOpt.ProjectName,'.xls');
    else
        fullSaveName = strcat(currFld,'\HMNL_results.xls');
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


% disp(' ')
% disp('MNL attributes');
% disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
% disp([char(EstimOpt.NamesA) ,blanks(EstimOpt.NVarA)',num2str(Results.DetailsA(:,1),'%8.4f'), star_sig(Results.DetailsA(:,3)), num2str(Results.DetailsA(:,2:3),'%8.4f %8.4f')])
% if EstimOpt.NVarM > 0
%     for i = 1:EstimOpt.NVarM
%         disp(' ');
%         disp(['Explanatory variable of random parameters'' means - ', char(EstimOpt.NamesM(i))]);
%         disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
%         disp([char(EstimOpt.NamesA),blanks(EstimOpt.NVarA)',num2str(Results.DetailsCM((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA,1),'%8.4f'), star_sig(Results.DetailsCM((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA,3)), num2str(Results.DetailsCM((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA,2:3),'%8.4f %8.4f')])
%     end
% end
%
% l = 0;
% for i = 1:EstimOpt.NLatent
%     disp(' ')
%     disp(num2str(i,'Interactions with Latent Variable %1.0f'));
%     disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
%     disp([char(EstimOpt.NamesA) ,blanks(EstimOpt.NVarA)',num2str(Results.DetailsL(l+1:l+EstimOpt.NVarA,1),'%8.4f'), star_sig(Results.DetailsL(l+1:l+EstimOpt.NVarA,3)), num2str(Results.DetailsL(l+1:l+EstimOpt.NVarA,2:3),'%8.4f %8.4f')])
%     l = l + EstimOpt.NVarA;
% end
%
% if EstimOpt.NVarS > 0
%     disp(' ')
%     disp('Covariates of scale')
%     disp(['var.',blanks(size(char(EstimOpt.NamesS),2)-2) ,'coef.      st.err.  p-value'])
%     disp([char(EstimOpt.NamesS) ,blanks(EstimOpt.NVarS)',num2str(Results.DetailsScale(:,1),'%8.4f'), star_sig(Results.DetailsScale(:,3)), num2str(Results.DetailsScale(:,2:3),'%8.4f %8.4f')])
%
% end
% l = 0;
% for i = 1:EstimOpt.NLatent
%     disp(' ')
%     disp(num2str(i,'Structural equation of Latent Variable %1.0f'));
%     disp(['var.', blanks(size(char(EstimOpt.NamesStr),2)-2) ,'coef.      st.err.  p-value'])
%     disp([char(EstimOpt.NamesStr) ,blanks(EstimOpt.NVarStr)',num2str(Results.DetailsS(l+1:l+EstimOpt.NVarStr,1),'%8.4f'), star_sig(Results.DetailsS(l+1:l+EstimOpt.NVarStr,3)), num2str(Results.DetailsS(l+1:l+EstimOpt.NVarStr,2:3),'%8.4f %8.4f')])
%     l = l + EstimOpt.NVarStr;
% end
%
% l = 0;
% for i = 1:size(INPUT.Xmea,2)
%     tmp = EstimOpt.NVarMeaExp*(EstimOpt.MeaExpMatrix(i) ~=0);
%     disp(' ')
%     disp([num2str('Measurment equation for '), char(EstimOpt.NamesMea(i))]);
%     if EstimOpt.MeaSpecMatrix(i) == 0
%         disp('Estimated using OLS')
%         disp('var.   coef.     st.err.  p-value')
%         Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1:3) = [exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)), Results.DetailsM(l+tmp+sum(EstimOpt.MeaMatrix(:,i))+2,2)*exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)),pv(exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)), Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,2)*exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)))];
%
%         disp([char(EstimOpt.NamesLV(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%7.4f %8.4f')])
%         l = l+sum(EstimOpt.MeaMatrix(:,i))+2+tmp;
%     elseif EstimOpt.MeaSpecMatrix(i) == 1
%         disp('Estimated using MNL')
%         for n = 1:(length(unique(INPUT.Xmea(:,i)))-1)
%             disp(num2str(n+1, 'Parameters for %1.0f alternative'))
%             disp('var.   coef.     st.err.  p-value')
%             disp([char(EstimOpt.NamesLV(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1), '%11.4f'), star_sig(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%8.4f %8.4f')])
%             l = l + 1 + sum(EstimOpt.MeaMatrix(:,i),1)+tmp;
%         end
%     elseif EstimOpt.MeaSpecMatrix(i) == 2
%         disp('Estimated using Ordered Probit')
%         disp('var.       coef.     st.err.  p-value')
%         disp([char(EstimOpt.NamesLV(l+1:l+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,repmat(blanks(sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',1,6),num2str(Results.DetailsM(l+1:l+sum(EstimOpt.MeaMatrix(:,i))+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+sum(EstimOpt.MeaMatrix(:,i))+tmp,3)), num2str(Results.DetailsM(l+1:l+sum(EstimOpt.MeaMatrix(:,i))+tmp,2:3),'%8.4f %8.4f')])
%         l = l+sum(EstimOpt.MeaMatrix(:,i))+tmp;
%         disp([num2str([1, Results.DetailsM(l+1,1)],'Cutoff %1.0f %7.4f'), star_sig(Results.DetailsM(l+1,3)), num2str(Results.DetailsM(l+1,2:3),'%8.4f %8.4f')])
%         if length(unique(INPUT.Xmea(:,i))) > 2 % if attitude is not binary
%             g = [Results.DetailsM(l+1,1) ; exp(Results.DetailsM(l+2:l+length(unique(INPUT.Xmea(:,i)))-1,1))];
%             for n = 2:length(unique(INPUT.Xmea(:,i)))-1
%                 stdx = sqrt(g(1:n)'*Results.ihess((EstimOpt.NVarStr)*EstimOpt.NLatent+l+1:(EstimOpt.NVarStr)*EstimOpt.NLatent+l+n,(EstimOpt.NVarStr)*EstimOpt.NLatent+ l+1:(EstimOpt.NVarStr)*EstimOpt.NLatent+l+n)*g(1:n));
%                 Results.DetailsM(l+n,1:3) = [sum(g(1:n),1), stdx, pv(sum(g(1:n),1), stdx)];
%                 disp([num2str([n, Results.DetailsM(l+n,1)],'Cutoff %1.0f %7.4f'), star_sig(Results.DetailsM(l+n,3)), num2str(Results.DetailsM(l+n,2:3),'%8.4f %8.4f')])
%             end
%         end
%         l = l+length(unique(INPUT.Xmea(:,i)))-1;
%     elseif EstimOpt.MeaSpecMatrix(i) == 3
%         disp('Estimated using Poisson regression')
%         disp('var.   coef.     st.err.  p-value')
%
%         disp([char(EstimOpt.NamesLV(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%7.4f %8.4f')])
%         l = l+sum(EstimOpt.MeaMatrix(:,i))+1+tmp;
%     elseif EstimOpt.MeaSpecMatrix(i) == 4
%         disp('Estimated using Negative Binomial regression')
%         disp('var.   coef.     st.err.  p-value')
%         Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1:3) = [exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)), Results.DetailsM(l+tmp+sum(EstimOpt.MeaMatrix(:,i))+2,2)*exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)),pv(exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)), Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,2)*exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)))];
%
%         disp([char(EstimOpt.NamesLV(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%7.4f %8.4f')])
%         l = l+sum(EstimOpt.MeaMatrix(:,i))+2+tmp;
%     elseif EstimOpt.MeaSpecMatrix(i) == 5
%         disp('Estimated using Zero Inflated Poisson regression')
%         disp('Probability of Non-participation (logit)')
%         disp('var.   coef.     st.err.  p-value')
%
%         disp([char(EstimOpt.NamesLV(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%7.4f %8.4f')])
%         l = l+sum(EstimOpt.MeaMatrix(:,i))+1+tmp;
%
%         disp('Poisson model')
%         disp('var.   coef.     st.err.  p-value')
%         disp([char(EstimOpt.NamesLV(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%7.4f %8.4f')])
%         l = l+sum(EstimOpt.MeaMatrix(:,i))+1+tmp;
%     elseif EstimOpt.MeaSpecMatrix(i) == 6
%         disp('Estimated using Zero Inflated Negative Binomial regression')
%         disp('Probability of Non-participation (logit)')
%         disp('var.   coef.     st.err.  p-value')
%         disp([char(EstimOpt.NamesLV(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+1+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%7.4f %8.4f')])
%         l = l+sum(EstimOpt.MeaMatrix(:,i))+1+tmp;
%
%         disp('Negative binomial model')
%         disp('var.   coef.     st.err.  p-value')
%         Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1:3) = [exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)), Results.DetailsM(l+tmp+sum(EstimOpt.MeaMatrix(:,i))+2,2)*exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)),pv(exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)), Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,2)*exp(Results.DetailsM(l+sum(EstimOpt.MeaMatrix(:,i))+tmp+2,1)))];
%         disp([char(EstimOpt.NamesLV(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)) ,blanks(2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp)',num2str(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,1),'%11.4f'), star_sig(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,3)), num2str(Results.DetailsM(l+1:l+2+sum(EstimOpt.MeaMatrix(:,i),1)+tmp,2:3),'%7.4f %8.4f')])
%         l = l+sum(EstimOpt.MeaMatrix(:,i))+2+tmp;
%     end
% end
%
%
% disp(' ')
% disp(['LL at convergence: ',num2str(Results.LL,'%8.4f')])
% disp(' ');
% Results.clocknote = clock;
% Results.tocnote = toc;
% [~,DayName] = weekday(now,'long');
% disp(['Estimation completed on ' DayName ', ' num2str(Results.clocknote(1)) '-' sprintf('%02.0f',Results.clocknote(2)) '-' sprintf('%02.0f',Results.clocknote(3)) ' at ' sprintf('%02.0f',Results.clocknote(4)) ':' sprintf('%02.0f',Results.clocknote(5)) ':' sprintf('%02.0f',Results.clocknote(6))])
% disp(['Estimation took ' num2str(Results.tocnote) ' seconds ('  num2str(floor(Results.tocnote/(60*60))) ' hours ' num2str(floor(rem(Results.tocnote,60*60)/60)) ' minutes ' num2str(rem(Results.tocnote,60)) ' seconds).']);
% disp(' ');
end

