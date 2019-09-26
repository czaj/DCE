function Results = LML(INPUT,Results_old,EstimOpt,OptimOpt)

% save tmp_LML
% return

global B_backup

tic

Results.bhat = [];
Results.R = [];
Results.R_out = {};
Results.stats = [];

NVarA = EstimOpt.NVarA;
NSdSim = EstimOpt.NSdSim;
NRep = EstimOpt.NRep;
NP = EstimOpt.NP;
NAlt = EstimOpt.NAlt;
NCT = EstimOpt.NCT;
NCTMiss = EstimOpt.NCTMiss;
NAltMiss = EstimOpt.NAltMiss;


%% Check data and inputs


if nargin < 3 % check no. of inputs
    error('Too few input arguments for LML(INPUT,EstimOpt,OptimOpt)')
end

disp(' ');
disp('__________________________________________________________________________________________________________________');
disp(' ');

warning off MATLAB:mir_warning_maybe_uninitialized_temporary

format shortG;
format compact;

if any(INPUT.W ~= 1)
    cprintf('Black','Estimating '); cprintf('*Black','weighted '); cprintf('Black','LML model...\n');
else
    disp('Estimating LML model ...')
end

if isfield(EstimOpt,'FullCov') == 0
    EstimOpt.FullCov = 0;
end
if ~isfield(EstimOpt,'WTP_space')
    EstimOpt.WTP_space = 0;
    EstimOpt.WTP_matrix = [];
elseif EstimOpt.WTP_space == 0
    EstimOpt.WTP_matrix = [];
end

if EstimOpt.FullCov == 0
    disp('with non-correlated random parameters ...')
    if EstimOpt.WTP_space > 0
        disp('in WTP-space ...')
    else
        disp('in preference-space ...')
    end
else
    disp('with correlated random parameters ...')
    if EstimOpt.WTP_space > 0
        disp('in WTP-space ...')
    else
        disp('in preference-space ...')
    end
end

if isfield(EstimOpt,'Dist') == 0 || isempty(EstimOpt.Dist)
    EstimOpt.Dist = zeros(1,NVarA);
    if EstimOpt.WTP_space == 0
        cprintf(rgb('DarkOrange'), 'WARNING: distributions for random parameters not specified - assuming approximate normality \n')
    else
        cprintf(rgb('DarkOrange'), 'WARNING: distributions for random parameters not specified - assuming approximate normality (monetary parameter(s) assumed approximate log-normal) \n')
        EstimOpt.Dist(end-EstimOpt.WTP_space+1:end) = 1; % cost in WTP-space models log-normally distributed
    end
else
    if ~isvector(EstimOpt.Dist)
        error('EstimOpt.Dist must be a vector')
    elseif length(EstimOpt.Dist) == 1
        EstimOpt.Dist = EstimOpt.Dist.*ones(1,NVarA);
    elseif length(EstimOpt.Dist) == NVarA
        EstimOpt.Dist = EstimOpt.Dist(:)';
    else
        error('Incorrect no. of random parameters'' distributions provided')
    end
end

if ~isfield(EstimOpt,'NOrder')
    EstimOpt.NOrder = 3;
end

disp(['Random parameters distributions: ', num2str(EstimOpt.Dist),' (0 - approximate normal, 1 - approximate lognormal, 2 - Legendre polynomial (normal), 3 - Legendre polynomial (log-normal), 4 - Step function, 5 - Linear Spline, 6 - Cubic Spline, 7 - Piece-wise Cubic Spline, 8 - Piece-wise Cubic Hermite Interpolating Spline'])
if any(EstimOpt.Dist == 0 | EstimOpt.Dist == 1)
    cprintf('Order of approximation for normal / lognormal distribution(s): ');
    cprintf('*blue',[num2str(EstimOpt.NOrder) ' ']);
    cprintf(' \n');
end
if any(EstimOpt.Dist == 2 | EstimOpt.Dist == 3)
    cprintf('Order of Legendre polynomial(s): ');
    cprintf('*blue',[num2str(EstimOpt.NOrder) ' ']);
    cprintf(' \n');
end
if any(EstimOpt.Dist == 4)
    cprintf('Number of step function segments: ');
    cprintf('*blue',[num2str(EstimOpt.NOrder) ' ']);
    cprintf(' \n');
end

if any(EstimOpt.Dist == 5 | EstimOpt.Dist == 6 | EstimOpt.Dist == 7 | EstimOpt.Dist == 8)
    cprintf('Number of spline knots (including bounds): ');
    cprintf('*blue',[num2str(EstimOpt.NOrder+2) ' ']);
    cprintf(' \n');
end

if EstimOpt.WTP_space > 0 && sum(EstimOpt.Dist(end-EstimOpt.WTP_space+1:end)==1 | EstimOpt.Dist(end-EstimOpt.WTP_space+1:end)==3) > 0 && any(mean(INPUT.Xa(:,end-EstimOpt.WTP_space+1:end)) >= 0)
    cprintf(rgb('DarkOrange'), 'WARNING: Cost attributes with log-normally distributed parameters should enter utility function with a ''-'' sign \n')
end

if EstimOpt.WTP_space > 0
    if isfield(EstimOpt, 'WTP_matrix') == 0
        WTP_att = (NVarA-EstimOpt.WTP_space)/EstimOpt.WTP_space;
        if rem(WTP_att,1) ~= 0
            error('EstimOpt.WTP_matrix associating attributes with cost parameters not provided')
        else
            if EstimOpt.WTP_space > 1
                disp(['EstimOpt.WTP_matrix associating attributes with cost parameters not provided - assuming equal shares for each of the ',num2str(EstimOpt.WTP_space),' monetary attributes'])
            end
            EstimOpt.WTP_matrix = NVarA - EstimOpt.WTP_space + kron(1:EstimOpt.WTP_space,ones(1,WTP_att));
            %         tic; EstimOpt.WTP_matrix = 1:EstimOpt.WTP_space;...
            %         EstimOpt.WTP_matrix = EstimOpt.WTP_matrix(floor((0:size(EstimOpt.WTP_matrix,2)*WTP_att-1)/WTP_att)+1); toc
        end
        %     elseif ~isequal(size(EstimOpt.WTP_matrix),[NVarA-EstimOpt.WTP_space,EstimOpt.WTP_space])
    elseif size(EstimOpt.WTP_matrix,2) ~= NVarA - EstimOpt.WTP_space
        error('Dimensions of EstimOpt.WTP_matrix not correct - for each non-monetary attribute provide no. of attribute to multiply it with')
    else
        EstimOpt.WTP_matrix = EstimOpt.WTP_matrix(:)';
    end
end

if isfield(EstimOpt,'NGrid') == 0 
    NGrid = 1000; % Train uses 1000
    EstimOpt.NGrid = NGrid;
else
    NGrid = EstimOpt.NGrid;
end

if isfield(EstimOpt, 'StepFun') == 1
    EstimOpt.StepVar = size(EstimOpt.StepFun(ones(EstimOpt.NVarA,1)),1);
    cprintf(rgb('DarkOrange'), 'Adding step function defined by user \n')
else
    EstimOpt.StepVar = 0;
end

if isfield(EstimOpt, 'PlotIndx') == 0
    EstimOpt.PlotIndx = 0; % Do not draw a plot
end

if isfield(EstimOpt, 'NoOutput') == 0
    EstimOpt.NoOutput = 0; % Do not draw a plot
end
if isfield(EstimOpt, 'FitCens') == 0
    EstimOpt.FitCens = 0; % Do not censor fitted values in estimation process
end

if isfield(EstimOpt,'NamesA') == 0 || isempty(EstimOpt.NamesA) || length(EstimOpt.NamesA) ~= NVarA
    EstimOpt.NamesA = (1:NVarA)';
    EstimOpt.NamesA = cellstr(num2str(EstimOpt.NamesA));
elseif size(EstimOpt.NamesA,1) ~= NVarA
    EstimOpt.NamesA = EstimOpt.NamesA';
end

% gcp; % start paralell pool (commented out - we don't use it for now)


%% Starting values

% save tmp2
% return

NVar = sum((EstimOpt.Dist == 0 | EstimOpt.Dist == 1)*EstimOpt.NOrder + ...
    (EstimOpt.Dist == 2 | EstimOpt.Dist == 3)*EstimOpt.NOrder + ...
    (EstimOpt.Dist == 4)*(EstimOpt.NOrder-1) + ...
    (EstimOpt.Dist == 5 | EstimOpt.Dist == 6 | EstimOpt.Dist == 7 | EstimOpt.Dist == 8)*(EstimOpt.NOrder+1),2) + ...
    EstimOpt.StepVar;

if exist('B_backup','var') && ~isempty(B_backup) && isvector(B_backup)
    B_backup = B_backup(:);
end

if EstimOpt.FullCov == 0
    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == NVar
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'LML_d') && isfield(Results_old.LML_d,'b0') % starting values provided
        Results_old.LML_d.b0_old = Results_old.LML_d.b0(:);
        Results_old.LML_d = rmfield(Results_old.LML_d,'b0');
        if length(Results_old.LML_d.b0_old) ~= NVar
            cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.LML_d = rmfield(Results_old.LML_d,'b0_old');
        else
            b0 = Results_old.LML_d.b0_old(:);
        end
    end
    if  ~exist('b0','var')
        b0 = zeros(NVar,1);
    end
else
    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == NVar + NVarA*(NVarA-1)/2
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'LML') && isfield(Results_old.LML,'b0') % starting values provided
        Results_old.LML.b0_old = Results_old.LML.b0(:);
        Results_old.LML = rmfield(Results_old.LML,'b0');
        if length(Results_old.LML.b0_old) ~= NVar + NVarA*(NVarA-1)/2
            cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.LML = rmfield(Results_old.LML,'b0_old');
        else
            b0 = Results_old.LML.b0_old(:);
        end
    end
    if  ~exist('b0','var')
        if isfield(Results_old,'LML_d') && isfield(Results_old.LML_d,'bhat') % starting values provided
            b0 = [Results_old.LML_d.bhat; zeros(NVarA*(NVarA-1)/2,1)];
        else
            b0 = zeros(NVar + NVarA*(NVarA-1)/2,1);
        end
    end
end


%% Optimization Options


if  isfield(EstimOpt,'BActive')
    EstimOpt.BActive = EstimOpt.BActive(:)';
    if size(EstimOpt.BActive,2) ~= NVar
        cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of constraints - ignoring \n')
        EstimOpt.BActive = ones(1,length(b0));
    end
else
    EstimOpt.BActive = ones(1,length(b0));
end

if EstimOpt.ConstVarActive == 1
    if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive) || sum(EstimOpt.BActive == 0) == 0
        error('Are there any constraints on model parameters (EstimOpt.ConstVarActive)? Constraints not provided (EstimOpt.BActive).')
    elseif length(b0) ~= length(EstimOpt.BActive)
        error('Check no. of constraints')
    end
    disp(['Starting values: ' mat2str(b0',2)])
    disp(['Parameters with zeros are constrained to their initial values: ' mat2str(EstimOpt.BActive')])
else
    if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive) || sum(EstimOpt.BActive == 0) == 0
        EstimOpt.BActive = ones(1,length(b0));
        disp(['Starting values: ' mat2str(b0',2)])
    else
        if length(b0) ~= length(EstimOpt.BActive)
            error('Check no. of constraints')
        else
            disp(['Starting values: ' mat2str(b0',2)])
            disp(['Parameters with zeros are constrained to their initial values: ' mat2str(EstimOpt.BActive')])
        end
    end
end


%% Bounds


if ~(isfield(EstimOpt,'Bounds') == 0 || isempty(EstimOpt.Bounds))
    if size(EstimOpt.Bounds,1) == 1 && size(EstimOpt.Bounds,2) == 2
        EstimOpt.Bounds = EstimOpt.Bounds(ones(1,NVarA),:);
    elseif size(EstimOpt.Bounds,1) == 2 && size(EstimOpt.Bounds,2) == 1
        EstimOpt.Bounds = EstimOpt.Bounds(:,ones(1,NVarA))';
    elseif size(EstimOpt.Bounds,1) == 2 && size(EstimOpt.Bounds,2) == NVarA
        EstimOpt.Bounds = EstimOpt.Bounds';
    elseif size(EstimOpt.Bounds,1) == NVarA && size(EstimOpt.Bounds,2) == 2
        % Size ok. Do nothing. Test bounds later.
    else
        error('Incorrect no. of Bounds provided')
    end
else
    if isfield(Results_old,'MXL') && isfield(Results_old.MXL,'bhat') && ~isempty(Results_old.MXL.bhat) && ... % MXL exists
            (size(Results_old.MXL.bhat(:),1) == ((NVarA + sum(1:NVarA)) + Results_old.MXL.EstimOpt.NVarS)) &&  ... % MXL has correct no. of parameters
            all(Results_old.MXL.EstimOpt.Dist == 0 | Results_old.MXL.EstimOpt.Dist == 1) % all parameters were normally or log-normally distributed
        VC = tril(ones(NVarA));
        VC(VC == 1) = Results_old.MXL.bhat(NVarA+1:NVarA+sum(1:NVarA));
        VC = VC*VC';
        EstimOpt.Bounds = [Results_old.MXL.bhat(1:NVarA) - 2*sqrt(diag(VC)),Results_old.MXL.bhat(1:NVarA) + 2*sqrt(diag(VC))];
        EstimOpt.Bounds(Results_old.MXL.EstimOpt.Dist == 1) = exp(EstimOpt.Bounds(Results_old.MXL.EstimOpt.Dist == 1));
    elseif isfield(Results_old,'MXL_d') && isfield(Results_old.MXL_d,'bhat') && ~isempty(Results_old.MXL_d.bhat) && ... % MXL exists
            (size(Results_old.MXL_d.bhat(:),1) == (NVarA*2 + Results_old.MXL_d.EstimOpt.NVarS)) && ... % MXL_d has correct no. of parameters
            all(Results_old.MXL_d.EstimOpt.Dist == 0 | Results_old.MXL_d.EstimOpt.Dist == 1) % all parameters were normally or log-normally distributed
        EstimOpt.Bounds = [Results_old.MXL_d.bhat(1:NVarA) - 2*abs(Results_oldMXL_dMXL.bhat(NVarA+1:NVarA*2)),Results_old.MXL_d.bhat(1:NVarA) + 2*abs(Results_old.MXL_d.bhat(NVarA+1:NVarA*2))];
        EstimOpt.Bounds(Results_old.MXL_d.EstimOpt.Dist == 1) = exp(EstimOpt.Bounds(Results_old.MXL_d.EstimOpt.Dist == 1)); %  median, not mean
    else % run quick MXL_d and use mean +/- 2*s.d.
        disp('Bounds not provided - using a quick MXL_d model to generate')
        EstimOpt_tmp = EstimOpt;
        %         EstimOpt_tmp.Display = 0;
        EstimOpt_tmp.NumGrad = 0;
        EstimOpt_tmp.NRep = 1e2;
        EstimOpt_tmp.Dist = [];
        EstimOpt_tmp.HessEstFix = 1;
        OptimOpt_tmp = optimoptions('fminunc');
        OptimOpt_tmp.Algorithm = 'quasi-newton';
        OptimOpt_tmp.GradObj = 'on';
        OptimOpt_tmp.Hessian = 'off';
        OptimOpt_tmp.Display = 'off';
        OptimOpt_tmp.FunValCheck= 'off';
        OptimOpt_tmp.Diagnostics = 'off';
        OptimOpt_tmp.OptimalityTolerance = 1e-3;
        OptimOpt_tmp.StepTolerance = 1e-3;
        Results_old.MXL_d = MXL(INPUT,Results_old,EstimOpt_tmp,OptimOpt_tmp);
        EstimOpt.Bounds = [Results_old.MXL_d.bhat(1:NVarA) - 2*abs(Results_old.MXL_d.bhat(NVarA+1:NVarA*2)),Results_old.MXL_d.bhat(1:NVarA) + 2*abs(Results_old.MXL_d.bhat(NVarA+1:NVarA*2))];
        disp(' ');
        disp('__________________________________________________________________________________________________________________');
        disp(' ');
    end
end

if any(EstimOpt.Bounds(EstimOpt.Dist == 1 | EstimOpt.Dist == 3,1) <= 0)
    cprintf(rgb('DarkOrange'),'WARNING: Lower bound of approximate log-normally distributed parameters must be  > 0. Adjusting offendig lower Bound(s) to realmin. \n')
    EstimOpt.Bounds((EstimOpt.Dist == 1 | EstimOpt.Dist == 3),1) = max(realmin,EstimOpt.Bounds((EstimOpt.Dist == 1 | EstimOpt.Dist == 3),1));
end

if any(EstimOpt.Bounds(:,1) >= EstimOpt.Bounds(:,2))
    error('Lower Bound(s) greater than upper Bound(s).')
end


%% Generate pseudo-random draws


if isfield(EstimOpt,'Seed1') == 1
    rng(EstimOpt.Seed1);
end
cprintf('Simulation with ');
cprintf('*blue',[num2str(NRep) ' ']);

if EstimOpt.Draws == 1
    cprintf('*blue','Pseudo-random '); cprintf('draws \n');
    err_mtx = rand(NP*NRep,NVarA);
elseif EstimOpt.Draws == 2 % LHS
    cprintf('*blue','Latin Hypercube Sampling '); cprintf('draws \n');
    err_mtx = lhsnorm(zeros((NVarA)*NP,1),diag(ones((NVarA)*NP,1)),NRep);
    err_mtx = reshape(err_mtx, NRep*NP, NVarA);
elseif EstimOpt.Draws >= 3 % Quasi random draws
    if EstimOpt.Draws == 3
        cprintf('*blue','Halton '); cprintf('draws (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = haltonset(NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
    elseif EstimOpt.Draws == 4 % apply reverse-radix scrambling
        cprintf('*blue','Halton '); cprintf('draws with reverse radix scrambling (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = haltonset(NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
        hm1 = scramble(hm1,'RR2');
    elseif EstimOpt.Draws == 5
        cprintf('*blue','Sobol '); cprintf('draws (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = sobolset(NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
    elseif EstimOpt.Draws == 6
        cprintf('*blue','Sobol '); cprintf('draws with random linear scramble and random digital shift (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = sobolset(NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
        hm1 = scramble(hm1,'MatousekAffineOwen');
    end
    err_mtx = net(hm1,NP*NRep); % this takes every point:
    clear hm1;
end

err_mtx = floor((NGrid).*err_mtx)' + 1;
GridMat = zeros(NVarA,NGrid);
for i = 1:NVarA
    GridMat(i,:) = EstimOpt.Bounds(i,1):((EstimOpt.Bounds(i,2) - EstimOpt.Bounds(i,1))/(NGrid-1)):EstimOpt.Bounds(i,2);
end

% if EstimOpt.WTP_space > 0
%     GridMat(1:end-EstimOpt.WTP_space,:) = GridMat(1:end-EstimOpt.WTP_space,:).*GridMat(EstimOpt.WTP_matrix,:);
% end

b_GridMat = B_lml(GridMat,EstimOpt); % NV x NGrid
% b_mtx = zeros(NVarA,NP*NRep);
% for i = 1:size(b_GridMat,1)
%     mod(i,NVarA)
% %     b_mtx(i,:) = b_GridMat(i,err_mtx(mod(i,NVarA),:)); % NV x NP*NRep 
% end

for i = 1:NVarA
    err_mtx(i,:) = GridMat(i,err_mtx(i,:));
end
b_mtx = B_lml(err_mtx,EstimOpt); % NV x NP*NRep % this is very fast too, so using elements from b_GridMat does not help much
% if EstimOpt.WTP_space > 0
%     err_mtx(1:end-EstimOpt.WTP_space,:) = err_mtx(1:end-EstimOpt.WTP_space,:).*err_mtx(EstimOpt.WTP_matrix,:);
% end
if isfield(EstimOpt, 'Drawskeep') && ~isempty(EstimOpt.Drawskeep) && EstimOpt.Drawskeep == 1
    Results.err = err_mtx;
end

if EstimOpt.StepVar > 0
    b_mtx = [b_mtx; EstimOpt.StepFun(err_mtx)];
end


%% Display Options
if EstimOpt.NoOutput == 1
    if EstimOpt.HessEstFix ~= 1
        cprintf(rgb('DarkOrange'),'WARNING: Setting HessEstFix to 1, output will not be generated anyway (EstimOpt.NoOutput = 1) \n')
    end
    EstimOpt.HessEstFix = 1;
end

if ((isfield(EstimOpt,'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    OptimOpt.GradObj = 'on';
end

if (isfield(EstimOpt,'ConstVarActive') == 0 || EstimOpt.ConstVarActive == 0) && isequal(OptimOpt.Algorithm,'quasi-newton') && isequal(OptimOpt.Hessian,'user-supplied')
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied Hessian off - quasi-newton algorithm does not use it anyway \n')
    OptimOpt.Hessian = 'off';
end

if EstimOpt.NumGrad == 1 && EstimOpt.ApproxHess == 0
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian only available if analythical gradient on \n')
    EstimOpt.ApproxHess = 1;
end

if EstimOpt.RobustStd == 1 && (EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2)
    EstimOpt.RobustStd = 0;
    cprintf(rgb('DarkOrange'),'WARNING: Setting off robust standard errors, they do not matter for BHHH aproximation of hessian \n')
end

fprintf('\n')
cprintf('Optimization algorithm: '); cprintf('*Black',[OptimOpt.Algorithm '\n'])

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
fprintf('\n')


%% Rescructure data


INPUT.XXa = reshape(INPUT.Xa,[NAlt*NCT,NP,NVarA]);
INPUT.XXa = permute(INPUT.XXa, [1,3,2]);
INPUT.YY = reshape(INPUT.Y,[NAlt*NCT,NP]);

cprintf('Conducting pre-estimation calculations for ');
cprintf('*blue',[num2str(NGrid) ' ']);
cprintf('grid points. \n');
tocnote_00 = toc;

b_gird = reshape(err_mtx,[NVarA,NRep,NP]);
if EstimOpt.WTP_space > 0
    b_gird(1:end-EstimOpt.WTP_space,:,:) = b_gird(1:end-EstimOpt.WTP_space,:,:).*b_gird(EstimOpt.WTP_matrix,:,:);
end

GridProbs = zeros([NP,NRep]);
XXa = INPUT.XXa;
% parfor n = 1:NP    
if any(isnan(XXa(:))) == 0 % faster version for complete dataset
    YYy = INPUT.YY==1;
    for n = 1:NP % switch parfor off for now and run Matlab in paralell processes instead
        U = reshape(XXa(:,:,n)*b_gird(:,:,n),[NAlt,NCT,NRep]);    
        U = exp(U - max(U,[],1)); % rescale utility to avoid exploding
        U_sum = reshape(sum(U,1),[NCT,NRep]);
        YYy_n = YYy(:,n);
        U_selected = reshape(U(YYy_n(:,ones(NRep,1))),[NCT,NRep]);
        GridProbs(n,:) = prod(U_selected./U_sum,1);
    end
else   
    for n = 1:NP % switch parfor off for now and run Matlab in paralell processes instead
        YnanInd = ~isnan(INPUT.YY(:,n));
        XXa_n = XXa(:,:,n);
        U = reshape(XXa_n(YnanInd,:)*b_gird(:,:,n),[NAltMiss(n),NCTMiss(n),NRep]);
        U = exp(U - max(U,[],1)); % rescale utility to avoid exploding
        U_sum = reshape(sum(U,1),[NCTMiss(n),NRep]);
        YYy_n = INPUT.YY(:,n)==1;
        U_selected = reshape(U(YYy_n(YnanInd,ones(NRep,1))),[NCTMiss(n),NRep]);
        GridProbs(n,:) = prod(U_selected./U_sum,1);
    end
end
tocnote_01 = toc-tocnote_00;
cprintf(['Pre-estimation completed in ' num2str(tocnote_01) ' seconds ('  num2str(floor(tocnote_01/(60*60))) ' hours ' num2str(floor(rem(tocnote_01,60*60)/60)) ' minutes ' num2str(rem(tocnote_01,60)) ' seconds).\n\n']);


%% Estimation


LLfun = @(B) LL_lml_MATlike(GridProbs,b_mtx,INPUT.W,EstimOpt,OptimOpt,B);

if EstimOpt.ConstVarActive == 0
    if EstimOpt.HessEstFix == 0
        [Results.bhat,LL,Results.exitf,Results.output,Results.g,Results.hess] = fminunc(LLfun,b0,OptimOpt);
    else
        [Results.bhat,LL,Results.exitf,Results.output,Results.g] = fminunc(LLfun,b0,OptimOpt);
    end
    %     options_tmp = optimset('MaxFunEvals',1e100,'MaxIter',1e3,'TolFun',1e-6,'TolX',1e-6,'OutputFcn',@outputf);
    %     [Results.beta,LL,Results.exitf,Results.output] = fminsearch(LLfun,b0,options_tmp);
elseif EstimOpt.ConstVarActive == 1 % equality constraints
    EstimOpt.CONS1 = diag(1 - EstimOpt.BActive);
    EstimOpt.CONS1(sum(EstimOpt.CONS1,1)==0,:) = [];
    EstimOpt.CONS2 = zeros(size(EstimOpt.CONS1,1),1);
    %     EstimOpt.CONS1 = sparse(EstimOpt.CONS1);
    %     EstimOpt.CONS2 = sparse(EstimOpt.CONS2);
    if EstimOpt.HessEstFix == 0
        [Results.bhat,LL,Results.exitf,Results.output,Results.lambda,Results.g,Results.hess] = fmincon(LLfun,b0,[],[],EstimOpt.CONS1,EstimOpt.CONS2,[],[],[],OptimOpt);
    else
        [Results.bhat,LL,Results.exitf,Results.output,Results.lambda,Results.g] = fmincon(LLfun,b0,[],[],EstimOpt.CONS1,EstimOpt.CONS2,[],[],[],OptimOpt);
    end
end


%% Hessian calculations

if EstimOpt.NoOutput == 0
    
    LLfun2 = @(B) LL_lml(GridProbs,b_mtx,EstimOpt,B);

    if EstimOpt.HessEstFix == 0 % this will fail if there is no gradient available!
        try
            [Results.LLdetailed,Results.jacobian] = LLfun2(Results.bhat);
        catch % theErrorInfo
            Results.LLdetailed = LLfun2(Results.bhat);
            Results.jacobian = numdiff(@(B) INPUT.W.*LLfun2(B),Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
            Results.jacobian = Results.jacobian.*INPUT.W(:,ones(1,size(Results.jacobian,2)));
        end
    elseif EstimOpt.HessEstFix == 1
        if isequal(OptimOpt.GradObj,'on') && EstimOpt.NumGrad == 0
            [Results.LLdetailed,Results.jacobian] = LLfun2(Results.bhat);
            Results.jacobian = Results.jacobian.*INPUT.W(:,ones(1,size(Results.jacobian,2)));
        else
            Results.LLdetailed = LLfun2(Results.bhat);
            Results.jacobian = numdiff(@(B) INPUT.W.*LLfun2(B),Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
            Results.jacobian = Results.jacobian.*INPUT.W(:,ones(1,size(Results.jacobian,2)));
        end
    elseif EstimOpt.HessEstFix == 2
        Results.jacobian = jacobianest(@(B) INPUT.W.*LLfun2(B),Results.bhat);
    elseif EstimOpt.HessEstFix == 3
        Results.LLdetailed = LLfun2(Results.bhat);
        Results.hess = hessian(@(B) sum(INPUT.W.*LLfun2(B)),Results.bhat);
    elseif EstimOpt.HessEstFix == 4

    end
    R2 = mean(exp(-Results.LLdetailed/EstimOpt.NCT),1);
    Results.LLdetailed = Results.LLdetailed.*INPUT.W;

    if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
        Results.hess = Results.jacobian'*Results.jacobian;
    end
    EstimOpt.BLimit = (sum(Results.hess) == 0 & EstimOpt.BActive == 1);
    EstimOpt.BActive(EstimOpt.BLimit == 1) = 0;
    Results.hess = Results.hess(EstimOpt.BActive == 1,EstimOpt.BActive == 1);
    Results.ihess = inv(Results.hess);
    Results.ihess = direcXpnd(Results.ihess,EstimOpt.BActive);
    Results.ihess = direcXpnd(Results.ihess',EstimOpt.BActive);


    %% Output


    % save tmp11
    % return

    Results.LL = -LL;
    Results.b_mtx = b_mtx;
    Results.GridMat = GridMat;

    % Results.P = P_lml(Results.bhat,b_GridMat);
    err_tmp = unique(err_mtx','rows')';
    b_tmp = B_lml(err_tmp,EstimOpt);
    if EstimOpt.StepVar > 0
        b_tmp = [b_tmp; EstimOpt.StepFun(err_tmp)];
    end

    [Results.P,Results.DistStats,Results.GridPlot] = P_lml(b_tmp,Results.bhat,err_tmp,Results.ihess,EstimOpt);

    EstimOpt.params = length(b0) - sum(EstimOpt.BActive == 0) + sum(EstimOpt.BLimit == 1);
    Results.stats = [Results.LL;Results_old.MNL0.LL;1-Results.LL/Results_old.MNL0.LL;R2;((2*EstimOpt.params-2*Results.LL))/EstimOpt.NObs;((log(EstimOpt.NObs)*EstimOpt.params-2*Results.LL))/EstimOpt.NObs;EstimOpt.NObs;EstimOpt.NP;EstimOpt.params];

    %File Output
    Results.EstimOpt = EstimOpt;
    Results.OptimOpt = OptimOpt;
    Results.INPUT = INPUT;
    Results.Dist = EstimOpt.Dist';

    % return

    if EstimOpt.PlotIndx > 0
        EstimOpt.Plot = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:NVarA
    %         Grid_i = mean(reshape(GridMat(i,:), [10,NGrid/10]),1); 
    %         P_tmp = sum(reshape(Results.P,[10,NGrid/10]),1);
            if rem(NVarA,2) == 0
                subplot(NVarA/2,2,i);
            else
                subplot(NVarA,1,i);
            end
    %         bar(Grid_i,P_tmp)
    %         tmp = sortrows([GridMat(i,:)',Results.P'])';
    %         plot(tmp(1,:),tmp(2,:))
            %plot(Results.GridPlot{i}',Results.P{i}')
            plot(Results.GridPlot(i,~isnan(Results.GridPlot(i,:)))',Results.P(i,~isnan(Results.P(i,:)))')
            title(EstimOpt.NamesA(i))
        end
        EstimOpt.Plot2 = gcf;
    end


    % save tmp2

    %% Output

    Template1 = {'Details'};
    Template2 = {'Details'};

    ST = {};
    Results.std = sqrt(diag(Results.ihess));
    % for i=1:length(Results.bhat)/EstimOpt.NVarA
    %     Results.(strcat('Details',num2str(i)))(1:NVarA,1) = Results.bhat(1+(i-1)*EstimOpt.NVarA:i*EstimOpt.NVarA);
    %     Results.(strcat('Details',num2str(i)))(1:NVarA,3:4) = [Results.std(1+(i-1)*EstimOpt.NVarA:i*EstimOpt.NVarA),pv(Results.bhat(1+(i-1)*EstimOpt.NVarA:i*EstimOpt.NVarA),Results.std(1+(i-1)*EstimOpt.NVarA:i*EstimOpt.NVarA))];
    %     Template1 = [Template1,{strcat('Details',num2str(i))}];
    %     Template2 = [Template2,{strcat('Details',num2str(i))}];
    %     Names.(strcat('Details',num2str(i))) = EstimOpt.NamesA;
    %     if i == 1 && length(Results.bhat)/EstimOpt.NVarA ~= 1 
    %         Heads.(strcat('Details',num2str(i))) = {strcat('Segment ',num2str(i));'tc'};
    %     elseif i<length(Results.bhat)/EstimOpt.NVarA && length(Results.bhat)/EstimOpt.NVarA ~= 1 
    %         Heads.(strcat('Details',num2str(i))) = {strcat('Segment ',num2str(i));'lc'};
    %     elseif i == length(Results.bhat)/EstimOpt.NVarA && length(Results.bhat)/EstimOpt.NVarA ~= 1 
    %         Heads.(strcat('Details',num2str(i))) = {strcat('Segment ',num2str(i));'lb'};
    %     else
    %         Heads.(strcat('Details',num2str(i))) = {strcat('Segment ',num2str(i));'tb'};
    %     end
    % end
    Heads.Details = {};
    if any(EstimOpt.Dist == 0 | EstimOpt.Dist == 1)
        type = ' (dist. approx.)';
        CorNOrder = EstimOpt.NOrder;
    elseif any(EstimOpt.Dist == 2 | EstimOpt.Dist == 3)
        type = ' (Legendre poly)';
        CorNOrder = EstimOpt.NOrder;
    elseif any(EstimOpt.Dist == 4)
        type = ' (step function)';
        CorNOrder = EstimOpt.NOrder-1;
    elseif any(EstimOpt.Dist == 5 | EstimOpt.Dist == 6 | EstimOpt.Dist == 7 | EstimOpt.Dist == 8)
        type = ' (spline knots)';
        CorNOrder = EstimOpt.NOrder+1;
    end
        
    for i=1:length(Results.bhat)/EstimOpt.NVarA
        Results.Details(1:NVarA,i*4-3) = Results.bhat(1+(i-1)*EstimOpt.NVarA:i*EstimOpt.NVarA);
        Results.Details(1:NVarA,i*4-1:i*4) = [Results.std(1+(i-1)*EstimOpt.NVarA:i*EstimOpt.NVarA),pv(Results.bhat(1+(i-1)*EstimOpt.NVarA:i*EstimOpt.NVarA),Results.std(1+(i-1)*EstimOpt.NVarA:i*EstimOpt.NVarA))];
        if i <= CorNOrder
            Heads.Details = [Heads.Details;{strcat('Segment ',num2str(i),type)}];
        elseif isfield(EstimOpt,'StepFun') == 1
            Heads.Details = [Heads.Details;{strcat('Segment ',num2str(i),' (user''s step function)')}];
        elseif EstimOpt.FullCov == 1
            Heads.Details = [Heads.Details;{strcat('Segment ',num2str(i),' (correlation)')}];
        end
    end
    Heads.Details = [Heads.Details;{'tb'}];

    Names.Details = EstimOpt.NamesA;

    %% Tworzenie naglowka

    Head = cell(1,2);
    if EstimOpt.FullCov == 0
        Head(1,1) = {'LML_d'};
    else
        Head(1,1) = {'LML'};
    end

    if EstimOpt.WTP_space > 0
        Head(1,2) = {'in WTP-space'};
    else
        Head(1,2) = {'in preference-space'};
    end

    %% Tworzenie stopki

    Tail = cell(18,2);
    Tail(2,1) = {'Model diagnostics'};
    Tail(3:18,1) = {'LL at convergence';'LL at constant(s) only';strcat('McFadden''s pseudo-R',char(178));strcat('Ben-Akiva-Lerman''s pseudo-R',char(178));'AIC/n';'BIC/n';'n (observations)';'r (respondents)';'k (parameters)';'';'Estimation method';'Simulation with';'Grid points';'Optimization method';'Gradient';'Hessian'};

    if isfield(Results_old,'MNL0') && isfield(Results_old.MNL0,'LL')
        Tail(3:11,2) = num2cell(Results.stats);
    end

    if any(INPUT.W ~= 1)
        Tail(13,2) = {'weighted simulated maximum likelihood'};
    else
        Tail(13,2) = {'simulated maximum likelihood'};
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

    Tail(15,2) = {num2str(NGrid)};
    Tail(16,2) = {OptimOpt.Algorithm};

    if strcmp(OptimOpt.GradObj,'on')
        if EstimOpt.NumGrad == 0
            Tail(17,2) = {'user-supplied, analytical'};
        else
            Tail(17,2) = {['user-supplied, numerical ',num2str(OptimOpt.FinDiffType)]};
        end
    else
        Tail(17,2) = {['built-in, ',num2str(OptimOpt.FinDiffType)]};
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
    Tail(18,2) = {outHessian};

    %% Statistics

    Statistics.Bounds = EstimOpt.Bounds;
    Statistics.Quantile = Results.DistStats(:,3:end,1);
    Statistics.Mean(:,1) = Results.DistStats(:,1,1);
    Statistics.Mean(:,3:4) = [Results.DistStats(:,1,2),pv(Results.DistStats(:,1,1),Results.DistStats(:,1,2))];
    Statistics.Sd(:,1) = Results.DistStats(:,2,1);
    Statistics.Sd(:,3:4) = [Results.DistStats(:,2,2),pv(Results.DistStats(:,2,1),Results.DistStats(:,2,2))];

    Names.Statistics = EstimOpt.NamesA;
    Heads.StatisticsQ = {'q0.025','q0.1','q0.25','q0.5','q0.75','q0.9','q0.975'};
    EstimOpt.StatTypes = ['b','b','m','m','m','m','m','m','m','m','q','q','q','q','q','q','q'];

    if EstimOpt.Display~=0

        Results.Dist = EstimOpt.Dist;

        Results.R_out = genOutput(EstimOpt, Results, Head, Tail, Names, Template1, Template2, Heads, ST, 'lml', Statistics);
    end
else
    Results.LL = -LL;
    Results.b_mtx = b_mtx;
    Results.GridMat = GridMat;
    EstimOpt.params = length(b0) - sum(EstimOpt.BActive == 0); % This may need a correction (Does not work when BLimit is ~= 0)
    Results.stats = [Results.LL;Results_old.MNL0.LL;1-Results.LL/Results_old.MNL0.LL;NaN;((2*EstimOpt.params-2*Results.LL))/EstimOpt.NObs;((log(EstimOpt.NObs)*EstimOpt.params-2*Results.LL))/EstimOpt.NObs;EstimOpt.NObs;EstimOpt.NP;EstimOpt.params];

end
