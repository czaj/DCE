function Results = LMXL(INPUT,Results_old,EstimOpt,OptimOpt)


% save tmp_MXL
% return

global B_backup

tic

Results.bhat = [];
Results.R = [];
Results.R_out = {};
Results.stats = [];

NVarA = EstimOpt.NVarA;
NSdSim = EstimOpt.NSdSim;


%% Check data and inputs


if nargin < 3 % check no. of inputs
    error('Too few input arguments for MXL(INPUT,EstimOpt,OptimOpt)')
end

disp(' ');
disp('__________________________________________________________________________________________________________________');
disp(' ');

warning off MATLAB:mir_warning_maybe_uninitialized_temporary

format shortG;
format compact;

if any(INPUT.W ~= 1)
    cprintf('Black','Estimating '); cprintf('*Black','weighted '); cprintf('Black','Logit-MXL model...\n');
else
    disp('Estimating Logit-MXL model ...')
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
    %EstimOpt.Dist(1) = 1; % scale distributed log-normally (does not matter for MXL)
    if EstimOpt.WTP_space == 0
        cprintf(rgb('DarkOrange'), 'WARNING: distributions for random parameters not specified - assuming normality \n')
    else
        cprintf(rgb('DarkOrange'), 'WARNING: distributions for random parameters not specified - assuming normality (monetary parameters assumed log-normal) \n')
        EstimOpt.Dist(end-EstimOpt.WTP_space+1:end) = 1; % cost in WTP-space models log-normally distributed
    end
else
    if length(EstimOpt.Dist) == 1
        EstimOpt.Dist = EstimOpt.Dist.*ones(1,NVarA);
    elseif length(EstimOpt.Dist) == NVarA
        EstimOpt.Dist = EstimOpt.Dist(:)';
    else
        error('Incorrect no. of random parameters'' distributions provided')
    end
end


if isfield(EstimOpt,'Bounds') == 0 || isempty(EstimOpt.Bounds) 
    error('Specify bounds for distributions')
else
    if size(EstimOpt.Bounds,1) == 1 && size(EstimOpt.Bounds,2) == 2
        EstimOpt.Bounds = EstimOpt.Bounds(ones(1,NVarA),:);
    elseif size(EstimOpt.Bounds,1) == NVarA && size(EstimOpt.Bounds,2) == 2
     
    else
        error('Incorrect no. of Bounds provided')
    end
    if any(EstimOpt.Bounds(EstimOpt.Dist == 1 | EstimOpt.Dist == 3,1) <= 0)
       error('For log normally distributed variables lower bound needs to be bigger than 0') 
    end
end

if isfield(EstimOpt,'Grid') == 0 || isempty(EstimOpt.Grid) 
    EstimOpt.Grid = 1000;
end

disp(['Random parameters distributions: ', num2str(EstimOpt.Dist),' (0 - normal (approx), 1 - lognormal (approx), 2 - Legendre (norm), 3 - Legendre (lnorm)'])

if EstimOpt.WTP_space > 0 && sum(EstimOpt.Dist(end-EstimOpt.WTP_space+1:end)==1) > 0 && any(mean(INPUT.Xa(:,end-EstimOpt.WTP_space+1:end)) >= 0)
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

if ~isfield(EstimOpt, 'Order')
    EstimOpt.Order = 3;
end


if isfield(EstimOpt,'NamesA') == 0 || isempty(EstimOpt.NamesA) || length(EstimOpt.NamesA) ~= NVarA
    EstimOpt.NamesA = (1:NVarA)';
    EstimOpt.NamesA = cellstr(num2str(EstimOpt.NamesA));
elseif size(EstimOpt.NamesA,1) ~= NVarA
    EstimOpt.NamesA = EstimOpt.NamesA';
end

%% Starting values
Params = sum((EstimOpt.Dist == 0 | EstimOpt.Dist == 1)*2 + (EstimOpt.Dist == 2 | EstimOpt.Dist == 3)*EstimOpt.Order,2);

if EstimOpt.FullCov == 0
    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == Params
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'LMXL_d') && isfield(Results_old.LMXL_d,'b0') % starting values provided
        Results_old.LMXL_d.b0_old = Results_old.LMXL_d.b0(:);
        Results_old.LMXL_d = rmfield(Results_old.LMXL_d,'b0');
        if length(Results_old.LMXL_d.b0_old) ~= Params
            cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.LMXL_d = rmfield(Results_old.LMXL_d,'b0_old');
        else
            b0 = Results_old.LMXL_d.b0_old(:);
        end
    end
    if  ~exist('b0','var')
        b0 = zeros(Params,1);
    end
else
    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == Params + NVarA*(NVarA-1)/2
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'LMXL') && isfield(Results_old.LMXL,'b0') % starting values provided
        Results_old.LMXL.b0_old = Results_old.LMXL.b0(:);
        Results_old.LMXL = rmfield(Results_old.LMXL,'b0');
        if length(Results_old.LMXL.b0_old) ~= Params + NVarA*(NVarA-1)/2
            cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.LMXL = rmfield(Results_old.LMXL,'b0_old');
        else
            b0 = Results_old.LMXL.b0_old(:);
        end
    end
    if  ~exist('b0','var')
        if isfield(Results_old,'LMXL_d') && isfield(Results_old.LMXL_d,'bhat') % starting values provided
            b0 = [Results_old.LMXL_d.bhat; zeros(NVarA*(NVarA-1)/2,1)];
        else
            b0 = zeros(Params+ NVarA*(NVarA-1)/2,1);
        end
    end
end

%% Optimization Options


if  isfield(EstimOpt,'BActive')
    EstimOpt.BActive = EstimOpt.BActive(:)';
else
    EstimOpt.BActive = ones(1,length(b0));
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
% cprintf('Simulation with ');
% cprintf('*blue',[num2str(EstimOpt.NRep) ' ']);
% 
% if EstimOpt.Draws == 1
%     cprintf('*blue','Pseudo-random '); cprintf('draws \n');
%     err_mtx = randn(EstimOpt.NP*EstimOpt.NRep, NVarA);
% elseif EstimOpt.Draws == 2 % LHS
%     cprintf('*blue','Latin Hypercube Sampling '); cprintf('draws \n');
%     err_mtx=lhsnorm(zeros((NVarA)*EstimOpt.NP,1),diag(ones((NVarA)*EstimOpt.NP,1)),EstimOpt.NRep);
%     err_mtx = reshape(err_mtx, EstimOpt.NRep*EstimOpt.NP, NVarA);
% elseif EstimOpt.Draws >= 3 % Quasi random draws
%     if EstimOpt.Draws == 3
%         cprintf('*blue','Halton '); cprintf('draws (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
%         hm1 = haltonset(NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
%     elseif EstimOpt.Draws == 4 % apply reverse-radix scrambling
%         cprintf('*blue','Halton '); cprintf('draws with reverse radix scrambling (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
%         hm1 = haltonset(NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
%         hm1 = scramble(hm1,'RR2');
%     elseif EstimOpt.Draws == 5
%         cprintf('*blue','Sobol '); cprintf('draws (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
%         hm1 = sobolset(NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
%     elseif EstimOpt.Draws == 6
%         cprintf('*blue','Sobol '); cprintf('draws with random linear scramble and random digital shift (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
%         hm1 = sobolset(NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
%         hm1 = scramble(hm1,'MatousekAffineOwen');
%     end
% 
%     err_mtx = net(hm1,EstimOpt.NP*EstimOpt.NRep); % this takes every point:
%     clear hm1;
%    
% end
% err_mtx = err_mtx';
GridMat = zeros(NVarA, EstimOpt.Grid);
err_mtx = zeros(NVarA, EstimOpt.NRep*EstimOpt.NP);
for i = 1:NVarA
    GridMat(i,:) = EstimOpt.Bounds(i,1):((EstimOpt.Bounds(i,2) - EstimOpt.Bounds(i,1))/(EstimOpt.Grid-1)):EstimOpt.Bounds(i,2);
    err_mtx(i,:) = randsample(GridMat(i,:),EstimOpt.NRep*EstimOpt.NP, true);
   % err_mtx(i,:) = EstimOpt.Bounds(i,1) + err_mtx(i,:)*(EstimOpt.Bounds(i,2) - EstimOpt.Bounds(i,1));
end

%% Display Options


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

INPUT.XXa = reshape(INPUT.Xa,EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP, NVarA);
INPUT.XXa = permute(INPUT.XXa, [1 3 2]);
INPUT.YY = reshape(INPUT.Y,EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP);


if isfield(EstimOpt, 'Drawskeep') && ~isempty(EstimOpt.Drawskeep) && EstimOpt.Drawskeep == 1
    Results.err = err_mtx;
end

%% Estimation

Q = LMXL_cond(INPUT.YY, INPUT.XXa, err_mtx, EstimOpt); % NP x NRep
Z = LMXL_vars(err_mtx, EstimOpt);

LLfun = @(B) LL_lmxl_MATlike(Q,Z,INPUT.W,EstimOpt,OptimOpt,B);

if EstimOpt.HessEstFix == 0
    [Results.bhat, LL, Results.exitf, Results.output, Results.g, Results.hess] = fminunc(LLfun, b0, OptimOpt);
else
    [Results.bhat, LL, Results.exitf, Results.output, Results.g] = fminunc(LLfun, b0, OptimOpt);
end


%% Output

% save tmp1
% return

Results.LL = -LL;

Results.P = evalProbs(Z, EstimOpt, Results.bhat);
Tmp =  reshape(err_mtx,NVarA,EstimOpt.NRep,EstimOpt.NP);
Results.B = Tmp(:,:,1);

Results.P_sort = zeros(NVarA, EstimOpt.NRep);
Results.B_sort = zeros(NVarA, EstimOpt.NRep);
for i = 1:NVarA
   [Results.B_sort(i,:),I] = sort(Results.B(i,:)); 
   Results.P_sort(i,:) = Results.P(I);
end
Results.P2_sort = sum(reshape(Results.P_sort, [NVarA, 10, EstimOpt.NRep/10] ),2);
Results.P2_sort = reshape(Results.P2_sort, [NVarA, EstimOpt.NRep/10]);

Results.B2_sort = mean(reshape(Results.B_sort, [NVarA, 10, EstimOpt.NRep/10] ),2);
Results.B2_sort = reshape(Results.B2_sort, [NVarA, EstimOpt.NRep/10]);


Results.Means = sum(Results.P(:,ones(NVarA,1))'.*Results.B,2);
Results.Stds = sqrt(sum(Results.P(:,ones(NVarA,1))'.*Results.B.^2,2) - Results.Means.^2);

disp(' ')
disp(['LL at convergence: ',num2str(Results.LL,'%8.4f')])
disp(' ');
clocknote = clock;
tocnote = toc;
[~,DayName] = weekday(now,'long');
disp(['Estimation completed on ' DayName ', ' num2str(clocknote(1)) '-' sprintf('%02.0f',clocknote(2)) '-' sprintf('%02.0f',clocknote(3)) ' at ' sprintf('%02.0f',clocknote(4)) ':' sprintf('%02.0f',clocknote(5)) ':' sprintf('%02.0f',clocknote(6))])
disp(['Estimation took ' num2str(tocnote) ' seconds ('  num2str(floor(tocnote/(60*60))) ' hours ' num2str(floor(rem(tocnote,60*60)/60)) ' minutes ' num2str(rem(tocnote,60)) ' seconds).']);
disp(' ');
Results.clocknote = clocknote;
Results.tocnote = clocknote;
