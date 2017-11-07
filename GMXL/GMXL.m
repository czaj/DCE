function Results = GMXL(INPUT,Results_old,EstimOpt,OptimOpt)

global B_backup;

% save tmp_GMXL
% return

% ModelOpt, EstimOpt.FullCov
% G_MNL_0:       b = [b_1,...,b_NVarA, l_1,...,l_NVarA, m_1,...,m_NVarA*NVarM, s_1,...,s_NVarS, t_1,...,t_NVarT, tau, gamma]
% G_MNL_1:       b = [b1,...,bNVarA, l_11,l_21,...,l_2NVarA,l_22,l_23,...,l_2NVarA,... ... ,l_NVarANVarA,  m_1,...,m_NVarA*NVarM, s_1,...,s_NVarS, t_1,...,t_NVarT, tau, gamma]

tic

Results.bhat = [];
Results.R = [];
Results.R_out = {};
Results.stats = [];


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
    cprintf('Black','Estimating '); cprintf('*Black','weighted '); cprintf('Black','GMXL model...\n');
else
    disp('Estimating GMXL model ...')
end
if isfield(EstimOpt,'FullCov') == 0
    EstimOpt.FullCov = 0;
end
if isfield(EstimOpt,'WTP_space') == 0
    EstimOpt.WTP_space = 0;
    EstimOpt.WTP_matrix = [];
elseif EstimOpt.WTP_space == 0
    EstimOpt.WTP_matrix = [];
end

if EstimOpt.WTP_space > 0
    error('Are you sure you want to estimate GMXL in WTP-space? We don''t think it makes sense ...')
end

if EstimOpt.FullCov == 0
    disp('with non-correlated random parameters in preference-space.')
else
    disp('with correlated random parameters in preference-space.')
end

if isfield(EstimOpt,'Dist') == 0 || isempty(EstimOpt.Dist)
    EstimOpt.Dist = zeros(1,EstimOpt.NVarA);
    EstimOpt.DistS = 1; % scale distributed log-normally
    disp('WARNING: distributions for random parameters not specified - assuming normality (scale assumed log-normal)')
else
    if length(EstimOpt.Dist) == 1
        EstimOpt.Dist = EstimOpt.Dist.*ones(1,EstimOpt.NVarA);
    elseif length(EstimOpt.Dist) == EstimOpt.NVarA
        EstimOpt.Dist = EstimOpt.Dist(:)';
    else
        error('Incorrect no. of random parameters'' distributions provided')
    end
end

if sum(EstimOpt.Dist == 3) > 0 || sum(EstimOpt.Dist == 4) > 0 || sum(EstimOpt.Dist > 5) > 0
    error('Incorrect distributions provided')
end

if isfield(EstimOpt,'DistS') == 0 || isempty(EstimOpt.DistS)
    EstimOpt.DistS = 1; % scale distributed log-normally
end

disp(['Random parameters distributions: ', num2str(EstimOpt.Dist),' (-1 - constant, 0 - normal, 1 - lognormal, 2 - spike, 5 - Weibull)'])

if EstimOpt.DistS == 1
    disp(['Scale parameter distribution: ', num2str(EstimOpt.DistS),' (-1 - constant, 0 - normal, 1 - lognormal, 2 - spike, 5 - Weibull)'])
else
    EstimOpt.DistS = 1;
    disp('WARNING: Distribution of scale set re-set to log-normal');
end

if isfield(INPUT,'Xs') == 0
    INPUT.Xs = zeros(size(INPUT.Y,1),0);
end

EstimOpt.NVarS = size(INPUT.Xs,2); % Number of covariates of scale

if isfield(INPUT,'Xm') == 0
    INPUT.Xm = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarM = size(INPUT.Xm,2); % Number of covariates of means of random parameters

if isfield(INPUT,'Xv') == 0 % This does not currently work
    INPUT.Xv = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarV = size(INPUT.Xv,2); % Number of covariates of variances of random parameters

if isfield(INPUT,'Xt') == 0
    INPUT.Xt = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarT = size(INPUT.Xt,2); % Number of covariates of tau

if isfield(EstimOpt,'tau0') == 0
    EstimOpt.tau0 = 1;
end

if isfield(EstimOpt,'Gamma0') == 0
    EstimOpt.Gamma0 = 0;
end

% if EstimOpt.WTP_space > 0
%     if isfield(EstimOpt,'WTP_matrix') == 0
%         WTP_att = (EstimOpt.NVarA-EstimOpt.WTP_space)/EstimOpt.WTP_space;
%         if rem(WTP_att,1) ~= 0
%             error('EstimOpt.WTP_matrix associating attributes with cost parameters not provided')
%         else
%             if EstimOpt.WTP_space > 1
%                 disp(['EstimOpt.WTP_matrix associating attributes with cost parameters not provided - assuming equal shares for each of the ',num2str(EstimOpt.WTP_space),' monetary attributes'])
%             end
%             EstimOpt.WTP_matrix = EstimOpt.NVarA - EstimOpt.WTP_space + kron(1:EstimOpt.WTP_space,ones(1,WTP_att));
%             %         tic; EstimOpt.WTP_matrix = 1:EstimOpt.WTP_space;...
%             %         EstimOpt.WTP_matrix = EstimOpt.WTP_matrix(floor((0:size(EstimOpt.WTP_matrix,2)*WTP_att-1)/WTP_att)+1); toc
%         end
%     elseif size(EstimOpt.WTP_matrix,2) ~= EstimOpt.NVarA - EstimOpt.WTP_space
%         error('Dimensions of EstimOpt.WTP_matrix not correct - for each non-monetary attribute provide no. of attribute to multiply it with')
%     else
%         EstimOpt.WTP_matrix = EstimOpt.WTP_matrix(:)';
%     end
% end

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
if EstimOpt.NVarS > 0
    if isfield(EstimOpt,'NamesS') == 0 || isempty(EstimOpt.NamesS) || length(EstimOpt.NamesS) ~= EstimOpt.NVarS
        EstimOpt.NamesS = (1:EstimOpt.NVarS)';
        EstimOpt.NamesS = cellstr(num2str(EstimOpt.NamesS));
    elseif size(EstimOpt.NamesS,1) ~= EstimOpt.NVarS
        EstimOpt.NamesS = EstimOpt.NamesS';
    end
end

if EstimOpt.NVarT > 0
    if isfield(EstimOpt,'NamesT') == 0 || isempty(EstimOpt.NamesT) || length(EstimOpt.NamesT) ~= EstimOpt.NVarT
        EstimOpt.NamesT = (1:EstimOpt.NVarT)';
        EstimOpt.NamesT = cellstr(num2str(EstimOpt.NamesT));
    elseif size(EstimOpt.NamesT,1) ~= EstimOpt.NVarT
        EstimOpt.NamesT = EstimOpt.NamesT';
    end
end

%% Starting values

if ~isfield(EstimOpt,'Gamma0')
    if  EstimOpt.FullCov == 0 && (~isfield(Results_old,'GMXL_d') || ~isfield(Results_old.GMXL_d,'b0'))
        EstimOpt.Gamma0 = 0;
    elseif EstimOpt.FullCov == 1 && (~isfield(Results_old,'GMXL') || ~isfield(Results_old.GMXL,'b0'))
        EstimOpt.Gamma0 = 0;
    end
end

if EstimOpt.FullCov == 0
    if exist('B_backup','var') && ~isempty(B_backup)
        if (EstimOpt.Gamma0 == 0 || EstimOpt.Gamma0 == 1)
            if size(B_backup,1) == EstimOpt.NVarA*2 + EstimOpt.NVarM*EstimOpt.NVarA + EstimOpt.NVarS + EstimOpt.NVarT + 1
                b0 = B_backup(:);
                disp('Using the starting values from Backup')
            end
        else
            if size(B_backup,1) == EstimOpt.NVarA*2 + EstimOpt.NVarM*EstimOpt.NVarA + EstimOpt.NVarS + EstimOpt.NVarT + 2
                b0 = B_backup(:);
                disp('Using the starting values from Backup')
            end
        end
    end
    if ~exist('b0','var') % There is no Backup
        if isfield(Results_old,'GMXL_d') && isfield(Results_old.GMXL_d,'b0') % starting values provided
            Results_old.GMXL_d.b0_old = Results_old.GMXL_d.b0(:);
            Results_old.GMXL_d = rmfield(Results_old.GMXL_d,'b0');
            if length(Results_old.GMXL_d.b0_old) ~=  EstimOpt.NVarA*2 + EstimOpt.NVarM*EstimOpt.NVarA + EstimOpt.NVarS + EstimOpt.NVarT + 2
                disp('WARNING: Incorrect no. of starting values or model specification')
                Results_old.GMXL_d = rmfield(Results_old.GMXL_d,'b0_old');
            else
                b0 = Results_old.GMXL_d.b0_old(:);
                if (EstimOpt.Gamma0 == 0 || EstimOpt.Gamma0 == 1)
                    b0 = b0(1:end-1);
                end
            end
        end
    end
    if ~exist('b0','var') % There is no Backup nor starting values provided
        if isfield(Results_old,'MXL_d') && isfield(Results_old.MXL_d,'bhat') %%isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat')
            disp('Using MXL_d results as starting values')
            Results_old.MXL_d.bhat = Results_old.MXL_d.bhat(:);
            b0 = [Results_old.MXL_d.bhat(1:end);zeros(EstimOpt.NVarT,1);sqrt(EstimOpt.tau0)];
        elseif isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat')
            disp('Using MNL results as starting values')
            b0 = [Results_old.MNL.bhat(1:EstimOpt.NVarA);max(1,sqrt(abs(Results_old.MNL.bhat(1:EstimOpt.NVarA))));0.1*ones(EstimOpt.NVarM.*EstimOpt.NVarA,1);Results_old.MNL.bhat(EstimOpt.NVarA+1:end);zeros(EstimOpt.NVarT,1);sqrt(EstimOpt.tau0)];
            if sum(EstimOpt.Dist == 1) > 0
                b0(EstimOpt.Dist == 1) = log(b0(EstimOpt.Dist == 1));
            end
        else
            error('No starting values available - run MNL or MXL_d first')
        end
        if (EstimOpt.Gamma0 ~= 0 && EstimOpt.Gamma0 ~= 1)
            b0 = [b0;log(EstimOpt.Gamma0./(1-EstimOpt.Gamma0))];
        end
    end

elseif EstimOpt.FullCov == 1

    if exist('B_backup','var') && ~isempty(B_backup) && ...
            ((EstimOpt.Gamma0 == 0 || EstimOpt.Gamma0 == 1) && (size(B_backup,1) == EstimOpt.NVarA + sum(1:EstimOpt.NVarA) + EstimOpt.NVarM*EstimOpt.NVarA + EstimOpt.NVarS + EstimOpt.NVarT + 1) || ...
            (EstimOpt.Gamma0 ~= 0 && EstimOpt.Gamma0 ~= 1 && size(B_backup,1) == EstimOpt.NVarA + sum(1:EstimOpt.NVarA) + EstimOpt.NVarM*EstimOpt.NVarA + EstimOpt.NVarS + EstimOpt.NVarT + 2))
        disp('Using the starting values from Backup')
        if (EstimOpt.Gamma0 == 0 || EstimOpt.Gamma0 == 1) && (size(B_backup,1) == EstimOpt.NVarA + sum(1:EstimOpt.NVarA) + EstimOpt.NVarM*EstimOpt.NVarA + EstimOpt.NVarS + EstimOpt.NVarT + 1)
            b0 = B_backup(:);
        elseif (EstimOpt.Gamma0 ~= 0 && EstimOpt.Gamma0 ~= 1) && (size(B_backup,1) == EstimOpt.NVarA + sum(1:EstimOpt.NVarA) + EstimOpt.NVarM*EstimOpt.NVarA + EstimOpt.NVarS + EstimOpt.NVarT + 2)
            b0 = B_backup(:);
        end
    elseif isfield(Results_old,'GMXL') && isfield(Results_old.GMXL,'b0') % starting values provided
        Results_old.GMXL.b0_old = Results_old.GMXL.b0(:);
        Results_old.GMXL = rmfield(Results_old.GMXL,'b0');
        if length(Results_old.GMXL.b0_old) ~=  EstimOpt.NVarA + sum(1:EstimOpt.NVarA) + EstimOpt.NVarM*EstimOpt.NVarA + EstimOpt.NVarS + EstimOpt.NVarT + 2
            disp('WARNING: Incorrect no. of starting values or model specification')
            Results_old.GMXL = rmfield(Results_old.GMXL,'b0_old');
        else
            b0 = Results_old.GMXL.b0_old(:);
            if (EstimOpt.Gamma0 == 0 || EstimOpt.Gamma0 == 1)
                b0 = b0(1:end-1);
            end
        end
    end
    if ~exist('b0','var') % There is no Backup nor starting values provided
        if isfield(Results_old,'MXL') && isfield(Results_old.MXL,'bhat')
            disp('Using MXL results as starting values')
            Results_old.MXL.bhat = Results_old.MXL.bhat(:);
            b0 = [Results_old.MXL.bhat;zeros(EstimOpt.NVarT,1);sqrt(EstimOpt.tau0)];
            if (EstimOpt.Gamma0 ~= 0 && EstimOpt.Gamma0 ~= 1)
                b0 = [b0;log(EstimOpt.Gamma0./(1-EstimOpt.Gamma0))];
            end
        elseif isfield(Results_old,'GMXL_d') && isfield(Results_old.GMXL_d,'bhat')
            disp('Using GMXL_d results as starting values')
            Results_old.GMXL_d.bhat = Results_old.GMXL_d.bhat(:);
            vc_tmp = diag(Results_old.GMXL_d.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2));
            b0 = [Results_old.GMXL_d.bhat(1:EstimOpt.NVarA);vc_tmp(tril(ones(size(vc_tmp))) == 1);Results_old.GMXL_d.bhat(EstimOpt.NVarA*2+1:end)];
        elseif isfield(Results_old,'MXL_d') && isfield(Results_old.MXL_d,'bhat')
            disp('Using MXL_d results as starting values')
            Results_old.MXL_d.bhat = Results_old.MXL_d.bhat(:);
            vc_tmp = diag(Results_old.MXL_d.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2));
            b0 = [Results_old.MXL_d.bhat(1:EstimOpt.NVarA);vc_tmp(tril(ones(size(vc_tmp))) == 1);Results_old.MXL_d.bhat(EstimOpt.NVarA*2+1:end);zeros(EstimOpt.NVarT,1);sqrt(EstimOpt.tau0)];
            if (EstimOpt.Gamma0 ~= 0 && EstimOpt.Gamma0 ~= 1)
                b0 = [b0;log(EstimOpt.Gamma0./(1-EstimOpt.Gamma0))];
            end
        elseif isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat')
            disp('Using MNL results as starting values')
            %             b0 = [Results_old.MNL.bhat(1:EstimOpt.NVarA);zeros(EstimOpt.NVarA*(EstimOpt.NVarM)+ncv(EstimOpt.NVarA),1);Results_old.MNL.bhat(EstimOpt.NVarA+1:end); zeros(EstimOpt.NVarT,1);EstimOpt.tau0; EstimOpt.Gamma0];
            b0 = [Results_old.MNL.bhat(1:EstimOpt.NVarA);0.1*ones(EstimOpt.NVarA*(EstimOpt.NVarM)+sum(1:EstimOpt.NVarA),1);0.1*ones(EstimOpt.NVarM.*EstimOpt.NVarA,1);Results_old.MNL.bhat(EstimOpt.NVarA+1:end);0.1*ones(EstimOpt.NVarT,1);sqrt(EstimOpt.tau0)];
            if sum(EstimOpt.Dist == 1) > 0
                b0(EstimOpt.Dist == 1) = log(b0(EstimOpt.Dist == 1));
            end
            if (EstimOpt.Gamma0 ~= 0 && EstimOpt.Gamma0 ~= 1)
                b0 = [b0;log(EstimOpt.Gamma0./(1-EstimOpt.Gamma0))];
            end
        else
            error('No starting values available')
        end
    end
end

if EstimOpt.Gamma0 == 0
    disp('GMXL type II - variance of residual taste heterogeneity fully scaled (gamma = 0)');
elseif EstimOpt.Gamma0 == 1
    disp('GMXL type I - variance of residual taste heterogeneity not scaled (gamma = 1)');
end

%if EstimOpt.Dist(1) == 5 && b0(end-1) == 0
%   disp('Tau cannot equal 0 for Weibull distributed scale - using tau=1 instead');
%   b0(end-1) = 1;
%end

%b0(end-1) = log(b0(end-1)); % tau entering the LL-function as exp(tau) later
%if isfield(EstimOpt,'Gamma0') == 1 && (EstimOpt.Gamma0 == 0 || EstimOpt.Gamma0 == 1)
%    b0 = b0(1:end-1);
%else
%    b0(end) = log(b0(end)./(1-b0(end)));
%end

%% Optimization Options

if  isfield(EstimOpt,'BActive')
    EstimOpt.BActive = EstimOpt.BActive(:)';
end

if sum(EstimOpt.Dist == -1) > 0
    if isfield(EstimOpt,'BActive') == 0 || isempty(EstimOpt.BActive)
        EstimOpt.BActive = ones(1,length(b0));
    end
    if EstimOpt.FullCov == 0
        EstimOpt.BActive(EstimOpt.NVarA+find(EstimOpt.Dist == -1)) = 0;
    elseif EstimOpt.FullCov == 1
        Vt = tril(ones(EstimOpt.NVarA));
        Vt(EstimOpt.Dist == -1,:) = 0;
        %         EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA)) = EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA)) .* (Vt(find(tril(ones(size(Vt)))))');
        EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA)) = EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA)).*(Vt(tril(ones(size(Vt))) ~= 0)');
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
    err_mtx = randn(EstimOpt.NP*EstimOpt.NRep,1+EstimOpt.NVarA); % to be cut down later
elseif EstimOpt.Draws == 2 % LHS
    cprintf('*blue','Latin Hypercube Sampling '); cprintf('draws \n');
    err_mtx = lhsnorm(zeros((1+EstimOpt.NVarA)*EstimOpt.NP,1),diag(ones((1+EstimOpt.NVarA)*EstimOpt.NP,1)),EstimOpt.NRep);
    err_mtx = reshape(err_mtx,[EstimOpt.NRep*EstimOpt.NP,1+EstimOpt.NVarA]);
elseif EstimOpt.Draws >= 3 % Quasi random draws
    if EstimOpt.Draws == 3
        cprintf('*blue','Halton '); cprintf(['draws (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n'])
        hm1 = haltonset(1+EstimOpt.NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
    elseif EstimOpt.Draws == 4 % apply reverse-radix scrambling
        cprintf('*blue','Halton '); cprintf(['draws with reverse radix scrambling (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n'])
        hm1 = haltonset(1+EstimOpt.NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
        hm1 = scramble(hm1,'RR2');
    elseif EstimOpt.Draws == 5
        cprintf('*blue','Sobol '); cprintf(['draws (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n'])
        hm1 = sobolset(1+EstimOpt.NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
    elseif EstimOpt.Draws == 6
        cprintf('*blue','Sobol '); cprintf(['draws with random linear scramble and random digital shift (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n'])
        hm1 = sobolset(1 + EstimOpt.NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
        hm1 = scramble(hm1,'MatousekAffineOwen');
    end

    err_mtx = net(hm1,EstimOpt.NP*EstimOpt.NRep); % this takes every point:
    clear hm1;

    if EstimOpt.NP*EstimOpt.NRep < 3e+7
        err_mtx(:,EstimOpt.Dist < 5) = icdf('Normal',err_mtx(:,EstimOpt.Dist < 5),0,1); % to be cut down later
        err_mtx(:,EstimOpt.Dist == 5) = -log(1-err_mtx(:,EstimOpt.Dist == 5));

    else % this is for very large number of draws * variables
        for i = 1:EstimOpt.NVarA
            if EstimOpt.Dist(i) < 5
                err_mtx(:,i) = icdf('Normal',err_mtx(:,i),0,1); %to be cut down later
            else
                err_mtx(:,i) = -log(1-err_mtx(:,i)); %to be cut down later
            end
        end
    end

    err_mtx(:,EstimOpt.Dist == -1) = 0;

end

if isfield(EstimOpt,'Drawskeep') && ~isempty(EstimOpt.Drawskeep) && EstimOpt.Drawskeep == 1
    Results.err = err_mtx;
end


%% Display Options

% if EstimOpt.NumGrad == 0
%    EstimOpt.NumGrad = 1;
%    OptimOpt.GradObj = 'off';
%    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - GMXL not supported by analytical gradient \n')
% end

% if EstimOpt.NumGrad == 0 && EstimOpt.WTP_space > 0
%     EstimOpt.NumGrad = 1;
%     OptimOpt.GradObj = 'off';
%     cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient to numerical - WTP-space not supported by analytical gradient \n')
% end
if EstimOpt.Gamma0 ~= 0 && EstimOpt.NumGrad == 0 && any(EstimOpt.Dist > 0)
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient to numerical - Non-normal random parameters not supported by analytical gradient for GMXL and GMXL of type I \n')
end

if EstimOpt.Gamma0 == 0 && EstimOpt.NumGrad == 0 && any(EstimOpt.Dist > 1)
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient to numerical - Non-normal or log-normal random parameters not supported by analytical gradient for GMXL of type II \n')
end

if EstimOpt.NumGrad == 0 && EstimOpt.DistS > 1
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient to numerical - Weibull Scale not supported by analytical gradient \n')
end

% if EstimOpt.NumGrad == 0 && (any(EstimOpt.NCTMiss ~= EstimOpt.NCT) || any(EstimOpt.NAltMiss ~= EstimOpt.NAlt))
%     EstimOpt.NumGrad = 1;
%     cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - Missing observations not supported by analytical gradient \n')
% end

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

%% Rescructure data

INPUT.YYY = reshape(INPUT.Y',[EstimOpt.NAlt,EstimOpt.NCT,EstimOpt.NP]);
idx = sum(reshape(INPUT.MissingInd,[EstimOpt.NAlt,EstimOpt.NCT,EstimOpt.NP])) == EstimOpt.NAlt;
INPUT.YYY(idx(ones(EstimOpt.NAlt,1),:,:)) = NaN; % replace YYY in missing choice-tasks with NaN
INPUT.YY = reshape(INPUT.YYY,[EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP]);
INPUT.Xa(INPUT.MissingInd == 1,:) = NaN;

XXa_tmp = reshape(INPUT.Xa',[EstimOpt.NVarA, EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP]);

INPUT.XXm = reshape(INPUT.Xm',[EstimOpt.NVarM, EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP]);
INPUT.XXm = squeeze(INPUT.XXm(:,1,:));
if EstimOpt.NVarM == 1
    INPUT.XXm = INPUT.XXm'; %NP x NVarM
end

INPUT.XXt = INPUT.Xt(1:EstimOpt.NCT*EstimOpt.NAlt:end,:);

for i = 1:EstimOpt.NP
    INPUT.XXa(:,:,i) = XXa_tmp(:,:,i)';
end

err_sliced = err_mtx';

EstimOpt.indx1 = [];
EstimOpt.indx2 = [];
if EstimOpt.NumGrad == 0 && EstimOpt.FullCov == 1
    for i = 1:EstimOpt.NVarA
        EstimOpt.indx1 = [EstimOpt.indx1,i:EstimOpt.NVarA];
        EstimOpt.indx2 = [EstimOpt.indx2,i*ones(1,EstimOpt.NVarA+1-i)];
    end
end

EstimOpt.indx3 = [];
if EstimOpt.NVarT > 0
    for i = 1:EstimOpt.NVarT
        EstimOpt.indx3 = [EstimOpt.indx3,1:EstimOpt.NVarA];
    end
end

%% Estimation

LLfun = @(B) LL_gmxl_MATlike(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,INPUT.XXt,err_sliced,INPUT.W,EstimOpt,OptimOpt,B);

if EstimOpt.ConstVarActive == 0
    if EstimOpt.HessEstFix == 0
        [Results.bhat,LL,Results.exitf,Results.output,Results.g,Results.hess] = fminunc(LLfun,b0,OptimOpt);
    else
        [Results.bhat,LL,Results.exitf,Results.output,Results.g] = fminunc(LLfun,b0,OptimOpt);
    end
elseif EstimOpt.ConstVarActive == 1 % equality constraints

    EstimOpt.CONS1 = diag(1-EstimOpt.BActive);
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
Results.b0 = b0;

Results.LLdetailed = LL_gmxl(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,INPUT.XXt,err_sliced,EstimOpt,Results.bhat);
Results.LLdetailed = Results.LLdetailed.*INPUT.W;
if any(INPUT.MissingInd == 1) % In case of some missing data
    idx = sum(reshape(INPUT.MissingInd,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP)) == EstimOpt.NAlt;
    idx = sum(reshape(idx,EstimOpt.NCT,EstimOpt.NP),1)'; % no. of missing NCT for every respondent
    idx = EstimOpt.NCT - idx;
    R2 = mean(exp(-Results.LLdetailed./idx),1);
else
    R2 = mean(exp(-Results.LLdetailed/EstimOpt.NCT),1);
end

if EstimOpt.HessEstFix == 1
    f = LL_gmxl(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,INPUT.XXt,err_sliced,EstimOpt,Results.bhat);
    Results.jacobian = numdiff(@(B) INPUT.W.*LL_gmxl(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,INPUT.XXt,err_sliced,EstimOpt,B),INPUT.W.*f,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) INPUT.W.*LL_gmxl(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,INPUT.XXt,err_sliced,EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.hess = hessian(@(B) sum(INPUT.W.*LL_gmxl(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,INPUT.XXt,err_sliced,EstimOpt,B),1),Results.bhat);
end

% if sum(EstimOpt.BActive == 0) > 0
%     if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
%         Results.jacobian = Results.jacobian(:, EstimOpt.BActive == 1);
%         Results.hess = Results.jacobian'*Results.jacobian;
%     elseif EstimOpt.HessEstFix == 0 || EstimOpt.HessEstFix == 3
%         Results.hess = Results.hess(EstimOpt.BActive == 1,EstimOpt.BActive == 1);
%     end
%     Results.ihess = inv(Results.hess);
%     Results.ihess = direcXpnd(Results.ihess,EstimOpt.BActive);
%     Results.ihess = direcXpnd(Results.ihess',EstimOpt.BActive);
%     Results.std = sqrt(diag(Results.ihess));
% 	Results.std(EstimOpt.BActive == 0) = NaN;
% else
% 	if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
%         Results.hess = Results.jacobian'*Results.jacobian;
%     end
%     Results.ihess = full(inv(sparse(Results.hess)));
%     Results.std = sqrt(diag(Results.ihess));
% end

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
        [~, Results.jacobian] = LL_gmxl(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,INPUT.XXt,err_sliced,EstimOpt,Results.bhat);
        Results.jacobian = Results.jacobian.*INPUT.W;
    else
        Results.jacobian = numdiff(@(B) INPUT.W.*LL_gmxl(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,INPUT.XXt,err_sliced,EstimOpt,B),Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
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
    Results.DetailsA(:,3:4) = [Results.std(1:EstimOpt.NVarA),pv(Results.bhat(1:EstimOpt.NVarA),Results.std(1:EstimOpt.NVarA))];
    Results.DetailsV(:,1) = abs(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2));
    Results.DetailsV(:,3:4) = [Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*2),pv(abs(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2)),Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*2))];
    Results.R = [Results.DetailsA,Results.DetailsV];
    if EstimOpt.NVarM > 0
        Results.DetailsM = [];
        for i=1:EstimOpt.NVarM
            Results.DetailsM(:,4*i-3) = Results.bhat(EstimOpt.NVarA*(2+i-1)+1:EstimOpt.NVarA*(2+i));
            Results.DetailsM(:,4*i-1:4*i) = [Results.std(EstimOpt.NVarA*(2+i-1)+1:EstimOpt.NVarA*(2+i)),pv(Results.bhat(EstimOpt.NVarA*(2+i-1)+1:EstimOpt.NVarA*(2+i)),Results.std(EstimOpt.NVarA*(2+i-1)+1:EstimOpt.NVarA*(2+i)))];
        end
        Results.R = [Results.R,Results.DetailsM];
    end
    if EstimOpt.NVarS > 0
        Results.DetailsS = [];
        Results.DetailsS(:,1) = Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NVarM)+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS);
        Results.DetailsS(:,3:4) = [Results.std(EstimOpt.NVarA*(2+EstimOpt.NVarM)+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS),pv(Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NVarM)+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS),Results.std(EstimOpt.NVarA*(2+EstimOpt.NVarM)+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS))];
    end
    if EstimOpt.NVarT > 0
        Results.DetailsT = [];
        Results.DetailsT(:,1) = Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS+EstimOpt.NVarT);
        Results.DetailsT(:,3:4) = [Results.std(EstimOpt.NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS+EstimOpt.NVarT),pv(Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS+EstimOpt.NVarT),Results.std(EstimOpt.NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS+1:EstimOpt.NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS+EstimOpt.NVarT))];
    end
    if isfield(EstimOpt,'Gamma0') && (EstimOpt.Gamma0 == 0 || EstimOpt.Gamma0 == 1) % transfer gamma* back to gamma
        Results.bhat = [Results.bhat;EstimOpt.Gamma0];
        Results.std = [Results.std;NaN];
        Results.DetailsGamma(1,1) = Results.bhat(end);
        Results.DetailsGamma(1,3:4) = [Results.std(end),NaN];
    else
        Results.DetailsGamma(1,1) = exp(Results.bhat(end))./(1+exp(Results.bhat(end)));
        Results.DetailsGamma(1,3:4) = [Results.std(end).*exp(Results.bhat(end))./(1+exp(Results.bhat(end))).^2,pv(exp(Results.bhat(end))./(1+exp(Results.bhat(end))),Results.std(end).*exp(Results.bhat(end))./(1+exp(Results.bhat(end))).^2)];
    end
    Results.DetailsTau(1,1) = exp(Results.bhat(end-1));
    Results.DetailsTau(1,3:4) = [Results.std(end-1).*exp(Results.bhat(end-1)),pv(exp(Results.bhat(end-1)),Results.std(end-1).*exp(Results.bhat(end-1)))];
    Results.R = [Results.R,NaN(size(Results.R,1),4)];
    Results.R(1:2,end-3:end) = [Results.DetailsTau;Results.DetailsGamma];
    % Hessian is calculated for the underlying form of bhat

else % => EstimOpt.FullCov == 1
    Results.DetailsA(:,1) = Results.bhat(1:EstimOpt.NVarA);
    Results.DetailsA(:,3:4) = [Results.std(1:EstimOpt.NVarA),pv(Results.bhat(1:EstimOpt.NVarA),Results.std(1:EstimOpt.NVarA))];
    Results.DetailsV = sdtri(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2),Results.ihess(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2,EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2),EstimOpt);
    Results.DetailsV = [Results.DetailsV(:,1),zeros(EstimOpt.NVarA,1),Results.DetailsV(:,2:3)];
    Results.R = [Results.DetailsA,Results.DetailsV];
    if EstimOpt.NVarM > 0
        Results.DetailsM = [];
        for i=1:EstimOpt.NVarM
            Results.DetailsM(:,4*i-3) = Results.bhat(EstimOpt.NVarA*(EstimOpt.NVarA/2+0.5+i)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+i));
            Results.DetailsM(:,4*i-1:4*i) = [Results.std(EstimOpt.NVarA*(EstimOpt.NVarA/2+0.5+i)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+i)),pv(Results.bhat(EstimOpt.NVarA*(EstimOpt.NVarA/2+0.5+i)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+i)),Results.std(EstimOpt.NVarA*(EstimOpt.NVarA/2+0.5+i)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+i)))];
        end
        Results.R = [Results.R, Results.DetailsM];
    end
    if EstimOpt.NVarS > 0
        Results.DetailsS = [];
        Results.DetailsS(:,1) = Results.bhat(EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+EstimOpt.NVarS);
        Results.DetailsS(:,3:4) = [Results.std(EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+EstimOpt.NVarS),pv(Results.bhat(EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+EstimOpt.NVarS),Results.std(EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+EstimOpt.NVarS))];
    end
    if EstimOpt.NVarT > 0
        Results.DetailsT = [];
        Results.DetailsT(:,1) = Results.bhat(EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+EstimOpt.NVarS+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+EstimOpt.NVarS+EstimOpt.NVarT);
        Results.DetailsT(:,3:4) = [Results.std(EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+EstimOpt.NVarS+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+EstimOpt.NVarS+EstimOpt.NVarT),pv(Results.bhat(EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+EstimOpt.NVarS+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+EstimOpt.NVarS+EstimOpt.NVarT),Results.std(EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+EstimOpt.NVarS+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+EstimOpt.NVarS+EstimOpt.NVarT))];
    end
    Results.chol = [Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5)),Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5)),pv(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5)),Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5)))];
    if isfield(EstimOpt,'Gamma0') && (EstimOpt.Gamma0 == 0 || EstimOpt.Gamma0 == 1) % transfer gamma* back to gamma
        Results.bhat = [Results.bhat;EstimOpt.Gamma0];
        Results.std = [Results.std;NaN];
        Results.DetailsGamma(1,1) = Results.bhat(end);
        Results.DetailsGamma(1,3:4) = [Results.std(end),NaN];
    else
        Results.DetailsGamma(1,1) = exp(Results.bhat(end))./(1+exp(Results.bhat(end)));
        Results.DetailsGamma(1,3:4) = [Results.std(end).*exp(Results.bhat(end))./(1+exp(Results.bhat(end))).^2,pv(exp(Results.bhat(end))./(1+exp(Results.bhat(end))),Results.std(end).*exp(Results.bhat(end))./(1+exp(Results.bhat(end))).^2)];
    end
    Results.DetailsTau(1,1) = exp(Results.bhat(end-1));
    Results.DetailsTau(1,3:4) = [Results.std(end-1).*exp(Results.bhat(end-1)),pv(exp(Results.bhat(end-1)),Results.std(end-1).*exp(Results.bhat(end-1)))];
    Results.R = [Results.R,NaN(size(Results.R,1),4)];
    Results.R(1:2,end-3:end) = [Results.DetailsTau;Results.DetailsGamma];
end

EstimOpt.params = length(b0);
if isfield(EstimOpt,'BActive')
    EstimOpt.params = EstimOpt.params - sum(EstimOpt.BActive == 0);
end

Results.stats = [Results.LL;Results_old.MNL0.LL;1 - Results.LL/Results_old.MNL0.LL;R2;((2*EstimOpt.params - 2*Results.LL))/EstimOpt.NObs;((log(EstimOpt.NObs)*EstimOpt.params - 2*Results.LL))/EstimOpt.NObs;EstimOpt.NObs;EstimOpt.NP;EstimOpt.params];

Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;
Results.INPUT = INPUT;

Template1 = {'DetailsA','DetailsV'};
Template2 = {'DetailsA','DetailsV'};
Names.DetailsA = EstimOpt.NamesA;
Heads.DetailsA = {'Means';'tc'};
Heads.DetailsV = {'Standard Deviations';'lb'};
ST = {};

if EstimOpt.NVarM > 0
    Template1 = [Template1,'DetailsM'];
    Temp = cell(1,size(Template2,2));
    Temp(1,1) = {'DetailsM'};
    Template2 = [Template2;Temp];
    Heads.DetailsM(:,2) = [EstimOpt.NamesM;{'lb'}];
    Heads.DetailsM(1:2,1) = {'Interactions of means';'lc'};
end

if EstimOpt.NVarS > 0
    Temp = cell(1,size(Template1,2));
    Temp(1,1) = {'DetailsS'};
    Template1 = [Template1;Temp];
    Temp = cell(1,size(Template2,2));
    Temp(1,1) = {'DetailsS'};
    Template2 = [Template2;Temp];
    Names.DetailsS = EstimOpt.NamesS;
    Heads.DetailsS = {'Covariates of Scale';'tb'};
    ST = [ST,{'DetailsS'}];
end

if EstimOpt.NVarT > 0
    Temp = cell(1,size(Template1,2));
    Temp(1,1) = {'DetailsT'};
    Template1 = [Template1;Temp];
    Temp = cell(1,size(Template2,2));
    Temp(1,1) = {'DetailsT'};
    Template2 = [Template2;Temp];
    Names.DetailsT = EstimOpt.NamesT;
    Heads.DetailsT= {'Covariates of tau';'tb'};
    ST = [ST,{'DetailsT'}];
end

Results.DetailsGT = [Results.DetailsGamma;Results.DetailsTau];
Temp = cell(1,size(Template1,2));
Temp(1,1) = {'DetailsGT'};
Template1 = [Template1;Temp];
Temp = cell(1,size(Template2,2));
Temp(1,1) = {'DetailsGT'};
Template2 = [Template2;Temp];
Names.DetailsGT = {'gamma';'tau'};
Heads.DetailsGT = {'GMXL parameters';'tb'};
ST = [ST,{'DetailsGT'}];

Head = cell(1,2);
if EstimOpt.FullCov == 0
    Head(1,1) = {'GMXL_d'};
else
    Head(1,1) = {'GMXL'};
end

% if EstimOpt.WTP_space > 0
%     Head(1,2) = {'in WTP-space'};
% else
%     Head(1,2) = {'in preference-space'};
% end
Head(1,2) = {'in preference-space'};
Tail = cell(17,2);
Tail(2,1) = {'Model diagnostics'};
Tail(3:17,1) = {'LL at convergence';'LL at constant(s) only';strcat('McFadden''s pseudo-R',char(178));strcat('Ben-Akiva-Lerman''s pseudo-R',char(178));'AIC/n';'BIC/n';'n (observations)';'r (respondents)';'k (parameters)';' ';'Estimation method';'Simulation with';'Optimization method';'Gradient';'Hessian'};

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

if EstimOpt.Display ~= 0
    Results.Dist = EstimOpt.Dist;
    Results.R_out = genOutput(EstimOpt,Results,Head,Tail,Names,Template1,Template2,Heads,ST);
end

end
