function Results = MXL2(INPUT,Results_old,EstimOpt,OptimOpt)


% save tmp_MXL
% return

global B_backup

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
    cprintf('Black','Estimating '); cprintf('*Black','weighted '); cprintf('Black','MXL model...\n');
else
    disp('Estimating MXL model ...')
end

cprintf(rgb('DarkOrange'), 'Model is estimated on differences of attributes values \n')

if isfield(EstimOpt,'FullCov') == 0;
    EstimOpt.FullCov = 0;
end
if ~isfield(EstimOpt,'WTP_space')
    EstimOpt.WTP_space = 0;
    EstimOpt.WTP_matrix = [];
elseif EstimOpt.WTP_space == 0;
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
        EstimOpt.Dist = zeros(1,EstimOpt.NVarA+1);
        EstimOpt.Dist(1) = 1; % scale distributed log-normally (does not matter for MXL)
    if EstimOpt.WTP_space == 0
        cprintf(rgb('DarkOrange'), 'WARNING: distributions for random parameters not specified - assuming normality \n')
    else
        cprintf(rgb('DarkOrange'), 'WARNING: distributions for random parameters not specified - assuming normality (monetary parameters assumed log-normal) \n')
        EstimOpt.Dist(end-EstimOpt.WTP_space+1:end) = 1; % cost in WTP-space models log-normally distributed
    end
else    
    if length(EstimOpt.Dist) == 1
        EstimOpt.Dist = EstimOpt.Dist.*ones(1,EstimOpt.NVarA+1);
    elseif length(EstimOpt.Dist) == 1 + EstimOpt.NVarA
        EstimOpt.Dist = EstimOpt.Dist(:)';
    else
        error('Incorrect no. of random parameters'' distributions provided')
    end
end
if isfield(EstimOpt, 'Triang') == 0 || length(EstimOpt.Triang) ~= sum(EstimOpt.Dist(2:end) == 3,2) % Needed only if any parameter has triangular distribution
    EstimOpt.Triang = zeros(1, sum(EstimOpt.Dist(2:end) == 3,2));
elseif length(EstimOpt.Triang) == 1
    EstimOpt.Triang = EstimOpt.Triang*ones(1, sum(EstimOpt.Dist(2:end) == 3,2));
else
    EstimOpt.Triang = EstimOpt.Triang(:)';
end

EstimOpt.Johnson = sum(EstimOpt.Dist(2:end) >= 5);

if (sum(EstimOpt.Dist(2:end) >=3 & EstimOpt.Dist(2:end) <=5) > 0 && any(find(EstimOpt.Dist(2:end) >=3 & EstimOpt.Dist(2:end) <=5) > sum(EstimOpt.Dist(2:end) >=3 & EstimOpt.Dist(2:end) <=5))) && EstimOpt.FullCov == 1
    cprintf(rgb('DarkOrange'), 'WARNING: It is recommended to put variables with random parameters with Triangular/Weibull/Sinh-Arcsinh distribution first \n')
end

disp(['Random parameters distributions: ', num2str(EstimOpt.Dist(2:end)),' (-1 - constant, 0 - normal, 1 - lognormal, 2 - spike, 3 - Triangular, 4  - Weibull, 5 - Sinh-Arcsinh, 6 - Johnson Sb, 7 - Johnson Su)'])

if EstimOpt.WTP_space > 0 && sum(EstimOpt.Dist(end-EstimOpt.WTP_space+1:end)==1) > 0 && any(mean(INPUT.Xa(:,end-EstimOpt.WTP_space+1:end)) >= 0)
    cprintf(rgb('DarkOrange'), 'WARNING: Cost attributes with log-normally distributed parameters should enter utility function with a ''-'' sign \n')
end

if isfield(INPUT, 'Xs') == 0
    INPUT.Xs = zeros(size(INPUT.Y,1),0);
end 
EstimOpt.NVarS = size(INPUT.Xs,2); % Number of covariates of scale
if isfield(INPUT, 'Xm') == 0
    INPUT.Xm = zeros(size(INPUT.Y,1),0);
end 
EstimOpt.NVarM = size(INPUT.Xm,2); % Number of covariates of means of random parameters

% This does not currently work:
if isfield(INPUT, 'Xv') == 0 
    INPUT.Xv = zeros(size(INPUT.Y,1),0);
end 
EstimOpt.NVarV = size(INPUT.Xv,2); % Number of covariates of variances of random parameters

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
%         tic; EstimOpt.WTP_matrix = 1:EstimOpt.WTP_space;...
%         EstimOpt.WTP_matrix = EstimOpt.WTP_matrix(floor((0:size(EstimOpt.WTP_matrix,2)*WTP_att-1)/WTP_att)+1); toc
        end
%     elseif ~isequal(size(EstimOpt.WTP_matrix),[EstimOpt.NVarA-EstimOpt.WTP_space,EstimOpt.WTP_space])
	elseif size(EstimOpt.WTP_matrix,2) ~= EstimOpt.NVarA - EstimOpt.WTP_space
        error('Dimensions of EstimOpt.WTP_matrix not correct - for each non-monetary attribute provide no. of attribute to multiply it with')
    else
        EstimOpt.WTP_matrix = EstimOpt.WTP_matrix(:)';
	end
end

if isfield(EstimOpt,'Scores') == 0 || isempty(EstimOpt.Scores)
    EstimOpt.Scores = 0;
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
if EstimOpt.NVarS > 0
    if isfield(EstimOpt,'NamesS') == 0 || isempty(EstimOpt.NamesS) || length(EstimOpt.NamesS) ~= EstimOpt.NVarS
        EstimOpt.NamesS = (1:EstimOpt.NVarS)';
        EstimOpt.NamesS = cellstr(num2str(EstimOpt.NamesS));
	elseif size(EstimOpt.NamesS,1) ~= EstimOpt.NVarS
        EstimOpt.NamesS = EstimOpt.NamesS';
    end
end


%% Starting values


if EstimOpt.FullCov == 0
	if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVarA*2 + EstimOpt.NVarM*EstimOpt.NVarA + EstimOpt.NVarS + 2*EstimOpt.Johnson
            b0 = B_backup(:);
            disp('Using the starting values from Backup')
    elseif isfield(Results_old,'MXL_d') && isfield(Results_old.MXL_d,'b0') % starting values provided
        Results_old.MXL_d.b0_old = Results_old.MXL_d.b0(:);
        Results_old.MXL_d = rmfield(Results_old.MXL_d,'b0');
        if length(Results_old.MXL_d.b0_old) ~= EstimOpt.NVarA*2 + EstimOpt.NVarM*EstimOpt.NVarA + EstimOpt.NVarS + 2*EstimOpt.Johnson
            cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.MXL_d = rmfield(Results_old.MXL_d,'b0_old');
        else
            b0 = Results_old.MXL_d.b0_old(:);
        end
	end
    if  ~exist('b0','var')
        if isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat') && length(Results_old.MNL.bhat) == (EstimOpt.NVarA*(1+ EstimOpt.NVarM) + EstimOpt.NVarS + 2*EstimOpt.Johnson)
            disp('Using MNL results as starting values')
            Results_old.MNL.bhat = Results_old.MNL.bhat(:);
%             b0 = [Results_old.MNL.bhat(1:EstimOpt.NVarA);max(1,sqrt(abs(Results_old.MNL.bhat(1:EstimOpt.NVarA))));0.1*ones(EstimOpt.NVarM.*EstimOpt.NVarA,1);Results_old.MNL.bhat(EstimOpt.NVarA+1:end)];
            b0 = [Results_old.MNL.bhat(1:EstimOpt.NVarA);max(1,abs(Results_old.MNL.bhat(1:EstimOpt.NVarA)));Results_old.MNL.bhat(EstimOpt.NVarA+1:end)];
            if sum(EstimOpt.Dist(2:end)==1) > 0
                b0(EstimOpt.Dist(2:EstimOpt.NVarA+1) == 1) = log(b0(EstimOpt.Dist(2:EstimOpt.NVarA+1) == 1));
            end
            if sum(EstimOpt.Dist(2:end) == 3) > 0 % Triangular
               indx = find( EstimOpt.Dist(2:end) == 3);
               b0([indx; indx+EstimOpt.NVarA]) = [log(b0(indx)- EstimOpt.Triang'); log(b0(indx)- EstimOpt.Triang')];
            end
            if sum(EstimOpt.Dist(2:end) == 4) > 0 % Weibull
               indx = find( EstimOpt.Dist(2:end) == 4);
               b0([indx; indx+EstimOpt.NVarA]) = [log(b0(indx)); zeros(length(indx),1)];
            end
            if sum(EstimOpt.Dist(2:end) >=5) > 0 % Johnson
               indx = find( EstimOpt.Dist(2:end) >= 5);
               tmp = [b0(indx); log(b0(indx+EstimOpt.NVarA))];
               b0([indx; indx+EstimOpt.NVarA]) = [zeros(length(indx),1), ones(length(indx),1)];
               b0 = [b0; tmp];
            end
        else
            error('No starting values available - run MNL first')
        end
    end    
else % EstimOpt.FullCov == 1    
	if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVarA*(1+EstimOpt.NVarM) + sum(1:EstimOpt.NVarA) + EstimOpt.NVarS + 2*EstimOpt.Johnson
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'MXL') && isfield(Results_old.MXL,'b0') % starting values provided
        Results_old.MXL.b0_old = Results_old.MXL.b0(:);
        Results_old.MXL = rmfield(Results_old.MXL,'b0');
        if length(Results_old.MXL.b0_old) ~= EstimOpt.NVarA*(1+EstimOpt.NVarM) + sum(1:EstimOpt.NVarA) + EstimOpt.NVarS + 2*EstimOpt.Johnson
            cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.MXL = rmfield(Results_old.MXL,'b0_old');
        else
            b0 = Results_old.MXL.b0_old;            
        end
	end
	if  ~exist('b0','var')
        if isfield(Results_old,'MXL_d') && isfield(Results_old.MXL_d,'bhat') && length(Results_old.MXL_d.bhat) == ((2+EstimOpt.NVarM)*EstimOpt.NVarA + EstimOpt.NVarS + 2*EstimOpt.Johnson)
            disp('Using MXL_d results as starting values')
            Results_old.MXL_d.bhat = Results_old.MXL_d.bhat(:);
            if sum(EstimOpt.Dist(2:end) >= 3) > 0
                vc_tmp = Results_old.MXL_d.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2);
                vc_tmp(EstimOpt.Dist(2:end) < 3) = vc_tmp(EstimOpt.Dist(2:end) < 3).^2;
                vc_tmp = diag(vc_tmp);
            else
                vc_tmp = (diag(Results_old.MXL_d.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2))).^2;
            end
            b0 = [Results_old.MXL_d.bhat(1:EstimOpt.NVarA); vc_tmp(tril(ones(size(vc_tmp)))==1);Results_old.MXL_d.bhat(EstimOpt.NVarA*2+1:end)];
        elseif isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat')
            disp('Using MNL results as starting values')
            Results_old.MNL.bhat = Results_old.MNL.bhat(:);
            b0 = [Results_old.MNL.bhat(1:EstimOpt.NVarA);zeros(sum(1:EstimOpt.NVarA),1); Results_old.MNL.bhat(EstimOpt.NVarA+1:end)];
            if sum(EstimOpt.Dist(2:end)==1) > 0
                b0(EstimOpt.Dist(2:EstimOpt.NVarA+1) == 1) = log(b0_old(EstimOpt.Dist(2:EstimOpt.NVarA+1) == 1));
            end
            if sum(EstimOpt.Dist(2:end) == 3) > 0 % Triangular
               b0(EstimOpt.Dist(2:end) == 3) = log(b0(EstimOpt.Dist(2:end) == 3)- EstimOpt.Triang'); 
            end
            if sum(EstimOpt.Dist(2:end) == 4) > 0 % Weibull
               b0(EstimOpt.Dist(2:end) == 4) = log(b0(EstimOpt.Dist(2:end) == 4));
            end
            if sum(EstimOpt.Dist(2:end) >=5) > 0 % Johnson
               indx = find( EstimOpt.Dist(2:end) >= 5);
               tmp = b0(indx); 
               b0(indx) = zeros(length(indx),1);
               b0 = [b0; tmp; zeros(length(indx),1)];
            end
        else
            error('No starting values available')
        end
	end
end


%% Optimization Options


if  isfield(EstimOpt,'BActive')
	EstimOpt.BActive = EstimOpt.BActive(:)';
else
    EstimOpt.BActive = ones(1,length(b0));
end

if sum(EstimOpt.Dist == -1) > 0
    if isfield(EstimOpt,'BActive') == 0 || isempty(EstimOpt.BActive)
        EstimOpt.BActive = ones(1,length(b0));
    end
    if EstimOpt.FullCov == 0
        EstimOpt.BActive(EstimOpt.NVarA+find(EstimOpt.Dist(2:end) == -1)) = 0;
    elseif EstimOpt.FullCov == 1
        Vt = tril(ones(EstimOpt.NVarA));
        Vt(EstimOpt.Dist(2:end)==-1,:) = 0;
%         EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA)) = EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA)) .* (Vt(find(tril(ones(size(Vt)))))');
        EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA)) = EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA)) .* (Vt(tril(ones(size(Vt)))~=0)');
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
    err_mtx = randn(EstimOpt.NP*EstimOpt.NRep, EstimOpt.NVarA+1); %to be cut down later   
elseif EstimOpt.Draws == 2 % LHS
    cprintf('*blue','Latin Hypercube Sampling '); cprintf('draws \n');
    err_mtx=lhsnorm(zeros((EstimOpt.NVarA+1)*EstimOpt.NP,1),diag(ones((EstimOpt.NVarA+1)*EstimOpt.NP,1)),EstimOpt.NRep); 
    err_mtx = reshape(err_mtx, EstimOpt.NRep*EstimOpt.NP, EstimOpt.NVarA+1);
elseif EstimOpt.Draws >= 3 % Quasi random draws 
    if EstimOpt.Draws == 3
        cprintf('*blue','Halton '); cprintf('draws (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = haltonset(EstimOpt.NVarA+1,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); % 
    elseif EstimOpt.Draws == 4 % apply reverse-radix scrambling
        cprintf('*blue','Halton '); cprintf('draws with reverse radix scrambling (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = haltonset(EstimOpt.NVarA+1,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); % 
        hm1 = scramble(hm1,'RR2');
    elseif EstimOpt.Draws == 5
        cprintf('*blue','Sobol '); cprintf('draws (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = sobolset(EstimOpt.NVarA+1,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); 
    elseif EstimOpt.Draws == 6
        cprintf('*blue','Sobol '); cprintf('draws with random linear scramble and random digital shift (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = sobolset(EstimOpt.NVarA+1,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); 
        hm1 = scramble(hm1,'MatousekAffineOwen');
    end
    
    err_mtx = net(hm1,EstimOpt.NP*EstimOpt.NRep); % this takes every point:
    clear hm1;
    err_mtx = err_mtx(:,2:EstimOpt.NVarA+1); 
    
    if EstimOpt.NP*EstimOpt.NRep < 3e+7   
        err_mtx = icdf('Normal',err_mtx,0,1); %to be cut down later  
    else % this is for very large number of draws * variables
        for i = 1:EstimOpt.NVarA
            err_mtx(:,i) = icdf('Normal',err_mtx(:,i),0,1); %to be cut down later
        end        
    end
    
    err_mtx(:,EstimOpt.Dist(2:end) == -1) = 0;

end


%% Display Options

    
if ((isfield(EstimOpt, 'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    OptimOpt.GradObj = 'on';
end

if EstimOpt.NVarS > 0 && EstimOpt.NumGrad == 0 && any(isnan(INPUT.Xa(:)))
	EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
	cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient off - covariates of scale not supported by analytical gradient \n')
end

if any(EstimOpt.Dist(2:EstimOpt.NVarA+1) > 1) && EstimOpt.NumGrad == 0
	EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
	cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient off - analytical gradient available for normally or lognormally distributed parameters only \n')
end


if ((isfield(EstimOpt, 'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    OptimOpt.GradObj = 'on';
end

if (isfield(EstimOpt, 'ConstVarActive') == 0 || EstimOpt.ConstVarActive == 0) && isequal(OptimOpt.Algorithm,'quasi-newton') && isequal(OptimOpt.Hessian,'user-supplied')
	cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied Hessian off - quasi-newton algorithm does not use it anyway \n')
    OptimOpt.Hessian = 'off';
end

if  EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4
	cprintf(rgb('DarkOrange'), 'WARNING: Exact Hessian not available - switch to MXL (instead of MXL2) \n')
    EstimOpt.ApproxHess = 1;
    EstimOpt.HessEstFix = 0;
end

if EstimOpt.RobustStd == 1 && (EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2)
    EstimOpt.RobustStd = 0; 
    cprintf(rgb('DarkOrange'), 'WARNING: Setting off robust standard errors, they do not matter for BHHH aproximation of hessian \n')
end

if  any(EstimOpt.Dist(2:end)>= 3 & EstimOpt.Dist(2:end) <= 5) && EstimOpt.NVarM ~= 0
    error('Covariates of means do not work with triangular/weibull/sinh-arcsinh distributions')
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


Ytmp = prod(isnan(reshape(INPUT.Y, EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP))); % which CT are missing
INPUT.Y = reshape(INPUT.Y, EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP);
INPUT.Y(1, Ytmp == 1) = 1;
INPUT.Y = reshape(INPUT.Y, EstimOpt.NAlt*EstimOpt.NCT*EstimOpt.NP,1);
XaChosen = INPUT.Xa(INPUT.Y == 1 ,:);
XaNChosen = INPUT.Xa(INPUT.Y ~= 1,:);
XaChosen = reshape(XaChosen, 1, EstimOpt.NCT*EstimOpt.NP, EstimOpt.NVarA);
XaChosen = reshape(XaChosen(ones(EstimOpt.NAlt-1,1),:,:), (EstimOpt.NAlt-1)*EstimOpt.NCT*EstimOpt.NP, EstimOpt.NVarA);
INPUT.Xa = XaNChosen - XaChosen;
INPUT.XXa = reshape(INPUT.Xa,(EstimOpt.NAlt-1)*EstimOpt.NCT,EstimOpt.NP, EstimOpt.NVarA);
INPUT.XXa = permute(INPUT.XXa, [1 3 2]);
INPUT.Y = INPUT.Y(INPUT.Y~=1);
INPUT.YY = reshape(INPUT.Y,(EstimOpt.NAlt-1)*EstimOpt.NCT,EstimOpt.NP);



% idx = sum(reshape(INPUT.MissingInd,[EstimOpt.NAlt,EstimOpt.NCT,EstimOpt.NP])) == EstimOpt.NAlt;
% INPUT.YYY(idx(ones(EstimOpt.NAlt,1),:,:)) = NaN; % replace YYY in missing choice-tasks with NaN
% INPUT.YY = reshape(INPUT.YYY,EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP)==1;
%INPUT.YY = reshape(INPUT.YYY,EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP);


INPUT.XXm = reshape(INPUT.Xm',EstimOpt.NVarM, EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP);
INPUT.XXm = squeeze(INPUT.XXm(:,1,:));
if EstimOpt.NVarM == 1
    INPUT.XXm = INPUT.XXm';
end
    

err_mtx = err_mtx';
% change err_mtx from NRep*NP x NVarA to NP*NRep x NVarA (incrasing the no. of draws only adds new draws for each respondent, does not change all draws per individual)
% err_mtx = reshape(permute(reshape(err_mtx,EstimOpt.NP,EstimOpt.NRep,EstimOpt.NVarA),[2,1,3]),EstimOpt.NP*EstimOpt.NRep,EstimOpt.NVarA)';
% problem - look at the first NRep draws for NVarA=1... all are positive... 
if isfield(EstimOpt, 'Drawskeep') && ~isempty(EstimOpt.Drawskeep) && EstimOpt.Drawskeep == 1
    Results.err = err_mtx;
end

VC = tril(ones(EstimOpt.NVarA));
VC(VC == 1) = (1:(EstimOpt.NVarA*(EstimOpt.NVarA-1)/2+EstimOpt.NVarA))';
EstimOpt.DiagIndex = diag(VC);

% Creating indices for analitical gradient

EstimOpt.indx1 = [];
EstimOpt.indx2 = [];
if EstimOpt.NumGrad == 0 && EstimOpt.FullCov == 1
   for i = 1:EstimOpt.NVarA
      EstimOpt.indx1 = [EstimOpt.indx1, i:EstimOpt.NVarA];
      EstimOpt.indx2 = [EstimOpt.indx2, i*ones(1,EstimOpt.NVarA+1-i)];
   end
end

% save tmp2
% return

% if EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4; %calculations needed for analitical Hessian
%     EstimOpt.XXX = permute(mmx('square',permute(INPUT.XXa,[2,4,1,3]),[]),[3,1,2,4])
%     EstimOpt.XXX = zeros(EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NVarA,EstimOpt.NVarA, EstimOpt.NP);
%     if EstimOpt.FullCov == 0
%         EstimOpt.VCx = zeros(EstimOpt.NVarA,EstimOpt.NVarA,EstimOpt.NRep,EstimOpt.NP);
%     else
%         EstimOpt.VCx = zeros(EstimOpt.NVarA*(EstimOpt.NVarA-1)/2+EstimOpt.NVarA,EstimOpt.NVarA*(EstimOpt.NVarA-1)/2+EstimOpt.NVarA,EstimOpt.NRep,EstimOpt.NP);
%     end
%     err_tmp = reshape(err_mtx,EstimOpt.NVarA,EstimOpt.NRep,EstimOpt.NP); 
%     for i = 1:EstimOpt.NP
%         for j = 1:EstimOpt.NAlt*EstimOpt.NCT
%             EstimOpt.XXX(j,:,:,i) = (INPUT.XXa(j,:,i)')*INPUT.XXa(j,:,i);
%         end
%         for j = 1:EstimOpt.NRep
%             if EstimOpt.FullCov == 0
%                 EstimOpt.VCx(:,:,j,i) = err_tmp(:,j,i)*err_tmp(:,j,i)';
%             else
%                 EstimOpt.VCx(:,:,j,i) = err_tmp(EstimOpt.indx2,j,i)*err_tmp(EstimOpt.indx2,j,i)';
%             end
%         end
%     end

% end


%% Estimation


LLfun = @(B) LL_mxl2_MATlike(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,err_mtx,INPUT.W,EstimOpt,OptimOpt,B);
     
if EstimOpt.ConstVarActive == 0  
    
    if EstimOpt.HessEstFix == 0
        [Results.bhat, LL, Results.exitf, Results.output, Results.g, Results.hess] = fminunc(LLfun, b0, OptimOpt);
    else
        [Results.bhat, LL, Results.exitf, Results.output, Results.g] = fminunc(LLfun, b0, OptimOpt);
    end      

%     [x,fval,exitflag,output,lambda,grad,hessian] = knitromatlab(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,extendedFeatures,options,KNITROOptions)
%     [Results.bhat,LL,Results.exitf,Results.output,Results.lambda,Results.g,Results.hess] = knitromatlab(LLfun,b0,[],[],[],[],[],[],[],[],[],'knitro.opt'); %
%     [Results.bhat,LL,Results.exitf,Results.output,Results.lambda,Results.g,Results.hess] = knitromatlab(LLfun,b0,[],[],[],[],[],[],[],[],OptimOpt,'knitro.opt'); %    
    
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


% save tmp1
% return

Results.LL = -LL;
Results.b0_old = b0;

LLfun2 =  @(B) LL_mxl2(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,err_mtx,EstimOpt,B);

if EstimOpt.HessEstFix == 4
    [Results.LLdetailed,~,Results.hess] = LLfun2(Results.bhat);
else
    Results.LLdetailed = LLfun2(Results.bhat);
end
Results.LLdetailed = Results.LLdetailed.*INPUT.W;
if any(INPUT.MissingInd == 1) % In case of some missing data
   idx = sum(reshape(INPUT.MissingInd,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP)) == EstimOpt.NAlt; ...
   idx = sum(reshape(idx, EstimOpt.NCT, EstimOpt.NP),1)'; % no. of missing NCT for every respondent
   idx = EstimOpt.NCT - idx;
   R2 = mean(exp(-Results.LLdetailed./idx),1);
else
    R2 = mean(exp(-Results.LLdetailed/EstimOpt.NCT),1);
end

if EstimOpt.Scores ~= 0
   Results.Scores =  BayesScores(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,err_mtx,EstimOpt,Results.bhat);
end

if EstimOpt.HessEstFix == 1
% 	f = LL_mxl(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,err_mtx,EstimOpt,Results.bhat);
%     Results.jacobian = numdiff(@(B) LL_mxl(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,err_mtx,EstimOpt,B),f,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
    Results.jacobian = numdiff(@(B) INPUT.W.*LLfun2(B),Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) INPUT.W.*LLfun2(B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.hess = hessian(@(B) sum(INPUT.W.*LLfun2(B)), Results.bhat);
elseif EstimOpt.HessEstFix == 4

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
           [~, Results.jacobian] = LLfun2(Results.bhat);
           Results.jacobian = Results.jacobian.*INPUT.W(:, ones(1,size(Results.jacobian,2)));
    else
        Results.jacobian = numdiff(@(B) INPUT.W.*LLfun2(B),Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
    end
    RobustHess = Results.jacobian'*Results.jacobian;
    Results.ihess = Results.ihess*RobustHess*Results.ihess;
end

Results.std = sqrt(diag(Results.ihess));
Results.std(EstimOpt.BActive == 0) = NaN;
Results.std(EstimOpt.BLimit == 1) = 0;
Results.std(imag(Results.std) ~= 0) = NaN;

% save out_MXL

if EstimOpt.FullCov == 0
    Results.DetailsA = [Results.bhat(1:EstimOpt.NVarA),Results.std(1:EstimOpt.NVarA),pv(Results.bhat(1:EstimOpt.NVarA),Results.std(1:EstimOpt.NVarA))];
    std_out = Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*2);...
    std_out(imag(Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*2)) ~= 0) = NaN; ...
%     Results.DetailsV = [Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2).^2,2.*std_out.*abs(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2)),pv((Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2)).^2,2.*std_out.*abs(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2)))];
    Results.DetailsV = [abs(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2)),std_out,pv(abs(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2)),std_out)];
    if sum(EstimOpt.Dist(2:end) == 3) > 0
        Results.DetailsA(EstimOpt.Dist(2:end) == 3,:) = [exp(Results.bhat(EstimOpt.Dist(2:end) == 3)) + EstimOpt.Triang', exp(Results.bhat(EstimOpt.Dist(2:end) == 3)).*Results.std(EstimOpt.Dist(2:end) == 3), pv(exp(Results.bhat(EstimOpt.Dist(2:end) == 3)) + EstimOpt.Triang', exp(Results.bhat(EstimOpt.Dist(2:end) == 3)).*Results.std(EstimOpt.Dist(2:end) == 3))];
        btmp = Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2);
        stdx = zeros(sum(EstimOpt.Dist(2:end) == 3),1);
        g = [exp(Results.bhat(EstimOpt.Dist(2:end) == 3)), exp(btmp(EstimOpt.Dist(2:end) == 3))];
        indx = find(EstimOpt.Dist(2:end) == 3);
        for i = 1:sum(EstimOpt.Dist(2:end) == 3)   
           stdx(i) = sqrt(g(i,:)*Results.ihess([indx(i), indx(i)+EstimOpt.NVarA], [indx(i), indx(i)+EstimOpt.NVarA])*g(i,:)');
        end
        Results.DetailsV(EstimOpt.Dist(2:end) == 3,:) = [exp(btmp(EstimOpt.Dist(2:end) == 3))+ exp(Results.bhat(EstimOpt.Dist(2:end) == 3)) + EstimOpt.Triang', stdx, pv(exp(btmp(EstimOpt.Dist(2:end) == 3))+ exp(Results.bhat(EstimOpt.Dist(2:end) == 3)) + EstimOpt.Triang', stdx) ];
    end
    if sum(EstimOpt.Dist(2:end) == 4) > 0
        Results.DetailsA(EstimOpt.Dist(2:end) == 4,:) = [exp(Results.bhat(EstimOpt.Dist(2:end) == 4)), exp(Results.bhat(EstimOpt.Dist(2:end) == 4)).*Results.std(EstimOpt.Dist(2:end) == 4), pv(exp(Results.bhat(EstimOpt.Dist(2:end) == 4)), exp(Results.bhat(EstimOpt.Dist(2:end) == 4)).*Results.std(EstimOpt.Dist(2:end) == 4))];
        btmp = Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2);
        stdx = exp(btmp).*Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*2);
        Results.DetailsV(EstimOpt.Dist(2:end) == 4,:) = [exp(btmp(EstimOpt.Dist(2:end) == 4)), stdx(EstimOpt.Dist(2:end) == 4), pv(exp(btmp(EstimOpt.Dist(2:end) == 4)), stdx(EstimOpt.Dist(2:end) == 4)) ];
    end
    Results.R = [Results.DetailsA, Results.DetailsV];   
    if EstimOpt.NVarM > 0
        Results.DetailsM = [];
        for i=1:EstimOpt.NVarM; Results.DetailsM = [Results.DetailsM, [Results.bhat(EstimOpt.NVarA*(2+i-1)+1:EstimOpt.NVarA*(2+i)),Results.std(EstimOpt.NVarA*(2+i-1)+1:EstimOpt.NVarA*(2+i)),pv(Results.bhat(EstimOpt.NVarA*(2+i-1)+1:EstimOpt.NVarA*(2+i)),Results.std(EstimOpt.NVarA*(2+i-1)+1:EstimOpt.NVarA*(2+i)))]];end
        Results.R = [Results.R, Results.DetailsM];
    end
    if EstimOpt.Johnson > 0
        ResultsJ = NaN(EstimOpt.NVarA, 6);
        % Location parameters
        Results.DetailsJL = [Results.bhat((end - 2*EstimOpt.Johnson+1):(end - EstimOpt.Johnson)),Results.std((end - 2*EstimOpt.Johnson+1):(end - EstimOpt.Johnson)),pv(Results.bhat((end - 2*EstimOpt.Johnson+1):(end - EstimOpt.Johnson)),Results.std((end - 2*EstimOpt.Johnson+1):(end - EstimOpt.Johnson)))];
        % Scale parameters
        Results.DetailsJS = [exp(Results.bhat((end - EstimOpt.Johnson+1):end)),exp(Results.bhat((end - EstimOpt.Johnson+1):end)).*Results.std((end - EstimOpt.Johnson+1):end),pv(exp(Results.bhat((end - EstimOpt.Johnson+1):end)),exp(Results.bhat((end - EstimOpt.Johnson+1):end)).*Results.std((end - EstimOpt.Johnson+1):end))];
        ResultsJ(EstimOpt.Dist(2:end) > 4 & EstimOpt.Dist(2:end) <= 7,1:3) =   Results.DetailsJL;
        ResultsJ(EstimOpt.Dist(2:end) > 4 & EstimOpt.Dist(2:end) <= 7,4:6) =   Results.DetailsJS;
        Results.R = [Results.R, ResultsJ];
    end
    if EstimOpt.NVarS > 0
        Results.DetailsS = [];
        for i=1:EstimOpt.NVarS; Results.DetailsS = [Results.DetailsS; [Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NVarM)+i),Results.std(EstimOpt.NVarA*(2+EstimOpt.NVarM)+i),pv(Results.bhat(EstimOpt.NVarA*(2+EstimOpt.NVarM)+i),Results.std(EstimOpt.NVarA*(2+EstimOpt.NVarM)+i))]];end
        DetailsS0 = NaN(EstimOpt.NVarA,3);
        DetailsS0(1:EstimOpt.NVarS,1:3) = Results.DetailsS;
        Results.R = [Results.R, DetailsS0]; % might not work if NVarS > NVarA
    end
    
elseif EstimOpt.FullCov == 1
    Results.DetailsA = [Results.bhat(1:EstimOpt.NVarA),Results.std(1:EstimOpt.NVarA),pv(Results.bhat(1:EstimOpt.NVarA),Results.std(1:EstimOpt.NVarA))];
    Results.DetailsV = sdtri(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2), Results.ihess(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2,EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA+3)/2),EstimOpt);
    if sum(EstimOpt.Dist(2:end) == 3) > 0
        Results.DetailsA(EstimOpt.Dist(2:end) == 3,:) = [exp(Results.bhat(EstimOpt.Dist(2:end) == 3)) + EstimOpt.Triang', exp(Results.bhat(EstimOpt.Dist(2:end) == 3)).*Results.std(EstimOpt.Dist(2:end) == 3), pv(exp(Results.bhat(EstimOpt.Dist(2:end) == 3)) + EstimOpt.Triang', exp(Results.bhat(EstimOpt.Dist(2:end) == 3)).*Results.std(EstimOpt.Dist(2:end) == 3))];
        btmp = Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA-1)/2+2*EstimOpt.NVarA);
        btmp = btmp(EstimOpt.DiagIndex);
        stdx = zeros(sum(EstimOpt.Dist(2:end) == 3),1);
        g = [exp(Results.bhat(EstimOpt.Dist(2:end) == 3)), exp(btmp(EstimOpt.Dist(2:end) == 3))];
        indx = find(EstimOpt.Dist(2:end) == 3);
        DiagIndex = EstimOpt.DiagIndex(EstimOpt.Dist(2:end) == 3);
        for i = 1:sum(EstimOpt.Dist(2:end) == 3)   
           stdx(i) = sqrt(g(i,:)*Results.ihess([indx(i), DiagIndex(i)+EstimOpt.NVarA], [indx(i), DiagIndex(i)+EstimOpt.NVarA])*g(i,:)');
        end
        Results.DetailsV(EstimOpt.Dist(2:end) == 3,:) = [exp(btmp(EstimOpt.Dist(2:end) == 3))+ exp(Results.bhat(EstimOpt.Dist(2:end) == 3)) + EstimOpt.Triang', stdx, pv(exp(btmp(EstimOpt.Dist(2:end) == 3))+ exp(Results.bhat(EstimOpt.Dist(2:end) == 3)) + EstimOpt.Triang', stdx) ];
    end
    if sum(EstimOpt.Dist(2:end) == 4) > 0
        Results.DetailsA(EstimOpt.Dist(2:end) == 4,:) = [exp(Results.bhat(EstimOpt.Dist(2:end) == 4)), exp(Results.bhat(EstimOpt.Dist(2:end) == 4)).*Results.std(EstimOpt.Dist(2:end) == 4), pv(exp(Results.bhat(EstimOpt.Dist(2:end) == 4)), exp(Results.bhat(EstimOpt.Dist(2:end) == 4)).*Results.std(EstimOpt.Dist(2:end) == 4))];
        btmp = Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA-1)/2+2*EstimOpt.NVarA);
        btmp = btmp(EstimOpt.DiagIndex);
        stdx = Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA-1)/2+2*EstimOpt.NVarA);
        stdx = stdx(EstimOpt.DiagIndex);
        stdx = exp(btmp).*stdx;
        Results.DetailsV(EstimOpt.Dist(2:end) == 4,:) = [exp(btmp(EstimOpt.Dist(2:end) == 4)), stdx(EstimOpt.Dist(2:end) == 4), pv(exp(btmp(EstimOpt.Dist(2:end) == 4)), stdx(EstimOpt.Dist(2:end) == 4)) ];
    end
    if sum(EstimOpt.Dist(2:end) == 5) > 0
        btmp = Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA-1)/2+2*EstimOpt.NVarA);
        btmp = btmp(EstimOpt.DiagIndex);
        stdtmp = Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA-1)/2+2*EstimOpt.NVarA);
        stdtmp = stdtmp(EstimOpt.DiagIndex);
        Results.DetailsV(EstimOpt.Dist(2:end) ==5,:) = [btmp(EstimOpt.Dist(2:end) == 5).^2, 2*btmp(EstimOpt.Dist(2:end) == 5).*stdtmp(EstimOpt.Dist(2:end) == 5), pv(btmp(EstimOpt.Dist(2:end) == 5).^2, 2*btmp(EstimOpt.Dist(2:end) == 5).*stdtmp(EstimOpt.Dist(2:end) == 5))];
    end
    Results.R = [Results.DetailsA, Results.DetailsV];        
    if EstimOpt.NVarM > 0
        Results.DetailsM = [];
        for i=1:EstimOpt.NVarM; Results.DetailsM = [Results.DetailsM, [Results.bhat(EstimOpt.NVarA*(EstimOpt.NVarA/2+0.5+i)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+i)),Results.std(EstimOpt.NVarA*(EstimOpt.NVarA/2+0.5+i)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+i)),pv(Results.bhat(EstimOpt.NVarA*(EstimOpt.NVarA/2+0.5+i)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+i)),Results.std(EstimOpt.NVarA*(EstimOpt.NVarA/2+0.5+i)+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+i)))]]; end
        Results.R = [Results.R, Results.DetailsM];
    end


    if EstimOpt.Johnson > 0        
        ResultsJ = NaN(EstimOpt.NVARA, 6);
        % Location parameters
        Results.DetailsJL = [Results.bhat((end - 2*EstimOpt.Johnson+1):(end - EstimOpt.Johnson)),Results.std((end - 2*EstimOpt.Johnson+1):(end - EstimOpt.Johnson)),pv(Results.bhat((end - 2*EstimOpt.Johnson+1):(end - EstimOpt.Johnson)),Results.std((end - 2*EstimOpt.Johnson+1):(end - EstimOpt.Johnson)))];
        % Scale parameters
        Results.DetailsJS = [exp(Results.bhat((end - EstimOpt.Johnson+1):end)),exp(Results.bhat((end - EstimOpt.Johnson+1):end)).*Results.std((end - EstimOpt.Johnson+1):end),pv(exp(Results.bhat((end - EstimOpt.Johnson+1):end)),exp(Results.bhat((end - EstimOpt.Johnson+1):end)).*Results.std((end - EstimOpt.Johnson+1):end))];
        ResultsJ(EstimOpt.Dist(2:end) > 4 & EstimOpt.Dist(2:end) <= 7,1:3) =   Results.DetailsJL;
        ResultsJ(EstimOpt.Dist(2:end) > 4 & EstimOpt.Dist(2:end) <= 7,4:6) =   Results.DetailsJS;
        Results.R = [Results.R, ResultsJ];
    end
	if EstimOpt.NVarS > 0
        Results.DetailsS = [];
        for i=1:EstimOpt.NVarS; Results.DetailsS = [Results.DetailsS; [Results.bhat(EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+i),Results.std(EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+i),pv(Results.bhat(EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+i),Results.std(EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5+EstimOpt.NVarM)+i))]];end
        DetailsS0 = NaN(EstimOpt.NVarA,3);
        DetailsS0(1:EstimOpt.NVarS,1:3) = Results.DetailsS;
        Results.R = [Results.R, DetailsS0]; % might not work of NVarS > NVarA
	end
        
    Results.chol = [Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5)),Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5)),pv(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5)),Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*(EstimOpt.NVarA/2+1.5)))]; %Results.R = (1:EstimOpt.NVarA*(EstimOpt.NVarA/2+0.5),size(Results.R,2)+1:size(Results.R,2)+3) = 
    Results.DetailsVcov = tril(ones(EstimOpt.NVarA)); ...
        choltmp = Results.chol(:,1);
    if sum(EstimOpt.Dist(2:end) >= 3 & EstimOpt.Dist(2:end) <= 5) > 0
        choltmp(EstimOpt.DiagIndex(EstimOpt.Dist(2:end) >= 3 & EstimOpt.Dist(2:end) <= 5)) = 1;
    end
    Results.DetailsVcov(Results.DetailsVcov == 1) = choltmp; ...
	if sum(EstimOpt.Dist(2:end) >= 3 & EstimOpt.Dist(2:end) <= 5) > 0
        choltmp = sqrt(sum(Results.DetailsVcov(EstimOpt.Dist(2:end) >= 3 & EstimOpt.Dist(2:end) <= 5,:).^2,2));
        Results.DetailsVcov(EstimOpt.Dist(2:end) >= 3 & EstimOpt.Dist(2:end) <= 5,:) = Results.DetailsVcov(EstimOpt.Dist(2:end) >= 3 & EstimOpt.Dist(2:end) <= 5,:)./choltmp(:, ones(1,EstimOpt.NVarA));
	end
    Results.DetailsVcov = Results.DetailsVcov*Results.DetailsVcov';
    Results.DetailsVcor = corrcov(Results.DetailsVcov);
end



EstimOpt.params = length(b0);
% if isfield(EstimOpt,'BActive')
% 	EstimOpt.params = EstimOpt.params - sum(EstimOpt.BActive == 0);
% end
EstimOpt.params = EstimOpt.params - sum(EstimOpt.BActive == 0) + sum(EstimOpt.BLimit == 1);

Results.stats = [Results_old.MNL0.LL; Results.LL; 1-Results.LL/Results_old.MNL0.LL;R2; ((2*EstimOpt.params-2*Results.LL) + 2*EstimOpt.params*(EstimOpt.params+1)/(EstimOpt.NObs-EstimOpt.params-1))/EstimOpt.NObs; EstimOpt.NObs; EstimOpt.params];

Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;
Results.INPUT = INPUT;

disp(' ');
disp('Means');
disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
disp([char(EstimOpt.NamesA) ,blanks(EstimOpt.NVarA)', num2str(Results.DetailsA(:,1),'%8.4f') star_sig(Results.DetailsA(:,3)) num2str(Results.DetailsA(:,2:3),'%8.4f %8.4f')])

disp(' ');
disp('Standard deviations');
disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
disp([char(EstimOpt.NamesA) ,blanks(EstimOpt.NVarA)', num2str(Results.DetailsV(:,1),'%8.4f') star_sig(Results.DetailsV(:,3)) num2str(Results.DetailsV(:,2:3),'%8.4f %8.4f')])

if EstimOpt.NVarM > 0 
    for i = 1:EstimOpt.NVarM
        disp(' ');
        disp(['Explanatory variable of random parameters'' means - ', char(EstimOpt.NamesM(i))]);
        disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
        disp([char(EstimOpt.NamesA),blanks(EstimOpt.NVarA)', num2str(Results.DetailsM(:,i*3-2),'%8.4f'), star_sig(Results.DetailsM(:,i*3)), num2str(Results.DetailsM(:,i*3-1:i*3),'%8.4f %8.4f')])
    end
%     headM = 'var.  coef.     st.err.  p-value';
%     formM = '%1.0f %8.4f %8.4f %8.4f';    
%     for i = 2:EstimOpt.NVarM
%         formM = [formM, ' %8.4f %8.4f %8.4f'];
%         headM = [headM, '    coef.  sd.err.  p-value'];        
%     end
%     disp(headM);
%     disp(num2str([(1:EstimOpt.NVarA)', Results.DetailsM],formM))

end


if EstimOpt.Johnson > 0
    disp(' ');
    disp('Johnson location parameters');
    disp(['var.', blanks(size(char(EstimOpt.NamesA(EstimOpt.Dist(2:end)>= 5)),2)-2) ,'coef.      st.err.  p-value'])
    disp([char(EstimOpt.NamesA(EstimOpt.Dist(2:end)>= 5)), blanks(EstimOpt.Johnson)', num2str(Results.DetailsJL(:,1),'%8.4f') star_sig(Results.DetailsJL(:,3)) num2str(Results.DetailsJL(:,2:3),'%8.4f %8.4f')])
    disp(' ');
    disp('Johnson scale parameters');
    disp(['var.', blanks(size(char(EstimOpt.NamesA(EstimOpt.Dist(2:end)>= 5)),2)-2) ,'coef.      st.err.  p-value'])
    disp([char(EstimOpt.NamesA(EstimOpt.Dist(2:end)>= 5)), blanks(EstimOpt.Johnson)', num2str(Results.DetailsJS(:,1),'%8.4f') star_sig(Results.DetailsJS(:,3)) num2str(Results.DetailsJS(:,2:3),'%8.4f %8.4f')])
end

if EstimOpt.NVarS > 0
    disp(' ');
    disp('Covariates of scale');
    disp(['var.',blanks(size(char(EstimOpt.NamesS),2)-2) ,'coef.      st.err.  p-value'])
    disp([char(EstimOpt.NamesS) ,blanks(EstimOpt.NVarS)', num2str(Results.DetailsS(:,1),'%8.4f') star_sig(Results.DetailsS(:,3)) num2str(Results.DetailsS(:,2:3),'%8.4f %8.4f')])
end

head = {'var.' , 'coef.', 'st.err.' , 'p-value'};
headx = [head, repmat(head(1,2:4),1,1+EstimOpt.NVarM+(EstimOpt.Johnson>0)*2)];

Results.R_out = cell(3+EstimOpt.NVarA + EstimOpt.NVarS + (EstimOpt.NVarS>0)*2 + 2 + 7, 1 + 3*(2+EstimOpt.NVarM+(EstimOpt.Johnson>0)*2));
   
if EstimOpt.FullCov == 0
    Results.R_out(1,1) = {'MXL_d'};
else
    Results.R_out(1,1) = {'MXL'};
end
Results.R_out(2,[2,5]) = {'Means', 'Standard Deviations'  };
Results.R_out(3,:) = headx;
Results.R_out(4:(EstimOpt.NVarA+3),1:7) = [EstimOpt.NamesA, num2cell(Results.DetailsA), num2cell(Results.DetailsV)];
if EstimOpt.NVarM > 0
    Results.R_out(4:(EstimOpt.NVarA+3),8:7+EstimOpt.NVarM*3) = num2cell(Results.DetailsM);
    Results.R_out(2,8:3:(5+EstimOpt.NVarM*3)) = EstimOpt.NamesM;
end

if EstimOpt.Johnson > 0
    Results.R_out(2,[8+EstimOpt.NVarM*3,8+EstimOpt.NVarM*3+3]) = {'Johnson location parameters','Johnson scale parameters'};
	Results.R_out(4:(EstimOpt.NVarA+3),8+EstimOpt.NVarM*3:13+EstimOpt.NVarM*3) = num2cell(ResultsJ);    
end

if EstimOpt.NVarS > 0
    Results.R_out(EstimOpt.NVarA+4,1) = {'Covariates of Scale'};
    Results.R_out(EstimOpt.NVarA+5,1:4) = head;
    Results.R_out(EstimOpt.NVarA+6:EstimOpt.NVarA + EstimOpt.NVarS + 5,1:4) = [EstimOpt.NamesS,num2cell(Results.DetailsS)];
end

Results.R_out(EstimOpt.NVarA + EstimOpt.NVarS + (EstimOpt.NVarS>0)*2 + 5,1) = {'Model characteristics'};
Results.R_out(EstimOpt.NVarA + EstimOpt.NVarS + (EstimOpt.NVarS>0)*2 + 6:end,1) = {'LL0'; 'LL' ; 'McFadden R2';'Ben-Akiva R2' ;'AIC/n' ; 'n'; 'k'};
Results.R_out(EstimOpt.NVarA + EstimOpt.NVarS + (EstimOpt.NVarS>0)*2 + 6:end,2) = num2cell(Results.stats);

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

% save(EstimOpt.fnameout, 'Results')

end

