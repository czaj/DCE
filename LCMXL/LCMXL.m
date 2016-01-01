function Results = LCMXL(INPUT,Results_old,EstimOpt,OptimOpt)


% save tmp_LCMXL
% return

global B_backup

tic

Results.bhat = [];
Results.R = [];
Results.R_out = {};
Results.stats = [];


%% Check data and inputs


if nargin < 3 % check no. of inputs
    error('Too few input arguments for LCMXL(INPUT,EstimOpt,OptimOpt)')
end

if isfield(EstimOpt,'NClass') == 0; 
    EstimOpt.NClass = 2; 
end

disp(' ');
disp('__________________________________________________________________________________________________________________');
disp(' ');

warning off MATLAB:mir_warning_maybe_uninitialized_temporary

format shortG;
format compact;

disp(num2str(EstimOpt.NClass,'Estimating LCMXL model with %1.0f classes ...'))

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
else
    disp('with correlated random parameters ...')            
    if EstimOpt.WTP_space > 0
        disp('in WTP-space.')
	else
        disp('in preference-space.')
    end
end

if isfield(EstimOpt,'Dist') == 0 || isempty(EstimOpt.Dist)
	EstimOpt.Dist = zeros(1,EstimOpt.NClass*EstimOpt.NVarA);
    if EstimOpt.WTP_space == 0
        cprintf(rgb('DarkOrange'), 'WARNING: distributions for random parameters not specified - assuming normality \n')
    else
        cprintf(rgb('DarkOrange'), 'WARNING: distributions for random parameters not specified - assuming normality (monetary parameters assumed log-normal) \n')
        EstimOpt.Dist = reshape(EstimOpt.Dist, EstimOpt.NVarA, EstimOpt.NClass);
        EstimOpt.Dist(end-EstimOpt.WTP_space+1:end,:) = 1; % cost in WTP-space models log-normally distributed
        EstimOpt.Dist = reshape(EstimOpt.Dist, 1, EstimOpt.NVarA*EstimOpt.NClass);
    end
elseif numel(EstimOpt.Dist) == 1
    EstimOpt.Dist = EstimOpt.Dist*ones(1,EstimOpt.NClass*EstimOpt.NVarA);
elseif length(EstimOpt.Dist) == EstimOpt.NVarA
    EstimOpt.Dist = repmat(EstimOpt.Dist(:)', 1,EstimOpt.NClass);
elseif length(EstimOpt.Dist) ~= EstimOpt.NClass*EstimOpt.NVarA
	EstimOpt.Dist = zeros(1,EstimOpt.NClass*EstimOpt.NVarA);
    if EstimOpt.WTP_space == 0
        cprintf(rgb('DarkOrange'), 'WARNING: incorrect number of distributions for random pararameters specified - overriding user setting and assuming normality \n')
    else
        cprintf(rgb('DarkOrange'), 'WARNING: incorrect number of distributions for random pararameters specified - overriding user setting and assuming normality (monetary parameters assumed log-normal) \n')
        EstimOpt.Dist = reshape(EstimOpt.Dist, EstimOpt.NVarA, EstimOpt.NClass);
        EstimOpt.Dist(end-EstimOpt.WTP_space+1:end,:) = 1; % cost in WTP-space models log-normally distributed
        EstimOpt.Dist = reshape(EstimOpt.Dist, 1, EstimOpt.NVarA*EstimOpt.NClass);
    end
else
    EstimOpt.Dist = EstimOpt.Dist(:)';
end
disp(['Random parameters distributions: ', num2str(EstimOpt.Dist),' (-1 - constant, 0 - normal, 1 - lognormal, 2 - spike, 5 - Weibull)'])

if EstimOpt.WTP_space > 0 && sum(EstimOpt.Dist(end-EstimOpt.WTP_space+1:end)==1) > 0 && any(mean(INPUT.Xa(:,end-EstimOpt.WTP_space+1:end)) >= 0)
    cprintf(rgb('DarkOrange'), 'WARNING: Cost attributes with log-normally distributed parameters should enter utility function with a ''-'' sign \n')
end

if isfield(INPUT, 'Xc') == 0 || numel(INPUT.Xc) == 0
    INPUT.Xc = ones(size(INPUT.Y,1),1);
else
    INPUT.Xc = [ones(size(INPUT.Y,1),1), INPUT.Xc];
end

if det(INPUT.Xc'*INPUT.Xc) == 0
    error('Xc matrix cointains collinear variables')
end

EstimOpt.NVarC = size(INPUT.Xc,2); % no. of variables explaining class probabilities

if isfield(EstimOpt, 'Scores') == 0
   EstimOpt.Scores = 0; 
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
%         tic; EstimOpt.WTP_matrix = 1:EstimOpt.WTP_space;...
%         EstimOpt.WTP_matrix = EstimOpt.WTP_matrix(floor((0:size(EstimOpt.WTP_matrix,2)*WTP_att-1)/WTP_att)+1); toc
        end
	elseif size(EstimOpt.WTP_matrix,2) ~= EstimOpt.NVarA - EstimOpt.WTP_space
        error('Dimensions of EstimOpt.WTP_matrix not correct - for each non-monetary attribute provide no. of attribute to multiply it with')
    else
        EstimOpt.WTP_matrix = EstimOpt.WTP_matrix(:)';
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
    elseif size(EstimOpt.NamesC,1) ~= EstimOpt.NVarC - 1
        EstimOpt.NamesC = EstimOpt.NamesC';
    end
    EstimOpt.NamesC = [{'Cons'};EstimOpt.NamesC];
else
    EstimOpt.NamesC = {'Cons'};
end


%% Starting values


EstimOpt.jitter1 = 0.8; % Jittering parameter (relative) for MNL or MXL starting values (attributes)
EstimOpt.jitter2 = 0.3; % Jittering parameter (absolute) for class probabilities starting values

if EstimOpt.FullCov == 0
	if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == 2*EstimOpt.NVarA*EstimOpt.NClass + EstimOpt.NVarC*(EstimOpt.NClass-1)
            b0 = B_backup(:);
            disp('Using the starting values from Backup')
    elseif isfield(Results_old,'LCMXL_d') && isfield(Results_old.LCMXL_d,'b0') % starting values provided
        Results_old.LCMXL_d.b0_old = Results_old.LCMXL_d.b0(:);
        Results_old.LCMXL_d = rmfield(Results_old.LCMXL_d,'b0');
        if length(Results_old.LCMXL_d.b0_old) ~= 2*EstimOpt.NVarA*EstimOpt.NClass + EstimOpt.NVarC*(EstimOpt.NClass-1)
            cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.LCMXL_d = rmfield(Results_old.LCMXL_d,'b0_old');
        else
            b0 = Results_old.LCMXL_d.b0_old(:);
        end
	end
    if  ~exist('b0','var')
        if isfield(Results_old,'LC') && isfield(Results_old.LC,'bhat')
            disp('Using LC coefficients for starting values')
            Results_old.LC.bhat = Results_old.LC.bhat(:);
            tmp = Results_old.LC.bhat(1:EstimOpt.NVarA*EstimOpt.NClass);
            if EstimOpt.WTP_space > 0
                tmp = reshape(tmp, EstimOpt.NVarA, EstimOpt.NClass);
                tmp(end - EstimOpt.WTP_space+1:end,:) = log(abs(tmp(end - EstimOpt.WTP_space+1:end,:)));
                tmp = reshape(tmp,EstimOpt.NVarA*EstimOpt.NClass,1); 
            end
            b0 = [tmp;0.3*ones(EstimOpt.NVarA*EstimOpt.NClass,1);Results_old.LC.bhat(EstimOpt.NClass*EstimOpt.NVarA+1:end) ];                        
        elseif isfield(Results_old,'MXL_d') && isfield(Results_old.MXL_d,'bhat')
            disp('Using MXL_d coefficients for starting values')
            Results_old.MXL_d.bhat = Results_old.MXL_d.bhat(:);
            b0 = [Results_old.MXL_d.bhat(repmat(reshape(1:2*EstimOpt.NVarA,EstimOpt.NVarA,2),EstimOpt.NClass,1),:) .* ...
                (EstimOpt.jitter1 .* unifrnd(0,ones(2*EstimOpt.NVarA.*EstimOpt.NClass,1))) ; ...
                EstimOpt.jitter2 + unifrnd(-1,ones(EstimOpt.NVarC.*(EstimOpt.NClass-1),1))];
        elseif isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat')
            disp('Using MNL results as starting values')
            Results_old.MNL.bhat = Results_old.MNL.bhat(:);
            if EstimOpt.WTP_space > 0
                Results_old.MNL.bhat(end - EstimOpt.WTP_space+1:end) = log(Results_old.MNL.bhat((end - EstimOpt.WTP_space+1:end)));
            end
            b0 = [Results_old.MNL.bhat((1:size(Results_old.MNL.bhat,1))'*ones(1,EstimOpt.NClass),:) .* ...
                (EstimOpt.jitter1 .* unifrnd(0,ones(EstimOpt.NVarA.*EstimOpt.NClass,1))) ; ...
                0.3*ones(EstimOpt.NVarA*EstimOpt.NClass,1); EstimOpt.jitter2 + unifrnd(-1,ones(EstimOpt.NVarC.*(EstimOpt.NClass-1),1))];
        else
            error('No starting values available - run MNL, LC or MXL_d first')
        end
    end    
    
else % EstimOpt.FullCov == 1
	if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NClass*(EstimOpt.NVarA + sum(1:EstimOpt.NVarA)) +EstimOpt.NVarC*(EstimOpt.NClass-1)
            b0 = B_backup(:);
            disp('Using the starting values from Backup')
    elseif isfield(Results_old,'LCMXL') && isfield(Results_old.LCMXL,'b0') % starting values provided
        Results_old.LCMXL.b0_old = Results_old.LCMXL.b0(:);
        Results_old.LCMXL = rmfield(Results_old.LCMXL,'b0');
        if length(Results_old.LCMXL.b0_old) ~= EstimOpt.NClass*(EstimOpt.NVarA + sum(1:EstimOpt.NVarA)) +EstimOpt.NVarC*(EstimOpt.NClass-1)
            cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.LCMXL = rmfield(Results_old.LCMXL,'b0_old');
        else
            b0 = Results_old.LCMXL.b0_old(:);
        end
	end
    if  ~exist('b0','var')
        if isfield(Results_old,'LCMXL_d') && isfield(Results_old.LCMXL_d,'bhat')
            disp('Using LCMXL_d coefficients for starting values')
            Results_old.LCMXL_d.bhat = Results_old.LCMXL_d.bhat(:);
            b0 = zeros(2*(EstimOpt.NVarA + sum(1:EstimOpt.NVarA))+EstimOpt.NVarC*(EstimOpt.NClass-1),1);
            b0(1:EstimOpt.NVarA*EstimOpt.NClass) = Results_old.LCMXL_d.bhat(1:EstimOpt.NVarA*EstimOpt.NClass);
            b0(EstimOpt.NClass*(EstimOpt.NVarA + sum(1:EstimOpt.NVarA))+1:end) = Results_old.LCMXL_d.bhat(2*EstimOpt.NVarA*EstimOpt.NClass+1:end);
            for i = 1: EstimOpt.NClass
                vc_tmp = diag(Results_old.LCMXL_d.bhat(EstimOpt.NVarA*EstimOpt.NClass +1 + (i-1)*EstimOpt.NVarA:EstimOpt.NVarA*EstimOpt.NClass +i*EstimOpt.NVarA ));
                b0(EstimOpt.NVarA*EstimOpt.NClass+1 + (i-1)*sum(1+EstimOpt.NVarA):EstimOpt.NVarA*EstimOpt.NClass+ i*sum(1:EstimOpt.NVarA)) = vc_tmp(tril(ones(size(vc_tmp)))==1);
            end
        elseif isfield(Results_old,'LC') && isfield(Results_old.LC,'bhat')
            disp('Using LC coefficients for starting values')
            Results_old.LC.bhat = Results_old.LC.bhat(:);
            b0 = [Results_old.LC.bhat(1:EstimOpt.NVarA*EstimOpt.NClass);zeros(sum(1:EstimOpt.NVarA)*EstimOpt.NClass,1);Results_old.LC.bhat(EstimOpt.NClass*EstimOpt.NVarA+1:end)];            
        elseif isfield(Results_old,'MXL') && isfield(Results_old.MXL,'bhat')
            disp('Using MXL coefficients for starting values')
            Results_old.MXL.bhat = Results_old.MXL.bhat(:);
            b0_tmp = [Results_old.MXL.bhat(repmat((1:EstimOpt.NVarA)',EstimOpt.NClass,1)) ;Results_old.MXL.bhat(repmat(((EstimOpt.NVarA+1):(EstimOpt.NVarA+sum(1:EstimOpt.NVarA)))',EstimOpt.NClass,1))];
            b0 = [b0_tmp.* ...
                (EstimOpt.jitter1 .* unifrnd(0,ones((EstimOpt.NVarA+sum(1:EstimOpt.NVarA))*EstimOpt.NClass,1))) ; ...
                EstimOpt.jitter2 + unifrnd(-1,ones(EstimOpt.NVarC.*(EstimOpt.NClass-1),1))];
        elseif isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat')
            disp('Using MNL results as starting values')
            Results_old.MNL.bhat = Results_old.MNL.bhat(:);
            b0 = [Results_old.MNL.bhat((1:size(Results_old.MNL.bhat,1))'*ones(1,EstimOpt.NClass),:) .* ...
                 (EstimOpt.jitter1 .* unifrnd(0,ones(EstimOpt.NVarA.*EstimOpt.NClass,1))) ; ...
                 zeros(sum(1:EstimOpt.NVarA)*EstimOpt.NClass,1); EstimOpt.jitter2 + unifrnd(-1,ones(EstimOpt.NVarC.*(EstimOpt.NClass-1),1))];
        else
            error('No starting values available - run MNL, LC or MXL first')
        end
    end        
end


%% Optimization Options

if  isfield(EstimOpt,'BActive')
	EstimOpt.BActive = EstimOpt.BActive(:)';
end

if sum(EstimOpt.Dist == -1) > 0
    if isfield(EstimOpt,'BActive') == 0 || isempty(EstimOpt.BActive)
        EstimOpt.BActive = ones(1,length(b0));
    end
    if EstimOpt.FullCov == 0
        EstimOpt.BActive(EstimOpt.NVarA*EstimOpt.NClass+find(EstimOpt.Dist == -1)) = 0;
    elseif EstimOpt.FullCov == 1
        for i = 1:EstimOpt.NClass
            Vt = tril(ones(EstimOpt.NVarA));
            Vt(EstimOpt.Dist((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA)==-1,:) = 0;
%             EstimOpt.BActive(EstimOpt.NVarA*EstimOpt.NClass +(i-1)*sum(1:EstimOpt.NVarA)+ 1:EstimOpt.NVarA*EstimOpt.NClass +i*sum(1:EstimOpt.NVarA)) = EstimOpt.BActive(EstimOpt.NVarA*EstimOpt.NClass +(i-1)*sum(1:EstimOpt.NVarA)+ 1:EstimOpt.NVarA*EstimOpt.NClass +i*sum(1:EstimOpt.NVarA)) .* Vt(find(tril(ones(size(Vt)))))';
            EstimOpt.BActive(EstimOpt.NVarA*EstimOpt.NClass +(i-1)*sum(1:EstimOpt.NVarA)+ 1:EstimOpt.NVarA*EstimOpt.NClass +i*sum(1:EstimOpt.NVarA)) = EstimOpt.BActive(EstimOpt.NVarA*EstimOpt.NClass +(i-1)*sum(1:EstimOpt.NVarA)+ 1:EstimOpt.NVarA*EstimOpt.NClass +i*sum(1:EstimOpt.NVarA)) .* Vt(tril(ones(size(Vt)))~=0)';
        end
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

if ~isfield(EstimOpt,'BActiveClass') || isempty(EstimOpt.BActiveClass) || sum(EstimOpt.BActiveClass == 0) == 0
    EstimOpt.BActiveClass = ones(EstimOpt.NVarA,1);
end


%% Generate pseudo-random Draws


if isfield(EstimOpt,'Seed1') == 1
    rng(EstimOpt.Seed1);
end
cprintf('Simulation with ');
cprintf('*blue',[num2str(EstimOpt.NRep) ' ']); 

if EstimOpt.Draws == 1
    cprintf('*blue','Pseudo-random '); cprintf('Draws \n');
    err_mtx = randn(EstimOpt.NP*EstimOpt.NRep, EstimOpt.NClass*EstimOpt.NVarA); %to be cut down later   
elseif EstimOpt.Draws == 2 % LHS
    cprintf('*blue','Latin Hypercube Sampling '); cprintf('Draws \n');
    err_mtx=lhsnorm(zeros((EstimOpt.NClass*EstimOpt.NVarA)*EstimOpt.NP,1),diag(ones((EstimOpt.NClass*EstimOpt.NVarA)*EstimOpt.NP,1)),EstimOpt.NRep); 
    err_mtx = reshape(err_mtx, EstimOpt.NRep*EstimOpt.NP, EstimOpt.NVarA+1);
elseif EstimOpt.Draws >= 3 % Quasi random Draws 
    if EstimOpt.Draws == 3
        cprintf('*blue','Halton '); cprintf('Draws (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = haltonset(EstimOpt.NClass*EstimOpt.NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); % 
    elseif EstimOpt.Draws == 4 % apply reverse-radix scrambling
        cprintf('*blue','Halton '); cprintf('Draws with reverse radix scrambling (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = haltonset(EstimOpt.NClass*EstimOpt.NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); % 
        hm1 = scramble(hm1,'RR2');
    elseif EstimOpt.Draws == 5
        cprintf('*blue','Sobol '); cprintf('Draws (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = sobolset(EstimOpt.NClass*EstimOpt.NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); 
    elseif EstimOpt.Draws == 6
        cprintf('*blue','Sobol '); cprintf('Draws with random linear scramble and random digital shift (skip = '); cprintf(num2str(EstimOpt.HaltonSkip)); cprintf('; leap = '); cprintf(num2str(EstimOpt.HaltonLeap)); cprintf(') \n')
        hm1 = sobolset(EstimOpt.NClass*EstimOpt.NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap); 
        hm1 = scramble(hm1,'MatousekAffineOwen');
    end
    
    err_mtx = net(hm1,EstimOpt.NP*EstimOpt.NRep); % this takes every point:
    clear hm1;
    
    if EstimOpt.NP*EstimOpt.NRep < 3e+7   
        err_mtx(:,EstimOpt.Dist < 3) = icdf('Normal',err_mtx(:,EstimOpt.Dist < 2),0,1); %to be cut down later
 
    else % this is for very large number of Draws * variables
        for i=1:EstimOpt.NClass*EstimOpt.NVarA
            if EstimOpt.Dist(i) < 3
                err_mtx(:,i) = icdf('Normal',err_mtx(:,i),0,1); %to be cut down later
            else
               % err_mtx(:,i) = -log(1-err_mtx(:,i)); %to be cut down later
                error('weibull does not work in LCMXL')
            end
        end        
    end  
    err_mtx(:,EstimOpt.Dist == -1) = 0;
end


%% Display Options


if ((isfield(EstimOpt, 'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    OptimOpt.GradObj = 'on';
end

% if any(EstimOpt.MissingNCT(:) == 1) && EstimOpt.NumGrad == 0
%    EstimOpt.NumGrad = 1;
%    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - missing choice tasks not supported by analytical gradient \n')
% end

if any(EstimOpt.MissingAlt(:) == 1) && EstimOpt.NumGrad == 0
   EstimOpt.NumGrad = 1;
   cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - missing alternatives not supported by analytical gradient \n')
end

% if any(EstimOpt.BActiveClass == 1) && EstimOpt.NumGrad == 0
%    EstimOpt.NumGrad = 1;
%    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - parameters constrained between classes not supported by analytical gradient \n')
% end
if any(EstimOpt.Dist > 1) && EstimOpt.NumGrad == 0
   EstimOpt.NumGrad = 1;
   cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - analytical gradient available for normally or lognormally distributed parameters only \n')
end

if EstimOpt.WTP_space > 1 && EstimOpt.NumGrad == 0
   EstimOpt.NumGrad = 1;
   cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - WTP_space > 1 not supported by analytical gradient \n')
end

cond = sum(reshape(EstimOpt.Dist', EstimOpt.NVarA, EstimOpt.NClass),2);
cond = sum((cond >0).*(cond < EstimOpt.NClass));
if cond > 0 && EstimOpt.NumGrad == 0
   EstimOpt.NumGrad = 1;
   disp('Setting user-supplied gradient to numerical - particular random parameter combinations not supported by analytical gradient \n')
end

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


%% Rescructure data


INPUT.XXc = INPUT.Xc(1:EstimOpt.NCT*EstimOpt.NAlt:end,:); % NP x NVarC

%INPUT.YY = reshape(INPUT.Y,EstimOpt.NAlt,EstimOpt.NCT,EstimOpt.NP);...
%INPUT.YY = INPUT.YY(:,:,:,ones(EstimOpt.NRep,1));
%INPUT.YY = permute(INPUT.YY, [1 2 4 3]); %NAlt x NCT x NRep  x NP

INPUT.YYY = reshape(INPUT.Y,EstimOpt.NAlt,EstimOpt.NCT,EstimOpt.NP);
idx = sum(reshape(INPUT.MissingInd,EstimOpt.NAlt,EstimOpt.NCT,EstimOpt.NP)) == EstimOpt.NAlt; ...
INPUT.YYY(idx(ones(EstimOpt.NAlt,1),:,:) == 1) = NaN; % replace YYY in missing choice-tasks with NaN
INPUT.YY = reshape(INPUT.YYY,EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP);

INPUT.Xa(INPUT.MissingInd == 1,:) = NaN;
INPUT.XXa  = reshape(INPUT.Xa, EstimOpt.NAlt*EstimOpt.NCT, EstimOpt.NP, EstimOpt.NVarA);
INPUT.XXa = permute(INPUT.XXa, [1 3 2]); % NAlt*NCT x NVarA  x NP

EstimOpt.Dist = reshape(EstimOpt.Dist, EstimOpt.NVarA, EstimOpt.NClass);

err_sliced = err_mtx'; % NVarA*class x NRep * NP
if isfield(EstimOpt, 'Drawskeep') && ~isempty(EstimOpt.Drawskeep) && EstimOpt.Drawskeep == 1
    Results.err = err_sliced;
end

EstimOpt.indx1 = [];
EstimOpt.indx2 = [];
if EstimOpt.NumGrad == 0 && EstimOpt.FullCov == 1
   for i = 1:EstimOpt.NVarA
      EstimOpt.indx1 = [EstimOpt.indx1, i:EstimOpt.NVarA];
      EstimOpt.indx2 = [EstimOpt.indx2, i*ones(1,EstimOpt.NVarA+1-i)];
   end
   tmp1 = EstimOpt.indx1;
   tmp2 = EstimOpt.indx2;
   for j =1: (EstimOpt.NClass-1)
       tmp1 = tmp1+EstimOpt.NVarA;
       tmp2 = tmp2+EstimOpt.NVarA;
       EstimOpt.indx1 = [EstimOpt.indx1,tmp1];
       EstimOpt.indx2 = [EstimOpt.indx2,tmp2];
   end
end


%% Estimation


LLfun = @(B) LL_lcmxl_MATlike(INPUT.YY,INPUT.XXa,INPUT.XXc,err_sliced,EstimOpt,OptimOpt,B);
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
Results.LLdetailed = LL_lcmxl(INPUT.YY,INPUT.XXa,INPUT.XXc,err_sliced,EstimOpt,Results.bhat);
if any(INPUT.MissingInd == 1) % In case of some missing data
   idx = sum(reshape(INPUT.MissingInd,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP)) == EstimOpt.NAlt; ...
   idx = sum(reshape(idx, EstimOpt.NCT, EstimOpt.NP),1)'; % no. of missing NCT for every respondent
   idx = EstimOpt.NCT - idx;
   R2 = mean(exp(-Results.LLdetailed./idx),1);
else
    R2 = mean(exp(-Results.LLdetailed/EstimOpt.NCT),1);
end

% save out_lcmxl

Results.LL = -LL;
if EstimOpt.Scores ~= 0 
    [Results.PScores, Results.BScores] = BayesScoresLCMXL(INPUT.YY,INPUT.XXa,INPUT.XXc,err_sliced,EstimOpt,Results.bhat);
end

if EstimOpt.HessEstFix == 1
	f = LL_lcmxl(INPUT.YY,INPUT.XXa,INPUT.XXc,err_sliced,EstimOpt,Results.bhat);
    Results.jacobian = numdiff(@(B) LL_lcmxl(INPUT.YY,INPUT.XXa,INPUT.XXc,err_sliced,EstimOpt,B),f, Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B)  LL_lcmxl(INPUT.YY,INPUT.XXa,INPUT.XXc,err_sliced,EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.jacobian = hessian(@(B)  sum(LL_lcmxl(INPUT.YY,INPUT.XXa,INPUT.XXc,err_sliced,EstimOpt,B),1), Results.bhat);
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
Results.std(imag(Results.std) ~= 0) = NaN;

if EstimOpt.FullCov == 0
    Results.DetailsA = [Results.bhat(1:EstimOpt.NClass*EstimOpt.NVarA), Results.std(1:EstimOpt.NClass*EstimOpt.NVarA), pv(Results.bhat(1:EstimOpt.NClass*EstimOpt.NVarA), Results.std(1:EstimOpt.NClass*EstimOpt.NVarA))];
    l = EstimOpt.NClass*EstimOpt.NVarA;
    Results.DetailsV = [Results.bhat(l+1:l+EstimOpt.NClass*EstimOpt.NVarA ).^2, 2*Results.std(l+1:l+EstimOpt.NClass*EstimOpt.NVarA).*abs(Results.bhat(l+1:l+EstimOpt.NClass*EstimOpt.NVarA )), pv(Results.bhat(l+1:l+EstimOpt.NClass*EstimOpt.NVarA), Results.std(l+1:l+EstimOpt.NClass*EstimOpt.NVarA))];
    l = 2*EstimOpt.NClass*EstimOpt.NVarA;
    Results.DetailsC = [Results.bhat(l+1:l+EstimOpt.NVarC*(EstimOpt.NClass-1)), Results.std(l+1:l+EstimOpt.NVarC*(EstimOpt.NClass-1)), pv(Results.bhat(l+1:l+EstimOpt.NVarC*(EstimOpt.NClass-1)), Results.std(l+1:l+EstimOpt.NVarC*(EstimOpt.NClass-1)))];
else
    Results.DetailsA = [Results.bhat(1:EstimOpt.NClass*EstimOpt.NVarA), Results.std(1:EstimOpt.NClass*EstimOpt.NVarA), pv(Results.bhat(1:EstimOpt.NClass*EstimOpt.NVarA), Results.std(1:EstimOpt.NClass*EstimOpt.NVarA))];
    l = EstimOpt.NClass*EstimOpt.NVarA;
    Results.DetailsV = [];
    for i = 1:EstimOpt.NClass
        Results.DetailsV = [Results.DetailsV ; sdtri(Results.bhat(l+1:l+sum(1:EstimOpt.NVarA)), Results.ihess(l+1:l+sum(1:EstimOpt.NVarA),l+1:l+sum(1:EstimOpt.NVarA)),EstimOpt)];
        l = l + sum(1:EstimOpt.NVarA);
    end
    l = EstimOpt.NClass*(EstimOpt.NVarA+sum(1:EstimOpt.NVarA));
    Results.DetailsC = [Results.bhat(l+1:l+EstimOpt.NVarC*(EstimOpt.NClass-1)), Results.std(l+1:l+EstimOpt.NVarC*(EstimOpt.NClass-1)), pv(Results.bhat(l+1:l+EstimOpt.NVarC*(EstimOpt.NClass-1)), Results.std(l+1:l+EstimOpt.NVarC*(EstimOpt.NClass-1)))];
end

Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;
Results.INPUT = INPUT;

bclass = reshape([Results.bhat(end - EstimOpt.NVarC*(EstimOpt.NClass-1) +1 :end); zeros(EstimOpt.NVarC,1)], EstimOpt.NVarC, EstimOpt.NClass);
V = exp(INPUT.XXc*bclass);% NP x NClass
Vsum = sum(V,2);
Results.PNClass = mean(V./Vsum(:,ones(EstimOpt.NClass,1)),1); %1 x NClass
clear bclass V Vsum

EstimOpt.params = length(b0);
if isfield(EstimOpt,'BActive')
	EstimOpt.params = EstimOpt.params - sum(EstimOpt.BActive == 0);
end
Results.R = [Results.DetailsA; Results.DetailsV; Results.DetailsC; [Results.PNClass',zeros(size(Results.PNClass')),zeros(size(Results.PNClass'))]];

Results.stats = [Results_old.MNL0.LL; Results.LL; 1-Results.LL/Results_old.MNL0.LL; R2; ((2*EstimOpt.params-2*Results.LL) + 2*EstimOpt.params*(EstimOpt.params+1)/(EstimOpt.NObs-EstimOpt.params-1))/EstimOpt.NObs; EstimOpt.NObs; EstimOpt.params];

Results.R_out = cell(4+EstimOpt.NVarA + 2+EstimOpt.NVarC +3+ 2 + 7, 1+6*EstimOpt.NClass);  
if EstimOpt.FullCov == 0
    Results.R_out(1,1) = {'LCMXL_d'};
else
    Results.R_out(1,1) = {'LCMXL'};
end
head = {'var.' , 'coef.', 'st.err.' , 'p-value'};
NClasses = {'NClass 1','NClass 2', 'NClass 3', 'NClass 4', 'NClass 5', 'NClass 6', 'NClass 7', 'NClass 8', 'NClass 9','NClass 10'};
headx = [head, repmat(head(1,2:4),1,2*EstimOpt.NClass-1)];
if EstimOpt.NClass<=10
    Results.R_out(2,2:6:(2+(EstimOpt.NClass-1)*6)) = NClasses(1,1:EstimOpt.NClass);
end
Results.R_out(3,2:3:(2+(EstimOpt.NClass)*6-3)) = repmat([{'Means'},{'Std. Deviations'}], 1,EstimOpt.NClass);
Results.R_out(4,:) = headx;
Results.R_out(5:(EstimOpt.NVarA+4),1) = EstimOpt.NamesA;
for i = 1: EstimOpt.NClass
    Results.R_out(5:(EstimOpt.NVarA+4),2 + 6*(i-1):1+6*i) = [num2cell(Results.DetailsA((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA,:)),num2cell(Results.DetailsV((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA,:))];
end
Results.R_out(EstimOpt.NVarA+5,1) = {'Latent class probability model'};
Results.R_out(EstimOpt.NVarA+6,1:1+3*(EstimOpt.NClass-1)) = headx(1:1+3*(EstimOpt.NClass-1));
Results.R_out(EstimOpt.NVarA+7:EstimOpt.NVarA+6+EstimOpt.NVarC,1) = EstimOpt.NamesC;
for i = 1: EstimOpt.NClass-1
    Results.R_out(EstimOpt.NVarA+7:EstimOpt.NVarA+6+EstimOpt.NVarC,2+ 3*(i-1):1+3*i) = num2cell(Results.DetailsC((i-1)*EstimOpt.NVarC+1:i*EstimOpt.NVarC,:));
end
Results.R_out(EstimOpt.NVarA+6+EstimOpt.NVarC+2,1) = {'Average class probabilities'};
Results.R_out(EstimOpt.NVarA+6+EstimOpt.NVarC+3,1:EstimOpt.NClass) = num2cell(Results.PNClass);
Results.R_out(EstimOpt.NVarA+6+EstimOpt.NVarC+5,1) = {'Model characteristics'};
Results.R_out(EstimOpt.NVarA+6+EstimOpt.NVarC+6:end,1) = {'LL0'; 'LL' ; 'McFadden R2';'Ben-Akiva R2' ;'AIC/n' ; 'n'; 'k'};
Results.R_out(EstimOpt.NVarA+6+EstimOpt.NVarC+6:end,2) = num2cell(Results.stats);

for j = 1:EstimOpt.NClass
	disp(' ')
    disp(num2str(j,'Latent class means %1.0f'));
    disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
    disp([char(EstimOpt.NamesA), blanks(EstimOpt.NVarA)',num2str(Results.DetailsA(((j-1)*EstimOpt.NVarA+1):j*EstimOpt.NVarA,1),'%8.4f'), star_sig(Results.DetailsA(((j-1)*EstimOpt.NVarA+1):j*EstimOpt.NVarA,3)), num2str(Results.DetailsA(((j-1)*EstimOpt.NVarA+1):j*EstimOpt.NVarA,2:3),'%8.4f %8.4f')])
    disp(' ')
    disp(num2str(j,'Latent class standard deviations %1.0f'));
    disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
    disp([char(EstimOpt.NamesA), blanks(EstimOpt.NVarA)',num2str(Results.DetailsV(((j-1)*EstimOpt.NVarA+1):j*EstimOpt.NVarA,1),'%8.4f'), star_sig(Results.DetailsV(((j-1)*EstimOpt.NVarA+1):j*EstimOpt.NVarA,3)), num2str(Results.DetailsV(((j-1)*EstimOpt.NVarA+1):j*EstimOpt.NVarA,2:3),'%8.4f %8.4f')])
end

for j = 1:EstimOpt.NClass -1 
    disp(' ')
    disp(num2str(j,'Latent class probability model %1.0f'))
    disp(['var.', blanks(size(char(EstimOpt.NamesC),2)-2) ,'coef.      st.err.  p-value'])
    disp([char(EstimOpt.NamesC) ,blanks(EstimOpt.NVarC)',num2str(Results.DetailsC(((j-1)*EstimOpt.NVarC+1):j*EstimOpt.NVarC,1),'%8.4f'), star_sig(Results.DetailsC(((j-1)*EstimOpt.NVarC+1):j*EstimOpt.NVarC,3)), num2str(Results.DetailsC(((j-1)*EstimOpt.NVarC+1):j*EstimOpt.NVarC,2:3),'%8.4f %8.4f')])
end

disp(' ')
disp('Avarage class probabilities')
disp('class   prob.')
disp(num2str([(1:EstimOpt.NClass)',Results.PNClass']))

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
Results.tocnote = tocnote;   

