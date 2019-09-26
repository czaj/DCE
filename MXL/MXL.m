function Results = MXL(INPUT,Results_old,EstimOpt,OptimOpt)

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
    cprintf('Black','Estimating '); cprintf('*Black','weighted '); cprintf('Black','MXL model...\n');
else
    disp('Estimating MXL model ...')
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

if isfield(EstimOpt,'NLTVariables') && ~isempty(EstimOpt.NLTVariables)
    EstimOpt.NLTVariables = EstimOpt.NLTVariables(:);
    EstimOpt.NVarNLT = length(unique(EstimOpt.NLTVariables));
    if ~ismember(unique(EstimOpt.NLTVariables),1:NVarA)
        error('Incorrect non-linear variable(s) specification')
    end
    if isfield(EstimOpt,'NLTType') == 0
        cprintf(rgb('DarkOrange'),'WARNING: Assuming Box-Cox transformation \n')
        EstimOpt.NLTType = 1;
    elseif EstimOpt.NLTType == 1
        disp('with Box-Cox transformed variable(s).')
    elseif EstimOpt.NLTType == 2
        disp('with Yeo-Johnson transformed variable(s)')
    else
        error('Incorrect transformation type')
    end
    if EstimOpt.NLTType == 1
        if any(INPUT.Xa(:,EstimOpt.NLTVariables) < 0)
            cprintf(rgb('DarkOrange'),'WARNING: Values of Box-Cox transformed variables < 0 \n')
        elseif any(INPUT.Xa(:,EstimOpt.NLTVariables) == 0) % not sure if this is stil necessary
            cprintf(rgb('DarkOrange'),'WARNING: Values of Box-Cox transformed variables including zeros shifted by 0.00001 \n')
            for i = 1:EstimOpt.NVarNLT
                if any(INPUT.Xa(:,EstimOpt.NLTVariables(i)) == 0)
                    INPUT.Xa(:,EstimOpt.NLTVariables(i)) = INPUT.Xa(:,EstimOpt.NLTVariables(i)) + 0.00001;
                end
            end
        end
    end
else
    EstimOpt.NVarNLT = 0;
    EstimOpt.NLTVariables = [];
    EstimOpt.NLTType = [];
end

if isfield(EstimOpt,'Dist') == 0 || isempty(EstimOpt.Dist)
    EstimOpt.Dist = zeros(1,NVarA);
    if EstimOpt.WTP_space == 0
        cprintf(rgb('DarkOrange'),'WARNING: distributions for random parameters not specified - assuming normality \n')
    else
        cprintf(rgb('DarkOrange'),'WARNING: distributions for random parameters not specified - assuming normality (monetary parameters assumed log-normal) \n')
        EstimOpt.Dist(end-EstimOpt.WTP_space+1:end) = 1; % cost in WTP-space models log-normally distributed
    end
else
    if length(EstimOpt.Dist) == 1
        EstimOpt.Dist = EstimOpt.Dist.*ones(1,NVarA); %needed?
    elseif length(EstimOpt.Dist) == NVarA
        EstimOpt.Dist = EstimOpt.Dist(:)';
    else
        error('Incorrect no. of random parameters'' distributions provided')
    end
end
if isfield(EstimOpt,'Triang') == 0 || length(EstimOpt.Triang) ~= sum(EstimOpt.Dist == 3,2) % Needed only if any parameter has triangular distribution
    EstimOpt.Triang = zeros(1,sum(EstimOpt.Dist == 3,2));
elseif length(EstimOpt.Triang) == 1
    EstimOpt.Triang = EstimOpt.Triang*ones(1,sum(EstimOpt.Dist == 3,2));
else
    EstimOpt.Triang = EstimOpt.Triang(:)';
end

EstimOpt.Johnson = sum(EstimOpt.Dist >= 5);

if (sum(EstimOpt.Dist >= 3 & EstimOpt.Dist <= 5) > 0 && any(find(EstimOpt.Dist >= 3 & EstimOpt.Dist <= 5) > sum(EstimOpt.Dist >= 3 & EstimOpt.Dist <= 5))) && EstimOpt.FullCov == 1
    cprintf(rgb('DarkOrange'),'WARNING: It is recommended to put variables with random parameters with Triangular/Weibull/Sinh-Arcsinh distribution first \n')
end

disp(['Random parameters distributions: ',num2str(EstimOpt.Dist),' (-1 - constant, 0 - normal, 1 - lognormal, 2 - spike, 3 - Triangular, 4  - Weibull, 5 - Sinh-Arcsinh, 6 - Johnson Sb, 7 - Johnson Su)'])

if EstimOpt.WTP_space > 0 && sum(EstimOpt.Dist(end-EstimOpt.WTP_space+1:end)==1) > 0 && any(mean(INPUT.Xa(:,end-EstimOpt.WTP_space+1:end)) >= 0)
    cprintf(rgb('DarkOrange'),'WARNING: Cost attributes with log-normally distributed parameters should enter utility function with a ''-'' sign \n')
end

if isfield(INPUT,'Xs') == 0
    INPUT.Xs = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarS = size(INPUT.Xs,2); % Number of covariates of scale
if isfield(INPUT,'Xm') == 0
    INPUT.Xm = zeros(size(INPUT.Y,1),0);
end
EstimOpt.NVarM = size(INPUT.Xm,2); % Number of covariates of means of random parameters

% This does not currently work:
% if isfield(INPUT, 'Xv') == 0
%     INPUT.Xv = zeros(size(INPUT.Y,1),0);
% end
% EstimOpt.NVarV = size(INPUT.Xv,2); % Number of covariates of variances of random parameters

if EstimOpt.WTP_space > 0
    if isfield(EstimOpt,'WTP_matrix') == 0
        WTP_att = (NVarA-EstimOpt.WTP_space)/EstimOpt.WTP_space;
        if rem(WTP_att,1) ~= 0
            error('EstimOpt.WTP_matrix associating attributes with cost parameters not provided')
        else
            if EstimOpt.WTP_space > 1
                disp(['EstimOpt.WTP_matrix associating attributes with cost parameters not provided - assuming equal shares for each of the ',num2str(EstimOpt.WTP_space),' monetary attributes'])
            end
            EstimOpt.WTP_matrix = NVarA - EstimOpt.WTP_space + kron(1:EstimOpt.WTP_space,ones(1,WTP_att));
        end
    elseif size(EstimOpt.WTP_matrix,2) ~= NVarA - EstimOpt.WTP_space
        error('Dimensions of EstimOpt.WTP_matrix not correct - for each non-monetary attribute provide no. of attribute to multiply it with')
    else
        EstimOpt.WTP_matrix = EstimOpt.WTP_matrix(:)';
    end
end

if isfield(EstimOpt,'Scores') == 0 || isempty(EstimOpt.Scores)
    EstimOpt.Scores = 0;
end

if isfield(EstimOpt,'NamesA') == 0 || isempty(EstimOpt.NamesA) || length(EstimOpt.NamesA) ~= NVarA
    EstimOpt.NamesA = (1:NVarA)';
    EstimOpt.NamesA = cellstr(num2str(EstimOpt.NamesA));
elseif size(EstimOpt.NamesA,1) ~= NVarA
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

NVar = (~EstimOpt.FullCov)*(NVarA*2 + EstimOpt.NVarM*NVarA + EstimOpt.NVarS + EstimOpt.NVarNLT + 2*EstimOpt.Johnson)+EstimOpt.FullCov*(NVarA*(1+EstimOpt.NVarM) + sum(1:NVarA) + EstimOpt.NVarS + EstimOpt.NVarNLT + 2*EstimOpt.Johnson);
if ~isfield(EstimOpt,'ExpB')
    EstimOpt.ExpB = [];
elseif ~isempty(EstimOpt.ExpB)
    EstimOpt.ExpB = EstimOpt.ExpB(:);
    if size(EstimOpt.ExpB,1) == 1
        EstimOpt.ExpB = EstimOpt.ExpB.*ones(NVar,1);
    elseif size(EstimOpt.ExpB,1) ~= NVar
        error('Dimensions of ExpB not correct - provide ExpB indicator for each parameter in B')
    elseif any((EstimOpt.ExpB ~= 0) & (EstimOpt.ExpB ~= 1))
        error('ExpB must only include logical (0 or 1) values')
    end
    EstimOpt.ExpB = (1:NVar)'.* EstimOpt.ExpB;
    EstimOpt.ExpB(EstimOpt.ExpB == 0) = [];
end

%% Starting values

if EstimOpt.FullCov == 0
    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == NVarA*(2+EstimOpt.NVarM) + EstimOpt.NVarS + EstimOpt.NVarNLT + 2*EstimOpt.Johnson
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'MXL_d') && isfield(Results_old.MXL_d,'b0') % starting values provided
        Results_old.MXL_d.b0_old = Results_old.MXL_d.b0(:);
        Results_old.MXL_d = rmfield(Results_old.MXL_d,'b0');
        if length(Results_old.MXL_d.b0_old) ~= NVarA*2 + EstimOpt.NVarM*NVarA + EstimOpt.NVarS + EstimOpt.NVarNLT + 2*EstimOpt.Johnson
            cprintf(rgb('DarkOrange'),'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.MXL_d = rmfield(Results_old.MXL_d,'b0_old');
        else
            b0 = Results_old.MXL_d.b0_old(:);
        end
    end
    if ~exist('b0','var')
        disp('Using MNL results as starting values')
        if ~(isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat') && length(Results_old.MNL.bhat) == (NVarA*(1+EstimOpt.NVarM) + EstimOpt.NVarS + EstimOpt.NVarNLT))
            EstimOpt_tmp = EstimOpt;
            EstimOpt_tmp.Display = 0;
            OptimOpt_tmp = optimoptions('fminunc');
            OptimOpt_tmp.Algorithm = 'quasi-newton';
            OptimOpt_tmp.GradObj = 'off';
            OptimOpt_tmp.Hessian = 'off';
            OptimOpt_tmp.Display = 'off';
            OptimOpt_tmp.FunValCheck= 'off';
            OptimOpt_tmp.Diagnostics = 'off';
            Results_old.MNL = MNL(INPUT,Results_old,EstimOpt_tmp,OptimOpt_tmp);            
        end        
        Results_old.MNL.bhat = Results_old.MNL.bhat(:);
        b0 = [Results_old.MNL.bhat(1:NVarA);max(1,abs(Results_old.MNL.bhat(1:NVarA)));Results_old.MNL.bhat(NVarA+1:end)];
        if sum(EstimOpt.Dist == 1) > 0
            if any(b0(EstimOpt.Dist == 1) < 0)
                cprintf(rgb('DarkOrange'),'WARNING: MNL estimates of log-normally distributed parameters negative - using arbitrary starting values (this may not solve the problem - sign of the attribute may need to be reversed \n')
                %b0(EstimOpt.Dist == 1 & b0(EstimOpt.Dist == 1) < 0) = 1.01;
                b0(b0(1:size(EstimOpt.Dist,2)) < 0 & EstimOpt.Dist' == 1) = 1.01;
            end
            b0(EstimOpt.Dist == 1) = log(b0(EstimOpt.Dist == 1));
        end
        if sum(EstimOpt.Dist == -1) > 0 % Fixed
            indx = find(EstimOpt.Dist == -1);
            b0([indx+NVarA]) = 0;
        end
        if sum(EstimOpt.Dist == 3) > 0 % Triangular
            indx = find(EstimOpt.Dist == 3);
            b0([indx;indx+NVarA]) = [log(b0(indx) - EstimOpt.Triang');log(b0(indx) - EstimOpt.Triang')];
        end
        if sum(EstimOpt.Dist == 4) > 0 % Weibull
            indx = find(EstimOpt.Dist == 4);
            b0([indx;indx+NVarA]) = [log(b0(indx));zeros(length(indx),1)];
        end
        if sum(EstimOpt.Dist >= 5) > 0 % Johnson
            indx = find(EstimOpt.Dist >= 5);
            tmp = [b0(indx);log(b0(indx+NVarA))];
            b0([indx;indx+NVarA]) = [zeros(length(indx),1),ones(length(indx),1)];
            b0 = [b0;tmp];
        end
%         else
%             error('No starting values available - run MNL first')
    end
    
else % EstimOpt.FullCov == 1

    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == NVarA*(1+EstimOpt.NVarM) + sum(1:NVarA) + EstimOpt.NVarS + EstimOpt.NVarNLT + 2*EstimOpt.Johnson
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'MXL') && isfield(Results_old.MXL,'b0') % starting values provided
        Results_old.MXL.b0_old = Results_old.MXL.b0(:);
        Results_old.MXL = rmfield(Results_old.MXL,'b0');
        if length(Results_old.MXL.b0_old) ~= NVarA*(1+EstimOpt.NVarM) + sum(1:NVarA) + EstimOpt.NVarS + EstimOpt.NVarNLT + 2*EstimOpt.Johnson
            cprintf(rgb('DarkOrange'),'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.MXL = rmfield(Results_old.MXL,'b0_old');
        else
            b0 = Results_old.MXL.b0_old;
        end
    end
    if ~exist('b0','var')
        if isfield(Results_old,'MXL_d') && isfield(Results_old.MXL_d,'bhat') && length(Results_old.MXL_d.bhat) == ((2+EstimOpt.NVarM)*NVarA + EstimOpt.NVarS + EstimOpt.NVarNLT + 2*EstimOpt.Johnson)
            disp('Using MXL_d results as starting values')
            Results_old.MXL_d.bhat = Results_old.MXL_d.bhat(:);
            if sum(EstimOpt.Dist >= 3) > 0
                vc_tmp = Results_old.MXL_d.bhat(NVarA+1:NVarA*2);
                vc_tmp(EstimOpt.Dist < 3) = vc_tmp(EstimOpt.Dist < 3).^2;
                vc_tmp = diag(vc_tmp);
            else
                vc_tmp = abs(diag(Results_old.MXL_d.bhat(NVarA+1:NVarA*2)));
            end
            b0 = [Results_old.MXL_d.bhat(1:NVarA);vc_tmp(tril(ones(size(vc_tmp))) == 1);Results_old.MXL_d.bhat(NVarA*2+1:end)];
        else
            disp('Using MNL results as starting values')
            if ~(isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat') && length(Results_old.MNL.bhat) == (NVarA*(1+EstimOpt.NVarM) + EstimOpt.NVarS + EstimOpt.NVarNLT))
                EstimOpt_tmp = EstimOpt;
                EstimOpt_tmp.Display = 0;         
                OptimOpt_tmp = optimoptions('fminunc');
                OptimOpt_tmp.Algorithm = 'quasi-newton';
                OptimOpt_tmp.GradObj = 'off';
                OptimOpt_tmp.Hessian = 'off';
                OptimOpt_tmp.Display = 'off';
                OptimOpt_tmp.FunValCheck= 'off';
                OptimOpt_tmp.Diagnostics = 'off';
                Results_old.MNL = MNL(INPUT,Results_old,EstimOpt_tmp,OptimOpt_tmp);
            end
            Results_old.MNL.bhat = Results_old.MNL.bhat(:);
            b0 = [Results_old.MNL.bhat(1:NVarA);zeros(sum(1:NVarA),1);Results_old.MNL.bhat(NVarA+1:end)];
            if sum(EstimOpt.Dist == 1) > 0
                if any(b0(EstimOpt.Dist == 1) < 0)
                    cprintf(rgb('DarkOrange'),'WARNING: MNL estimates of log-normally distributed parameters negative - using arbitrary starting values (this may not solve the problem - sign of the attribute may need to be reversed \n')
                    %b0(EstimOpt.Dist == 1 & b0(EstimOpt.Dist == 1) < 0) = 1.01;
                    b0(b0(1:size(EstimOpt.Dist,2)) < 0 & EstimOpt.Dist' == 1) = 1.01;
                end
                b0(EstimOpt.Dist == 1) = log(b0(EstimOpt.Dist == 1));
            end
            if sum(EstimOpt.Dist == -1) > 0 % Fixed
                b0(EstimOpt.Dist == -1) = 0;
            end
            if sum(EstimOpt.Dist == 3) > 0 % Triangular
                b0(EstimOpt.Dist == 3) = log(b0(EstimOpt.Dist == 3) -EstimOpt.Triang');
            end
            if sum(EstimOpt.Dist == 4) > 0 % Weibull
                b0(EstimOpt.Dist == 4) = log(b0(EstimOpt.Dist == 4));
            end
            if sum(EstimOpt.Dist >= 5) > 0 % Johnson
                indx = find(EstimOpt.Dist >= 5);
                tmp = b0(indx);
                b0(indx) = zeros(length(indx),1);
                b0 = [b0;tmp;zeros(length(indx),1)];
            end
%         else
%             error('No starting values available')
        end
    end
end

%% Optimization Options

if isfield(EstimOpt,'BActive')
    EstimOpt.BActive = EstimOpt.BActive(:)';
else
    EstimOpt.BActive = ones(1,length(b0));
end

if sum(EstimOpt.Dist == -1) > 0
    if isfield(EstimOpt,'BActive') == 0 || isempty(EstimOpt.BActive)
        EstimOpt.BActive = ones(1,length(b0));
    end
    if EstimOpt.FullCov == 0
        EstimOpt.BActive(NVarA+find(EstimOpt.Dist == -1)) = 0;
    elseif EstimOpt.FullCov == 1
        Vt = tril(ones(NVarA));
        Vt(EstimOpt.Dist == -1,:) = 0;
        EstimOpt.BActive(NVarA+1:NVarA+sum(1:NVarA)) = EstimOpt.BActive(NVarA+1:NVarA+sum(1:NVarA)).*(Vt(tril(ones(size(Vt)))~=0)');
    end
end

if ~isempty(EstimOpt.ExpB) && EstimOpt.NumGrad == 0
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
    if EstimOpt.Display ~= 0
        cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient off - ExpB not supported by analytical gradient \n')
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
    err_mtx = randn(EstimOpt.NP*EstimOpt.NRep,NVarA);
elseif EstimOpt.Draws == 2 % LHS
    cprintf('*blue','Latin Hypercube Sampling '); cprintf('draws \n');
    err_mtx = lhsnorm(zeros((NVarA)*EstimOpt.NP,1),diag(ones((NVarA)*EstimOpt.NP,1)),EstimOpt.NRep);
    err_mtx = reshape(err_mtx,[EstimOpt.NRep*EstimOpt.NP,NVarA]);
elseif EstimOpt.Draws >= 3 % Quasi random draws
    if EstimOpt.Draws == 3
        cprintf('*blue','Halton '); cprintf(['draws (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n']);
        hm1 = haltonset(NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
    elseif EstimOpt.Draws == 4 % apply reverse-radix scrambling
        cprintf('*blue','Halton '); cprintf(['draws with reverse radix scrambling (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n']);
        hm1 = haltonset(NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
        hm1 = scramble(hm1,'RR2');
    elseif EstimOpt.Draws == 5
        cprintf('*blue','Sobol '); cprintf(['draws (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n']);
        hm1 = sobolset(NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
    elseif EstimOpt.Draws == 6
        cprintf('*blue','Sobol '); cprintf(['draws with random linear scramble and random digital shift (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n']);
        hm1 = sobolset(NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
        hm1 = scramble(hm1,'MatousekAffineOwen');
    end
    
    err_mtx = net(hm1,EstimOpt.NP*EstimOpt.NRep); % this takes every point:
    clear hm1;
    %     err_mtx = err_mtx(:,2:NVarA+1);
    
    if EstimOpt.NP*EstimOpt.NRep < 3e+7
        err_mtx = icdf('Normal',err_mtx,0,1); %to be cut down later
    else % this is for very large number of draws * variables
        for i = 1:NVarA
            err_mtx(:,i) = icdf('Normal',err_mtx(:,i),0,1); %to be cut down later
        end
    end
end

err_mtx(:,EstimOpt.Dist == -1) = 0;


%% Display Options


if ((isfield(EstimOpt,'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    OptimOpt.GradObj = 'on';
end

% if EstimOpt.NVarS > 0 && EstimOpt.NumGrad == 0 && any(isnan(INPUT.Xa(:)))
% 	EstimOpt.NumGrad = 1;
%     OptimOpt.GradObj = 'off';
% 	cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient off - covariates of scale not supported by analytical gradient \n')
% end

if any(EstimOpt.Dist > 1) && EstimOpt.NumGrad == 0
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient off - analytical gradient available for normally or lognormally distributed parameters only \n')
end


% This is not necessary any more (?)
% if any(var(EstimOpt.NAltMissInd)) ~= 0 && EstimOpt.NumGrad == 0 && EstimOpt.WTP_space > 1
%     EstimOpt.NumGrad = 1;
%     OptimOpt.GradObj = 'off';
%     cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient off - analytical gradient not available when the number of alternatives differs and WTP_space > 1  \n')
% end

if ((isfield(EstimOpt,'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'on';
end

% if EstimOpt.NVarNLT > 0 && EstimOpt.NLTType == 2 && EstimOpt.NumGrad == 0
% 	EstimOpt.NumGrad = 1;
% 	if EstimOpt.Display ~= 0
%         cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - Yeo-Johnston transformation not supported by analytical gradient \n')
% 	end
% end

if (isfield(EstimOpt,'ConstVarActive') == 0 || EstimOpt.ConstVarActive == 0) && isequal(OptimOpt.Algorithm,'quasi-newton') && isequal(OptimOpt.Hessian,'user-supplied')
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied Hessian off - quasi-newton algorithm does not use it anyway \n')
    OptimOpt.Hessian = 'off';
end

if EstimOpt.NumGrad == 1 && EstimOpt.ApproxHess == 0
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian only available if analythical gradient on \n')
    EstimOpt.ApproxHess = 1;
end

if EstimOpt.NVarS > 0 && (EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4)
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian not available for models with covariates of scale \n')
    EstimOpt.ApproxHess = 1;
    EstimOpt.HessEstFix = 0;
end

if EstimOpt.NVarM > 0 && (EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4)
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian not available for models with covariates of means \n')
    EstimOpt.ApproxHess = 1;
    EstimOpt.HessEstFix = 0;
end

if EstimOpt.WTP_space > 0 && (EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4)
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian not available for models in WTP-space \n')
    EstimOpt.ApproxHess = 1;
    EstimOpt.HessEstFix = 0;
end

if any(isnan(INPUT.Xa(:))) == 1 && (EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4)
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian not available with missing data \n')
    EstimOpt.ApproxHess = 1;
    EstimOpt.HessEstFix = 0;
end

if any(EstimOpt.Dist ~= 0) && (EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4)
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian available for models with normally distributed parameters only \n')
    EstimOpt.ApproxHess = 1;
    EstimOpt.HessEstFix = 0;
end

if EstimOpt.NVarNLT > 0 && (EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4)
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian not available for models with non-linear transformation(s) of variable(s) \n')
    EstimOpt.ApproxHess = 1;
    EstimOpt.HessEstFix = 0;
end

if any(INPUT.W ~= 1) && ((EstimOpt.ApproxHess == 0 && EstimOpt.NumGrad == 0) || EstimOpt.HessEstFix == 4)
    INPUT.W = ones(EstimOpt.NP,1);
    cprintf(rgb('DarkOrange'),'WARNING: Setting all weights to 1, they are not supported with analytical hessian \n')
end

if EstimOpt.RobustStd == 1 && (EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2)
    EstimOpt.RobustStd = 0;
    cprintf(rgb('DarkOrange'),'WARNING: Setting off robust standard errors, they do not matter for BHHH aproximation of hessian \n')
end

if  any(EstimOpt.Dist >= 3 & EstimOpt.Dist <= 5) && EstimOpt.NVarM ~= 0
    error('Covariates of means do not work with triangular/weibull/sinh-arcsinh distributions')
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


INPUT.XXa = reshape(INPUT.Xa,[EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP,NVarA]);
INPUT.XXa = permute(INPUT.XXa,[1 3 2]);
INPUT.YY = reshape(INPUT.Y,[EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP]);

% idx = sum(reshape(INPUT.MissingInd,[EstimOpt.NAlt,EstimOpt.NCT,EstimOpt.NP])) == EstimOpt.NAlt;
% INPUT.YYY(idx(ones(EstimOpt.NAlt,1),:,:)) = NaN; % replace YYY in missing choice-tasks with NaN
% INPUT.YY = reshape(INPUT.YYY,EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP)==1;
%INPUT.YY = reshape(INPUT.YYY,EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP);

INPUT.XXm = reshape(INPUT.Xm',[EstimOpt.NVarM,EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP]);
% INPUT.XXm = squeeze(INPUT.XXm(:,1,:));
% if EstimOpt.NVarM == 1
%     INPUT.XXm = INPUT.XXm';
% end
INPUT.XXm = reshape(INPUT.XXm(:,1,:),[EstimOpt.NVarM,EstimOpt.NP]);


err_mtx = err_mtx';
% change err_mtx from NRep*NP x NVarA to NP*NRep x NVarA (incrasing the no. of draws only adds new draws for each respondent, does not change all draws per individual)
% err_mtx = reshape(permute(reshape(err_mtx,EstimOpt.NP,EstimOpt.NRep,NVarA),[2,1,3]),EstimOpt.NP*EstimOpt.NRep,NVarA)';
% problem - look at the first NRep draws for NVarA=1... all are positive...
if isfield(EstimOpt,'Drawskeep') && ~isempty(EstimOpt.Drawskeep) && EstimOpt.Drawskeep == 1
    Results.err = err_mtx;
end

VC = tril(ones(NVarA));
VC(VC == 1) = (1:(NVarA*(NVarA-1)/2+NVarA))';
EstimOpt.DiagIndex = diag(VC);

% Creating indices for analitical gradient

EstimOpt.indx1 = [];
EstimOpt.indx2 = [];
if EstimOpt.NumGrad == 0 && EstimOpt.FullCov == 1
    for i = 1:NVarA
        EstimOpt.indx1 = [EstimOpt.indx1,i:NVarA];
        EstimOpt.indx2 = [EstimOpt.indx2,i*ones(1,NVarA+1-i)];
    end
end

% save tmp2
% return

% if EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4; %calculations needed for analitical Hessian
%     EstimOpt.XXX = permute(mmx('square',permute(INPUT.XXa,[2,4,1,3]),[]),[3,1,2,4])
%     EstimOpt.XXX = zeros(EstimOpt.NAlt*EstimOpt.NCT,NVarA,NVarA, EstimOpt.NP);
%     if EstimOpt.FullCov == 0
%         EstimOpt.VCx = zeros(NVarA,NVarA,EstimOpt.NRep,EstimOpt.NP);
%     else
%         EstimOpt.VCx = zeros(NVarA*(NVarA-1)/2+NVarA,NVarA*(NVarA-1)/2+NVarA,EstimOpt.NRep,EstimOpt.NP);
%     end
%     err_tmp = reshape(err_mtx,NVarA,EstimOpt.NRep,EstimOpt.NP);
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


LLfun = @(B) LL_mxl_MATlike(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,err_mtx,INPUT.W,EstimOpt,OptimOpt,B);

if EstimOpt.ConstVarActive == 0
    
    if EstimOpt.HessEstFix == 0
        [Results.bhat,LL,Results.exitf,Results.output,Results.g,Results.hess] = fminunc(LLfun,b0,OptimOpt);
    else
        [Results.bhat,LL,Results.exitf,Results.output,Results.g] = fminunc(LLfun,b0,OptimOpt);
    end

    %         options_tmp = optimset('MaxFunEvals',1e100,'MaxIter',1e3,'TolFun',1e-6,'TolX',1e-6,'OutputFcn',@outputf);
    %
    %         [Results.beta, LL,Results.exitf,Results.output] = fminsearch(LLfun, b0,options_tmp);
    %
    %         save tmp1
    
    %     [x,fval,exitflag,output,lambda,grad,hessian] = knitromatlab(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,extendedFeatures,options,KNITROOptions)
    %     [Results.bhat,LL,Results.exitf,Results.output,Results.lambda,Results.g,Results.hess] = knitromatlab(LLfun,b0,[],[],[],[],[],[],[],[],[],'knitro.opt'); %
    %     [Results.bhat,LL,Results.exitf,Results.output,Results.lambda,Results.g,Results.hess] = knitromatlab(LLfun,b0,[],[],[],[],[],[],[],[],OptimOpt,'knitro.opt'); %

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


%% Output


% save tmp1
% return

Results.LL = -LL;
Results.b0_old = b0;

LLfun2 = @(B) LL_mxl(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,err_mtx,EstimOpt,B);

if EstimOpt.HessEstFix == 0 % this will fail if there is no gradient available!
    try
        [Results.LLdetailed,Results.jacobian] = LLfun2(Results.bhat);
        Results.jacobian = Results.jacobian.*INPUT.W;
    catch % theErrorInfo
        Results.LLdetailed = LLfun2(Results.bhat);
        Results.jacobian2 = numdiff(@(B) INPUT.W.*LLfun2(B),Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
        Results.jacobian2 = Results.jacobian2.*INPUT.W;
    end
elseif EstimOpt.HessEstFix == 1
    if isequal(OptimOpt.GradObj,'on') && EstimOpt.NumGrad == 0
        [Results.LLdetailed,Results.jacobian] = LLfun2(Results.bhat);
        Results.jacobian = Results.jacobian.*INPUT.W;
    else
        Results.LLdetailed = LLfun2(Results.bhat);
        Results.jacobian = numdiff(@(B) INPUT.W.*LLfun2(B),Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
        Results.jacobian = Results.jacobian.*INPUT.W;
    end
elseif EstimOpt.HessEstFix == 2
    Results.LLdetailed = LLfun2(Results.bhat);
    Results.jacobian = jacobianest(@(B) INPUT.W.*LLfun2(B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.LLdetailed = LLfun2(Results.bhat);
    Results.hess = hessian(@(B) sum(INPUT.W.*LLfun2(B)),Results.bhat);
elseif EstimOpt.HessEstFix == 4
    [Results.LLdetailed,~,Results.hess] = LLfun2(Results.bhat);
    % no weighting?
end
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

if EstimOpt.RobustStd == 1
    if ~isfield(Results,'jacobian')
        if EstimOpt.NumGrad == 0
            [~,Results.jacobian] = LLfun2(Results.bhat);
            Results.jacobian = Results.jacobian.*INPUT.W(:,ones(1,size(Results.jacobian,2)));
        else
            Results.jacobian = numdiff(@(B) INPUT.W.*LLfun2(B),Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
        end
    end
    RobustHess = Results.jacobian'*Results.jacobian;
    Results.ihess = Results.ihess*RobustHess*Results.ihess;
end

Results.std = sqrt(diag(Results.ihess));
Results.std(EstimOpt.BActive == 0) = NaN;
Results.std(EstimOpt.BLimit == 1) = 0;
Results.std(imag(Results.std) ~= 0) = NaN;

if any(INPUT.MissingInd == 1) % In case of some missing data
    idx = sum(reshape(INPUT.MissingInd,[EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP])) == EstimOpt.NAlt;
    idx = sum(reshape(idx,[EstimOpt.NCT,EstimOpt.NP]),1)'; % no. of missing NCT for every respondent
    idx = EstimOpt.NCT - idx;
    R2 = mean(exp(-Results.LLdetailed./idx),1);
    Results.CrossEntropy = mean(Results.LLdetailed./idx,1);
else
    R2 = mean(exp(-Results.LLdetailed/EstimOpt.NCT),1);
    Results.CrossEntropy = mean(Results.LLdetailed./EstimOpt.NCT,1);
end

if EstimOpt.Scores ~= 0
    Results.Scores = BayesScoresMXL(INPUT.YY,INPUT.XXa,INPUT.XXm,INPUT.Xs,err_mtx,EstimOpt,Results.bhat);
end

% save out_MXL1
% return

if EstimOpt.FullCov == 0
    Results.DetailsA(1:NVarA,1) = Results.bhat(1:NVarA);
    Results.DetailsA(1:NVarA,3:4) = [Results.std(1:NVarA),pv(Results.bhat(1:NVarA),Results.std(1:NVarA))];
    Results.DetailsV(1:NVarA,1) = abs(Results.bhat(NVarA+1:NVarA*2));
    Results.DetailsV(1:NVarA,3:4) = [Results.std(NVarA+1:NVarA*2),pv(Results.bhat(NVarA+1:NVarA*2),Results.std(NVarA+1:NVarA*2))];
    if EstimOpt.NVarM > 0
        Results.DetailsM = [];
        for i=1:EstimOpt.NVarM
            Results.DetailsM(1:NVarA,4*i-3) = Results.bhat(NVarA*(2+i-1)+1:NVarA*(2+i));
            Results.DetailsM(1:NVarA,4*i-1:4*i) = [Results.std(NVarA*(2+i-1)+1:NVarA*(2+i)),pv(Results.bhat(NVarA*(2+i-1)+1:NVarA*(2+i)),Results.std(NVarA*(2+i-1)+1:NVarA*(2+i)))];
        end
    end
    
    %     if any(EstimOpt.Dist == 1) %     transform normal to lognormal for display (delta):
    %         log_idx = find(EstimOpt.Dist==1);
    %         Results.DetailsA(log_idx,1) = exp(Results.DetailsA(log_idx,1) + Results.DetailsV(log_idx,1).^2./2);
    %         Results.std_lognormal = Results.std;
    %         for i = 1:size(log_idx,2)
    %             Results.std_lognormal(log_idx(i)) = (exp(2*Results.bhat(log_idx(i)) + Results.DetailsV(log_idx(i),1).^2).*...
    %                 (Results.ihess(log_idx(i),log_idx(i)) + Results.DetailsV(log_idx(i),1).*(2.*Results.ihess(log_idx(i),NVarA + log_idx(i)) + Results.DetailsV(log_idx(i),1).*Results.ihess(NVarA + log_idx(i),NVarA + log_idx(i)))))^0.5;
    %             Results.std_lognormal(NVarA + log_idx(i)) = exp(-1 + 2*abs(Results.bhat(log_idx(i))) + 2*Results.DetailsV(log_idx(i),1).^2).*...
    %                 (Results.ihess(log_idx(i),log_idx(i)) + Results.DetailsV(log_idx(i),1).*(4.*Results.ihess(log_idx(i),NVarA + log_idx(i)) + 4.*Results.DetailsV(log_idx(i),1).*Results.ihess(NVarA + log_idx(i),NVarA + log_idx(i)))).^2;
    %         end
    %         Results.DetailsA(log_idx,3:4) = [Results.std_lognormal(log_idx),pv(Results.DetailsA(log_idx,1),Results.std_lognormal(log_idx))];
    %         Results.DetailsV(log_idx,1) = ((exp(Results.DetailsV(log_idx,1).^2) - 1).*exp(2.*Results.bhat(log_idx) + Results.DetailsV(log_idx,1).^2)).^0.5;
    %         Results.DetailsV(log_idx,3:4) = [Results.std_lognormal(NVarA + log_idx),pv(Results.DetailsV(log_idx,1),Results.std_lognormal(NVarA + log_idx))];
    %     end
    
    if isfield(EstimOpt,'EffectiveMoments') && EstimOpt.EffectiveMoments == 1
    if any(EstimOpt.Dist == 1) %     transform normal to lognormal for display (simulation):
        Results.DetailsA_underlying = Results.DetailsA;
        Results.DetailsV_underlying = Results.DetailsV;
        if EstimOpt.NVarM > 0
            Results.DetailsM_underlying = Results.DetailsM;
        end
        log_idx = find(EstimOpt.Dist==1);
        try % in case Results.ihess is not positive semidefinite (avoid mvnrnd error)
            bhat_sim = mvnrnd(Results.bhat(1:NVarA*(2+EstimOpt.NVarM)),Results.ihess(1:NVarA*(2+EstimOpt.NVarM),1:NVarA*(2+EstimOpt.NVarM)),NSdSim)';
            if EstimOpt.NVarM > 0
                BM_i = permute(reshape(bhat_sim(NVarA*2+1:NVarA*2+NVarA*EstimOpt.NVarM,:),[NVarA,EstimOpt.NVarM,NSdSim]),[3,1,2]);
                bhat_new = zeros(NSdSim,5*NVarA);
                bhatM_new = zeros(NSdSim,5*NVarA*EstimOpt.NVarM);
                parfor i = 1:NSdSim
                    bhat_i = bhat_sim(:,i);
                    B0_i = mvnrnd(bhat_i(1:NVarA)',diag(bhat_i(NVarA+1:NVarA*2).^2),NSdSim);
                    B0exp_i = B0_i;
                    B0exp_i(:,log_idx) = exp(B0exp_i(:,log_idx));
                    bhat_new(i,:) = [mean(B0exp_i),median(B0exp_i),std(B0exp_i),quantile(B0exp_i,0.025),quantile(B0exp_i,0.975)];
                    
                    B0M_i = BM_i;
                    B0M_i(:,log_idx) = exp(B0M_i(:,log_idx) + B0_i(:,log_idx)) - B0exp_i(:,log_idx);                    
                    bhatM_new(i,:) = [mean(B0M_i),median(B0M_i),std(B0M_i),quantile(B0M_i,0.025),quantile(B0M_i,0.975)];    
                end
            else
                bhat_new = zeros(NSdSim,5*NVarA);
                parfor i = 1:NSdSim
                    bhat_i = bhat_sim(:,i);
                    B0_i = mvnrnd(bhat_i(1:NVarA)',diag(bhat_i(NVarA+1:NVarA*2).^2),NSdSim);
                    %     b0v = bhat_i(NVarA+1:end)';
                    %     VC = tril(ones(NVarA));
                    %     VC(VC==1) = b0v;
                    %     VC = VC*VC';
                    %     B0_i = mvnrnd(bhat_i(1:NVarA)',VC,sim2);
                    B0_i(:,log_idx) = exp(B0_i(:,log_idx));
                    bhat_new(i,:) = [mean(B0_i),median(B0_i),std(B0_i),quantile(B0_i,0.025),quantile(B0_i,0.975)];
                end
            end
            Results.DistStats_Coef = reshape(median(bhat_new,1),[NVarA,5]);
            Results.DistStats_Std = reshape(std(bhat_new,[],1),[NVarA,5]);
            Results.DistStats_Q025 = reshape(quantile(bhat_new,0.025),[NVarA,5]);
            Results.DistStats_Q975 = reshape(quantile(bhat_new,0.975),[NVarA,5]);            
            Results.DetailsA(log_idx,1) = exp(Results.DetailsA(log_idx,1) + Results.DetailsV(log_idx,1).^2./2); % We use analytical formula instead of simulated moment
            Results.DetailsA(log_idx,3:4) = [Results.DistStats_Std(log_idx,1),pv(Results.DetailsA(log_idx,1),Results.DistStats_Std(log_idx,1))];
            Results.DetailsV(log_idx,1) = ((exp(Results.DetailsV(log_idx,1).^2) - 1).*exp(2.*Results.bhat(log_idx) + Results.DetailsV(log_idx,1).^2)).^0.5; % We use analytical formula instead of simulated moment
            Results.DetailsV(log_idx,3:4) = [Results.DistStats_Std(log_idx,3),pv(Results.DetailsV(log_idx,1),Results.DistStats_Std(log_idx,3))];
            if EstimOpt.NVarM > 0
                Results.DistStatsM_Coef = reshape(median(bhatM_new,1),[NVarA,5]);
                Results.DistStatsM_Std = reshape(std(bhatM_new,[],1),[NVarA,5]);
                Results.DistStatsM_Q025 = reshape(quantile(bhatM_new,0.025),[NVarA,5]);
                Results.DistStatsM_Q975 = reshape(quantile(bhatM_new,0.975),[NVarA,5]);
%                 Results.DetailsM(log_idx,1) = Results.DistStatsM_Coef(log_idx,1); 
                Results.DetailsM(log_idx,1) = exp(Results.DetailsA_underlying(log_idx,1) + Results.DetailsM(log_idx,1) + Results.DetailsV_underlying(log_idx,1).^2./2) - Results.DetailsA(log_idx,1); % We use analytical formula instead of simulated moments
                Results.DetailsM(log_idx,3:4) = [Results.DistStatsM_Std(log_idx,1),pv(Results.DetailsM(log_idx,1),Results.DistStatsM_Std(log_idx,1))];
            end
        catch % theErrorInfo
            Results.DetailsA(log_idx,[1,3:4]) = NaN;
            Results.DetailsV(log_idx,[1,3:4]) = NaN;
        end
    end
    end
    
    if sum(EstimOpt.Dist == 3) > 0
        Results.DetailsA(EstimOpt.Dist == 3,1) = exp(Results.bhat(EstimOpt.Dist == 3)) + EstimOpt.Triang';
        Results.DetailsA(EstimOpt.Dist == 3,3:4) = [exp(Results.bhat(EstimOpt.Dist == 3)).*Results.std(EstimOpt.Dist == 3),pv(exp(Results.bhat(EstimOpt.Dist == 3)) + EstimOpt.Triang',exp(Results.bhat(EstimOpt.Dist == 3)).*Results.std(EstimOpt.Dist == 3))];
        btmp = Results.bhat(NVarA+1:NVarA*2);
        stdx = zeros(sum(EstimOpt.Dist == 3),1);
        g = [exp(Results.bhat(EstimOpt.Dist == 3)),exp(btmp(EstimOpt.Dist == 3))];
        indx = find(EstimOpt.Dist == 3);
        for i = 1:sum(EstimOpt.Dist == 3)
            stdx(i) = sqrt(g(i,:)*Results.ihess([indx(i),indx(i)+NVarA],[indx(i),indx(i)+NVarA])*g(i,:)');
        end
        %Results.DetailsV(EstimOpt.Dist == 3,1) = exp(btmp(EstimOpt.Dist == 3));
        Results.DetailsV(EstimOpt.Dist == 3,1) = exp(btmp(EstimOpt.Dist == 3)) + exp(Results.bhat(EstimOpt.Dist == 3)) + EstimOpt.Triang';
        %Results.DetailsV(EstimOpt.Dist == 3,3:4) = [stdx(EstimOpt.Dist == 3),pv(exp(btmp(EstimOpt.Dist == 3)),stdx(EstimOpt.Dist == 3))];
        Results.DetailsV(EstimOpt.Dist == 3,3:4) = [stdx,pv(exp(btmp(EstimOpt.Dist == 3))+ exp(Results.bhat(EstimOpt.Dist == 3)) + EstimOpt.Triang',stdx)];
    end
    if sum(EstimOpt.Dist == 4) > 0
        Results.DetailsA(EstimOpt.Dist == 4,1) = exp(Results.bhat(EstimOpt.Dist == 4));
        Results.DetailsA(EstimOpt.Dist == 4,3:4) = [exp(Results.bhat(EstimOpt.Dist == 4)).*Results.std(EstimOpt.Dist == 4),pv(exp(Results.bhat(EstimOpt.Dist == 4)),exp(Results.bhat(EstimOpt.Dist == 4)).*Results.std(EstimOpt.Dist == 4))];
        btmp = Results.bhat(NVarA+1:NVarA*2);
        stdx = exp(btmp).*Results.std(NVarA+1:NVarA*2);
        Results.DetailsV(EstimOpt.Dist == 4,1) = exp(btmp(EstimOpt.Dist == 4));
        Results.DetailsV(EstimOpt.Dist == 4,3:4) = [stdx(EstimOpt.Dist == 4),pv(exp(btmp(EstimOpt.Dist == 4)),stdx(EstimOpt.Dist == 4))];
    end
    Results.R = [Results.DetailsA,Results.DetailsV];
    if EstimOpt.NVarM > 0
%         Results.DetailsM = [];
%         for i=1:EstimOpt.NVarM
%             Results.DetailsM(1:NVarA,4*i-3) = Results.bhat(NVarA*(2+i-1)+1:NVarA*(2+i));
%             Results.DetailsM(1:NVarA,4*i-1:4*i) = [Results.std(NVarA*(2+i-1)+1:NVarA*(2+i)),pv(Results.bhat(NVarA*(2+i-1)+1:NVarA*(2+i)),Results.std(NVarA*(2+i-1)+1:NVarA*(2+i)))];
%         end
        Results.R = [Results.R,Results.DetailsM];
    end
    if EstimOpt.NVarNLT > 0
        Results.DetailsNLT = [];
        for i=1:EstimOpt.NVarNLT
            Results.DetailsNLT(i,1) = Results.bhat(NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS+i);
            Results.DetailsNLT(i,3:4) = [Results.std(NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS+i),pv(Results.bhat(NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS+i),Results.std(NVarA*(2+EstimOpt.NVarM)+EstimOpt.NVarS+i))];
        end
        Results.DetailsNLT0 = NaN(NVarA,4);
        Results.DetailsNLT0(EstimOpt.NLTVariables,:) = Results.DetailsNLT;
        Results.R = [Results.R,Results.DetailsNLT0];
    end
    if EstimOpt.Johnson > 0
        Results.ResultsJ = NaN(NVarA,8);
        % Location parameters
        Results.DetailsJL(:,1) = Results.bhat((end-2*EstimOpt.Johnson+1):(end-EstimOpt.Johnson));
        Results.DetailsJL(:,3) = Results.std((end-2*EstimOpt.Johnson+1):(end-EstimOpt.Johnson));
        Results.DetailsJL(:,4) = pv(Results.bhat((end-2*EstimOpt.Johnson+1):(end-EstimOpt.Johnson)),Results.std((end - 2*EstimOpt.Johnson+1):(end - EstimOpt.Johnson)));
        % Scale parameters
        Results.DetailsJS(:,1) = exp(Results.bhat((end-EstimOpt.Johnson+1):end));
        Results.DetailsJS(:,3) = exp(Results.bhat((end-EstimOpt.Johnson+1):end)).*Results.std((end-EstimOpt.Johnson+1):end);
        Results.DetailsJS(:,4) = pv(exp(Results.bhat((end-EstimOpt.Johnson+1):end)),exp(Results.bhat((end-EstimOpt.Johnson+1):end)).*Results.std((end-EstimOpt.Johnson+1):end));
        Results.ResultsJ(EstimOpt.Dist > 4 & EstimOpt.Dist <= 7,1:4) = Results.DetailsJL;
        Results.ResultsJ(EstimOpt.Dist > 4 & EstimOpt.Dist <= 7,5:8) = Results.DetailsJS;
        Results.R = [Results.R,Results.ResultsJ];
    end
    if EstimOpt.NVarS > 0
        Results.DetailsS = [];
        for i=1:EstimOpt.NVarS
            Results.DetailsS(i,1) = Results.bhat(NVarA*(2+EstimOpt.NVarM)+i);
            Results.DetailsS(i,3:4) = [Results.std(NVarA*(2+EstimOpt.NVarM)+i),pv(Results.bhat(NVarA*(2+EstimOpt.NVarM)+i),Results.std(NVarA*(2+EstimOpt.NVarM)+i))];
        end
        DetailsS0 = NaN(EstimOpt.NVarS,4);
        DetailsS0(1:EstimOpt.NVarS,1:4) = Results.DetailsS;
        
        if EstimOpt.NVarS <= NVarA % will not work if NVarS > NVarA
            Results.R = [Results.R;[DetailsS0,NaN(size(DetailsS0,1),size(Results.R,2)-size(DetailsS0,2))]];
        end
    end
    
elseif EstimOpt.FullCov == 1
    Results.DetailsA(1:NVarA,1) = Results.bhat(1:NVarA);
    Results.DetailsA(1:NVarA,3:4) = [Results.std(1:NVarA),pv(Results.bhat(1:NVarA),Results.std(1:NVarA))];
    Results.DetailsV = sdtri(Results.bhat(NVarA+1:NVarA*(NVarA+3)/2),Results.ihess(NVarA+1:NVarA*(NVarA+3)/2,NVarA+1:NVarA*(NVarA+3)/2),EstimOpt);
    Results.DetailsV = [Results.DetailsV(:,1),zeros(NVarA,1),Results.DetailsV(:,2:3)];
    if EstimOpt.NVarM > 0
        Results.DetailsM = [];
        for i=1:EstimOpt.NVarM
            Results.DetailsM(1:NVarA,4*i-3) = Results.bhat(NVarA*(NVarA/2+0.5+i)+1:NVarA*(NVarA/2+1.5+i));
            Results.DetailsM(1:NVarA,4*i-1:4*i) = [Results.std(NVarA*(NVarA/2+0.5+i)+1:NVarA*(NVarA/2+1.5+i)),pv(Results.bhat(NVarA*(NVarA/2+0.5+i)+1:NVarA*(NVarA/2+1.5+i)),Results.std(NVarA*(NVarA/2+0.5+i)+1:NVarA*(NVarA/2+1.5+i)))];
        end    
    end
    
    if isfield(EstimOpt,'EffectiveMoments') && EstimOpt.EffectiveMoments == 1
    if any(EstimOpt.Dist == 1) %     transform normal to lognormal for display (simulation):
        Results.DetailsA_underlying = Results.DetailsA;
        Results.DetailsV_underlying = Results.DetailsV;
        if EstimOpt.NVarM > 0
            Results.DetailsM_underlying = Results.DetailsM;
        end
        log_idx = find(EstimOpt.Dist==1);
        try % in case Results.ihess is not positive semidefinite (avoid mvnrnd error)
            bhat_sim = mvnrnd(Results.bhat,Results.ihess,NSdSim)';
            if EstimOpt.NVarM > 0
                BM_i = permute(reshape(bhat_sim(NVarA+sum(1:NVarA)+1:NVarA+sum(1:NVarA)+NVarA*EstimOpt.NVarM,:),[NVarA,EstimOpt.NVarM,NSdSim]),[3,1,2]);                
                bhat_new = zeros(NSdSim,5*NVarA);
                bhatM_new = zeros(NSdSim,5*NVarA*EstimOpt.NVarM);
                parfor i = 1:NSdSim
                    bhat_i = bhat_sim(:,i);
                    b0v = bhat_i(NVarA+1:NVarA+sum(1:NVarA))';
                    VC = tril(ones(NVarA));
                    VC(VC == 1) = b0v;
                    VC = VC*VC';
                    B0_i = mvnrnd(bhat_i(1:NVarA)',VC,NSdSim);
                    B0exp_i = B0_i;
                    B0exp_i(:,log_idx) = exp(B0exp_i(:,log_idx));                                        
                    bhat_new(i,:) = [mean(B0exp_i),median(B0exp_i),std(B0exp_i),quantile(B0exp_i,0.025),quantile(B0exp_i,0.975)];
                    
                    B0M_i = BM_i;
                    B0M_i(:,log_idx) = exp(B0M_i(:,log_idx) + B0_i(:,log_idx)) - B0exp_i(:,log_idx);                    
                    bhatM_new(i,:) = [mean(B0M_i),median(B0M_i),std(B0M_i),quantile(B0M_i,0.025),quantile(B0M_i,0.975)];                    
                end
            else
                bhat_new = zeros(NSdSim,5*NVarA);
                parfor i = 1:NSdSim
                    bhat_i = bhat_sim(:,i);
                    %             B0_i = mvnrnd(bhat_i(1:NVarA)',diag(bhat_i(NVarA+1:NVarA*2).^2),NSdSim);
                    b0v = bhat_i(NVarA+1:NVarA+sum(1:NVarA))';
                    VC = tril(ones(NVarA));
                    VC(VC == 1) = b0v;
                    VC = VC*VC';
                    B0_i = mvnrnd(bhat_i(1:NVarA)',VC,NSdSim);
                    B0_i(:,log_idx) = exp(B0_i(:,log_idx));
                    bhat_new(i,:) = [mean(B0_i),median(B0_i),std(B0_i),quantile(B0_i,0.025),quantile(B0_i,0.975)];
                end
            end

            Results.DistStats_Coef = reshape(median(bhat_new,1),[NVarA,5]);
            Results.DistStats_Std = reshape(std(bhat_new,[],1),[NVarA,5]);
            Results.DistStats_Q025 = reshape(quantile(bhat_new,0.025),[NVarA,5]);
            Results.DistStats_Q975 = reshape(quantile(bhat_new,0.975),[NVarA,5]);
            Results.DetailsA(log_idx,1) = exp(Results.DetailsA(log_idx,1) + Results.DetailsV(log_idx,1).^2./2); % We use analytical formula instead of simulated moments
            Results.DetailsA(log_idx,3:4) = [Results.DistStats_Std(log_idx,1),pv(Results.DetailsA(log_idx,1),Results.DistStats_Std(log_idx,1))];
            Results.DetailsV(log_idx,1) = ((exp(Results.DetailsV(log_idx,1).^2) - 1).*exp(2.*Results.bhat(log_idx) + Results.DetailsV(log_idx,1).^2)).^0.5; % We use analytical formula instead of simulated moments
            Results.DetailsV(log_idx,3:4) = [Results.DistStats_Std(log_idx,3),pv(Results.DetailsV(log_idx,1),Results.DistStats_Std(log_idx,3))];
            if EstimOpt.NVarM > 0
                Results.DistStatsM_Coef = reshape(median(bhatM_new,1),[NVarA,5]);
                Results.DistStatsM_Std = reshape(std(bhatM_new,[],1),[NVarA,5]);
                Results.DistStatsM_Q025 = reshape(quantile(bhatM_new,0.025),[NVarA,5]);
                Results.DistStatsM_Q975 = reshape(quantile(bhatM_new,0.975),[NVarA,5]);
%                 Results.DetailsM(log_idx,1) = Results.DistStatsM_Coef(log_idx,1);
                Results.DetailsM(log_idx,1) = exp(Results.DetailsA_underlying(log_idx,1) + Results.DetailsM(log_idx,1) + Results.DetailsV_underlying(log_idx,1).^2./2) - Results.DetailsA(log_idx,1); % We use analytical formula instead of simulated moments
                Results.DetailsM(log_idx,3:4) = [Results.DistStatsM_Std(log_idx,1),pv(Results.DetailsM(log_idx,1),Results.DistStatsM_Std(log_idx,1))];
            end
            
        catch % theErrorInfo
            Results.DetailsA(log_idx,[1,3:4]) = NaN;
            Results.DetailsV(log_idx,[1,3:4]) = NaN;
        end
    end
    end
    if sum(EstimOpt.Dist == 3) > 0
        Results.DetailsA(EstimOpt.Dist == 3,1) = exp(Results.bhat(EstimOpt.Dist == 3)) + EstimOpt.Triang';
        Results.DetailsA(EstimOpt.Dist == 3,3:4) = [exp(Results.bhat(EstimOpt.Dist == 3)).*Results.std(EstimOpt.Dist == 3),pv(exp(Results.bhat(EstimOpt.Dist == 3)) + EstimOpt.Triang',exp(Results.bhat(EstimOpt.Dist == 3)).*Results.std(EstimOpt.Dist == 3))];
        btmp = Results.bhat(NVarA+1:NVarA*(NVarA-1)/2+2*NVarA);
        btmp = btmp(EstimOpt.DiagIndex);
        stdx = zeros(sum(EstimOpt.Dist == 3),1);
        g = [exp(Results.bhat(EstimOpt.Dist == 3)),exp(btmp(EstimOpt.Dist == 3))];
        indx = find(EstimOpt.Dist == 3);
        DiagIndex = EstimOpt.DiagIndex(EstimOpt.Dist == 3);
        for i = 1:sum(EstimOpt.Dist == 3)
            stdx(i) = sqrt(g(i,:)*Results.ihess([indx(i),DiagIndex(i)+NVarA],[indx(i),DiagIndex(i)+NVarA])*g(i,:)');
        end
        Results.DetailsV(EstimOpt.Dist == 3,1) = exp(btmp(EstimOpt.Dist == 3)) + exp(Results.bhat(EstimOpt.Dist == 3)) + EstimOpt.Triang';
        Results.DetailsV(EstimOpt.Dist == 3,3:4) = [stdx,pv(exp(btmp(EstimOpt.Dist == 3))+exp(Results.bhat(EstimOpt.Dist == 3))+EstimOpt.Triang',stdx)];
    end
    if sum(EstimOpt.Dist == 4) > 0
        Results.DetailsA(EstimOpt.Dist == 4,1) = exp(Results.bhat(EstimOpt.Dist == 4));
        Results.DetailsA(EstimOpt.Dist == 4,3:4) = [exp(Results.bhat(EstimOpt.Dist == 4)).*Results.std(EstimOpt.Dist == 4),pv(exp(Results.bhat(EstimOpt.Dist == 4)),exp(Results.bhat(EstimOpt.Dist == 4)).*Results.std(EstimOpt.Dist == 4))];
        btmp = Results.bhat(NVarA+1:NVarA*(NVarA-1)/2+2*NVarA);
        btmp = btmp(EstimOpt.DiagIndex);
        stdx = Results.std(NVarA+1:NVarA*(NVarA-1)/2+2*NVarA);
        stdx = stdx(EstimOpt.DiagIndex);
        stdx = exp(btmp).*stdx;
        Results.DetailsV(EstimOpt.Dist == 4,1) = exp(btmp(EstimOpt.Dist == 4));
        Results.DetailsV(EstimOpt.Dist == 4,3:4) = [stdx(EstimOpt.Dist == 4),pv(exp(btmp(EstimOpt.Dist == 4)),stdx(EstimOpt.Dist == 4))];
    end
    if sum(EstimOpt.Dist == 5) > 0
        btmp = Results.bhat(NVarA+1:NVarA*(NVarA-1)/2+2*NVarA);
        btmp = btmp(EstimOpt.DiagIndex);
        stdtmp = Results.std(NVarA+1:NVarA*(NVarA-1)/2+2*NVarA);
        stdtmp = stdtmp(EstimOpt.DiagIndex);
        Results.DetailsV(EstimOpt.Dist == 5,1) = btmp(EstimOpt.Dist == 5).^2;
        Results.DetailsV(EstimOpt.Dist == 5,3:4) = [2*btmp(EstimOpt.Dist == 5).*stdtmp(EstimOpt.Dist == 5),pv(btmp(EstimOpt.Dist == 5).^2,2*btmp(EstimOpt.Dist == 5).*stdtmp(EstimOpt.Dist == 5))];
    end
    Results.R = [Results.DetailsA,Results.DetailsV];
    if EstimOpt.NVarM > 0
%         Results.DetailsM = [];
%         for i=1:EstimOpt.NVarM
%             Results.DetailsM(1:NVarA,4*i-3) = Results.bhat(NVarA*(NVarA/2+0.5+i)+1:NVarA*(NVarA/2+1.5+i));
%             Results.DetailsM(1:NVarA,4*i-1:4*i) = [Results.std(NVarA*(NVarA/2+0.5+i)+1:NVarA*(NVarA/2+1.5+i)),pv(Results.bhat(NVarA*(NVarA/2+0.5+i)+1:NVarA*(NVarA/2+1.5+i)),Results.std(NVarA*(NVarA/2+0.5+i)+1:NVarA*(NVarA/2+1.5+i)))];
%         end
        Results.R = [Results.R,Results.DetailsM];
    end
    if EstimOpt.NVarNLT > 0
        Results.DetailsNLT = [];
        for i=1:EstimOpt.NVarNLT
            Results.DetailsNLT(i,1) = Results.bhat(NVarA+sum(1:NVarA)+EstimOpt.NVarM+EstimOpt.NVarS+i);
            Results.DetailsNLT(i,3:4) = [Results.std(NVarA+sum(1:NVarA)+EstimOpt.NVarM+EstimOpt.NVarS+i),pv(Results.bhat(NVarA+sum(1:NVarA)+EstimOpt.NVarM+EstimOpt.NVarS+i),Results.std(NVarA+sum(1:NVarA)+EstimOpt.NVarM+EstimOpt.NVarS+i))];
        end
        Results.DetailsNLT0 = NaN(NVarA,4);
        Results.DetailsNLT0(EstimOpt.NLTVariables,:) = Results.DetailsNLT;
        Results.R = [Results.R,Results.DetailsNLT0];
    end
    if EstimOpt.Johnson > 0
        Results.ResultsJ = NaN(NVarA,8);
        % Location parameters
        Results.DetailsJL(:,1) = Results.bhat((end-2*EstimOpt.Johnson+1):(end-EstimOpt.Johnson));
        Results.DetailsJL(:,3) = Results.std((end-2*EstimOpt.Johnson+1):(end-EstimOpt.Johnson));
        Results.DetailsJL(:,4) = pv(Results.bhat((end-2*EstimOpt.Johnson+1):(end-EstimOpt.Johnson)),Results.std((end-2*EstimOpt.Johnson+1):(end-EstimOpt.Johnson)));
        % Scale parameters
        Results.DetailsJS(:,1) = exp(Results.bhat((end-EstimOpt.Johnson+1):end));
        Results.DetailsJS(:,3) = exp(Results.bhat((end-EstimOpt.Johnson+1):end)).*Results.std((end-EstimOpt.Johnson+1):end);
        Results.DetailsJS(:,4) = pv(exp(Results.bhat((end - EstimOpt.Johnson+1):end)),exp(Results.bhat((end-EstimOpt.Johnson+1):end)).*Results.std((end-EstimOpt.Johnson+1):end));
        Results.ResultsJ(EstimOpt.Dist > 4 & EstimOpt.Dist <= 7,1:4) = Results.DetailsJL;
        Results.ResultsJ(EstimOpt.Dist > 4 & EstimOpt.Dist <= 7,5:8) = Results.DetailsJS;
        Results.R = [Results.R,Results.ResultsJ];
    end
    if EstimOpt.NVarS > 0
        Results.DetailsS = [];
        for i=1:EstimOpt.NVarS
            Results.DetailsS(i,1) = Results.bhat(NVarA*(NVarA/2+1.5+EstimOpt.NVarM)+i);
            Results.DetailsS(i,3:4) = [Results.std(NVarA*(NVarA/2+1.5+EstimOpt.NVarM)+i),pv(Results.bhat(NVarA*(NVarA/2+1.5+EstimOpt.NVarM)+i),Results.std(NVarA*(NVarA/2+1.5+EstimOpt.NVarM)+i))];
        end
        DetailsS0 = NaN(EstimOpt.NVarS,4);
        DetailsS0(1:EstimOpt.NVarS,1:4) = Results.DetailsS;
        if EstimOpt.NVarS == NVarA % will not work if NVarS > NVarA
            Results.R = [Results.R,DetailsS0];
        end
    end
    Results.chol = [Results.bhat(NVarA+1:NVarA*(NVarA/2+1.5)),Results.std(NVarA+1:NVarA*(NVarA/2+1.5)),pv(Results.bhat(NVarA+1:NVarA*(NVarA/2+1.5)),Results.std(NVarA+1:NVarA*(NVarA/2+1.5)))];
    Results.DetailsVcov = tril(ones(NVarA));
    choltmp = Results.chol(:,1);
    if sum(EstimOpt.Dist >= 3 & EstimOpt.Dist <= 5) > 0
        choltmp(EstimOpt.DiagIndex(EstimOpt.Dist >= 3 & EstimOpt.Dist <= 5)) = 1;
    end
    Results.DetailsVcov(Results.DetailsVcov == 1) = choltmp;
    if sum(EstimOpt.Dist >= 3 & EstimOpt.Dist <= 5) > 0
        choltmp = sqrt(sum(Results.DetailsVcov(EstimOpt.Dist >= 3 & EstimOpt.Dist <= 5,:).^2,2));
        Results.DetailsVcov(EstimOpt.Dist >= 3 & EstimOpt.Dist <= 5,:) = Results.DetailsVcov(EstimOpt.Dist >= 3 & EstimOpt.Dist <= 5,:)./choltmp(:,ones(1,NVarA));
    end
    Results.DetailsVcov = Results.DetailsVcov*Results.DetailsVcov';
    Results.DetailsVcor = corrcov(Results.DetailsVcov);
end

EstimOpt.params = length(b0) - sum(EstimOpt.BActive == 0) + sum(EstimOpt.BLimit == 1);
Results.stats = [Results.LL;Results_old.MNL0.LL;1-Results.LL/Results_old.MNL0.LL;R2;((2*EstimOpt.params-2*Results.LL))/EstimOpt.NObs;((log(EstimOpt.NObs)*EstimOpt.params-2*Results.LL))/EstimOpt.NObs;EstimOpt.NObs;EstimOpt.NP;EstimOpt.params];

%File Output
Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;
Results.INPUT = INPUT;
Results.Dist = transpose(EstimOpt.Dist);
EstimOpt.JSNVariables = find(EstimOpt.Dist > 4 & EstimOpt.Dist <= 7);


%% Tworzebnie templatek do printu


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

if EstimOpt.NVarNLT > 0
    Template1 = [Template1,'DetailsNLT0'];
    Temp = cell(1,size(Template2,2));
    Temp(1,1) = {'DetailsNLT0'};
    Template2 = [Template2;Temp];
    if EstimOpt.NLTType == 1
        Heads.DetailsNLT0 = {'Box-Cox transformation parameters';'tb'};
    elseif EstimOpt.NLTType == 2
        Heads.DetailsNLT0 = {'Yeo-Johnson transformation parameters';'tb'};
    end
end

if EstimOpt.Johnson > 0
    Heads.ResultsJ = {'Johnson location parameters';'Johnson scale parameters';'tb'}; %heads need to be written vertically
    Template1 = [Template1,'ResultsJ'];
    Temp = cell(1,size(Template2,2));
    Temp(1,1) = {'ResultsJ'};
    Template2 = [Template2;Temp];
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
    ST = {'DetailsS'};
end


%% Tworzenie naglowka


Head = cell(1,2);
if EstimOpt.FullCov == 0
    Head(1,1) = {'MXL_d'};
else
    Head(1,1) = {'MXL'};
end

if EstimOpt.WTP_space > 0
    Head(1,2) = {'in WTP-space'};
else
    Head(1,2) = {'in preference-space'};
end


%% Tworzenie stopki


Tail = cell(17,2);
Tail(2,1) = {'Model diagnostics'};
Tail(3:17,1) = {'LL at convergence';'LL at constant(s) only';strcat('McFadden''s pseudo-R',char(178));strcat('Ben-Akiva-Lerman''s pseudo-R',char(178));'AIC/n';'BIC/n';'n (observations)';'r (respondents)';'k (parameters)';'';'Estimation method';'Simulation with';'Optimization method';'Gradient';'Hessian'};

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


%% Tworzenie ResultsOut, drukowanie na ekran i do pliku .xls


if EstimOpt.Display~=0
    Results.R_out = genOutput(EstimOpt,Results,Head,Tail,Names,Template1,Template2,Heads,ST);
end


end