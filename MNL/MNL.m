function Results = MNL(INPUT,Results_old,EstimOpt,OptimOpt)


% save tmp_MNL
% return

global B_backup

tic

Results.bhat = [];
Results.R = [];
Results.R_out = {};
Results.stats = [];


%% Check data and inputs

if nargin < 3
    error('Too few input arguments for MNL(INPUT,EstimOpt,OptimOpt)')
end

warning off MATLAB:mir_warning_maybe_uninitialized_temporary

format shortG;
format compact;

if isfield(EstimOpt, 'Display') == 0
    EstimOpt.Display = 0;
end

if EstimOpt.Display ~= 0
    disp(' ');
    disp('__________________________________________________________________________________________________________________');
    disp(' ');
    if any(INPUT.W ~= 1)
        cprintf('Black','Estimating '); cprintf('*Black','weighted'); cprintf('Black',' MNL model...\n');
    else
        disp('Estimating MNL model ...')
    end
end

if isfield(EstimOpt, 'WTP_space') == 0
    EstimOpt.WTP_space = 0;
    EstimOpt.WTP_matrix = [];
elseif EstimOpt.WTP_space == 0
    EstimOpt.WTP_matrix = [];
end

if EstimOpt.Display ~= 0
    if EstimOpt.WTP_space > 0
        disp('in WTP-space ...')
    else
        disp('in preference-space ...')
    end
    if isfield(EstimOpt,'NLTVariables') && ~isempty(EstimOpt.NLTVariables)
        disp('with non-linear transformation(s) ... ')
    end
end

if isfield(EstimOpt,'NLTVariables')
    EstimOpt.NLTVariables = EstimOpt.NLTVariables(:);
    EstimOpt.NVarNLT = length(unique(EstimOpt.NLTVariables));
    if ~ismember(unique(EstimOpt.NLTVariables),1:EstimOpt.NVarA)
        error('Incorrect non-linear variable(s) specification')
    end
    if isfield(EstimOpt, 'NLTType') == 0
        cprintf(rgb('DarkOrange'), 'WARNING: Assuming Box-Cox transformation \n')
        EstimOpt.NLTType = 1;
    elseif EstimOpt.NLTType == 1
        disp('using Box-Cox transformation(s)')
    elseif EstimOpt.NLTType == 2
        disp('using Yeo-Johnson transformation(s)')
    else
        error('Incorrect transformation type')
    end
    if EstimOpt.NLTType == 1
        if any(INPUT.Xa(:, EstimOpt.NLTVariables) < 0)
            cprintf(rgb('DarkOrange'), 'WARNING: Values of Box-Cox transformed variables < 0 \n')
        elseif any(INPUT.Xa(:, EstimOpt.NLTVariables) == 0) % not sure if this is stil necessary
            cprintf(rgb('DarkOrange'), 'WARNING: Values of Box-Cox transformed variables including zeros shifted by 0.00001 \n')
            for i = 1:EstimOpt.NVarNLT
                if any(INPUT.Xa(:, EstimOpt.NLTVariables(i)) == 0)
                    INPUT.Xa(:, EstimOpt.NLTVariables(i)) = INPUT.Xa(:, EstimOpt.NLTVariables(i)) + 0.00001;
                end
            end
        end
    end
else
    EstimOpt.NVarNLT = 0;
end

if isfield(INPUT, 'Xs') == 0
    INPUT.Xs = zeros(size(INPUT.Y,1),0);
end
if isfield(EstimOpt,'SCEXP')==0
    EstimOpt.SCEXP = 1;
end
EstimOpt.NVarS = size(INPUT.Xs,2); % Number of covariates of scale

if isfield(INPUT, 'Xm') == 0 || size(INPUT.Xm,1) ~= size(INPUT.Xa,1)
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
            %         tic; EstimOpt.WTP_matrix = 1:EstimOpt.WTP_space;...
            %         EstimOpt.WTP_matrix = EstimOpt.WTP_matrix(floor((0:size(EstimOpt.WTP_matrix,2)*WTP_att-1)/WTP_att)+1); toc
        end
        %     elseif ~isequal(size(EstimOpt.WTP_matrix),[EstimOpt.NVarA-EstimOpt.WTP_space,EstimOpt.WTP_space])
    elseif size(EstimOpt.WTP_matrix,2) ~= EstimOpt.NVarA - EstimOpt.WTP_space
        error('Dimensions of EstimOpt.WTP_matrix not correct - for each non-monetary attribute provide no. of attribute to multiply it with')
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
if EstimOpt.NVarS > 0
    if isfield(EstimOpt,'NamesS') == 0 || isempty(EstimOpt.NamesS) || length(EstimOpt.NamesS) ~= EstimOpt.NVarS
        EstimOpt.NamesS = (1:EstimOpt.NVarS)';
        EstimOpt.NamesS = cellstr(num2str(EstimOpt.NamesS));
    elseif size(EstimOpt.NamesS,1) ~= EstimOpt.NVarS
        EstimOpt.NamesS = EstimOpt.NamesS';
    end
end

if ~isfield(EstimOpt,'ExpB')
    EstimOpt.ExpB = [];
elseif ~isempty(EstimOpt.ExpB)
    EstimOpt.ExpB = EstimOpt.ExpB(:);
    if size(EstimOpt.ExpB,1) == 1
        EstimOpt.ExpB = EstimOpt.ExpB.*ones(size(EstimOpt.NVarA*(1+EstimOpt.NVarM) + EstimOpt.NVarS + EstimOpt.NVarNLT,1));
    elseif size(EstimOpt.ExpB,1) ~= EstimOpt.NVarA*(1+EstimOpt.NVarM) + EstimOpt.NVarS + EstimOpt.NVarNLT
        error('Dimensions of ExpB not correct - provide ExpB indicator for each parameter in B')
    elseif any((EstimOpt.ExpB ~= 0) & (EstimOpt.ExpB ~= 1))
        error('ExpB must only include logical (0 or 1) values')
    end
    EstimOpt.ExpB = (1:EstimOpt.NVarA*(1+EstimOpt.NVarM) + EstimOpt.NVarS + EstimOpt.NVarNLT)' .* EstimOpt.ExpB;
    EstimOpt.ExpB(EstimOpt.ExpB == 0) = [];
end


%% Starting values


if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVarA*(1+EstimOpt.NVarM) + EstimOpt.NVarS + EstimOpt.NVarNLT
    b0 = B_backup(:);
    if EstimOpt.Display ~= 0
        disp('Using the starting values from Backup')
    end
elseif isfield(Results_old,'MNL') && isfield(Results_old.MNL,'b0') && (length(Results_old.MNL.b0) == EstimOpt.NVarA*(1+EstimOpt.NVarM) + EstimOpt.NVarS + EstimOpt.NVarNLT) % MNL starting values provided
    %     disp('Using MNL results as starting values')
    Results_old.MNL.b0_old = Results_old.MNL.b0;
    Results_old.MNL = rmfield(Results_old.MNL,'b0');
    if length(Results_old.MNL.b0_old) ~= EstimOpt.NVarA*(1+EstimOpt.NVarM) + EstimOpt.NVarS + EstimOpt.NVarNLT
        if EstimOpt.Display ~= 0
            cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
        end
        Results_old.MNL = rmfield(Results_old.MNL,'b0_old');
    else
        b0 = Results_old.MNL.b0_old(:);
    end
end
if ~exist('b0','var')
    if EstimOpt.Display ~= 0
        disp('Using linear regression estimates as starting values')
    end
    if EstimOpt.NVarS > 0
        b00 = zeros(EstimOpt.NVarS,1); ...
            Y = INPUT.Y(INPUT.MissingInd == 0);
        Xa = INPUT.Xa(INPUT.MissingInd == 0,:);
        if EstimOpt.NVarM > 0
            Xm = reshape(INPUT.Xm(INPUT.MissingInd == 0,:), size(Xa,1),1, EstimOpt.NVarM);
            Xm = reshape(Xm(:, ones(1, EstimOpt.NVarA),:), size(Xa,1),EstimOpt.NVarA*EstimOpt.NVarM);
            Xa2 = reshape(Xa(:,:, ones(1, EstimOpt.NVarM)), size(Xa,1),EstimOpt.NVarA*EstimOpt.NVarM);
            b0 = [regress(Y,[Xa, Xa2.*Xm]);b00;ones(EstimOpt.NVarNLT,1)];
        else
            b0 = [regress(Y,Xa);b00;ones(EstimOpt.NVarNLT,1)];
        end
        %         if EstimOpt.WTP_space > 0
        %             b0(1:EstimOpt.NVarA-EstimOpt.WTP_space) = b0(1:EstimOpt.NVarA-EstimOpt.WTP_space) .* b0(EstimOpt.WTP_matrix,:);
        %         end
    else
        if EstimOpt.NVarM > 0
            Xm = reshape(INPUT.Xm, size(INPUT.Xa,1),1, EstimOpt.NVarM);
            Xm = reshape(Xm(:, ones(1, EstimOpt.NVarA),:), size(INPUT.Xa,1),EstimOpt.NVarA*EstimOpt.NVarM);
            Xa2 = reshape(INPUT.Xa(:,:, ones(1, EstimOpt.NVarM)), size(INPUT.Xa,1),EstimOpt.NVarA*EstimOpt.NVarM);
            b0 = [regress(INPUT.Y,[INPUT.Xa, Xa2.*Xm]);ones(EstimOpt.NVarNLT,1)];
        else
            b0 = [regress(INPUT.Y,INPUT.Xa); ones(EstimOpt.NVarNLT,1)];
        end
        %         if EstimOpt.WTP_space > 0
        %             b0(1:EstimOpt.NVarA-EstimOpt.WTP_space) = b0(1:EstimOpt.NVarA-EstimOpt.WTP_space) .* b0(EstimOpt.WTP_matrix,:);
        %         end
    end
    if ~isempty(EstimOpt.ExpB)
        b0(EstimOpt.ExpB) = max(log(b0(EstimOpt.ExpB)),1);
    end
end


%% Optimization Options


% if any(EstimOpt.MissingAlt(:) == 1) && EstimOpt.NumGrad == 0
% % 	EstimOpt.NumGrad = 1;
%     OptimOpt.GradObj = 'off';
%     if EstimOpt.Display ~= 0
% %         cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - missing alternatives not supported by analytical gradient \n')
%         cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient off - missing alternatives not supported by analytical gradient \n')
%     end
% end

if isfield(EstimOpt,'BActive')
    EstimOpt.BActive = EstimOpt.BActive(:)';
end

if EstimOpt.ConstVarActive == 1
    if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive) || sum(EstimOpt.BActive == 0) == 0
        error ('Are there any constraints on model parameters (EstimOpt.ConstVarActive)? Constraints not provided (EstimOpt.BActive).')
    elseif length(b0) ~= length(EstimOpt.BActive)
        error('Check no. of constraints')
    end
    if EstimOpt.Display ~= 0
        disp(['Initial values: ' mat2str(b0',2)])
        disp(['Parameters with zeros are constrained to their initial values: ' mat2str(EstimOpt.BActive')])
    end
else
    if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive) || sum(EstimOpt.BActive == 0) == 0
        EstimOpt.BActive = ones(1,length(b0));
        if EstimOpt.Display ~= 0
            disp(['Initial values: ' mat2str(b0',2)])
        end
    else
        if length(b0) ~= length(EstimOpt.BActive)
            error('Check no. of constraints')
        else
            if EstimOpt.Display ~= 0
                disp(['Initial values: ' mat2str(b0',2)])
                disp(['Parameters with zeros are constrained to their initial values: ' mat2str(EstimOpt.BActive')])
            end
        end
    end
end

if ((isfield(EstimOpt, 'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    if EstimOpt.Display ~= 0
        cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    end
    OptimOpt.GradObj = 'on';
end

if EstimOpt.NVarS > 0 && EstimOpt.NumGrad == 0
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
    if EstimOpt.Display ~= 0
        %         cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - covariates of scale not supported by analytical gradient \n')
        cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient off - covariates of scale not supported by analytical gradient \n')
    end
end

% if EstimOpt.NVarNLT > 0 && EstimOpt.NLTType == 2 && EstimOpt.NumGrad == 0
% 	EstimOpt.NumGrad = 1;
% 	if EstimOpt.Display ~= 0
%         cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - Yeo-Johnston transformation not supported by analytical gradient \n')
% 	end
% end

if ~isempty(EstimOpt.ExpB) && EstimOpt.NumGrad == 0
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
    if EstimOpt.Display ~= 0
        cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient off - ExpB not supported by analytical gradient \n')
    end
end

if (isfield(EstimOpt, 'ConstVarActive') == 0 || EstimOpt.ConstVarActive == 0) && isequal(OptimOpt.Algorithm,'quasi-newton') && isequal(OptimOpt.Hessian,'user-supplied')
    if EstimOpt.Display ~= 0
        cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied Hessian off - quasi-newton algorithm does not use it anyway \n')
    end
    OptimOpt.Hessian = 'off';
end

if EstimOpt.Display ~= 0
    fprintf('\n')
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
    fprintf('\n')
end


%% Restucturing Data - gets rid of not completed choice tasks, but leaves missing alternatives


idx = sum(reshape(INPUT.MissingInd,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP)) == EstimOpt.NAlt; ...
    idx = reshape(idx(ones(EstimOpt.NAlt,1),:), EstimOpt.NAlt*EstimOpt.NCT*EstimOpt.NP,1);
INPUT.Y = INPUT.Y(idx == 0);

INPUT.Xa(INPUT.MissingInd == 1,:) = NaN;
INPUT.Xs(INPUT.MissingInd == 1,:) = NaN;
%Xa and Xs has NaNs where there is missing alternative
INPUT.Xa = INPUT.Xa(idx == 0,:);
if EstimOpt.NVarM > 0 && (EstimOpt.WTP_space > 0 || EstimOpt.NVarNLT > 0)
    INPUT.Xm = INPUT.Xm(1:EstimOpt.NCT*EstimOpt.NAlt:end,:);
    EstimOpt.XmIndx = zeros(sum(EstimOpt.NCTMiss),1);
    EstimOpt.XmIndx2 = zeros(sum(EstimOpt.NCTMiss)*EstimOpt.NAlt,1);
    for i = 1:EstimOpt.NP
        NCTno = sum(EstimOpt.NCTMiss(1:i-1));
        EstimOpt.XmIndx(NCTno+1:NCTno+EstimOpt.NCTMiss(i)) = i;
        EstimOpt.XmIndx2(NCTno*EstimOpt.NAlt+1:(NCTno+EstimOpt.NCTMiss(i))*EstimOpt.NAlt) = i;
    end
    NVarMOld = EstimOpt.NVarM;
elseif EstimOpt.NVarM > 0 && EstimOpt.WTP_space == 0 && EstimOpt.NVarNLT == 0
    INPUT.Xm = INPUT.Xm(idx == 0,:);
    Xm = reshape(INPUT.Xm, size(INPUT.Xa,1),1, EstimOpt.NVarM);
    Xm = reshape(Xm(:, ones(1, EstimOpt.NVarA),:), size(INPUT.Xa,1),EstimOpt.NVarA*EstimOpt.NVarM);
    Xa2 = reshape(INPUT.Xa(:,:, ones(1, EstimOpt.NVarM)), size(INPUT.Xa,1),EstimOpt.NVarA*EstimOpt.NVarM);
    INPUT.Xa = [INPUT.Xa, Xa2.*Xm];
    EstimOpt.NVarA = EstimOpt.NVarA*(1+EstimOpt.NVarM);
    NVarMOld = EstimOpt.NVarM;
    EstimOpt.NVarM = 0;
    INPUT.Xm = zeros(size(INPUT.Y,1),0);
else
    NVarMOld = EstimOpt.NVarM;
end
if EstimOpt.NVarS > 0
    INPUT.Xs = INPUT.Xs(idx == 0,:);
end
if any(INPUT.W ~= 1)
    Wtmp = zeros(sum(EstimOpt.NCTMiss),1);
    Wtmp(1:EstimOpt.NCTMiss(1)) = INPUT.W(1);
    for i = 2:EstimOpt.NP
        Wtmp(sum(EstimOpt.NCTMiss(1:(i-1)))+1:sum(EstimOpt.NCTMiss(1:i))) = INPUT.W(i);
    end
    INPUT.W = Wtmp;
    INPUT.W = sum(EstimOpt.NCTMiss)*INPUT.W/sum(INPUT.W);
else
    INPUT.W = ones(sum(EstimOpt.NCTMiss),1);
end


%% Estimation


LLfun = @(B) LL_mnl_MATlike(INPUT.Y, INPUT.Xa,INPUT.Xm, INPUT.Xs,INPUT.W, EstimOpt,OptimOpt,B);

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

% save tmp_MNL_output

Results.LL = -LL;
Results.b0_old = b0;

if isfield(EstimOpt,'R2type') == 0
    EstimOpt.R2type = 0;
end

Results.LLdetailed = LL_mnl(INPUT.Y,INPUT.Xa,INPUT.Xm, INPUT.Xs,EstimOpt,Results.bhat);
Results.LLdetailed = Results.LLdetailed.*INPUT.W;
if any(INPUT.MissingInd == 1) % In case of some missing data
    idx = sum(reshape(INPUT.MissingInd,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP)) == EstimOpt.NAlt; ...
        idx = sum(reshape(idx, EstimOpt.NCT, EstimOpt.NP),1)'; % no. of missing NCT for every respondent
    R2 = zeros(EstimOpt.NP,1);
    idx = EstimOpt.NCT - idx;
    l = 1;
    for i = 1:EstimOpt.NP
        R2(i) = prod(exp(Results.LLdetailed(l:l-1+idx(i)))).^(1/idx(i));
        l = l+idx(i);
    end
    R2 = mean(R2);
else
    R2 = mean(prod(reshape(exp(Results.LLdetailed), EstimOpt.NCT, EstimOpt.NP),1).^(1/EstimOpt.NCT),2);
end

if EstimOpt.HessEstFix == 1
    f = LL_mnl(INPUT.Y,INPUT.Xa,INPUT.Xm,INPUT.Xs,EstimOpt,Results.bhat);
    Results.jacobian = numdiff(@(B) -INPUT.W.*LL_mnl(INPUT.Y,INPUT.Xa,INPUT.Xm,INPUT.Xs,EstimOpt,B),-INPUT.W.*f,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) -INPUT.W.*LL_mnl(INPUT.Y,INPUT.Xa,INPUT.Xm,INPUT.Xs,EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.hess = hessian(@(B) -sum(INPUT.W.*LL_mnl(INPUT.Y,INPUT.Xa,INPUT.Xm,INPUT.Xs,EstimOpt,B),1), Results.bhat);
elseif EstimOpt.HessEstFix == 4 % analytical - missing
    Results.hess = hessian(@(B) -sum(INPUT.W.*LL_mnl(INPUT.Y,INPUT.Xa,INPUT.Xm,INPUT.Xs,EstimOpt,B),1), Results.bhat);
    %     EstimOpt_tmp = EstimOpt;
    %     EstimOpt_tmp.NumGrad = 0;
    %     EstimOpt_tmp.ApproxHess = 0;
    %     [~,~,Results.hess] = LL_mnl(INPUT.Y,INPUT.Xa,INPUT.Xs,EstimOpt_tmp,Results.bhat);
end

if sum(EstimOpt.BActive == 0) > 0
    if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
        Results.jacobian = Results.jacobian(:, EstimOpt.BActive == 1);
        Results.hess = Results.jacobian'*Results.jacobian;
    elseif EstimOpt.HessEstFix == 0 || EstimOpt.HessEstFix == 3
        Results.hess = Results.hess(EstimOpt.BActive == 1,EstimOpt.BActive == 1);
    end
    Results.ihess = inv(Results.hess);
    Results.ihess = direcXpnd(Results.ihess,EstimOpt.BActive);
    Results.ihess = direcXpnd(Results.ihess',EstimOpt.BActive);
    
else
    if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
        Results.hess = Results.jacobian'*Results.jacobian;
    end
    Results.ihess = inv(Results.hess);
end

if EstimOpt.RobustStd == 1
    if EstimOpt.NumGrad == 0
        [~, Results.jacobian] = LL_mnl(INPUT.Y,INPUT.Xa,INPUT.Xs,EstimOpt,Results.bhat);
        Results.jacobian = -Results.jacobian.*INPUT.W(:, ones(1,size(Results.jacobian,2)));
    else
        Results.jacobian = numdiff(@(B) -INPUT.W.*LL_mnl(INPUT.Y,INPUT.Xa,INPUT.Xs,EstimOpt,B),Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
    end
    RobJacob = zeros(EstimOpt.NP, size(Results.jacobian,2));
    RobJacob(1,:) = sum(Results.jacobian(1:EstimOpt.NCTMiss(1),:),1);
    for i = 2:EstimOpt.NP
        RobJacob(i,:) = sum(Results.jacobian(sum(EstimOpt.NCTMiss(1:(i-1)))+1:sum(EstimOpt.NCTMiss(1:i)),:),1);
    end
    RobustHess = RobJacob'*RobJacob;
    Results.ihess = Results.ihess*RobustHess*Results.ihess;
end
Results.std = sqrt(diag(Results.ihess));

if sum(EstimOpt.BActive == 0) > 0
    Results.std(EstimOpt.BActive == 0) = NaN;
end

Results.std(imag(Results.std) ~= 0) = NaN;
Results.R = [Results.bhat , Results.std , pv(Results.bhat , Results.std)];

EstimOpt.Params = length(b0);
if isfield(EstimOpt,'BActive')
    EstimOpt.Params = EstimOpt.Params - sum(EstimOpt.BActive == 0);
end

if isfield(Results_old,'MNL0') && isfield(Results_old.MNL0,'LL')
    Results.stats = [Results.LL; Results_old.MNL0.LL;  1-Results.LL/Results_old.MNL0.LL;R2; ((2*EstimOpt.Params-2*Results.LL))/EstimOpt.NObs; ((log(EstimOpt.NObs)*EstimOpt.Params-2*Results.LL))/EstimOpt.NObs ;EstimOpt.NObs; EstimOpt.NP; EstimOpt.Params];
end

if EstimOpt.WTP_space == 0
    % Results.WTP = IP_MNL_delta(Results.bhat(1:EstimOpt.NVarA) , Results.ihess(1:EstimOpt.NVarA,1:EstimOpt.NVarA),EstimOpt)';
else
    Results.WTP = [Results.bhat(1:end-EstimOpt.WTP_space) Results.std(1:end-EstimOpt.WTP_space) Results.bhat(1:end-EstimOpt.WTP_space)-norminv(0.975,0,1)*Results.std(1:end-EstimOpt.WTP_space) Results.bhat(1:end-EstimOpt.WTP_space)+norminv(0.975,0,1)*Results.std(1:end-EstimOpt.WTP_space)];
end

Results.INPUT = INPUT;
Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;


% clocknote = clock;
% tocnote = toc;
% [~,DayName] = weekday(now,'long');
if NVarMOld > 0 && EstimOpt.WTP_space == 0 && EstimOpt.NVarNLT == 0
    EstimOpt.NVarA = EstimOpt.NVarA/(1+NVarMOld);
end

Results.DetailsA(1:EstimOpt.NVarA,1) = Results.bhat(1:EstimOpt.NVarA);
Results.DetailsA(1:EstimOpt.NVarA,3:4) = [Results.std(1:EstimOpt.NVarA),pv(Results.bhat(1:EstimOpt.NVarA),Results.std(1:EstimOpt.NVarA))];

if NVarMOld > 0
    Results.DetailsM = [];
    for i = 1:NVarMOld
        Results.DetailsM = [Results.DetailsM, [Results.bhat(EstimOpt.NVarA*(i)+1:EstimOpt.NVarA*(i+1)),zeros(EstimOpt.NVarA,1),Results.std(EstimOpt.NVarA*(i)+1:EstimOpt.NVarA*(i+1)),pv(Results.bhat(EstimOpt.NVarA*(i)+1:EstimOpt.NVarA*(i+1)),Results.std(EstimOpt.NVarA*(i)+1:EstimOpt.NVarA*(i+1)))]];
    end
    %     Results.R(end-EstimOpt.NVarA*NVarMOld+1:end,:) = [];
    %     Results.R = [Results.R, Results.DetailsM];
end

if EstimOpt.NVarNLT > 0
    Results.DetailsNLT(:,1) = Results.bhat(EstimOpt.NVarA+EstimOpt.NVarS+EstimOpt.NVarA*NVarMOld+1:end);
    Results.DetailsNLT(:,3:4) = [Results.std(EstimOpt.NVarA+EstimOpt.NVarS+EstimOpt.NVarA*NVarMOld+1:end),pv(Results.bhat(EstimOpt.NVarA+EstimOpt.NVarS+EstimOpt.NVarA*NVarMOld+1:end),Results.std(EstimOpt.NVarA+EstimOpt.NVarS+EstimOpt.NVarA*NVarMOld+1:end))];
    Results.DetailsNLT0 = NaN(EstimOpt.NVarA,4);
    Results.DetailsNLT0(EstimOpt.NLTVariables,:) = Results.DetailsNLT;
end
if EstimOpt.NVarS > 0
    Results.DetailsS = [];
    Results.DetailsS(:,1) = Results.bhat(EstimOpt.NVarA*(1+NVarMOld)+1:EstimOpt.NVarA*(1+NVarMOld)+EstimOpt.NVarS);
    Results.DetailsS(:,3:4) = [Results.std(EstimOpt.NVarA*(1+NVarMOld)+1:EstimOpt.NVarA*(1+NVarMOld)+EstimOpt.NVarS), pv(Results.bhat(EstimOpt.NVarA*(1+NVarMOld)+1:EstimOpt.NVarA*(1+NVarMOld)+EstimOpt.NVarS),Results.std(EstimOpt.NVarA*(1+NVarMOld)+1:EstimOpt.NVarA*(1+NVarMOld)+EstimOpt.NVarS))];
end


%% Tworzebnie templatek do printu


Template1 = {'DetailsA'};
Template2 = {'DetailsA'};
Names.DetailsA = EstimOpt.NamesA;
Heads.DetailsA = {'Means'};
ST = {'DetailsA'};

if NVarMOld > 0
    Template1 = [Template1, 'DetailsM'];
    Template2 = [Template2; 'DetailsM'];
    Heads.DetailsM = EstimOpt.NamesM;
    ST = [ST, 'DetailsM'];
end

if EstimOpt.NVarNLT > 0
    Template1 = [Template1, 'DetailsNLT0'];
    Template2 = [Template2; 'DetailsNLT0'];
    ST = [ST, 'DetailsNLT0'];
    if EstimOpt.NLTType == 1
        Heads.DetailsNLT0 = {'Box-Cox transformation parameters'};
    elseif EstimOpt.NLTType == 2
        Heads.DetailsNLT0 = {'Yeo-Johnson transformation parameters'};
    end
end

if EstimOpt.NVarS > 0
    Temp = cell(1, size(Template1,2));
    Temp(1,1) = {'DetailsS'};
    Template1 = [Template1; Temp];
    Template2 = [Template2; 'DetailsS'];
    Names.DetailsS = EstimOpt.NamesS;
    Heads.DetailsS = {'Covariates of Scale'};
    ST = [ST, 'DetailsS'];
end

%% Tworzenie naglowka
Head = cell(1,2);
Head(1,1) = {'MNL'};
if EstimOpt.WTP_space > 0
    Head(1,2) = {'in WTP-space'};
else
    Head(1,2) = {'in preference-space'};
end
%% Tworzenie stopki
Tail = cell(16,2);
Tail(2,1) = {'Model diagnostics'};
Tail(3:16,1) = { 'LL at convergence' ; 'LL at constant(s) only'; strcat('McFadden''s pseudo-R',char(178));strcat('Ben-Akiva-Lerman''s pseudo-R',char(178))  ;'AIC/n' ;'BIC/n'; 'n (observations)'; 'r (respondents)';'k (parameters)';' ';'Estimation method';'Optimization method';'Gradient';'Hessian'};

if isfield(Results_old,'MNL0') && isfield(Results_old.MNL0,'LL')
    Tail(3:11,2) = num2cell(Results.stats);
end

if any(INPUT.W ~= 1)
    Tail(13,2) = {'weighted'};
else
    Tail(13,2) = {'maximum likelihood'};
end

Tail(14,2) = {OptimOpt.Algorithm;};

if strcmp(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        Tail(15,2) = {'user-supplied, analytical'};
    else
        Tail(15,2) = {['user-supplied, numerical ',num2str(OptimOpt.FinDiffType)]};
    end
else
    Tail(15,2) = {['built-in, ',num2str(OptimOpt.FinDiffType)]};
    
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

Tail(16,2) = {outHessian};
%% Tworzenie ResultsOut, drukowanie na ekran i do pliku .xls
EstimOpt.Dist = -ones(1,EstimOpt.NVarA+1);
if EstimOpt.Display~=0
    Results.Dist = transpose(EstimOpt.Dist(:,2:end));
    Results.R_out = genOutput(EstimOpt, Results, Head, Tail, Names, Template1, Template2, Heads, ST);
    fullOrgTemplate = which('template.xls');
    currFld = pwd;
    if isfield(EstimOpt,'ProjectName')
        fullSaveName = strcat(currFld,'\MNL_results_',EstimOpt.ProjectName,'.xls');
    else
        fullSaveName = strcat(currFld,'\MNL_results.xls');
    end
    
    copyfile(fullOrgTemplate, 'templateTMP.xls')
    fullTMPTemplate = which('templateTMP.xls');
    excel = actxserver('Excel.Application');
    %     try WbookCheck = excel.Workbooks('fullTMPTemplate');
    %         Close(WbookCheck)
    %         catch
    %              excelWorkbook = excel.Workbooks.Open(fullTMPTemplate);
    %     end
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
    i = 1;
    if isfield(EstimOpt,'xlsOverwrite') && EstimOpt.xlsOverwrite == 0
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
end