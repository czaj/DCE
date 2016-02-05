function Results = LC(INPUT,Results_old,EstimOpt,OptimOpt)


% save tmp_LC
% return

tic

global B_backup

Results.bhat = [];
Results.R = [];
Results.R_out = {};
Results.stats = [];


%% Check data and inputs


if nargin < 3 % check no. of inputs
    error('Too few input arguments for LC(INPUT,EstimOpt,OptimOpt)')
end

warning off MATLAB:mir_warning_maybe_uninitialized_temporary

format shortG;
format compact;

if isfield(EstimOpt, 'Display') == 0
    EstimOpt.Display = 1;
end

if isfield(EstimOpt,'NClass') == 0; 
    EstimOpt.NClass = 2; 
end

if EstimOpt.Display == 1
    disp(' ');
    disp('__________________________________________________________________________________________________________________');
    disp(' ');
end
if EstimOpt.Display == 1
    disp(num2str(EstimOpt.NClass,'Estimating LC Model with %1.0f classes ...'))
end

if isfield(EstimOpt, 'WTP_space') == 0
    EstimOpt.WTP_space = 0;
    EstimOpt.WTP_matrix = [];
elseif EstimOpt.WTP_space == 0;
	EstimOpt.WTP_matrix = [];
end

if EstimOpt.Display == 1
    if EstimOpt.WTP_space > 0
        disp('in WTP-space.')
    else
        disp('in preference-space.')
    end
end

if isfield(INPUT, 'Xc') == 0 || numel(INPUT.Xc) == 0
    INPUT.Xc = ones(size(INPUT.Y,1),1);
else
    INPUT.Xc = [ones(size(INPUT.Y,1),1), INPUT.Xc];
end
EstimOpt.NVarC = size(INPUT.Xc,2); % no. of variables explaining class probabilities

if isfield(EstimOpt, 'ClassScores') == 0
   EstimOpt.ClassScores = 0; 
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
%     elseif ~isequal(size(EstimOpt.WTP_matrix),[EstimOpt.NVarA-EstimOpt.WTP_space,EstimOpt.WTP_space])
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

if ~isfield(EstimOpt,'BActiveClass') || isempty(EstimOpt.BActiveClass) || sum(EstimOpt.BActiveClass == 0) == 0
    EstimOpt.BActiveClass = ones(EstimOpt.NVarA,1);
end

% for later:
% - covariates of scale (class-specific?)
% - allow for non-constant scale parameters for classes


%% Starting values


EstimOpt.jitter1 = 1.1; % Jittering parameter (relative) for MNL starting values (attributes)
EstimOpt.jitter2 = 0.1; % Jittering parameter (absolute) for class probabilities starting values

if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVarA*EstimOpt.NClass + (EstimOpt.NClass-1)*EstimOpt.NVarC
    b0 = B_backup(:);
    if EstimOpt.Display == 1
        disp('Using the starting values from Backup')
    end
elseif isfield(Results_old,'LC') && isfield(Results_old.LC,'b0') % starting values provided
    Results_old.LC.b0_old = Results_old.LC.b0(:);
    Results_old.LC = rmfield(Results_old.LC,'b0');
    if length(Results_old.LC.b0_old) ~= EstimOpt.NVarA*EstimOpt.NClass + (EstimOpt.NClass-1)*EstimOpt.NVarC
        if EstimOpt.Display == 1
            cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
        end
        Results_old.LC = rmfield(Results_old.LC,'b0_old');
    else
        b0 = Results_old.LC.b0_old(:);
    end
end
if  ~exist('b0','var')
    if isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat')
        if EstimOpt.Display == 1
            disp('Using MNL results as starting values')
        end
        b0 = [Results_old.MNL.bhat((1:size(Results_old.MNL.bhat,1))'*ones(1,EstimOpt.NClass),:) .* ...
            (EstimOpt.jitter1 .* unifrnd(0,ones(EstimOpt.NVarA.*EstimOpt.NClass,1))) ; ...
            EstimOpt.jitter2 + unifrnd(-1,ones(EstimOpt.NVarC.*(EstimOpt.NClass-1),1))];    
    else
            error('No starting values available - run MNL first')
    end
end


%% Optimization Options

if  isfield(EstimOpt,'BActive')
	EstimOpt.BActive = EstimOpt.BActive(:)';
end

if isfield(EstimOpt,'BActiveClass') && ~isempty(EstimOpt.BActiveClass) && sum(EstimOpt.BActiveClass == 0) > 0 && ...
    isfield(EstimOpt,'BActive') && ~isempty(EstimOpt.BActive) && sum(EstimOpt.BActive == 0) > 0
    disp('WARINING: The model does not currently support simultaneous BActive and BActiveClass constraints')
end

if isfield(EstimOpt,'BActiveClass')
	EstimOpt.BActiveClass = EstimOpt.BActiveClass(:);
end
if ~isfield(EstimOpt,'BActiveClass') || isempty(EstimOpt.BActiveClass) || sum(EstimOpt.BActiveClass == 0) == 0
    EstimOpt.BActiveClass = ones(EstimOpt.NVarA,1);
else
    if length(EstimOpt.BActiveClass) ~= EstimOpt.NVarA
        error('Check no. of constraints (BActiveClass)')
    else
        disp(['Parameters of the attributes with zeros constrained equal between classes: ' mat2str(EstimOpt.BActiveClass)])
    end
end

if EstimOpt.ConstVarActive == 1
    if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive) || sum(EstimOpt.BActive == 0) == 0
        error ('Are there any constraints on model parameters (EstimOpt.ConstVarActive)? Constraints not provided (EstimOpt.BActive).')
    elseif length(b0) ~= length(EstimOpt.BActive)
        error('Check no. of constraints')
    end
    if EstimOpt.Display == 1
        disp(['Initial values: ' mat2str(b0',2)])
        disp(['Parameters with zeros are constrained to their initial values: ' mat2str(EstimOpt.BActive')])
    end
else    
    if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive) || sum(EstimOpt.BActive == 0) == 0
        EstimOpt.BActive = ones(1,length(b0));
        if EstimOpt.Display == 1
            disp(['Initial values: ' mat2str(b0',2)])
        end
    else        
        if length(b0) ~= length(EstimOpt.BActive)
            error('Check no. of constraints')
        else
            if EstimOpt.Display == 1
                disp(['Initial values: ' mat2str(b0',2)])
                disp(['Parameters with zeros are constrained to their initial values: ' mat2str(EstimOpt.BActive')])
            end
        end
    end
end

if any(EstimOpt.MissingCT(:) == 1) && EstimOpt.NumGrad == 0
   EstimOpt.NumGrad = 1;
   OptimOpt.GradObj = 'off';
   cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - missing choice tasks not supported by analytical gradient \n')
end

if any(EstimOpt.MissingAlt(:) == 1) && EstimOpt.NumGrad == 0
   EstimOpt.NumGrad = 1;
   OptimOpt.GradObj = 'off';
   cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - missing alternatives not supported by analytical gradient \n')
end

if ((isfield(EstimOpt, 'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    OptimOpt.GradObj = 'on';
end

if any(EstimOpt.BActiveClass == 0) && EstimOpt.NumGrad == 0
   EstimOpt.NumGrad = 1;
   OptimOpt.GradObj = 'off';
   cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - parameters constrained between classes not supported by analytical gradient \n')
end

% if EstimOpt.WTP_space > 1 && EstimOpt.NumGrad == 0
%    EstimOpt.NumGrad = 1;
%    cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - WTP_space > 1 not supported by analytical gradient \n')
% end

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


%% RESTRUNCTURING DATA

INPUT.XXc = INPUT.Xc(1:EstimOpt.NCT*EstimOpt.NAlt:end,:); % NP x NVarC

INPUT.YY = reshape(INPUT.Y,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP);...
INPUT.YY = INPUT.YY(:,[1:size(INPUT.YY,2)]'*ones(1,EstimOpt.NClass)); %NAlt x NCT*NP*NClass

INPUT.MissingInd = INPUT.MissingInd([1:size(INPUT.MissingInd,1)]' * ones(1,EstimOpt.NClass),:);...
INPUT.MissingInd = reshape(INPUT.MissingInd,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP*EstimOpt.NClass); %NAlt x NCT*NP*NClass

if sum(EstimOpt.BActiveClass == 0,1) > 0
   bactive_1 = EstimOpt.BActive(1:EstimOpt.NClass*EstimOpt.NVarA);
   b0_1 = b0(1:EstimOpt.NClass*EstimOpt.NVarA);
   bactive_2 =  EstimOpt.BActive(EstimOpt.NClass*EstimOpt.NVarA+1:end);
   b0_2 = b0(EstimOpt.NClass*EstimOpt.NVarA+1:end);
   b0 = b0_1(1:EstimOpt.NVarA);
   EstimOpt.BActive = bactive_1(1:EstimOpt.NVarA);

   for i = 2:EstimOpt.NClass
       b0x = b0_1((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA);
       bactivex = bactive_1((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA);
       b0 = [b0; b0x(EstimOpt.BActiveClass ==1)];
       EstimOpt.BActive = [EstimOpt.BActive, bactivex(EstimOpt.BActiveClass ==1)];
   end
   
   b0 = [b0; b0_2];
   EstimOpt.BActive = [EstimOpt.BActive, bactive_2];
   clear bactive_1 b0_1 bactive_2 b0_2 b0x bactivex
end

%% ESTIMATION

LLfun = @(B) LL_lc_MATlike(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.MissingInd,EstimOpt,OptimOpt,B);
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

% save tmp_LC

%% OUTPUT

Results.LL = -LL;
Results.b0_old = b0;

Results.LLdetailed = LL_lc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.MissingInd,EstimOpt,Results.bhat);
if any(INPUT.MissingInd == 1) % In case of some missing data
   idx = sum(reshape(INPUT.MissingInd,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP)) == EstimOpt.NAlt; ...
   idx = sum(reshape(idx, EstimOpt.NCT, EstimOpt.NP),1)'; % no. of missing NCT for every respondent
   idx = EstimOpt.NCT - idx;
   R2 = mean(exp(-Results.LLdetailed./idx),1);
else
    R2 = mean(exp(-Results.LLdetailed/EstimOpt.NCT),1);
end


if EstimOpt.ClassScores ~= 0 
    Results.ClassScores = BayesProbs(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.MissingInd,EstimOpt,Results.bhat);
end

if EstimOpt.HessEstFix == 1
	f = LL_lc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.MissingInd,EstimOpt,Results.bhat);
    Results.jacobian = numdiff(@(B) LL_lc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.MissingInd,EstimOpt,B),f,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) LL_lc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.MissingInd,EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.hess = hessian(@(B) sum(LL_lc(INPUT.YY,INPUT.Xa,INPUT.XXc,INPUT.MissingInd,EstimOpt,B)), Results.bhat);
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
Results.std = sqrt(diag(Results.ihess));
Results.std(EstimOpt.BActive == 0) = NaN;
Results.std(EstimOpt.BLimit == 1) = 0;

Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;
Results.INPUT = INPUT;

EstimOpt.params = length(Results.bhat);
EstimOpt.params = EstimOpt.params - sum(EstimOpt.BActive == 0) + sum(EstimOpt.BLimit == 1);

if sum(EstimOpt.BActiveClass == 0,1) == 0
    Results.DetailsA = [Results.bhat(1:EstimOpt.NClass*EstimOpt.NVarA), Results.std(1:EstimOpt.NClass*EstimOpt.NVarA), pv(Results.bhat(1:EstimOpt.NClass*EstimOpt.NVarA), Results.std(1:EstimOpt.NClass*EstimOpt.NVarA))];
    l = EstimOpt.NClass*EstimOpt.NVarA;
else
    Bclass = Results.bhat(1:EstimOpt.NVarA)*ones(1,EstimOpt.NClass);
    stdclass = Results.std(1:EstimOpt.NVarA)*ones(1,EstimOpt.NClass);
   for i = 1: (EstimOpt.NClass - 1)
       Bclass(EstimOpt.BActiveClass == 1,i+1) = Results.bhat(EstimOpt.NVarA + (i-1)*sum(EstimOpt.BActiveClass,1)+1:EstimOpt.NVarA + i*sum(EstimOpt.BActiveClass,1));
       stdclass(EstimOpt.BActiveClass == 1,i+1) = Results.std(EstimOpt.NVarA + (i-1)*sum(EstimOpt.BActiveClass,1)+1:EstimOpt.NVarA + i*sum(EstimOpt.BActiveClass,1));
   end
   stdclass(EstimOpt.BActiveClass == 0,2:end) = NaN;
   Results.DetailsA = [reshape(Bclass, EstimOpt.NClass*EstimOpt.NVarA,1), reshape(stdclass, EstimOpt.NClass*EstimOpt.NVarA,1), pv(reshape(Bclass, EstimOpt.NClass*EstimOpt.NVarA,1), reshape(stdclass, EstimOpt.NClass*EstimOpt.NVarA,1))];
   l = EstimOpt.NVarA + (EstimOpt.NClass - 1)*sum(EstimOpt.BActiveClass,1);
end
Results.DetailsV = [Results.bhat(l+1:end), Results.std(l+1:end), pv(Results.bhat(l+1:end), Results.std(l+1:end))];

if sum(EstimOpt.BActiveClass == 0,1) == 0
    bclass = reshape([Results.bhat(EstimOpt.NVarA*EstimOpt.NClass+1:end); zeros(EstimOpt.NVarC,1)], EstimOpt.NVarC, EstimOpt.NClass);	
else
	bclass = reshape([Results.bhat((EstimOpt.NClass-1)*sum(EstimOpt.BActiveClass,1)+EstimOpt.NVarA+1:end); zeros(EstimOpt.NVarC,1)], EstimOpt.NVarC, EstimOpt.NClass);
    Results.bhat = [Results.DetailsA(:,1);Results.DetailsV(:,1)];
end

V = exp(INPUT.XXc*bclass); ...% NP x NClass
Vsum = sum(V,2); ...
Results.PClass = mean(V./Vsum(:,ones(EstimOpt.NClass,1)),1); %1 x NClass
clear bclass l stdclass i f V Vsum

Results.R = [Results.DetailsA; Results.DetailsV; [Results.PClass',zeros(size(Results.PClass')),zeros(size(Results.PClass'))]];

Results.stats = [Results_old.MNL0.LL; Results.LL; 1-Results.LL/Results_old.MNL0.LL;R2; ((2*EstimOpt.params-2*Results.LL) + 2*EstimOpt.params*(EstimOpt.params+1)/(EstimOpt.NObs-EstimOpt.params-1))/EstimOpt.NObs; EstimOpt.NObs; EstimOpt.params];

Results.R_out = cell(3+EstimOpt.NVarA + 2 + EstimOpt.NVarC + 3 + 2 + 7, 1 + 3*EstimOpt.NClass);  
Results.R_out(1,1) = {'LC'};
head = {'var.' , 'coef.', 'st.err.' , 'p-value'};
NClasses = {'Class 1','Class 2', 'Class 3', 'Class 4', 'Class 5', 'Class 6', 'Class 7', 'Class 8', 'Class 9','Class 10'};
headx = [head, repmat(head(1,2:4),1,EstimOpt.NClass-1)];
if EstimOpt.NClass<=10
    Results.R_out(2,2:3:(2+(EstimOpt.NClass-1)*3)) = NClasses(1,1:EstimOpt.NClass);
end
Results.R_out(3,:) = headx;
Results.R_out(4:(EstimOpt.NVarA+3),1) = EstimOpt.NamesA;
for i = 1: EstimOpt.NClass
    Results.R_out(4:(EstimOpt.NVarA+3),2 + 3*(i-1):1+3*i) = num2cell(Results.DetailsA((i-1)*EstimOpt.NVarA+1:i*EstimOpt.NVarA,:));
end
Results.R_out(EstimOpt.NVarA+4,1) = {'Latent class probability model'};
Results.R_out(EstimOpt.NVarA+5,1:end-3) = headx(1:end-3);
Results.R_out(EstimOpt.NVarA+6:EstimOpt.NVarA+5+EstimOpt.NVarC,1) = EstimOpt.NamesC;
for i = 1: EstimOpt.NClass-1
    Results.R_out(EstimOpt.NVarA+6:EstimOpt.NVarA+5+EstimOpt.NVarC,2+ 3*(i-1):1+3*i) = num2cell(Results.DetailsV((i-1)*EstimOpt.NVarC+1:i*EstimOpt.NVarC,:));
end
Results.R_out(EstimOpt.NVarA+5+EstimOpt.NVarC+2,1) = {'Average class probabilities'};
Results.R_out(EstimOpt.NVarA+5+EstimOpt.NVarC+3,1:EstimOpt.NClass) = num2cell(Results.PClass);
Results.R_out(EstimOpt.NVarA+5+EstimOpt.NVarC+5,1) = {'Model characteristics'};
Results.R_out(EstimOpt.NVarA+5+EstimOpt.NVarC+6:end,1) = {'LL0'; 'LL' ; 'McFadden R2';'Ben-Akiva R2' ;'AIC/n' ; 'n'; 'k'};
Results.R_out(EstimOpt.NVarA+5+EstimOpt.NVarC+6:end,2) = num2cell(Results.stats);



if EstimOpt.Display == 1
    for j = 1:EstimOpt.NClass
        disp(' ')
        disp(num2str(j,'Latent class parameters %1.0f'));
        disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
        disp([char(EstimOpt.NamesA),blanks(EstimOpt.NVarA)',num2str(Results.DetailsA(((j-1)*EstimOpt.NVarA+1):j*EstimOpt.NVarA,1),'%8.4f'), star_sig(Results.DetailsA(((j-1)*EstimOpt.NVarA+1):j*EstimOpt.NVarA,3)), num2str(Results.DetailsA(((j-1)*EstimOpt.NVarA+1):j*EstimOpt.NVarA,2:3),'%8.4f %8.4f')])
    end

    for j = 1:EstimOpt.NClass -1 
        disp(' ')
        disp(num2str(j,'Latent class probability model %1.0f'))
        disp(['var.', blanks(size(char(EstimOpt.NamesC),2)-2) ,'coef.      st.err.  p-value'])
        disp([char(EstimOpt.NamesC) ,blanks(EstimOpt.NVarC)',num2str(Results.DetailsV(((j-1)*EstimOpt.NVarC+1):j*EstimOpt.NVarC,1),'%8.4f'), star_sig(Results.DetailsV(((j-1)*EstimOpt.NVarC+1):j*EstimOpt.NVarC,3)), num2str(Results.DetailsV(((j-1)*EstimOpt.NVarC+1):j*EstimOpt.NVarC,2:3),'%8.4f %8.4f')])
    end
    
    disp(' ')
    disp(['Avarage class probabilities'])
    disp(['class   prob.'])
    disp(num2str([(1:EstimOpt.NClass)',Results.PClass']))
        
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
end
