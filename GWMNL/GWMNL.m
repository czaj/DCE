function Results = GWMNL(INPUT,Results_old,EstimOpt,OptimOpt)
% GWMNL creates Geographically Weighted Multinomial Logit Model.
%
% Syntax:   GWMNL(INPUT,EstimOpt,OptimOpt)
%           GWMNL(INPUT,Results_old,EstimOpt,OptimOpt)
%
% Inputs:
%    INPUT - clean, updated INPUT data from DataCleanDCE
%    EstimOpt - Estimation Options (check below)
%    OptimOpt - Optimizer Options define how algorithms converges to the final result. They are set by default based on provided EstimOpt in DataCleanDCE, however, they are subject to change.
%    Results_old - here one can provide old results to use as starting
%    values
%
% EstimOpt Options:
% Set them by e.g. Estimopt.DataFile = 'Project'
%
% GWMNL uses bandwidth:
% •	BandType = 1; by default uses global bandwidth, set to 2 for spatially varying bandwidth, set to 3 for square root of spatially varying bandwidth
% •	BandSearch = 0; model uses fixed bandwidth as provided in BandVal; set to 1 to find optimal bandwidth value
% •	BandVal = 3; default bandwidth value
% •	BandLower = 0; bandwidth lower bound
% •	BandUpper = 5; bandwidth upper bound
% •	Clustered = 1; uses clustered standard errors, set to 0 otherwise
% •	WeightSum = 1; no particular summing, set to 2 for summing to 1, set to 3 for summing to NP
% •	BStrap = 0; set to 1 for bootstrapping
% o	NBStrap = 100; No of bootstrap replications
% 
% 
% General basics:
% •	DataFile – path/name of the .mat data file
% •	Display – 1; shows output, set to 0 to hide it 
% •	ProjectName – Name of the project/model
% •	WTP_space – set to 1 for estimation in WTP space. If missing or set to 0, MNL uses Preference Space
% •	NCT - Number of choice tasks per person 
% •	NAlt - Number of alternatives
% •	NP – Number of respondents
% 
% 
% Variables options:
% •	NamesA – Names of variables in list e.g. {'-Opt out';’-Cost (EUR)'}
% •	ConstVarActive = 0; set to 1 to constrain model parameters to its initial values 
% 
% 
% Modelling options from DataCleanDCE:
% •	NumGrad = 0; uses analytical gradient in calculations, set to 1 for numerical gradient
% •	HessEstFix = 0; Options: 
% o	0 - use optimization Hessian, 
% o	1 - use jacobian-based (BHHH) Hessian, 
% o	2 - use high-precision jacobian-based (BHHH) Hessian,
% o	3 - use numerical Hessian, 
% o	4 - use analytical Hessian
% 
% Example: 
%    Results.GWMNL = GWMNL(INPUT,Results,EstimOpt,OptimOpt);
%
% Author: Mikolaj Czajkowski, Professor
% University of Warsaw, Faculty of Economic Sciences
% email address: mik@czaj.org 
% Website: http://czaj.org/#

% save tmp_MNL
% return

global B_backup

tic 

Results.bhat = [];
Results.R = [];
Results.stats = [];


%% Check data and inputs


if nargin < 3
    error('Too few input arguments for GWMNL(INPUT,EstimOpt,OptimOpt)')
end

warning off MATLAB:mir_warning_maybe_uninitialized_temporary

format shortG;
format compact;



disp(' ');
disp('__________________________________________________________________________________________________________________');
disp(' ');
disp('Estimating Geographically Weighted MNL model ...')
if isfield(EstimOpt, 'WTP_space') == 0
    EstimOpt.WTP_space = 0;
    EstimOpt.WTP_matrix = [];
elseif EstimOpt.WTP_space == 0;
	EstimOpt.WTP_matrix = [];
end
if EstimOpt.WTP_space > 0
    disp('in WTP-space ...')
else
    disp('in preference-space ...') 
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
	end
end


if isfield(EstimOpt,'HessEstFix') == 0 
    EstimOpt.HessEstFix = 1;
end

if isfield(EstimOpt,'NamesA') == 0 || isempty(EstimOpt.NamesA) || length(EstimOpt.NamesA) ~= EstimOpt.NVarA
    EstimOpt.NamesA = (1:EstimOpt.NVarA)';
    EstimOpt.NamesA = cellstr(num2str(EstimOpt.NamesA));
elseif size(EstimOpt.NamesA,1) ~= EstimOpt.NVarA
    EstimOpt.NamesA = EstimOpt.NamesA';
end


if isfield(EstimOpt, 'BandType') == 0
    EstimOpt.BandType =1; % global bandwidth
end

if isfield(EstimOpt, 'BandVal') == 0
    EstimOpt.BandVal =3; % global bandwidth value
end


if isfield(EstimOpt, 'BandSearch') == 0
    EstimOpt.BandSearch = 0; % Estimate model, without Bandwidth search
end
if EstimOpt.BandType == 1
    disp('with global bandwidth')
elseif EstimOpt.BandType == 2
    disp('with spatially varying bandwidth')
elseif EstimOpt.BandType == 3    
    disp('with square root of spatially varying bandwidth')
end

    
if isfield(EstimOpt, 'Clustered') == 0
    EstimOpt.Clustered =1; 
end
if isfield(EstimOpt, 'WeightSum') == 0
    EstimOpt.WeightSum =1; % no particular summing, 2 - summing to 1, 3- summing to NP
end
if isfield(EstimOpt, 'BandLower') == 0
    EstimOpt.BandLower =0; %
end
if isfield(EstimOpt, 'BandUpper') == 0
    EstimOpt.BandUpper = 5; % 
end

if EstimOpt.Clustered == 1
    disp('with clustered standard errors')
end

if EstimOpt.BandSearch == 0
    disp(['Model estimated with fixed bandwidth value = ' ,num2str(EstimOpt.BandVal,'%2.4f')])
else
    disp('Multiple models will be estimated to find optimal bandwidth value')
end

if isfield(INPUT, 'Crds') == 0 || size(INPUT.Crds,1) ~= EstimOpt.NP || size(INPUT.Crds,2) ~= 2
   error('Coordinates are misspecified - check INPUT.Crds') 
end

if isfield(INPUT, 'Crds2') == 0 || size(INPUT.Crds2,1) ~= EstimOpt.NP 
   disp('Only geographical weights included')
   EstimOpt.BandNo = 0;
else
   if isfield(EstimOpt, 'BandIndx') == 0 || length(EstimOpt.BandIndx)~= size(INPUT.Crds2,2)
       cprintf(rgb('DarkOrange'), 'WARNING: Incorrectly specified EstimOpt.BandIndx - assuming only one bandwidth parameter\n')
       EstimOpt.BandIndx = ones(1,size(INPUT.Crds2,2));
       EstimOpt.BandNo = length(unique(EstimOpt.BandIndx));
   else
       EstimOpt.BandIndx = EstimOpt.BandIndx(:)';
       EstimOpt.BandNo = length(unique(EstimOpt.BandIndx));
   end
   if isfield(EstimOpt, 'BandVal2') == 0 || length(EstimOpt.BandVal2)~= EstimOpt.BandNo
      error('Misspecified values of bandwidth parameters for additional weights') 
   else
       EstimOpt.BandVal2 = EstimOpt.BandVal2(:);
   end
end

if EstimOpt.BandSearch == 1 && EstimOpt.BandNo > 0
   error('BandSearch does not work with additional bandwidth parameters') 
end
if isfield(EstimOpt, 'BStrap') == 0
    EstimOpt.BStrap = 0; % No bootstraping  
end

if isfield(EstimOpt, 'NBStrap') == 0
    EstimOpt.NBStrap = 100; % No of bootstrap replications 
end
%% Starting values

if isfield(EstimOpt,'StartMatFull') && EstimOpt.StartMatFull == 1
    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVarA && size(B_backup,2) == EstimOpt.NP
        b0 = B_backup;
        disp('Using the starting values from Backup')
    else
        error('No starting values - use EstimOpt.StartMatFull = 0 and run MNL first');
    end

else
    EstimOpt.StartMatFull = 0;
    if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVarA && size(B_backup,2) == 1
        b0 = B_backup;
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat') && length(Results_old.MNL.bhat) == EstimOpt.NVarA  % MNL starting values provided
        disp('Using MNL results as starting values')
        Results_old.MNL.b0_old = Results_old.MNL.bhat;
        if length(Results_old.MNL.b0_old) ~= EstimOpt.NVarA
            cprintf(rgb('DarkOrange'), 'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.MNL = rmfield(Results_old.MNL,'b0_old');        
        else
            b0 = Results_old.MNL.b0_old(:);
        end
    else
        error('No starting values - run MNL first');
    end
    B_backup = b0(:, ones(EstimOpt.NP,1));
    
end

%% Optimization Options

if (isfield(EstimOpt, 'ConstVarActive') == 0 || EstimOpt.ConstVarActive == 0) && isequal(OptimOpt.Algorithm,'quasi-newton') && isequal(OptimOpt.Hessian,'user-supplied')
    if EstimOpt.Display ~= 0
        cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied Hessian off - quasi-newton algorithm does not use it anyway \n')
    end
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
            cprintf('*Black','ex-post using Dekker (2014) formula \n')
        
    end
else
    if strcmp(OptimOpt.Hessian,'user-supplied')
        cprintf('Hessian: '); cprintf('*Black','user-supplied, BHHH, ')
    else
        cprintf('Hessian: '); cprintf('*Black',['built-in, ' OptimOpt.HessUpdate ', '])
    end
    switch EstimOpt.HessEstFix
        case 0
            cprintf('*Black','retained from optimization \n')
        case 1
            cprintf('*Black','ex-post using Dekker (2014) formula \n')
    end
end 


%% Weights computations and BandWidth search
OptimOpt.Diagnostics = 'off'; 
OptimOpt.Display = 'off';
OptimOpt.OutputFcn = {};

% finding set of unique coordinates
if EstimOpt.BandNo == 0
    Crds = unique(INPUT.Crds, 'rows'); 
    NoReg = size(Crds,1); % no. of local regressions to run
    CrdsVec = pdist(INPUT.Crds); 
    CrdsMat = squareform(CrdsVec);
else
    Crds = unique([INPUT.Crds, INPUT.Crds2], 'rows'); 
    NoReg = size(Crds,1); % no. of local regressions to run
    CrdsVec = pdist(INPUT.Crds); 
    CrdsMat = squareform(CrdsVec);
    CrdsMat2 = zeros(EstimOpt.NP,EstimOpt.NP,EstimOpt.BandNo);
    BandIndx = unique(EstimOpt.BandIndx);
    for i = 1:EstimOpt.BandNo
        CrdsVec2 = pdist(INPUT.Crds2(:,EstimOpt.BandIndx==BandIndx(i))); 
        CrdsMat2(:,:, i) = squareform(CrdsVec2);
    end
end




if EstimOpt.BandSearch ~= 0 
    tic
    disp(' ')
    disp('Starting Bandwidth search... ')
    disp(['in interval: [',num2str(EstimOpt.BandLower,'%3.4f'),', ',num2str(EstimOpt.BandUpper,'%3.4f'), ']' ])
    EstimOpt.b0 = b0;
    FunCV = @(B) BandSearch(INPUT.Y, INPUT.Xa,INPUT.Crds, Crds, CrdsMat ,NoReg, EstimOpt,OptimOpt,B);
    %[Results.BandVal, Results.CV, Results.exitfCV] = fminsearch(FunCV, EstimOpt.BandVal, optimset('OutputFcn',@outputf_gw));
    [Results.BandVal, Results.CV, Results.exitfCV] = fminbnd(FunCV, EstimOpt.BandLower,EstimOpt.BandUpper,  optimset('OutputFcn',@outputf_gw));
    disp(['Bandwidth value found = ',num2str(Results.BandVal,'%3.4f')])
    Results.tocnoteCV = toc;
    disp(['BandWidth search took ' num2str(Results.tocnoteCV) ' seconds ('  num2str(floor(Results.tocnoteCV/(60*60))) ' hours ' num2str(floor(rem(Results.tocnoteCV,60*60)/60)) ' minutes ' num2str(rem(Results.tocnoteCV,60)) ' seconds).']);
    disp(' ')
    disp('Proceeding with model estimation:')
    EstimOpt.BandVal = Results.BandVal;
end



if EstimOpt.BandType == 1 % global kernel
    Weights = exp(-0.5*(CrdsMat/EstimOpt.BandVal).^2);
    
else % spatially varying kernel
    CrdsVec2 = pdist(Crds); 
    CrdsMat2 = squareform(CrdsVec2);
    Rij = zeros(EstimOpt.NP, EstimOpt.NP); 
    for i = 1:NoReg
        
       [B,I] = sort(CrdsMat2(:,i));
       Rtmp = (0:NoReg-1)';
       %Rtmp = Rtmp(I);
       CrdsX = Crds(I,:);
       FindIndx = find(INPUT.Crds(:,1) == Crds(i,1) & INPUT.Crds(:,2) == Crds(i,2));
       for j =1:NoReg
           FindIndx2 = find(INPUT.Crds(:,1) == CrdsX(j,1) & INPUT.Crds(:,2) == CrdsX(j,2));
           Rij(FindIndx2, FindIndx) = Rtmp(j);
       end
    end
    Results.Rij = Rij;
    if EstimOpt.BandType == 2
        Weights = exp(-Rij/EstimOpt.BandVal);
    else
        Weights = exp(-(Rij.^0.5)/EstimOpt.BandVal);
    end
end
if EstimOpt.BandNo > 0
     for i = 1:EstimOpt.BandNo
        Weights = Weights.*exp(-0.5*(CrdsMat2(:,:,i)/EstimOpt.BandVal2(i)).^2);
     end
end
%% Estimation

HX = zeros(NoReg,1);
LLX = zeros(NoReg,1); % sum of LL from given location
LLXX = zeros(NoReg,1); % sum of all LL for every model
LLXXX = zeros(NoReg,1); % sum of all LL for every model, but weighted
%Ytmp = reshape(Y, EstimOpt.NAlt*EstimOpt.NCT, EstimOpt.NP);
%XXa = reshape(Xa, EstimOpt.NAlt*EstimOpt.NCT, EstimOpt.NP, EstimOpt.NVarA);
Xtmp = reshape(INPUT.Xa, EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP, EstimOpt.NVarA);
Xstar = Xtmp(2:end,:,:)- Xtmp(ones(EstimOpt.NAlt-1,1),:,:);
Xstar = reshape(Xstar, (EstimOpt.NAlt-1)*EstimOpt.NCT*EstimOpt.NP, EstimOpt.NVarA);

LocNoX = zeros(NoReg,1); % no. of locations

tic
% OptimOpt.Algorithm = 'trust-region';
% OptimOpt.Hessian = 'off';

if EstimOpt.BStrap == 0
    Results.bhat = zeros(EstimOpt.NVarA, EstimOpt.NP);
    Results.std = zeros(EstimOpt.NVarA, EstimOpt.NP);
    Results.ihess = zeros(EstimOpt.NVarA,EstimOpt.NVarA, EstimOpt.NP);
    Results.LL = zeros(EstimOpt.NP,1);
    Results.exitf = zeros(EstimOpt.NP,1);
    
    for n = 1:NoReg
        disp(['Estimating ', num2str(n,'%1.0f'), '/', num2str(NoReg, '%1.0f'), ' model'])
        if EstimOpt.BandNo == 0
            FindIndx = find(INPUT.Crds(:,1) == Crds(n,1) & INPUT.Crds(:,2) == Crds(n,2));
        else
            TmpX = INPUT.Crds(:,1) == Crds(n,1) & INPUT.Crds(:,2) == Crds(n,2);
            for aa = 1:size(INPUT.Crds2,2)
                TmpX = TmpX & INPUT.Crds2(:,aa) == Crds(n,2+aa);
            end
            FindIndx = find(TmpX);
        end
        LocNo = length(FindIndx);
        LocNoX(n) = LocNo;
        if EstimOpt.WeightSum == 1
            Weights_n = Weights(:,FindIndx(1));
        elseif EstimOpt.WeightSum == 2
            Weights_n = Weights(:,FindIndx(1))/sum(Weights(:,FindIndx(1)),1);
        elseif EstimOpt.WeightSum == 3
            Weights_n = Weights(:,FindIndx(1))/sum(Weights(:,FindIndx(1)),1)*EstimOpt.NP;
        end
        
        %Weights_n = ones(EstimOpt.NP,1);
       % Weights_n = ones(EstimOpt.NP,1)/4.34766037059130;
        LLfun = @(B) LL_gwmnl_MATlike(INPUT.Y, INPUT.Xa, Weights_n, EstimOpt,OptimOpt,B);
        if EstimOpt.StartMatFull == 1
            %[bhat, LL, exitf, output, g, hess] = fminunc(LLfun, b0(:,FindIndx(1)), OptimOpt);
            if EstimOpt.HessEstFix == 1
                [bhat, LL, exitf] = fminunc(LLfun, b0(:,FindIndx(1)), OptimOpt);
            else
                [bhat, LL, exitf, ~, ~, hess] = fminunc(LLfun, b0(:,FindIndx(1)), OptimOpt);
            end
        else
            %[bhat, LL, exitf, output, g, hess] = fminunc(LLfun, b0, OptimOpt);
            if EstimOpt.HessEstFix == 1
                [bhat, LL, exitf] = fminunc(LLfun, b0, OptimOpt);
            else
                [bhat, LL, exitf, ~, ~, hess] = fminunc(LLfun, b0, OptimOpt);
            end
            b0 = bhat;
        end
        Results.bhat(:, FindIndx) =bhat(:, ones(LocNo,1));
        [LL,J, probs] = LL_gwmnl_bs(INPUT.Y, INPUT.Xa,EstimOpt,bhat);
        LL = reshape(LL, EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP);
        LL = LL(2:end,:);
        Weights_tmp = Weights_n;
        Weights_n = reshape(Weights_n(:, ones(EstimOpt.NCT*(EstimOpt.NAlt-1),1))',EstimOpt.NCT*(EstimOpt.NAlt-1)*EstimOpt.NP,1) ;
        XXstar = Xstar.*sqrt(Weights_n(:, ones(EstimOpt.NVarA,1)));
        % To doda³em
        VXl = zeros(EstimOpt.NVarA, (EstimOpt.NAlt-1)*EstimOpt.NCT*EstimOpt.NP);
        for j = 1:EstimOpt.NCT*EstimOpt.NP
            Vloc = -LL(:,j)*LL(:,j)'+diag(LL(:,j));
            VXl(:, (j-1)*(EstimOpt.NAlt-1)+1:j*(EstimOpt.NAlt-1)) = XXstar((j-1)*(EstimOpt.NAlt-1)+1:j*(EstimOpt.NAlt-1),:)'*Vloc;
        end
        if EstimOpt.HessEstFix == 1 % use Dekker formula
            Omega = VXl*XXstar;
        else
           Omega = hess; 
        end
        Indx1 = [];
        Indx2 = [];
        for j = 1:length(FindIndx)
            Indx1 = [Indx1, (FindIndx(j)-1)*EstimOpt.NCT*(EstimOpt.NAlt-1)+1:FindIndx(j)*EstimOpt.NCT*(EstimOpt.NAlt-1)];
            Indx2 = [Indx2, (FindIndx(j)-1)*EstimOpt.NCT+1:FindIndx(j)*EstimOpt.NCT];
        end
        
        if EstimOpt.Clustered == 1
            InvOmega = inv(Omega);
            J = J.*Weights_n(1:EstimOpt.NAlt-1:end,ones(EstimOpt.NVarA,1));
            J = reshape(J, EstimOpt.NCT, EstimOpt.NP, EstimOpt.NVarA);
            J = squeeze(sum(J,1));
            U = J'*J;
            InvOmega = InvOmega*U*InvOmega;
        else
            InvOmega = inv(Omega);
        end
        
        %InvOm = Omega\(XXstar');
        stdx = sqrt(diag(InvOmega));
        Results.ihess(:,:,FindIndx) = InvOmega(:,:, ones(LocNo,1));
        Results.std(:, FindIndx) =stdx(:, ones(LocNo,1));
        HX(n) = trace(XXstar(Indx1,:)*InvOmega*VXl(:,Indx1));
        LLX(n) = sum(probs(Indx2));
        LLXX(n) = sum(probs);
        probs = reshape(probs, EstimOpt.NCT, EstimOpt.NP);
        LLXXX(n) = sum(sum(probs'.*Weights_tmp(:, ones(1,EstimOpt.NCT))));
        %Results.LL(FindIndx) = -LL;
        Results.exitf(FindIndx) = exitf;

        B_backup(:, FindIndx) = bhat(:, ones(LocNo,1));
    end
    Results.AIC = -2*sum(LLX)/(EstimOpt.NCT*EstimOpt.NP)+(2*sum(HX)+1)/(EstimOpt.NCT*EstimOpt.NP-sum(HX)-2);
    Results.LL = LLX;
    Results.LL2 = LLXX;
    Results.LL3 = LLXXX;
    Results.LocNo = LocNoX;
    Results.Htrace = HX;
else
    Results.bhat = zeros(EstimOpt.NVarA, EstimOpt.NP);
    Results.bhatBS = zeros(EstimOpt.NVarA, EstimOpt.NP, EstimOpt.NBStrap);
    Results.std = zeros(EstimOpt.NVarA, EstimOpt.NP);
    %Results.ihess = zeros(EstimOpt.NVarA,EstimOpt.NVarA, EstimOpt.NP);
    Results.LL = zeros(EstimOpt.NP,1);
    Results.exitf = zeros(EstimOpt.NP,1);
    Results.Samples = zeros(EstimOpt.NP, EstimOpt.NBStrap);
    vec = (1:EstimOpt.NP)';
    for i = 1:EstimOpt.NBStrap
        Results.Samples(:,i) = randsample(vec, EstimOpt.NP, true);
    end
    Y = reshape(INPUT.Y, EstimOpt.NAlt*EstimOpt.NCT, EstimOpt.NP);
    Xa = reshape(INPUT.Xa, EstimOpt.NAlt*EstimOpt.NCT, EstimOpt.NP, EstimOpt.NVarA);
    disp(['Bootstrapped standard errors with ', num2str(EstimOpt.NBStrap,'%1.0f'), ' replications'])
    disp(' ')
    for n = 1:NoReg
        disp(['Estimating ', num2str(n,'%1.0f'), '/', num2str(NoReg, '%1.0f'), ' model'])

        FindIndx = find(INPUT.Crds(:,1) == Crds(n,1) & INPUT.Crds(:,2) == Crds(n,2));
        LocNo = length(FindIndx);
        Weights_n = Weights(:,FindIndx(1))/sum(Weights(:,FindIndx(1)),1)*EstimOpt.NP;

        LLfun = @(B) LL_gwmnl_MATlike(INPUT.Y, INPUT.Xa, Weights_n, EstimOpt,OptimOpt,B);
        if EstimOpt.StartMatFull == 1
            [bhat, LL, exitf] = fminunc(LLfun, b0(:,FindIndx(1)), OptimOpt);
        else
            [bhat, LL, exitf] = fminunc(LLfun, b0, OptimOpt);
            b0 = bhat;
        end
        for j = 1:EstimOpt.NBStrap
            cprintf('.')
            WeightsBS = Weights_n(Results.Samples(:,j));
            YBS = Y(:, Results.Samples(:,j));
            XaBS = Xa(:,Results.Samples(:,j),:);
            YBS = reshape(YBS, EstimOpt.NAlt*EstimOpt.NCT*EstimOpt.NP,1);
            XaBS = reshape(XaBS, EstimOpt.NAlt*EstimOpt.NCT*EstimOpt.NP,EstimOpt.NVarA);
            LLfunBS = @(B) LL_gwmnl_MATlike(YBS, XaBS, WeightsBS, EstimOpt,OptimOpt,B);
            if EstimOpt.StartMatFull == 1
                bhatBS = fminunc(LLfunBS, b0(:,FindIndx(1)), OptimOpt);
            else
                bhatBS = fminunc(LLfunBS, b0, OptimOpt);
            end
            Results.bhatBS(:, FindIndx, j) = bhatBS(:, ones(LocNo,1));
        end
        cprintf('\n')
        Results.bhat(:, FindIndx) =bhat(:, ones(LocNo,1));
        stdx = squeeze(Results.bhatBS(:, FindIndx(1), :));
        stdx = std(stdx,0,2);
        Results.std(:, FindIndx) =stdx(:, ones(LocNo,1));
        Results.LL(FindIndx) = -LL;
        Results.exitf(FindIndx) = exitf;

        B_backup(:, FindIndx) = bhat(:, ones(LocNo,1));
    end
end


Results.NoReg = NoReg;
Results.Weights = Weights;
Results.Crds = Crds;
Results.CrdsMat = CrdsMat;
%% Output

Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;

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
