function [B_out,Summary,Res] = LML_search(INPUT,Results,EstimOpt,OptimOpt,varargin)

% save LML_search_out

% INPUT:
% EstimOpt.LMLSearchNOrder - maximum NOrder
% EstimOpt.LMLSearchNTrials - number of search trials
% optionally provide starting B matrix (NOrder x NParam x NDist x 2 (FullCov = [0,1])) - can be smaller or larger, use NaN for missing parameter values

% OUTPUT:
% B0 - B matrix of the best specification found (EstimOpt.LMLSearchNOrder-1 x NParam x NDist x 2)
% Summary - summary of the best specifications found ((EstimOpt.LMLSearchNOrder-1) * NDist x [NOrder,Dist,FullCov,LL,NParam,AIC/n,BIC/n]*2+1)
% Res - Results of all trials (6 (NStats) x EstimOpt.LMLSearchNOrder-1 x NDist x 2 x EstimOpt.LMLSearchNTrials)

% ISSUES / TODO
% FullCov = 1 may not work correctly with EstimOpt.StepVar > 0

if ~isfield(EstimOpt,'WTP_space')
    EstimOpt.WTP_space = 0;
end

try
    
    global B_backup
    
    if isfield(EstimOpt,'LMLSearchNOrder') == 0
        EstimOpt.LMLSearchNOrder = 10;
    elseif EstimOpt.LMLSearchNOrder < 2
        error('EstimOpt.LMLSearchNOrder cannot be < 2')
    end
    
    if isfield(EstimOpt,'LMLSearchNTrials') == 0
        EstimOpt.LMLSearchNTrials = 10;
    elseif mod(EstimOpt.LMLSearchNTrials,1) ~= 0
        error('EstimOpt.LMLSearchNTrials must be integer')
    end
    
    % Prepare B0 matrix:
    if isfield(EstimOpt, 'StepFun') == 1
        EstimOpt.StepVar = size(EstimOpt.StepFun(ones(EstimOpt.NVarA,1)),1);
    else
        EstimOpt.StepVar = 0;
    end
    B0 = NaN(EstimOpt.LMLSearchNOrder,EstimOpt.NVarA.*(EstimOpt.LMLSearchNOrder+1) + EstimOpt.StepVar + EstimOpt.NVarA*(EstimOpt.NVarA-1)/2,7,2);
    Summary = NaN((EstimOpt.LMLSearchNOrder-1)*7,15);
    SimIdx = true(EstimOpt.LMLSearchNOrder-1,7,2);
    
    if nargin > 4
        B0_in = varargin{1}; % NOrder x NParams x NDist x 2 (FullCov = [0,1])
%             save tmp1
        %     EstimOpt.LMLSearchNOrder = EstimOpt.LMLSearchNOrder + 1; % an additional evaluation for the starting valuaes provided
        %     % test B0 size
        %     if ndims(B0_in) ~= 4
        %         cprintf(rgb('DarkOrange'), 'WARNING: incorrect number of dimensions of starting values matrix - must be 4'); % could be modified to work with uncorrelated LML only
        %     end
        if (ndims(B0_in) ~= ndims(B0)) || (any(size(B0_in) ~= size(B0)))
            B0_tmp = NaN(size(B0));
            B0_tmp(1:min(size(B0_in,1),size(B0,1)),1:min(size(B0_in,2),size(B0,2)),1:min(size(B0_in,3),size(B0,3)),1:min(size(B0_in,4),size(B0,4))) = B0_in(1:min(size(B0_in,1),size(B0,1)),1:min(size(B0_in,2),size(B0,2)),1:min(size(B0_in,3),size(B0,3)),1:min(size(B0_in,4),size(B0,4)));
            B0_in = B0_tmp; % Extend B0_in to match B0
        end
        %     if any(size(B0_in) > size(B0))
        %         B0_in = B0_in(1:size(B0,1),1:size(B0,2),1:size(B0,3),1:size(B0,4));
        %     end
        B0(~isnan(B0_in)) = B0_in(~isnan(B0_in));
        if nargin > 5
            Summary_in = varargin{2}; % NOrder x (2*NParams + 1)
            if any(size(Summary_in) ~= size(Summary))
                Summary_tmp = NaN(size(Summary));
                Summary_tmp(1:min(size(Summary_in,1),size(Summary,1)),1:min(size(Summary_in,2),size(Summary,2)),1:min(size(Summary_in,3),size(Summary,3)),1:min(size(Summary_in,4),size(Summary,4))) = Summary_in(1:min(size(Summary_in,1),size(Summary,1)),1:min(size(Summary_in,2),size(Summary,2)),1:min(size(Summary_in,3),size(Summary,3)),1:min(size(Summary_in,4),size(Summary,4)));
                Summary_in = Summary_tmp; % Extend B0_in to match B0
            end
            Summary(~isnan(Summary_in)) = Summary_in(~isnan(Summary_in));
            if nargin > 6
                Res_in = varargin{3};
                if nargin > 7
                    SimIdx_in = varargin{4}; % EstimOpt.LMLSearchNOrder-1 x 7 (Dist) x 2
                    if any(size(SimIdx_in) ~= [EstimOpt.LMLSearchNOrder-1,7,2])
                        error('Incorrect size of the matrix providing indices for simulations to carry out (needs to be EstimOpt.LMLSearchNOrder-1 x 7 (NDist) x 2 (FullCov))')
                    else
                        SimIdx = SimIdx_in;
                    end
                end
            end
        end
    end
    
    EstimOpt.PlotIndx = 0; % No plots
    EstimOpt.NoOutput = 1; % No output
    
    % save tmp1
    % return
    
    B0_start = unifrnd(-5,5,[size(B0),EstimOpt.LMLSearchNTrials]); % Random starting values    
%     Res = cell(3,EstimOpt.LMLSearchNOrder-1,7,2,EstimOpt.LMLSearchNTrials); % To save results    
    B_out = NaN(size(B0));
    LL = NaN(EstimOpt.LMLSearchNOrder-1,7,2);
    Res.bhat = NaN(EstimOpt.LMLSearchNOrder-1, 7, 2, EstimOpt.LMLSearchNTrials, EstimOpt.NVarA.*(EstimOpt.LMLSearchNOrder+1) + EstimOpt.StepVar + EstimOpt.NVarA*(EstimOpt.NVarA-1)/2);
    Res.LL = NaN(EstimOpt.LMLSearchNOrder-1, 7, 2, EstimOpt.LMLSearchNTrials);
    Res.stats = NaN(EstimOpt.LMLSearchNOrder-1, 7, 2, EstimOpt.LMLSearchNTrials,9);
    
%     save tmp1
    
    for j = 1:7 % Loop over Dist
        if j == 1
            Dist = 0;
            EstimOpt.Dist = [Dist * ones(size(INPUT.Xa,2)-1,1); Dist + (EstimOpt.WTP_space==1)]; % for WTP-space models cost (the last attribute) is always log-normal
        elseif j == 2
            Dist = j;
            EstimOpt.Dist = [Dist * ones(size(INPUT.Xa,2)-1,1); Dist + (EstimOpt.WTP_space==1)]; % for WTP-space models cost (the last attribute) is always log-normal
        elseif j > 2
            Dist = j+1;
            EstimOpt.Dist = [Dist * ones(size(INPUT.Xa,2)-1,1); Dist];
        end
        
        for i = 2:EstimOpt.LMLSearchNOrder % Loop over NOrder
            
            EstimOpt.NOrder = i;
            
            if i == 2 % NOrder = 2
                
                if SimIdx(i-1,j,1) == 1
                    
                    disp(['Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 0, Trial = ',num2str(1),'/',num2str(EstimOpt.LMLSearchNTrials)])
                    
                    EstimOpt.FullCov = 0;
                    
                    NVar = sum((EstimOpt.Dist == 0 | EstimOpt.Dist == 1)*EstimOpt.NOrder + ...
                        (EstimOpt.Dist == 2 | EstimOpt.Dist == 3)*EstimOpt.NOrder + ...
                        (EstimOpt.Dist == 4)*(EstimOpt.NOrder-1) + ...
                        (EstimOpt.Dist == 5 | EstimOpt.Dist == 6 | EstimOpt.Dist == 7 | EstimOpt.Dist == 8)*(EstimOpt.NOrder+1),1) + ...
                        EstimOpt.StepVar;
                    
                    B_backup = B0(i-1,~isnan(B0(i-1,:,j,1)),j,1); % Start from B0_in or 0
                    if size(B_backup(:),1) ~= NVar
                        B_backup = zeros(1,NVar); % If B_in not provided use vector of 0s as starting values in Trial 1
                    else
                        B0_start(i-1,1:NVar,j,1,2) = zeros(1,NVar); % If B_in provided use vector of 0s as starting values in Trial 2
                    end
                    
                    try
                        Results.LML_d = LML(INPUT,Results,EstimOpt,OptimOpt);
                        B_out(i-1,1:size(Results.LML_d.bhat,1),j,1) = Results.LML_d.bhat;
                        LL(i-1,j,1) = Results.LML_d.LL;
                        Summary(j + (i-2)*7,1:7) = [EstimOpt.NOrder,Dist,0,Results.LML_d.LL,Results.LML_d.stats(9,1),Results.LML_d.stats(5:6,1)'];
                        Res.bhat(i-1,j,1,1,1:length(Results.LML_d.bhat)) = Results.LML_d.bhat;
                        Res.LL(i-1,j,1,1) = Results.LML_d.LL;
                        Res.stats(i-1,j,1,1,1:length(Results.LML_d.stats)) = Results.LML_d.stats;
                    catch theErrorInfo
                        cprintf(rgb('DarkOrange'), ['WARNING: Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 0, Trial = ',num2str(1),'/',num2str(EstimOpt.LMLSearchNTrials),' - resulted in estimation error\n']); 
%                         rethrow(theErrorInfo)
                    end
                    
                    for k = 1:EstimOpt.LMLSearchNTrials-1
                        
                        disp(['Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 0, Trial = ',num2str(k+1),'/',num2str(EstimOpt.LMLSearchNTrials)])
                        
                        B_backup = B0_start(i-1,1:NVar,j,1,k);
                        
                        try
                            Results.LML_d = LML(INPUT,Results,EstimOpt,OptimOpt);
                            if isnan(LL(i-1,j,1)) || (Results.LML_d.LL > LL(i-1,j,1)) % Update LL and B_out to keep the best results
                                LL(i-1,j,1) = Results.LML_d.LL;
                                B_out(i-1,1:size(Results.LML_d.bhat,1),j,1) = Results.LML_d.bhat;
                                Summary(j + (i-2)*7,1:7) = [EstimOpt.NOrder,Dist,0,Results.LML_d.LL,Results.LML_d.stats(9,1),Results.LML_d.stats(5:6,1)'];
                            end
                            Res.bhat(i-1,j,1,k+1,1:length(Results.LML_d.bhat)) = Results.LML_d.bhat;
                            Res.LL(i-1,j,1,k+1) = Results.LML_d.LL;
                            Res.stats(i-1,j,1,k+1,1:length(Results.LML_d.stats)) = Results.LML_d.stats;
                        catch theErrorInfo
                            cprintf(rgb('DarkOrange'), ['WARNING: Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 0, Trial = ',num2str(k+1),'/',num2str(EstimOpt.LMLSearchNTrials),' - resulted in estimation error\n'])
%                             rethrow(theErrorInfo)
                        end
                        
                    end
                end
                
                if SimIdx(i-1,j,2) == 1
                    
                    EstimOpt.FullCov = 1;
                    
                    disp(['Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 1, Trial = ',num2str(1),'/',num2str(EstimOpt.LMLSearchNTrials)])
                    
                    B_backup = B0(i-1,~isnan(B0(i-1,:,j,2)),j,2); % Start from B0_in or 0
                    if size(B_backup(:),1) ~= (NVar + EstimOpt.NVarA*(EstimOpt.NVarA-1)/2)
                        %B_backup = zeros(1,NVar + EstimOpt.NVarA*(EstimOpt.NVarA-1)/2); % If B_in not provided use vector of 0s as starting values in Trial 1
                        B_backup = [B_out(i-1,1:NVar,j,1)';zeros(EstimOpt.NVarA*(EstimOpt.NVarA-1)/2,1)]; % Use values from model without correlations + zeros
                    else
                        B0_start(i-1,1:NVar + EstimOpt.NVarA*(EstimOpt.NVarA-1)/2,j,2,2) = zeros(1,NVar + EstimOpt.NVarA*(EstimOpt.NVarA-1)/2); % If B_in provided use vector of 0s as starting values in Trial 2
                    end
                    
                    try
                        Results.LML = LML(INPUT,Results,EstimOpt,OptimOpt);
                        B_out(i-1,1:size(Results.LML.bhat,1),j,2) = Results.LML.bhat;
                        LL(i-1,j,2) = Results.LML_d.LL;
                        Summary(j + (i-2)*7,9:15) = [EstimOpt.NOrder,Dist,1,Results.LML.LL,Results.LML.stats(9,1),Results.LML.stats(5:6,1)'];
                        Res.bhat(i-1,j,2,1,1:length(Results.LML.bhat)) = Results.LML.bhat;
                        Res.LL(i-1,j,2,1) = Results.LML.LL;
                        Res.stats(i-1,j,2,1,1:length(Results.LML.stats)) = Results.LML.stats;
                    catch theErrorInfo
                        cprintf(rgb('DarkOrange'), ['WARNING: Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 1, Trial = ',num2str(1),'/',num2str(EstimOpt.LMLSearchNTrials),' - resulted in estimation error\n'])
%                         rethrow(theErrorInfo)
                    end
                    
                    for k = 1:EstimOpt.LMLSearchNTrials-1
                        disp(['Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 1, Trial = ',num2str(k+1),'/',num2str(EstimOpt.LMLSearchNTrials)])
                        
                        B_backup = B0_start(i-1,1:NVar + EstimOpt.NVarA*(EstimOpt.NVarA-1)/2,j,2,k);
                        
                        try
                            Results.LML = LML(INPUT,Results,EstimOpt,OptimOpt);
                            if isnan(LL(i-1,j,2)) || (Results.LML.LL > LL(i-1,j,2)) % Update LL and B_out to keep the best results
                                LL(i-1,j,2) = Results.LML.LL;
                                B_out(i-1,1:size(Results.LML.bhat,1),j,2) = Results.LML.bhat;
                                Summary(j + (i-2)*7,9:15) = [EstimOpt.NOrder,Dist,1,Results.LML.LL,Results.LML.stats(9,1),Results.LML.stats(5:6,1)'];
                            end
                            Res.bhat(i-1,j,2,k+1,1:length(Results.LML.bhat)) = Results.LML.bhat;
                            Res.LL(i-1,j,2,k+1) = Results.LML.LL;
                            Res.stats(i-1,j,2,k+1,1:length(Results.LML.stats)) = Results.LML.stats;
                        catch theErrorInfo
                            cprintf(rgb('DarkOrange'), ['WARNING: Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 1, Trial = ',num2str(k+1),'/',num2str(EstimOpt.LMLSearchNTrials),' - resulted in estimation error\n'])
%                             rethrow(theErrorInfo)
                        end
                    end
                end
                
            else % NOrder > 2
                
                if SimIdx(i-1,j,1) == 1
                    
                    disp(['Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 0, Trial = ',num2str(1),'/',num2str(EstimOpt.LMLSearchNTrials)])
                    
                    EstimOpt.FullCov = 0;
                    
                    NVar = sum((EstimOpt.Dist == 0 | EstimOpt.Dist == 1)*EstimOpt.NOrder + ...
                        (EstimOpt.Dist == 2 | EstimOpt.Dist == 3)*EstimOpt.NOrder + ...
                        (EstimOpt.Dist == 4)*(EstimOpt.NOrder-1) + ...
                        (EstimOpt.Dist == 5 | EstimOpt.Dist == 6 | EstimOpt.Dist == 7 | EstimOpt.Dist == 8)*(EstimOpt.NOrder+1),1) + ...
                        EstimOpt.StepVar;
                    
                    B_backup = B0(i-1,~isnan(B0(i-1,:,j,1)),j,1); % Start from B0_in or 0
                    if size(B_backup(:),1) ~= NVar
                        if (j == 1) || (j == 2)
                            B_backup = B_out(i-2,1:NVar,j,1);
                            B_backup(isnan(B_backup)) = 0; % If B_in not provided use vector of Polynomial of the i-1'th order and 0s as starting values in Trial 1
                        else % j > 2
                            B_backup = zeros(1,NVar); % If B_in not provided use vector of 0s as starting values in Trial 1
                        end
                    else
                        B0_start(i-1,1:NVar,j,1,2) = zeros(1,NVar); % If B_in provided use vector of 0s as starting values in Trial 2
                    end
                    
                    try
                        Results.LML_d = LML(INPUT,Results,EstimOpt,OptimOpt);
                        B_out(i-1,1:size(Results.LML_d.bhat,1),j,1) = Results.LML_d.bhat;
                        LL(i-1,j,1) = Results.LML_d.LL;
                        Summary(j + (i-2)*7,1:7) = [EstimOpt.NOrder,Dist,0,Results.LML_d.LL,Results.LML_d.stats(9,1),Results.LML_d.stats(5:6,1)'];
                        Res.bhat(i-1,j,1,1,1:length(Results.LML_d.bhat)) = Results.LML_d.bhat;
                        Res.LL(i-1,j,1,1) = Results.LML_d.LL;
                        Res.stats(i-1,j,1,1,1:length(Results.LML_d.stats)) = Results.LML_d.stats;
                    catch theErrorInfo
                        cprintf(rgb('DarkOrange'), ['WARNING: Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 0, Trial = ',num2str(1),'/',num2str(EstimOpt.LMLSearchNTrials),' - resulted in estimation error\n'])
%                         rethrow(theErrorInfo)
                    end
                    
                    for k = 1:EstimOpt.LMLSearchNTrials-1
                        
                        disp(['Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 0, Trial = ',num2str(k+1),'/',num2str(EstimOpt.LMLSearchNTrials)])
                        
                        if (j == 1) || (j == 2)
                            B_backup = B_out(i-2,1:NVar,j,1); % This time use B_out from the same NOrder
                            B_backup(isnan(B_out(i-2,1:NVar,j,1))) = B0_start(i-1,isnan(B_out(i-2,1:NVar,j,1)),j,1,k); % but only keep elements for which of NOdrer-1 is not NaN, use B0_start for others
                        else % j > 2
                            B_backup = B0_start(i-1,isnan(B_out(i-2,1:NVar,j,1)),j,1,k); % If B_in not provided use vector of 0s as starting values in Trial 1
                        end
                        
                        try
                            Results.LML_d = LML(INPUT,Results,EstimOpt,OptimOpt);
                            if isnan(LL(i-1,j,1)) || (Results.LML_d.LL > LL(i-1,j,1)) % Update LL and B_out to keep the best results
                                LL(i-1,j,1) = Results.LML_d.LL;
                                B_out(i-1,1:size(Results.LML_d.bhat,1),j,1) = Results.LML_d.bhat;
                                Summary(j + (i-2)*7,1:7) = [EstimOpt.NOrder,Dist,0,Results.LML_d.LL,Results.LML_d.stats(9,1),Results.LML_d.stats(5:6,1)'];
                            end
                            Res.bhat(i-1,j,1,k+1,1:length(Results.LML_d.bhat)) = Results.LML_d.bhat;
                            Res.LL(i-1,j,1,k+1) = Results.LML_d.LL;
                            Res.stats(i-1,j,1,k+1,1:length(Results.LML_d.stats)) = Results.LML_d.stats;
                        catch theErrorInfo
                            cprintf(rgb('DarkOrange'), ['WARNING: Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 0, Trial = ',num2str(k+1),'/',num2str(EstimOpt.LMLSearchNTrials),' - resulted in estimation error\n'])
%                             rethrow(theErrorInfo)
                        end
                        
                    end
                end
                
                if SimIdx(i-1,j,2) == 1
                    
                    EstimOpt.FullCov = 1;
                    
                    disp(['Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 1, Trial = ',num2str(1),'/',num2str(EstimOpt.LMLSearchNTrials)])
                    
                    B_backup = B0(i-1,~isnan(B0(i-1,:,j,1)),j,2); % Start from B0_in or 0
                    if size(B_backup(:),1) ~= NVar + EstimOpt.NVarA*(EstimOpt.NVarA-1)/2
                        if (j == 1) || (j == 2)
                            B_backup = [B_out(i-2,1:NVar - EstimOpt.NVarA,j,2),zeros(1,EstimOpt.NVarA,1),B_out(i-2,NVar - EstimOpt.NVarA + 1:NVar - EstimOpt.NVarA + (EstimOpt.NVarA)*(EstimOpt.NVarA - 1)/2,j,2)]; % If B_in not provided use vector of Polynomial of the i-1'th order and 0s as starting values in Trial 1
                        else % j > 2
                            B_backup = zeros(1,NVar + EstimOpt.NVarA*(EstimOpt.NVarA-1)/2); % If B_in not provided use vector of 0s as starting values in Trial 1
                        end
                    else
                        B0_start(i-1,1:NVar + EstimOpt.NVarA*(EstimOpt.NVarA-1)/2,j,2,2) = zeros(1,NVar + EstimOpt.NVarA*(EstimOpt.NVarA-1)/2); % If B_in provided use vector of 0s as starting values in Trial 2
                    end
                    
                    try
                        Results.LML = LML(INPUT,Results,EstimOpt,OptimOpt);
                        B_out(i-1,1:size(Results.LML.bhat,1),j,2) = Results.LML.bhat;
                        LL(i-1,j,2) = Results.LML.LL;
                        Summary(j + (i-2)*7,9:15) = [EstimOpt.NOrder,Dist,1,Results.LML.LL,Results.LML.stats(9,1),Results.LML.stats(5:6,1)'];
                        Res.bhat(i-1,j,2,1,1:length(Results.LML.bhat)) = Results.LML.bhat;
                        Res.LL(i-1,j,2,1) = Results.LML.LL;
                        Res.stats(i-1,j,2,1,1:length(Results.LML.stats)) = Results.LML.stats;
                    catch theErrorInfo
                        cprintf(rgb('DarkOrange'), ['WARNING: Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 1, Trial = ',num2str(1),'/',num2str(EstimOpt.LMLSearchNTrials),' - resulted in estimation error\n'])
%                         rethrow(theErrorInfo)
                    end
                    
                    for k = 1:EstimOpt.LMLSearchNTrials-1
                        
                        disp(['Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 1, Trial = ',num2str(k+1),'/',num2str(EstimOpt.LMLSearchNTrials)])
                        
                        if (j == 1) || (j == 2)
                            B_backup = [B_out(i-2,1:NVar - EstimOpt.NVarA,j,2),B0_start(i-1,NVar - EstimOpt.NVarA + 1:NVar,j,2,k),B_out(i-2,NVar - EstimOpt.NVarA + 1:NVar - EstimOpt.NVarA + (EstimOpt.NVarA)*(EstimOpt.NVarA - 1)/2,j,2)]; % This time use B_out from the same NOrder, and for new elements (in the middle) use values from B0_start
                        else % j > 2
                            B_backup = B0_start(i-1,isnan(B_out(i-2,1:NVar + EstimOpt.NVarA*(EstimOpt.NVarA-1)/2,j,1)),j,2,k); % If B_in not provided use vector of 0s as starting values in Trial 1
                        end      
                        
                        try
                            Results.LML = LML(INPUT,Results,EstimOpt,OptimOpt);
                            if isnan(LL(i-1,j,2)) || (Results.LML.LL > LL(i-1,j,2)) % Update LL and B_out to keep the best results
                                LL(i-1,j,2) = Results.LML.LL;
                                B_out(i-1,1:size(Results.LML.bhat,1),j,2) = Results.LML.bhat;
                                Summary(j + (i-2)*7,9:15) = [EstimOpt.NOrder,Dist,1,Results.LML.LL,Results.LML.stats(9,1),Results.LML.stats(5:6,1)'];
                            end
                            Res.bhat(i-1,j,2,k+1,1:length(Results.LML.bhat)) = Results.LML.bhat;
                            Res.LL(i-1,j,2,k+1) = Results.LML.LL;
                            Res.stats(i-1,j,2,k+1,1:length(Results.LML.stats)) = Results.LML.stats;
                        catch theErrorInfo
                            cprintf(rgb('DarkOrange'), ['WARNING: Dist = ',num2str(j),'/',num2str(7),', NOrder = ',num2str(i),'/',num2str(EstimOpt.LMLSearchNOrder),', FullCov = 1, Trial = ',num2str(k+1),'/',num2str(EstimOpt.LMLSearchNTrials),' - resulted in estimation error\n'])
%                             rethrow(theErrorInfo)
                        end
                    end
                end
            end
        end
    end
    
%     save tmp1

    if exist('Res_in','var') == 1 && ~isempty(Res_in)
        Res.bhat = cat(4,Res.bhat,Res_in.bhat);
        Res.LL = cat(4,Res.LL,Res_in.LL);
        Res.stats = cat(4,Res.stats,Res_in.stats);
        %         {:,:,:,:,size(Res,5)+1:size(Res,5)+size(Res_in,5)} = cat(5,Res,Res_in);
        %         Res = cat(max(ndims(Res_in),ndims(Res))+1,Res_in,Res);
    end
    
catch theErrorInfo
    
%     save LML_search_results % save current state in case of error
%     rethrow(theErrorInfo)
    
end


% Summary1 = zeros(EstimOpt.LMLSearchNOrder - 1,6,7); % for FullCov = 0
% Summary2 = zeros(EstimOpt.LMLSearchNOrder - 1,6,7); % for FullCov = 1
%
% for i = 2:EstimOpt.LMLSearchNOrder
%     for j = 1:7
%         Tmp = Res{3,i-1,j,1};
%         Summary1(i-1,:,j) = [i;0;Tmp([1 9 5 6])]; % NOrder FullCov LL Params  AIC BIC
%
%         Tmp = Res{3,i-1,j,2};
%         Summary2(i-1,:,j) = [i;1;Tmp([1 9 5 6])]; % NOrder FullCov LL Params  AIC BIC
%     end
% end

%         R((j-2)*7+i,1:15) = [EstimOpt.NOrder, 0, Results.LML_d.EstimOpt.Dist(1), Results.LML_d.LL, Results.LML_d.stats(9), Results.LML_d.stats(5), Results.LML_d.stats(6), ...
%             NaN, EstimOpt.NOrder, 1, Results.LML.EstimOpt.Dist(1), Results.LML.LL, Results.LML.stats(9), Results.LML.stats(5), Results.LML.stats(6)];

% Summary = cat(4,Summary1,Summary2); % NOrder x stats x dist
