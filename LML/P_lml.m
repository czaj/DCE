function [P_t,DistStats,Grid_t] = P_lml(b_GridMat,bhat,varargin)

P = mnlquick(b_GridMat,bhat); % Calculate probabilities 1xNGrid

if nargout > 1
    if nargin >= 5
        GridMat = varargin{1};
        iHess = varargin{2};
        EstimOpt = varargin{3};
%         KR_idx = 0;
        if nargin == 6
            KR_idx = varargin{4};
        else
            KR_idx = 0;
        end
    else
        error('Simulating distribution statistics requires providing GridMat, iHess, EstimOpt as inputs');
    end
    M_info = memory;
    if 0.6*M_info.MaxPossibleArrayBytes/8 < size(b_GridMat,1)*size(GridMat,1)*size(GridMat,2)
         cprintf(rgb('DarkOrange'), 'WARNING: Probably not enough memory, switching to numerical approximation (slower, but using less memory) \n')
         h = @(b) sum(mnlquick(b_GridMat,b).*GridMat,2); % Calculates Mean
         %H = jacobianest(h,bhat);
         H = numdiff(h,h(bhat),bhat,1,[]);
         h = @(b) sqrt(sum(mnlquick(b_GridMat,b).*(GridMat.^2),2) - sum(mnlquick(b_GridMat,b).*GridMat,2).^2); % Calculates Std. Dev
         %H2 = jacobianest(h,bhat);
         H2 = numdiff(h,h(bhat),bhat,1,[]);
    else
        [H,H2] = jacobian1(b_GridMat,GridMat,bhat);
    end
    h = @(b) sum(mnlquick(b_GridMat,b).*GridMat,2); % Calculates Mean
    %H = jacobianest(h,bhat);
    M.Mean = [h(bhat),zeros(EstimOpt.NVarA,1),sqrt(diag(H*iHess*H')),pv(h(bhat),sqrt(diag(H*iHess*H')))];
    h = @(b) sqrt(sum(mnlquick(b_GridMat,b).*(GridMat.^2),2) - sum(mnlquick(b_GridMat,b).*GridMat,2).^2); % Calculates Std. Dev
    %H = jacobianest(h,bhat);
    M.Std = [h(bhat),zeros(EstimOpt.NVarA,1),sqrt(diag(H2*iHess*H2')),pv(h(bhat),sqrt(diag(H2*iHess*H2')))];
    [P_t,Grid_t] = P_transform(P,GridMat,EstimOpt.NVarA, EstimOpt.NGrid);
    M.Quantile = zeros(EstimOpt.NVarA,7);
    Quantiles = [0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975];
    for i = 1:EstimOpt.NVarA
        %P_cumsum = cumsum(P_t{i},2);
        P_cumsum = cumsum(P_t(i,~isnan(P_t(i,:))),2);
        %Grid_i = Grid_t{i};
        Grid_i = Grid_t(i,~isnan(Grid_t(i,:)));
        for j = 1:7
            M.Quantile(i,j) = Grid_i(:,find(P_cumsum > Quantiles(j),1));
        end
    end
    if KR_idx == 1 % using K&R (simulates s.e. and 95% c.i.)
        bhat_mtx = mvnrnd(bhat,iHess,EstimOpt.NSdSim);
        Stats_mtx = zeros(EstimOpt.NVarA,7,EstimOpt.NSdSim);
        for i = 1:EstimOpt.NSdSim
            P_cumsum_i = cumsum(mnlquick(b_GridMat,bhat_mtx(i,:)'),2);
            GridMat_cumsum_i = sum(mnlquick(b_GridMat,bhat_mtx(i,:)').*GridMat,2);
            GridMat_std_i = sqrt(sum(mnlquick(b_GridMat,bhat_mtx(i,:)').*GridMat.^2,2) - sum(mnlquick(b_GridMat,bhat_mtx(i,:)').*GridMat,2).^2);
            Stats_mtx(:,:,i) = [GridMat_cumsum_i(:,end),GridMat_std_i,GridMat(:,find(P_cumsum_i > 0.1,1)),GridMat(:,find(P_cumsum_i > 0.25,1)),GridMat(:,find(P_cumsum_i > 0.5,1)),GridMat(:,find(P_cumsum_i > 0.75,1)),GridMat(:,find(P_cumsum_i > 0.9,1))];
        end
        DistStats = cat(3,median(Stats_mtx,3),std(Stats_mtx,[],3),quantile(Stats_mtx,0.025,3),quantile(Stats_mtx,0.975,3)); % NVarA x (mean, std, q0.1, q0.25, q0.5, q0.75, q0.9) x (point,s.e.,l.b.95%,u.b.95%)
        % [M.Mean(:,1),Stats(:,1,1)]
        % [M.Std(:,1),Stats(:,2,1)]
        % [M.Quantile,Stats(:,3:end,1)]
        % replace simulated values for point estimates with values based on bhat:
    end
    
    DistStats(:,1,1) = M.Mean(:,1); % point estimates
    DistStats(:,1,2) = M.Mean(:,3); % s.e.
    DistStats(:,2,1) = M.Std(:,1); % point estimates
    DistStats(:,2,2) = M.Std(:,3); % s.e.
    DistStats(:,3:9,1) = M.Quantile; % 
    
end

end


%% supplementary functions


function PX = mnlquick(b_GridMat,bhat)
    Fit = b_GridMat'*bhat; % NGrid x 1
    Fit(Fit > prctile(Fit,99)) = prctile(Fit,99); % censoring to account for extreme values
    Fit = exp(Fit - max(Fit));
    Fit_sum = sum(Fit);
    PX = (Fit./Fit_sum)';
end

function [P_t,Grid_t] = P_transform(P,Grid,NVarA, NGrid)
%     P_t = cell(NVarA,1);
%     Grid_t = cell(NVarA,1);
%     for i = 1:NVarA
%        Grid_i = Grid(i,:);
%        U = unique(Grid_i);
%        P_i = zeros(length(U),1);
%        for j = 1:length(U)
%           % Indx = find(Grid == U(j));
%            P_i(j) = sum(P(Grid_i == U(j)));
%        end
%        P_t{i} = P_i';
%        Grid_t{i} = U;
%     end
    
    P_t = NaN(NVarA,NGrid);
    Grid_t = NaN(NVarA,NGrid);
    for i = 1:NVarA
       Grid_i = Grid(i,:);
%        U = unique(Grid_i);
%        P_i = zeros(length(U),1);
%        for j = 1:length(U)
%           % Indx = find(Grid == U(j));
%            P_i(j) = sum(P(Grid_i == U(j)));
%        end
       [U,~,c] = unique(Grid_i);
       P_i = accumarray(c,P);
       P_t(i,:) = P_i';
       Grid_t(i,:) = U;
    end
end

function [J1, J2] = jacobian1(b_GridMat,GridMat,bhat) % calculates jacobian for means and variance
    P = mnlquick(b_GridMat,bhat); % 1 x NP*NRep
    PBG = P.*b_GridMat; % Var x NP*NRep
    sumPBG = sum(PBG,2);
    Tmp = PBG - P.*sumPBG; % Var x NP*NRep
    Tmp = reshape(Tmp, [size(Tmp,1), 1, size(Tmp,2)]);
    Grid = reshape(GridMat, [1, size(GridMat,1),size(GridMat,2)]);    
    J1 = sum(Tmp.*Grid,3)'; % NVarA x Var
    
    
    Meanx = sum(P.*GridMat,2);
    Stdx = sqrt(sum(P.*(GridMat.^2),2) - Meanx.^2);
    J2_1 = sum(Tmp.*(Grid.^2),3)'; % NVarA x Var
    J2_2 = 2*Meanx.*J1;
    J2 = (J2_1 - J2_2)./(2*Stdx);
end
