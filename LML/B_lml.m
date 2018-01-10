function b_mtx = B_lml(GridMat,EstimOpt)

% save tmp_B_lml
% return

% TODO: make NOrder NVarA-specific

NP = EstimOpt.NP;
NRep = EstimOpt.NRep;
NVarA = EstimOpt.NVarA;
Dist = EstimOpt.Dist;
NOrder = EstimOpt.NOrder;
FullCov = EstimOpt.FullCov;
Bounds = EstimOpt.Bounds;

NVarAApprox = sum(Dist == 0 | Dist == 1);
NVarAPoly = sum(Dist == 2 | Dist == 3);
NVarAStep = sum(Dist == 4);
% NVarASpline = sum(Dist == 5 | Dist == 6 | Dist == 7 | Dist == 8);
NVarALSpline = sum(Dist == 5); % Linear (bounds extension does not matter)
NVarACSpline = sum(Dist == 6); % Cubic (with bounds extension)
NVarAPCSpline = sum(Dist == 7); % Piece-wise Cubic (with bounds extension)
NVarAPCHISpline = sum(Dist == 8); % Piece-wise Cubic Hermite Interpolating (with bounds extension)

if FullCov == 1 && any(Dist == 3)
    GridMat_old = GridMat;
end

b_mtx = NaN([NVarA,size(GridMat,2),NOrder+1]); % to be cut later

% Approximate normal and lognormal:
if NVarAApprox > 0
    b_mtx(Dist == 0,:,1) = GridMat(Dist == 0,:);
    b_mtx(Dist == 1,:,1) = log(GridMat(Dist == 1,:));
    b_mtx(Dist == 0,:,2) = GridMat(Dist == 0,:).^2;
    b_mtx(Dist == 1,:,2) = log(GridMat(Dist == 1,:)).^2;
end

% Polynomials:
if NVarAPoly > 0
    GridMat(Dist == 2,:) = (GridMat(Dist == 2,:) - Bounds(Dist == 2,1))./(Bounds(Dist == 2,2) - Bounds(Dist == 2,1));
    GridMat(Dist == 3,:) = log(GridMat(Dist == 3,:));
    GridMat(Dist == 3,:) = (GridMat(Dist == 3,:) - Bounds(Dist == 3,1))./(Bounds(Dist == 3,2) - Bounds(Dist == 3,1));
    b_mtx(Dist == 2 | Dist == 3,:,1) = GridMat(Dist == 2 | Dist == 3,:);
    b_mtx(Dist == 2 | Dist == 3,:,2) = ((2*2-1)/2)*GridMat(Dist == 2 | Dist == 3,:).*b_mtx(Dist == 2 | Dist == 3,:,1)-(2-1)/2;
    if NOrder >= 3
        b_mtx(Dist == 2 | Dist == 3,:,3) = ((2*3-1)/3)*GridMat(Dist == 2 | Dist == 3,:).*b_mtx(Dist == 2 | Dist == 3,:,2)-((3-1)/3)*b_mtx(Dist == 2 | Dist == 3,:,1);
    end
    if NOrder > 3
        for i = 4:NOrder
            b_mtx(Dist == 2 | Dist == 3,:,i) = ((2*i-1)/i)*GridMat(Dist == 2 | Dist == 3,:).*b_mtx(Dist == 2 | Dist == 3,:,i-1)-((i-1)/i)*b_mtx(Dist == 2 | Dist == 3,:,i-2);
        end
    end
end

% Step functions:
if NVarAStep > 0
    BoundsStep = EstimOpt.Bounds(Dist == 4,:);
    Thresholds = zeros(NVarAStep,NOrder+1);
    for i = 1:NVarAStep
        Thresholds(i,:) = BoundsStep(i,1) : (BoundsStep(i,2) - BoundsStep(i,1))/(NOrder) : BoundsStep(i,2); % this is marginally different than train's
    end
    for i = 1:NOrder-1 % the last level is the reference
        b_mtx(Dist == 4,:,i) = (GridMat(Dist == 4,:) >= Thresholds(:,i)) & (GridMat(Dist == 4,:) < Thresholds(:,i+1)); % NVarAStep x NRep*NP x NOrder
    end
end

% Splines:
if NVarALSpline > 0 % Linear
    BoundsKnots = EstimOpt.Bounds(Dist == 5,:);
    Knots = zeros(NVarALSpline,NOrder+2);
    for i = 1:NVarALSpline
        Knots(i,:) = BoundsKnots(i,1) : (BoundsKnots(i,2) - BoundsKnots(i,1))/(NOrder+1) : BoundsKnots(i,2); % this is marginally different than train's
    end
    b_mtx(Dist == 5,:,1:NOrder+1) = 0;
    for i = 1:NOrder+1 % the last level is the reference
        segment_idx = (GridMat(Dist == 5,:) >= Knots(:,i)) & (GridMat(Dist == 5,:) < Knots(:,i+1)); % NVarAStep x NRep*NP
        b_mtx(Dist == 5,:,i) = b_mtx(Dist == 5,:,i) + segment_idx.*(1 - (GridMat(Dist == 5,:) - Knots(:,i))./(Knots(:,i+1)-Knots(:,i))); % NVarASpline x NRep*NP x NOrder
        if i < NOrder+1
            b_mtx(Dist == 5,:,i+1) = segment_idx.*((GridMat(Dist == 5,:) - Knots(:,i))./(Knots(:,i+1)-Knots(:,i))); % NVarASpline x NRep*NP x NOrder
        end
    end
end

if NVarACSpline > 0 % Cubic
    BoundsKnots = EstimOpt.Bounds(Dist == 6,:);
    Knots = zeros(NVarACSpline,NOrder+2);
    for i = 1:NVarACSpline
        Knots(i,:) = BoundsKnots(i,1) : (BoundsKnots(i,2) - BoundsKnots(i,1))/(NOrder+1) : BoundsKnots(i,2); % this is marginally different than train's
    end
    segment_idx = NaN(size(b_mtx(Dist == 6,:,1:NOrder+1)));
    for i = 1:NOrder+1
        segment_idx(:,:,i) = (GridMat(Dist == 6,:) >= Knots(:,i)) & (GridMat(Dist == 6,:) < Knots(:,i+1)); % NVarAStep x NRep*NP
    end
    b_mtx_tmp = zeros(size(b_mtx(Dist == 6,:,1:NOrder+1)));
    err_mtx_tmp = GridMat(Dist == 6,:);
    for i = 1:NVarACSpline
        b_mtx_tmp(i,segment_idx(i,:,1) == 1,1) = spline([2*Knots(i,1) - Knots(i,2),Knots(i,:),2*Knots(i,end) - Knots(i,end-1)],[0,1,zeros(1,size(Knots(i,:),2))],err_mtx_tmp(i,segment_idx(i,:,1) == 1));
        for j = 2:NOrder+1
            b_mtx_tmp(i,segment_idx(i,:,j-1) == 1,j) = spline([2*Knots(i,1) - Knots(i,2),Knots(i,:),2*Knots(i,end) - Knots(i,end-1)],[zeros(1,j),1,zeros(1,size(Knots(i,:),2)-j+1)],err_mtx_tmp(i,segment_idx(i,:,j-1) == 1));
            b_mtx_tmp(i,segment_idx(i,:,j) == 1,j) = spline([2*Knots(i,1) - Knots(i,2),Knots(i,:),2*Knots(i,end) - Knots(i,end-1)],[zeros(1,j),1,zeros(1,size(Knots(i,:),2)-j+1)],err_mtx_tmp(i,segment_idx(i,:,j) == 1));
        end
    end
    b_mtx(Dist == 6,:,1:NOrder+1) = b_mtx_tmp;
end

if NVarAPCSpline > 0 % Piece-wise Cubic
    BoundsKnots = EstimOpt.Bounds(Dist == 7,:);
    Knots = zeros(NVarAPCSpline,NOrder+2);
    for i = 1:NVarAPCSpline
        Knots(i,:) = BoundsKnots(i,1) : (BoundsKnots(i,2) - BoundsKnots(i,1))/(NOrder+1) : BoundsKnots(i,2); % this is marginally different than train's
    end
    segment_idx = NaN(size(b_mtx(Dist == 7,:,1:NOrder+1)));
    for i = 1:NOrder+1
        segment_idx(:,:,i) = (GridMat(Dist == 7,:) >= Knots(:,i)) & (GridMat(Dist == 7,:) < Knots(:,i+1)); % NVarAStep x NRep*NP
    end
    b_mtx_tmp = zeros(size(b_mtx(Dist == 7,:,1:NOrder+1)));
    err_mtx_tmp = GridMat(Dist == 7,:);
    for i = 1:NVarAPCSpline
        b_mtx_tmp(i,segment_idx(i,:,1) == 1,1) = spline([2*Knots(i,1) - Knots(i,2),Knots(i,1),Knots(i,2)],[0,1,0],err_mtx_tmp(i,segment_idx(i,:,1) == 1));
        for j = 2:NOrder+1
            b_mtx_tmp(i,segment_idx(i,:,j-1) == 1,j) = spline([Knots(i,j-1),Knots(i,j),Knots(i,j+1)],[0,1,0],err_mtx_tmp(i,segment_idx(i,:,j-1) == 1));
            b_mtx_tmp(i,segment_idx(i,:,j) == 1,j) = spline([Knots(i,j-1),Knots(i,j),Knots(i,j+1)],[0,1,0],err_mtx_tmp(i,segment_idx(i,:,j) == 1));
        end
        %         b_mtx_tmp(i,segment_idx(i,:,end) == 1,end) = spline([Knots(i,end-1),Knots(i,end),2*Knots(i,end) - Knots(i,end-1)],[0,1,0],err_mtx_tmp(i,segment_idx(i,:,end) == 1)); reference (0)
    end
    b_mtx(Dist == 7,:,1:NOrder+1) = b_mtx_tmp;
end

if NVarAPCHISpline > 0 % Piece-wise Cubic Hermite Interpolating (with bounds extension)
    BoundsKnots = EstimOpt.Bounds(Dist == 8,:);
    Knots = zeros(NVarAPCHISpline,NOrder+2);
    for i = 1:NVarAPCHISpline
        Knots(i,:) = BoundsKnots(i,1) : (BoundsKnots(i,2) - BoundsKnots(i,1))/(NOrder+1) : BoundsKnots(i,2); % this is marginally different than train's
    end
    segment_idx = NaN(size(b_mtx(Dist == 8,:,1:NOrder+1)));
    for i = 1:NOrder+1
        segment_idx(:,:,i) = (GridMat(Dist == 8,:) >= Knots(:,i)) & (GridMat(Dist == 8,:) < Knots(:,i+1)); % NVarAStep x NRep*NP
    end
    b_mtx_tmp = zeros(size(b_mtx(Dist == 8,:,1:NOrder+1)));
    err_mtx_tmp = GridMat(Dist == 8,:);
    for i = 1:NVarAPCHISpline
        b_mtx_tmp(i,segment_idx(i,:,1) == 1,1) = pchip([2*Knots(i,1) - Knots(i,2),Knots(i,:),2*Knots(i,end) - Knots(i,end-1)],[0,1,zeros(1,size(Knots(i,:),2))],err_mtx_tmp(i,segment_idx(i,:,1) == 1));
        for j = 2:NOrder+1
            b_mtx_tmp(i,segment_idx(i,:,j-1) == 1,j) = pchip([2*Knots(i,1) - Knots(i,2),Knots(i,:),2*Knots(i,end) - Knots(i,end-1)],[zeros(1,j),1,zeros(1,size(Knots(i,:),2)-j+1)],err_mtx_tmp(i,segment_idx(i,:,j-1) == 1));
            b_mtx_tmp(i,segment_idx(i,:,j) == 1,j) = pchip([2*Knots(i,1) - Knots(i,2),Knots(i,:),2*Knots(i,end) - Knots(i,end-1)],[zeros(1,j),1,zeros(1,size(Knots(i,:),2)-j+1)],err_mtx_tmp(i,segment_idx(i,:,j) == 1));
        end
    end
    b_mtx(Dist == 8,:,1:NOrder+1) = b_mtx_tmp;
end


% Cubic splines:
% if NVarACubicSpline > 0
%     BoundsKnots = EstimOpt.Bounds(Dist == 5,:);
%     Knots = zeros(NVarASpline,NOrder+2);
%     for i = 1:NVarASpline
%         Knots(i,:) = BoundsKnots(i,1) : (BoundsKnots(i,2) - BoundsKnots(i,1))/(NOrder+1) : BoundsKnots(i,2); % this is marginally different than train's
%     end
%     segment_idx = NaN(size(b_mtx(Dist == 5,:,1:NOrder+1)));
%     for i = 1:NOrder+1
%         segment_idx(:,:,i) = (err_mtx(Dist == 5,:) >= Knots(:,i)) & (err_mtx(Dist == 5,:) < Knots(:,i+1)); % NVarAStep x NRep*NP
%     end
%     b_mtx_tmp = zeros(size(b_mtx(Dist == 5,:,1:NOrder+1)));
%     err_mtx_tmp = err_mtx(Dist == 5,:);
%     for i = 1:NVarACubicSpline
%         %         b_mtx_tmp(i,segment_idx(i,:,1) == 1,1) = spline([2*Knots(i,1) - Knots(i,2),Knots(i,1),Knots(i,2)],[0,1,0],err_mtx_tmp(i,segment_idx(i,:,1) == 1));
%         %         b_mtx_tmp(i,segment_idx(i,:,1) == 1,1) = spline([2*Knots(i,1) - Knots(i,2),Knots(i,:),2*Knots(i,end) - Knots(i,end-1)],[0,1,zeros(1,size(Knots(i,:),2))],err_mtx_tmp(i,segment_idx(i,:,1) == 1));
%         b_mtx_tmp(i,segment_idx(i,:,1) == 1,1) = pchip([2*Knots(i,1) - Knots(i,2),Knots(i,:),2*Knots(i,end) - Knots(i,end-1)],[0,1,zeros(1,size(Knots(i,:),2))],err_mtx_tmp(i,segment_idx(i,:,1) == 1));
%         for j = 2:NOrder+1
%             %             b_mtx_tmp(i,segment_idx(i,:,j-1) == 1,j) = spline([Knots(i,j-1),Knots(i,j),Knots(i,j+1)],[0,1,0],err_mtx_tmp(i,segment_idx(i,:,j-1) == 1));
%             %             b_mtx_tmp(i,segment_idx(i,:,j-1) == 1,j) = spline([2*Knots(i,1) - Knots(i,2),Knots(i,:),2*Knots(i,end) - Knots(i,end-1)],[zeros(1,j),1,zeros(1,size(Knots(i,:),2)-j+1)],err_mtx_tmp(i,segment_idx(i,:,j-1) == 1));
%             b_mtx_tmp(i,segment_idx(i,:,j-1) == 1,j) = pchip([2*Knots(i,1) - Knots(i,2),Knots(i,:),2*Knots(i,end) - Knots(i,end-1)],[zeros(1,j),1,zeros(1,size(Knots(i,:),2)-j+1)],err_mtx_tmp(i,segment_idx(i,:,j-1) == 1));
%             %             b_mtx_tmp(i,segment_idx(i,:,j) == 1,j) = spline([Knots(i,j-1),Knots(i,j),Knots(i,j+1)],[0,1,0],err_mtx_tmp(i,segment_idx(i,:,j) == 1));
%             %             b_mtx_tmp(i,segment_idx(i,:,j) == 1,j) = spline([2*Knots(i,1) - Knots(i,2),Knots(i,:),2*Knots(i,end) - Knots(i,end-1)],[zeros(1,j),1,zeros(1,size(Knots(i,:),2)-j+1)],err_mtx_tmp(i,segment_idx(i,:,j) == 1));
%             b_mtx_tmp(i,segment_idx(i,:,j) == 1,j) = pchip([2*Knots(i,1) - Knots(i,2),Knots(i,:),2*Knots(i,end) - Knots(i,end-1)],[zeros(1,j),1,zeros(1,size(Knots(i,:),2)-j+1)],err_mtx_tmp(i,segment_idx(i,:,j) == 1));
%         end
%         %         b_mtx_tmp(i,segment_idx(i,:,end) == 1,end) = spline([Knots(i,end-1),Knots(i,end),2*Knots(i,end) - Knots(i,end-1)],[0,1,0],err_mtx_tmp(i,segment_idx(i,:,end) == 1)); reference (0)
%     end
% end
%     plot(b_mtx(2,:,2))
% plot(b_mtx1(5,:,2))
% plot(b_mtx2(5,:,2))
% plot(b_mtx3(5,:,2))
%     plot(b_mtx(5,:,3))
%     plot(b_mtx(5,:,4))
% end

b_mtx = reshape(permute(b_mtx,[1,3,2]),[size(b_mtx,1)*size(b_mtx,3),size(b_mtx,2)]); % NVarA * NOrder+1, NRep*NP
% b_mtx = reshape(permute(b_mtx,[3,1,2]),[size(b_mtx,1)*size(b_mtx,3),size(b_mtx,2)]); % train's Z is ordered by NOrder first and NV later, like this

% Correlations
if FullCov == 1
    if any(Dist == 3)
        GridMat = GridMat_old;
    end
    indx1 = tril(repmat(1:NVarA,[NVarA,1])',-1);
    indx1 = indx1(indx1~=0);
    indx2 = tril(repmat(1:NVarA,[NVarA,1]),-1);
    indx2 = indx2(indx2~=0);
    b_mtx = [b_mtx;GridMat(indx1,:).*GridMat(indx2,:)];
end

b_mtx(all(isnan(b_mtx),2),:) = [];