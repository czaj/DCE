function b_mtx = B_lml(err_mtx,EstimOpt)

save tmp_B_lml
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
NVarASpline = sum(Dist == 5);

if FullCov == 1 && any(Dist == 3)
    err_mtx_old = err_mtx;
end

b_mtx = NaN([NVarA,NP*NRep,NOrder+1]); % to be cut later

% Approximate normal and lognormal:
if NVarAApprox > 0
    b_mtx(Dist == 0,:,1) = err_mtx(Dist == 0,:);
    b_mtx(Dist == 1,:,1) = log(err_mtx(Dist == 1,:));
    b_mtx(Dist == 0,:,2) = err_mtx(Dist == 0,:).^2;
    b_mtx(Dist == 1,:,2) = log(err_mtx(Dist == 1,:)).^2;
end

% Polynomials:
if NVarAPoly > 0
    err_mtx(Dist == 2,:) = (err_mtx(Dist == 2,:) - Bounds(Dist == 2,1))./(Bounds(Dist == 2,2) - Bounds(Dist == 2,1));
    err_mtx(Dist == 3,:) = log(err_mtx(Dist == 3,:));
    err_mtx(Dist == 3,:) = (err_mtx(Dist == 3,:) - Bounds(Dist == 3,1))./(Bounds(Dist == 3,2) - Bounds(Dist == 3,1));
    b_mtx(Dist == 2 | Dist == 3,:,1) = err_mtx(Dist == 2 | Dist == 3,:);
    b_mtx(Dist == 2 | Dist == 3,:,2) = ((2*2-1)/2)*err_mtx(Dist == 2 | Dist == 3,:).*b_mtx(Dist == 2 | Dist == 3,:,1)-(2-1)/2;
    b_mtx(Dist == 2 | Dist == 3,:,3) = ((2*3-1)/3)*err_mtx(Dist == 2 | Dist == 3,:).*b_mtx(Dist == 2 | Dist == 3,:,2)-((3-1)/3)*b_mtx(Dist == 2 | Dist == 3,:,1);
    if NOrder > 3
        for i = 4:NOrder
            b_mtx(Dist == 2 | Dist == 3,:,i) = ((2*i-1)/i)*err_mtx(Dist == 2 | Dist == 3,:).*b_mtx(Dist == 2 | Dist == 3,:,i-1)-((i-1)/i)*b_mtx(Dist == 2 | Dist == 3,:,i-2);
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
        b_mtx(Dist == 4,:,i) = (err_mtx(Dist == 4,:) >= Thresholds(:,i)) & (err_mtx(Dist == 4,:) < Thresholds(:,i+1)); % NVarAStep x NRep*NP x NOrder
    end
end

% Splines:
if NVarASpline > 0
    BoundsKnots = EstimOpt.Bounds(Dist == 5,:);
    Knots = zeros(NVarASpline,NOrder+2);
    for i = 1:NVarASpline
        Knots(i,:) = BoundsKnots(i,1) : (BoundsKnots(i,2) - BoundsKnots(i,1))/(NOrder+1) : BoundsKnots(i,2); % this is marginally different than train's
    end
    b_mtx(Dist == 5,:,1:NOrder+1) = 0;
    for i = 1:NOrder+1 % the last level is the reference
        segment_idx = (err_mtx(Dist == 5,:) >= Knots(:,i)) & (err_mtx(Dist == 5,:) < Knots(:,i+1)); % NVarAStep x NRep*NP
        b_mtx(Dist == 5,:,i) = b_mtx(Dist == 5,:,i) + segment_idx.*(1 - (err_mtx(Dist == 5,:) - Knots(:,i))./(Knots(:,i+1)-Knots(:,i))); % NVarASpline x NRep*NP x NOrder         
        if i < NOrder+1
            b_mtx(Dist == 5,:,i+1) = segment_idx.*((err_mtx(Dist == 5,:) - Knots(:,i))./(Knots(:,i+1)-Knots(:,i))); % NVarASpline x NRep*NP x NOrder
        end
    end
end

b_mtx = reshape(permute(b_mtx,[1,3,2]),[size(b_mtx,1)*size(b_mtx,3),size(b_mtx,2)]); 
% b_mtx = reshape(permute(b_mtx,[3,1,2]),[size(b_mtx,1)*size(b_mtx,3),size(b_mtx,2)]); % train's Z is ordered by NOrder first and NV later, like this

% Correlations
if FullCov == 1
    if any(Dist == 3)
        err_mtx = err_mtx_old;
    end
    indx1 = tril(repmat(1:NVarA,[NVarA,1])',-1);
    indx1 = indx1(indx1~=0);
    indx2 = tril(repmat(1:NVarA,[NVarA,1]),-1);
    indx2 = indx2(indx2~=0);
    b_mtx = [b_mtx;err_mtx(indx1,:).*err_mtx(indx2,:)];
end

b_mtx(all(isnan(b_mtx),2),:) = [];