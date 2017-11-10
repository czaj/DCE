function Z = LMXL_vars(err_mtx, EstimOpt)

NP = EstimOpt.NP;
NRep = EstimOpt.NRep;
NVarA = EstimOpt.NVarA;
Dist = EstimOpt.Dist;
Order = EstimOpt.Order;
FullCov = EstimOpt.FullCov;
Bounds = EstimOpt.Bounds;

Leg = sum(Dist == 2 | Dist == 3);

if FullCov == 1 && any(Dist == 3)
    err_mtx_old = err_mtx;
end
% err_mtx(Dist == 2,:) = (err_mtx(Dist == 2,:) -min(err_mtx(Dist == 2,:),[],2))./(max(err_mtx(Dist == 2,:),[],2) - min(err_mtx(Dist == 2,:),[],2));
err_mtx(Dist == 2,:) = (err_mtx(Dist == 2,:) -Bounds(Dist == 2,1))./(Bounds(Dist == 2,2) - Bounds(Dist == 2,1));
err_mtx(Dist == 3,:) = log(err_mtx(Dist == 3,:));
err_mtx(Dist == 3,:) = (err_mtx(Dist == 3,:) -Bounds(Dist == 3,1))./(Bounds(Dist == 3,2) - Bounds(Dist == 3,1));
%err_mtx(Dist == 3,:) = (err_mtx(Dist == 3,:) -min(err_mtx(Dist == 3,:),[],2))./(max(err_mtx(Dist == 3,:),[],2) - min(err_mtx(Dist == 3,:),[],2));

Z = zeros(NVarA, NRep*NP);
Z(Dist == 0,:) = err_mtx(Dist == 0,:);
Z(Dist == 1,:) = log(err_mtx(Dist == 1,:));

Z2 = zeros(NVarA, NRep*NP);
Z2(Dist == 0,:) = err_mtx(Dist == 0,:).^2;
Z2(Dist == 1,:) = log(err_mtx(Dist == 1,:)).^2;

if Leg > 0
    Z_tmp = LegP(Order, err_mtx(Dist == 2 | Dist == 3,:));
    Z(Dist == 2 | Dist == 3,:) = Z_tmp(:,:,1);
    Z2(Dist == 2 | Dist == 3,:) = Z_tmp(:,:,2);
%     Z(Dist == 2,:) = legendreP(1,err_mtx(Dist == 2,:));
%     Z(Dist == 3,:) = legendreP(1,err_mtx(Dist == 3,:));
%     Z2(Dist == 2,:) = legendreP(2,err_mtx(Dist == 2,:));
%     Z2(Dist == 3,:) = legendreP(2,err_mtx(Dist == 3,:));
end
Z = [Z;Z2];

if Leg > 0
    Z = [Z; reshape(permute(Z_tmp(:,:, 3:end), [1 3 2]), [(Order - 2)*sum(Dist == 2 | Dist == 3), NRep*NP])];
end

if FullCov == 1
    if any(Dist == 3)
        err_mtx = err_mtx_old;
    end
    indx1 = [];
    indx2 = [];
    for i = 1:NVarA
        indx1 = [indx1, (i+1):NVarA];
        indx2 = [indx2, i*ones(1,NVarA-i)];
    end
    
    Z = [Z; err_mtx(indx1,:).*err_mtx(indx2,:)];
end


