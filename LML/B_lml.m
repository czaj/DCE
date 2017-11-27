function b_mtx = B_lml(err_mtx, EstimOpt)

NP = EstimOpt.NP;
NRep = EstimOpt.NRep;
NVarA = EstimOpt.NVarA;
Dist = EstimOpt.Dist;
Order = EstimOpt.Order;
FullCov = EstimOpt.FullCov;
Bounds = EstimOpt.Bounds;

NVarAApprox = sum(Dist == 0 | Dist == 1);
NVarAPoly = sum(Dist == 2 | Dist == 3);

if FullCov == 1 && any(Dist == 3)
    err_mtx_old = err_mtx;
end

b_mtx_1 = zeros([NVarA,NP*NRep]);
b_mtx_2 = zeros([NVarA,NP*NRep]);

% Approximate normal and lognormal:
if NVarAApprox > 0
    b_mtx_1(Dist == 0,:) = err_mtx(Dist == 0,:);
    b_mtx_1(Dist == 1,:) = log(err_mtx(Dist == 1,:));
    b_mtx_2(Dist == 0,:) = err_mtx(Dist == 0,:).^2;
    b_mtx_2(Dist == 1,:) = log(err_mtx(Dist == 1,:)).^2;
end

% Polynomials:
if NVarAPoly > 0
    err_mtx(Dist == 2,:) = (err_mtx(Dist == 2,:) - Bounds(Dist == 2,1))./(Bounds(Dist == 2,2) - Bounds(Dist == 2,1));
    % err_mtx(Dist == 2,:) = (err_mtx(Dist == 2,:) -min(err_mtx(Dist == 2,:),[],2))./(max(err_mtx(Dist == 2,:),[],2) - min(err_mtx(Dist == 2,:),[],2));
    err_mtx(Dist == 3,:) = log(err_mtx(Dist == 3,:));
    err_mtx(Dist == 3,:) = (err_mtx(Dist == 3,:) - Bounds(Dist == 3,1))./(Bounds(Dist == 3,2) - Bounds(Dist == 3,1));
    % err_mtx(Dist == 3,:) = (err_mtx(Dist == 3,:) -min(err_mtx(Dist == 3,:),[],2))./(max(err_mtx(Dist == 3,:),[],2) - min(err_mtx(Dist == 3,:),[],2));

    %     b_mtx_tmp = LegP(Order,err_mtx(Dist == 2 | Dist == 3,:));
    b_mtx_tmp = zeros([size(err_mtx(Dist == 2 | Dist == 3,:)),Order]);
    b_mtx_tmp(:,:,1) = err_mtx(Dist == 2 | Dist == 3,:);
    b_mtx_tmp(:,:,2) = ((2*2-1)/2)*err_mtx(Dist == 2 | Dist == 3,:).*b_mtx_tmp(:,:,1)-(2-1)/2;
    b_mtx_tmp(:,:,3) = ((2*3-1)/3)*err_mtx(Dist == 2 | Dist == 3,:).*b_mtx_tmp(:,:,2)-((3-1)/3)*b_mtx_tmp(:,:,1);
    if Order > 3
        for i = 1:(Order-3)
            Tmp = 3+i;
            b_mtx_tmp(:,:,Tmp) = ((2*Tmp-1)/Tmp)*err_mtx(Dist == 2 | Dist == 3,:).*b_mtx_tmp(:,:,Tmp-1)-((Tmp-1)/Tmp)*b_mtx_tmp(:,:,Tmp-2);
        end
    end
    b_mtx_1(Dist == 2 | Dist == 3,:) = b_mtx_tmp(:,:,1);
    b_mtx_2(Dist == 2 | Dist == 3,:) = b_mtx_tmp(:,:,2);
    %     b_mtx(Dist == 2,:) = legendreP(1,err_mtx(Dist == 2,:));
    %     b_mtx(Dist == 3,:) = legendreP(1,err_mtx(Dist == 3,:));
    %     b_mtx_2(Dist == 2,:) = legendreP(2,err_mtx(Dist == 2,:));
    %     b_mtx_2(Dist == 3,:) = legendreP(2,err_mtx(Dist == 3,:));
end

b_mtx = [b_mtx_1;b_mtx_2];
if (NVarAPoly > 0) && (Order > 2)
    b_mtx = [b_mtx;reshape(permute(b_mtx_tmp(:,:,3:end),[1,3,2]),[(Order - 2)*NVarAPoly,NRep*NP])];
end

if FullCov == 1
    if any(Dist == 3)
        err_mtx = err_mtx_old;
    end
    indx1 = tril(repmat(1:NVarA,[NVarA,1])',-1);
    indx1 = indx1(indx1~=0);
    indx2 = tril(repmat(1:NVarA,[NVarA,1]),-1);
    indx2 = indx2(indx2~=0);
    %     indx1 = [];
    %     indx2 = [];
    %     for i = 1:NVarA
    %         indx1 = [indx1, (i+1):NVarA];
    %         indx2 = [indx2, i*ones(1,NVarA-i)];
    %     end
    b_mtx = [b_mtx;err_mtx(indx1,:).*err_mtx(indx2,:)];
end