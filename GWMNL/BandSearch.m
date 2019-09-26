function AIC = BandSearch(Y, Xa,CrdsIN, Crds, CrdsMat ,NoReg, EstimOpt, OptimOpt,B0)

%% Calculation of weights

if EstimOpt.BandType == 1 % global kernel
    Weights = exp(-0.5*(CrdsMat/B0).^2);
    
else % spatially varying kernel
    CrdsVec2 = pdist(Crds); 
    CrdsMat2 = squareform(CrdsVec2);
    Rij = zeros(EstimOpt.NP, EstimOpt.NP); 
    for i = 1:NoReg
        
       [B,I] = sort(CrdsMat2(:,i));
       Rtmp = (0:NoReg-1)';
       %Rtmp = Rtmp(I);
       CrdsX = Crds(I,:);
       FindIndx = find(CrdsIN(:,1) == Crds(i,1) & CrdsIN(:,2) == Crds(i,2));
       for j =1:NoReg
           FindIndx2 = find(CrdsIN(:,1) == CrdsX(j,1) & CrdsIN(:,2) == CrdsX(j,2));
           Rij(FindIndx2, FindIndx) = Rtmp(j);
       end
    end
    
    if EstimOpt.BandType == 2
        Weights = exp(-Rij/B0);
    else
        Weights = exp(-(Rij.^0.5)/B0);
    end
end

HX = zeros(NoReg,1);
LLX = zeros(NoReg,1);
%Ytmp = reshape(Y, EstimOpt.NAlt*EstimOpt.NCT, EstimOpt.NP);
%XXa = reshape(Xa, EstimOpt.NAlt*EstimOpt.NCT, EstimOpt.NP, EstimOpt.NVarA);
Xtmp = reshape(Xa, EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP, EstimOpt.NVarA);
Xstar = Xtmp(2:end,:,:)- Xtmp(ones(EstimOpt.NAlt-1,1),:,:);
Xstar = reshape(Xstar, (EstimOpt.NAlt-1)*EstimOpt.NCT*EstimOpt.NP, EstimOpt.NVarA);

for n = 1:NoReg

    FindIndx = find(CrdsIN(:,1) == Crds(n,1) & CrdsIN(:,2) == Crds(n,2));
    % LocNo = length(FindIndx);
    
 
    if EstimOpt.WeightSum == 1
        Weights_n = Weights(:,FindIndx(1));
    elseif EstimOpt.WeightSum == 2
        Weights_n = Weights(:,FindIndx(1))/sum(Weights(:,FindIndx(1)),1);
    elseif EstimOpt.WeightSum == 3
        Weights_n = Weights(:,FindIndx(1))/sum(Weights(:,FindIndx(1)),1)*EstimOpt.NP;
    end
    %Weights_n = Weights(:,FindIndx(1))/sum(Weights(:,FindIndx(1)),1);
    
    LLfun = @(B) LL_gwmnl_MATlike(Y, Xa, Weights_n, EstimOpt,OptimOpt,B);
    if EstimOpt.StartMatFull == 1
        [bhat] = fminunc(LLfun, EstimOpt.b0(:,FindIndx(1)), OptimOpt);
    else
        [bhat] = fminunc(LLfun, EstimOpt.b0, OptimOpt);
        EstimOpt.b0 = bhat;
    end
    
    [LL,J, probs] = LL_gwmnl_bs(Y,Xa,EstimOpt,bhat);
    LL = reshape(LL, EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP);
    LL = LL(2:end,:);
    Weights_n = reshape(Weights_n(:, ones(EstimOpt.NCT*(EstimOpt.NAlt-1),1))',EstimOpt.NCT*(EstimOpt.NAlt-1)*EstimOpt.NP,1) ;
    XXstar = Xstar.*sqrt(Weights_n(:, ones(EstimOpt.NVarA,1)));
    
    % To doda³em
    VXl = zeros(EstimOpt.NVarA, (EstimOpt.NAlt-1)*EstimOpt.NCT*EstimOpt.NP);
    for j = 1:EstimOpt.NCT*EstimOpt.NP
        Vloc = -LL(:,j)*LL(:,j)'+diag(LL(:,j));
        VXl(:, (j-1)*(EstimOpt.NAlt-1)+1:j*(EstimOpt.NAlt-1)) = XXstar((j-1)*(EstimOpt.NAlt-1)+1:j*(EstimOpt.NAlt-1),:)'*Vloc;
    end



    
    Indx1 = [];
    Indx2 = [];
    for j = 1:length(FindIndx)
       Indx1 = [Indx1, (FindIndx(j)-1)*EstimOpt.NCT*(EstimOpt.NAlt-1)+1:FindIndx(j)*EstimOpt.NCT*(EstimOpt.NAlt-1)];
       Indx2 = [Indx2, (FindIndx(j)-1)*EstimOpt.NCT+1:FindIndx(j)*EstimOpt.NCT];
    end
    %size(Omega)
    %size(XXstar)
    %InvOm = Omega\(XXstar');
    Omega = VXl*XXstar;
    if EstimOpt.Clustered  == 0
        InvOm = Omega\VXl(:,Indx1);
        HX(n) = trace(XXstar(Indx1,:)*InvOm);
    else
        InvOm = inv(Omega);
        J = J.*Weights_n(1:EstimOpt.NAlt-1:end,ones(EstimOpt.NVarA,1));
        J = reshape(J, EstimOpt.NCT, EstimOpt.NP, EstimOpt.NVarA);
        J = squeeze(sum(J,1));
        U = J'*J;
        InvOmega = InvOm*U*InvOm;
        InvOm = InvOmega*VXl(:,Indx1);
        HX(n) = trace(XXstar(Indx1,:)*InvOm);
    end
    %HX(n) = trace(XXstar(Indx1,:)*InvOm*Vl(:,Indx1));
    
    LLX(n) = sum(probs(Indx2));
end
% LLX
% pause
AIC = -2*sum(LLX)/(EstimOpt.NCT*EstimOpt.NP)+(2*sum(HX)+1)/(EstimOpt.NCT*EstimOpt.NP-sum(HX)-2);
%AIC = -2*sum(LLX)/(EstimOpt.NCT*EstimOpt.NP);

%% 