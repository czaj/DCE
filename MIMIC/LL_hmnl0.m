function f = LL_hmnl0(Xmea,EstimOpt,b)

l = 0;
L_mea = ones(EstimOpt.NP,1);
LogFact = 0;
% save tmp1
for i = 1:size(Xmea,2)
    if EstimOpt.MeaSpecMatrix(i) == 0 % OLS - might not work for MissingIndMea
        X = ones(EstimOpt.NP,1);
        bx = b(l+1:l+2);
        X_mea_n = Xmea(:,i);
        L = normpdf(X_mea_n,X*bx(1),exp(bx(2)));
        L_mea(EstimOpt.MissingIndMea(:,i) == 0,:) = L_mea(EstimOpt.MissingIndMea(:,i) == 0,:).*L;
        l = l + 2;
    elseif EstimOpt.MeaSpecMatrix(i) == 1 % MNL       
%         EstimOpt.NVarcut0 + length(unique(INPUT.Xmea(INPUT.MissingInd==0 & (EstimOpt.MissingIndMea(:,i) == 0),i))) - 1
        %         UniqueMea = unique(Xmea(:,i));
        UniqueMea = unique(Xmea(EstimOpt.MissingInd_tmp==0 & (EstimOpt.MissingIndMea(:,i) == 0),i));        
        k = length(UniqueMea) - 1;
        bx = [0,b(l+1:l+k)'];
        l = l + k;
        V = exp(ones(EstimOpt.NP,1)*bx); % NP x Unique values
        V = V./sum(V,2);
%         L = zeros(EstimOpt.NP,1);
        L = NaN(EstimOpt.NP,1);
        for j = 1:length(UniqueMea)
            L(Xmea(:,i) == UniqueMea(j)) = V(Xmea(:,i) == UniqueMea(j),j);
        end
%         L_mea = L_mea.*L;
        L_mea(~isnan(L)) = L_mea(~isnan(L)).*L(~isnan(L)); % this is a work-around for missing Ind values for MNL that do not matter anyway (because of INPUT.MissingInd)
    elseif EstimOpt.MeaSpecMatrix(i) == 2 % Ordered Probit
        UniqueMea = unique(Xmea(EstimOpt.MissingIndMea(:,i) == 0,i));
        k = length(UniqueMea) - 1;
%         save tmp1
        bx = b(l+1:l+k);        
        bx(2:end) = exp(bx(2:end));
        bx = cumsum(bx);
        L = zeros(sum(EstimOpt.MissingIndMea(:,i) == 0),1);
        L(Xmea(EstimOpt.MissingIndMea(:,i) == 0,i) == min(UniqueMea)) = normcdf(bx(1));
        L(Xmea(EstimOpt.MissingIndMea(:,i) == 0,i) == max(UniqueMea)) = 1 - normcdf(bx(end));
        for j = 2:k
            L(Xmea(EstimOpt.MissingIndMea(:,i) == 0,i) == UniqueMea(j)) = normcdf(bx(j)) - normcdf(bx(j-1));
        end
        L_mea(EstimOpt.MissingIndMea(:,i) == 0,:) = L_mea(EstimOpt.MissingIndMea(:,i) == 0,:).*L;
        l = l + k;
    elseif EstimOpt.MeaSpecMatrix(i) == 3 % Poisson
        bx = b(l+1);
        L = exp(bx*Xmea(:,i)-exp(bx));
        L_mea = L_mea.*L;
        LogFact = LogFact + gammaln(Xmea(:,i)+1);
        l = l + 1;
    elseif EstimOpt.MeaSpecMatrix(i) == 4 % NB
        bx = b(l+1);
        theta = exp(b(l+2));
        u = theta./(theta+exp(bx));
        L = min(gamma(theta+Xmea(:,i)),realmax)./(gamma(theta).*min(gamma(Xmea(:,i)+1),realmax));
        L = L.*(u.^theta).*((1-u).^Xmea(:,i));
        L_mea = L_mea.*L;
        l = l + 2;
    elseif EstimOpt.MeaSpecMatrix(i) == 5 % ZIP
        bx = b(l+1:l+2);
        p = exp(bx(1));
        p = p/(1+p);
        L = zeros(EstimOpt.NP,1);
        L(Xmea(:,i) == 0) = p + (1-p)*exp(-exp(bx(2)));
        L(Xmea(:,i) ~= 0) = (1-p)*exp(bx(2)*Xmea(Xmea(:,i) ~= 0,i) - exp(bx(2)))./min(gamma(Xmea(Xmea(:,i) ~= 0,i)+1),realmax);
        L_mea = L_mea.*L;
        l = l + 2;
    elseif EstimOpt.MeaSpecMatrix(i) == 6 % ZINB
        bx = b(l+1:l+2);
        theta = exp(b(l+3));
        p = exp(bx(1));
        p = p/(1+p);
        L = zeros(EstimOpt.NP,1);
        lam = exp(bx(2));
        u = theta./(theta+lam);
        L(Xmea(:,i) == 0) = p + (1-p)*u^theta;
        L(Xmea(:,i) ~= 0) = (1-p)*exp(gammaln(theta+Xmea(Xmea(:,i) ~= 0,i))-gammaln(Xmea(Xmea(:,i) ~= 0,i)+1)-gammaln(theta)).*(u.^theta).*((1-u).^Xmea(Xmea(:,i) ~= 0,i));
        L_mea = L_mea.*L;
        l = l + 3;
    end
end
f = -sum(log(L_mea)) + sum(LogFact);