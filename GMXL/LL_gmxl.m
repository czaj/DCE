function f = LL_gmxl(YY,XXa,XXm,XXs,XXt,err,EstimOpt,b0)

% save input1
% return
NAlt = EstimOpt.NAlt; 
NCT = EstimOpt.NCT; 
NP = EstimOpt.NP; 
NRep = EstimOpt.NRep; 
NVarA = EstimOpt.NVarA; 
NVarM = EstimOpt.NVarM; 
NVarS = EstimOpt.NVarS; 
NVarT = EstimOpt.NVarT; 
Dist = EstimOpt.Dist; 
WTP_space = EstimOpt.WTP_space; 
WTP_matrix = EstimOpt.WTP_matrix; 

b0a = b0(1:NVarA);

if EstimOpt.FullCov == 0
    b0v = (b0(NVarA+1:NVarA*2)).^2; 
%     b0v = b0(NVarA+1:NVarA*2);
    VC = diag(b0v);
    b0m = b0(NVarA*2+1:NVarA*(NVarM+2));
    b0m = reshape(b0m,NVarA, NVarM);
    b0s = b0(NVarA*(NVarM+2)+1:NVarA*(NVarM+2)+NVarS);
    b0t = b0(NVarA*(NVarM+2)+NVarS+1:NVarA*(NVarM+2)+NVarS+NVarT);
    tau = (b0(NVarA*(NVarM+2)+NVarS+NVarT+1)).^2;
%     tau = exp(b0(NVarA*(NVarM+2)+NVarS+NVarT+1));
    if isfield(EstimOpt,'gamma0') == 1 && (EstimOpt.gamma0 == 0 || EstimOpt.gamma0 == 1)
        gamma = EstimOpt.gamma0;
    else
        gamma = exp(b0(NVarA*(NVarM+2)+NVarS+NVarT+2)) ./ (1 + exp(b0(NVarA*(NVarM+2)+NVarS+NVarT+2)));
    end
else
	b0v = b0(NVarA+1:NVarA+sum(1:NVarA));
    VC = tril(ones(NVarA));
    VC(VC==1) = b0v;
    b0m = b0(NVarA*(NVarA/2+1.5)+1:NVarA*(NVarA/2+1.5+NVarM));
    b0m = reshape(b0m,NVarA, NVarM);
    b0s = b0(NVarA*(NVarA/2+1.5+NVarM)+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS);
    b0t = b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarT);   
    tau = exp(b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarT+1));
    if isfield(EstimOpt,'gamma0') == 1 && (EstimOpt.gamma0 == 0 || EstimOpt.gamma0 == 1)
        gamma = EstimOpt.gamma0;
    else
        gamma = exp(b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarT+2)) ./ (1 + exp(b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarT+2)));
    end
end
b0n = b0a(:,ones(NP,1)) + b0m*XXm;
b0n = reshape(b0n([1:size(b0n,1)]'*ones(1,NRep),[1:size(b0n,2)]'),NVarA,NRep*NP);
if size(XXt,2) == 0
    if EstimOpt.Dist(1) == 1
        errx = reshape(err(1,:), NRep, NP)'; % NP x NRep
        sigma_bar = - log(mean(exp(tau*errx),1)); %1 x NRep
        %sigma_bar = -log(mean(exp(tau*(squeeze(err_sliced(1,:,:))')),1));
        sigma = exp(sigma_bar(ones(NP,1),:) + tau* errx); % NP x NRep
	elseif EstimOpt.Dist(1) == 5
        shape_n = 1./tau; 
        errx = reshape(err(1,:), NRep, NP)';
        sigma_bar = mean(errx.^shape_n,1); %1 x NRep
        sigma = (errx.^shape_n)./ sigma_bar(ones(NP,1),:);
        %sigma_bar = mean(squeeze(err_sliced(1,:,:)).^shape_n(ones(EstimOpt.NRep,1),ones(1,EstimOpt.NP)),2)';
    end    
else    
    if EstimOpt.Dist(1) == 1
        cov_tau = exp(XXt*b0t);  % NP x 1
        cov_tau = cov_tau(:,ones(1,NRep));
        errx = reshape(err(1,:), NRep, NP)';
        sigma_bar = -log(mean(exp(tau*cov_tau.*errx),1)); % 1 x NRep
        %sigma_bar = -log(mean(exp(tau*cov_tau(:,ones(1,EstimOpt.NRep)).*(squeeze(err_sliced(1,:,:))')),1));
        sigma = exp(sigma_bar(ones(NP,1),:) + tau*cov_tau.*errx);
    elseif EstimOpt.Dist(1) == 5
        shape_n = 1./(tau*exp(XXt*b0t));  % NP x 1
        shape_n = shape_n(:,  ones(1, NRep));
        errx = reshape(err(1,:), NRep, NP)'; 
        sigma_bar = mean(errx.^shape_n,1); % 1 x NRep
        %sigma_bar = mean(squeeze(err_sliced(1,:,:)).^shape_n(ones(EstimOpt.NRep,1),:),2)';
        sigma = (errx.^shape_n)./ sigma_bar(ones(NP,1),:); % NP x NRep
    end
end

sigma = reshape(sigma', 1, NRep*NP);

if sum(Dist(2:end) > 0) == 0;  % Normal
    if gamma == 0
        b_mtx = sigma(ones(NVarA,1),:).*(b0n + VC*err(2:end,:));  % NVarA x NRep*NP
    elseif gamma == 1
        b_mtx = sigma(ones(NVarA,1),:).*b0n + VC*err(2:end,:);
    else
        b_mtx = sigma(ones(NVarA,1),:).*(b0n + (1-gamma)*VC*err(2:end,:)) + gamma*VC*err(2:end,:);
    end
elseif sum(Dist(2:end)==1) > 0;  % Log - normal
    if gamma == 0
        b_mtx  = b0n + VC*err(2:end,:); 
        b_mtx(Dist(2:end)==1,:) = exp(b_mtx(Dist(2:end)==1,:)); 
        b_mtx = sigma(ones(NVarA,1),:).*b_mtx;  % NVarA x NRep*NP
    elseif gamma == 1
        b_mtx = sigma(ones(NVarA,1),:).*b0n + VC*err(2:end,:);
        b_mtx(Dist(2:end)==1,:) = exp(b_mtx(Dist(2:end)==1,:)); 
    else
        b_mtx = sigma(ones(NVarA,1),:).*(b0n + (1-gamma)*VC*err(2:end,:)) + gamma*VC*err(2:end,:);
        b_mtx(Dist(2:end)==1,:) = exp(b_mtx(Dist(2:end)==1,:)); 
    end
    
elseif sum(Dist(2:end)==2) > 0;  % Spike       
    if gamma == 0
        b_mtx  = b0n + VC*err(2:end,:); 
        b_mtx(Dist(2:end)==1,:) = max(b_mtx(Dist(2:end)==1,:),0); 
        b_mtx = sigma(ones(NVarA,1),:).*b_mtx;  % NVarA x NRep*NP
    elseif gamma == 1
        b_mtx = sigma(ones(NVarA,1),:).*b0n + VC*err(2:end,:);
        b_mtx(Dist(2:end)==1,:) = max(b_mtx(Dist(2:end)==1,:),0); 
    else
        b_mtx = sigma(ones(NVarA,1),:).*(b0n + (1-gamma)*VC*err(2:end,:)) + gamma*VC*err(2:end,:);
        b_mtx(Dist(2:end)==1,:) = max(b_mtx(Dist(2:end)==1,:),0); 
    end
elseif sum(Dist(2:end)==5) > 0; 
    if sum(sum(VC.*(1-eye(size(b0n,1)))~=0))~=0; error ('Weibull distribution can only be used with non-correlated parameters'); end; 
    
    if gma ~= 0
        error ('Weibull distributed attriute parameters possible only with G-MNL Type II')
    else
        b_mtx = zeros(NVarA,NP*NRep); 
        err2 = err(2:end,:);
        b_mtx(Dist(2:end)==0,:) = b0n(Dist(2:end)==0,:) + VC(Dist(2:end)==0,Dist(2:end)==0)*err2(Dist(2:end)==0,:); 
        Wexp = (1./diag(VC(Dist(2:end)==5,Dist(2:end)==5)));
        b_mtx(Dist(2:end)==5,:) = b0n(Dist(2:end)==5).*err2(Dist(2:end)==5,:).^Wexp(:,ones(1,NRep)); 
        b_mtx = b_mtx.*sigma(ones(NVarA,1),:); 
    end
    
end; 

if WTP_space > 0
    b_mtx(1:end-WTP_space,:) = b_mtx(1:end-WTP_space,:).*b_mtx(WTP_matrix,:); 
end

cs = reshape(exp(XXs*b0s),NAlt*NCT,1,NP); 
XXa_n = XXa .* cs(:,ones(1,NVarA,1),:) ; 

p0 = zeros(NP,1);

%for n = 1:NP
%    U = reshape(exp(XXa_n(:,:,n)*b_mtx(:,((n-1)*NRep+1):n*NRep)),NAlt,NCT,NRep); 
%    U(isnan(U)) = 0;  % skip missing ALT
%    U_sum = sum(U,1); 
%    YY_n = YY(:,:,n); 
%    P = sum(YY_n(:,:,ones(1,1,NRep)).*U ./ U_sum(ones(NAlt,1),:,:),1); 
%    P(isnan(P)) = 1;  % skip missing NCT, alternative: A2(1,sum(MissingInd2(:,:,n),1)==NAlt,:) = 1; 
%    p0(n) = mean(prod(P,2),3); 

%end

if any(isnan(XXa(:))) == 0  % faster version for complete dataset
    for n = 1:NP 
        U = reshape(XXa_n(:,:,n)*b_mtx(:,((n-1)*NRep+1):n*NRep),NAlt,NCT,NRep); 
        U_max = max(U); 
        U = exp(U - U_max(ones(NAlt,1),:,:));  % rescale utility to avoid exploding
        U_sum = reshape(sum(U,1),NCT,NRep); 
        U_selected = reshape(U(YY(:,n*ones(NRep,1))==1),NCT,NRep);    
        p0(n) = mean(prod(U_selected ./ U_sum,1),2);
    end; 
else  % this works only if NAlt is constant for each respondent and if missing ALT is not the first in NCT
    for n = 1:NP 
%         U = reshape(XXa_n(~isnan(YY(:,n)),:,n)*b_mtx(:,((n-1)*NRep+1):n*NRep),NAlt,NCT-sum(isnan(YY(1:NAlt:end,n))),NRep); 
% from MXL: U = reshape(XXa_n(~isnan(YY(:,n)),:,n)*b_mtx(:,((n-1)*NRep+1):n*NRep),NAlt,NCT-sum(isnan(YY(1:NAlt:end,n))),NRep); % this would be faster if there are no ALT missing
        U = reshape(XXa_n(~isnan(YY(:,n)),:,n)*b_mtx(:,((n-1)*NRep+1):n*NRep),numel(YY(~isnan(YY(:,n)),n))./(NCT-sum(isnan(YY(1:NAlt:end,n)))),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep);         
        U_max = max(U); 
%         U = exp(U - U_max(ones(NAlt,1),:,:)); 
	    U = exp(U - U_max(ones(numel(YY(~isnan(YY(:,n)),n))./(NCT-sum(isnan(YY(1:NAlt:end,n)))),1),:,:)); 
        U_sum = reshape(nansum(U,1),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep); 
        U_selected = reshape(U(YY(~isnan(YY(:,n)),n*ones(NRep,1))==1),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep);    
        p0(n) = mean(prod(U_selected ./ U_sum,1),2);
    end; 
end 

f = -log(p0);
% f = -log(max(p0,realmin));

