function [f,g,h] = LL_mxl(YY,XXa,XXm,Xs,err,EstimOpt,b0)
% f is a vector of -log(P_ni), where P_ni is a simulated choice probability
% P for person 'n' and alternative 'i'
% It is summed in the calling function to get the SLL
%
% g is a matrix of size [NP,2*NVarA+NVarNLT+NVarS]
% It consists of gradients for each person.
% Each row is calculated as \frac{\partial ln(P*), \partial \THETA}
% where \THETA is a vector of parameters for each alternative.
% It is summed (in NP dimension) in the calling function to get the grad(SLL).
%
% err -- matrix of random draws from standard normal distribution
%
%     YY = gpuArray(YY);
%     XXa = gpuArray(XXa);
%     XXm = gpuArray(XXm);
%     Xs = gpuArray(Xs);
%     err = gpuArray(err);
%     b0 = gpuArray(b0);


% save LL_mxl
% return

% issues:
% - NLT, Johnson, nargout == 3 - still use squeeze, ones, bsxfun

NAlt = EstimOpt.NAlt;
NCT = EstimOpt.NCT;
NP = EstimOpt.NP;
NRep = EstimOpt.NRep;
NVarA = EstimOpt.NVarA;
NVarM = EstimOpt.NVarM;
NVarS = EstimOpt.NVarS;
Dist = EstimOpt.Dist;
WTP_space = EstimOpt.WTP_space;
WTP_matrix = EstimOpt.WTP_matrix;
FullCov = EstimOpt.FullCov;
DiagIndex = EstimOpt.DiagIndex;
Triang = EstimOpt.Triang;
NVarNLT = EstimOpt.NVarNLT;
NLTVariables = EstimOpt.NLTVariables;
NLTType = EstimOpt.NLTType;
Johnson = EstimOpt.Johnson;
NCTMiss = EstimOpt.NCTMiss;
NAltMiss = EstimOpt.NAltMiss;
NAltMissInd = EstimOpt.NAltMissInd;
NAltMissIndExp = EstimOpt.NAltMissIndExp;
%NCTMissIndExp = EstimOpt.NCTMissIndExp;
MissingCT = EstimOpt.MissingCT;
RealMin = EstimOpt.RealMin;
ExpB = EstimOpt.ExpB;

% if nargout == 3
%    XXX = permute(mmx('square',permute(XXa,[2,4,1,3]),[]),[3,1,2,4]);
% %    VCx = EstimOpt.VCx;
% end
if FullCov == 1 && nargout > 1
    indx1 = EstimOpt.indx1;
    indx2 = EstimOpt.indx2;
else
    indx1 = [];
    indx2 = [];
end

if ~isempty(ExpB)
    b0(ExpB) = exp(b0(ExpB));
end

%% First parameters of the distributions
b0a = b0(1:NVarA);  % first parameters of the distributions

% Additional variables for distributions other than normal
% Introduction of new variables requires zeroing the corresponding elements of the parameters vector
if any(Dist == 3)   % triangular distribution
    b0triang_c = exp(b0a(Dist == 3)) + Triang'; % exp of parameters' values
    b0a(Dist == 3) = 0; % Set the first parameter to mu of standard normal (0)
end
if any(Dist == 4)   % weibull distribution
    b0weibA = exp(b0a(Dist == 4));
    b0a(Dist == 4) = 0;
end
if any(Dist == 5)   % Sinh-Arcsinh distribution
    b0sinhA = b0a(Dist == 5);
    b0a(Dist == 5) = 0;
end

%% Second parameters of the distributions (std devs as default)
if FullCov == 0 % no covariance, only variances
    b0v = (b0(NVarA+1:NVarA*2));    % Second parameters of the distributions (std devs as default)
    if any(Dist == 3)   % tiangular
        b0triag_b = exp(b0v(Dist == 3)) + b0triang_c;
        b0v(Dist == 3) = 1; % Set the second parameter to sigma of standard normal (1)
    end
    if any(Dist == 4)   % weibull
        b0weibB = exp(-b0v(Dist == 4)); % b0weibB = 1/exp(b0v(Dist == 4))
        b0v(Dist == 4) = 1;
    else
        b0weibB = []; % initialize variable for parfor loop
    end
    if any(Dist == 5)   % sinh-arcsinh
        b0sinhB = b0v(Dist == 5).^2;
        b0v(Dist == 5) = 1;
    end
    %     b0v = b0v.^2;
    VC = diag(b0v); % no covariance elements (varcov matrix; here var only)
    b0m = b0(NVarA*2+1:NVarA*(NVarM+2));
    b0m = reshape(b0m,[NVarA,NVarM]);   % means
    b0s = b0(NVarA*(NVarM+2)+1:NVarA*(NVarM+2)+NVarS);  % scale
    b0t = b0(NVarA*(NVarM+2)+NVarS+1:NVarA*(NVarM+2)+NVarS+NVarNLT);    % nonlinear transformation
    b0j = b0(NVarA*(NVarM+2)+NVarS+NVarNLT+1:NVarA*(NVarM+2)+NVarS+NVarNLT+2*Johnson);  % Johnson distribution additional parameters
else    % with covariance elements
    b0v = b0(NVarA+1:NVarA+sum(1:NVarA));   % second parameters of the distributions (with covariance)
    tmp = b0v(DiagIndex);   % second parameters of the distributions (std devs as default)
    b0v(DiagIndex(Dist >=3 & Dist <=5)) = 1;    % triangular, weibull and sinh-arcsinh distributions
    if any(Dist == 3)
        b0triag_b = exp(tmp(Dist == 3)) + b0triang_c;
    end
    if any(Dist == 4)
        b0weibB = exp(-tmp(Dist == 4));
    else
        b0weibB = [];
    end
    if any(Dist == 5)
        b0sinhB = tmp(Dist == 5).^2;
    end
    VC = tril(ones(NVarA)); % with covariance elements
    VC(VC == 1) = b0v;
    if any(Dist >= 3 & Dist <= 5)
        tmp = sqrt(sum(VC(Dist >= 3 & Dist <= 5,:).^2,2));
        VC(Dist >= 3 & Dist <= 5,:) = VC(Dist >= 3 & Dist <= 5,:)./tmp;
    end
    b0m = b0(NVarA*(NVarA/2+1.5)+1:NVarA*(NVarA/2+1.5+NVarM));
    b0m = reshape(b0m,[NVarA,NVarM]);
    b0s = b0(NVarA*(NVarA/2+1.5+NVarM)+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS);
    b0t = b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarNLT);
    b0j = b0(NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarNLT+1:NVarA*(NVarA/2+1.5+NVarM)+NVarS+NVarNLT+2*Johnson);
end

%% Nonlinear transformations
if NVarNLT > 0
    % IndTransNon0 = (abs(bt) > 0.00001)';
    IndTransNon0 = (abs(b0t) > eps)';
    Xt = XXa(:,NLTVariables,:);	% nonlinearly transformed variables
    %     bt_tmp = permute(b0t(:, ones(NAlt*NCT,1), ones(NP,1)), [2 1 3]);
    bt_tmp = b0t(:,ones(NAlt*NCT,1))';
    bt_tmp = bt_tmp(:,:,ones(1,1,NP));
    
    % Transform variables with the chosen transformation type (Box-Cox, Yeo-Johnson)
    if NLTType == 1 % BC
        Xt(:,IndTransNon0,:) = -(Xt(:,IndTransNon0,:).^bt_tmp(:,IndTransNon0,:) - 1)./bt_tmp(:,IndTransNon0,:);
        Xt(:,~IndTransNon0,:) = -log(Xt(:,~IndTransNon0,:));
    elseif NLTType == 2 % YJ
        IndXtNon0 = (Xt >= 0);
        IndXtCase1 = IndXtNon0 & IndTransNon0; % X >= 0, lam ~= 0
        IndXtCase2 = IndXtNon0 & ~IndTransNon0; % X >= 0, lam = 0
        %     IndTransNon2 = (abs(bt - 2) < 0.00001)';
        IndTransNon2 = (abs(b0t - 2) > eps)';
        IndXtCase3 = ~IndXtNon0 & IndTransNon2;  % X < 0, lam ~= 2
        IndXtCase4 = ~IndXtNon0 & ~IndTransNon2; % X < 0, lam = 2
        %         bt_tmp = b0t(:,ones(size(XXa,1),1))';
        Xt(IndXtCase1) = ((Xt(IndXtCase1) + 1).^bt_tmp(IndXtCase1) - 1)./bt_tmp(IndXtCase1);
        Xt(IndXtCase2) = log(Xt(IndXtCase2) + 1);
        Xt(IndXtCase3) = -((-Xt(IndXtCase3) + 1).^(2 - bt_tmp(IndXtCase3)) - 1)./(2 - bt_tmp(IndXtCase3));
        Xt(IndXtCase4) = -log(-Xt(IndXtCase4) + 1);
    end
    
    if EstimOpt.NumGrad == 0 %
        if NLTType == 1 % BC
            XXt = XXa(:,NLTVariables,:);
            XXt(:,IndTransNon0,:) = -(XXt(:,IndTransNon0,:).^bt_tmp(:,IndTransNon0,:).*(bt_tmp(:,IndTransNon0,:).*log(XXt(:,IndTransNon0,:))-1)+1)./(bt_tmp(:,IndTransNon0,:).^2);
            XXt(:,IndTransNon0 == 0,:) = -0.5*log(XXt(:,IndTransNon0 == 0)).^2;
        elseif NLTType == 2 % YJ
            XXt = XXa(:,NLTVariables,:);
            XXt(IndXtCase1) = ((XXt(IndXtCase1)+1).^bt_tmp(IndXtCase1).*(bt_tmp(IndXtCase1).*log(XXt(IndXtCase1)+1)-1)+1)./(bt_tmp(IndXtCase1).^2);% X >= 0, lam ~= 0
            XXt(IndXtCase2) = 0.5*log(XXt(IndXtCase2)+1).^2;% X >= 0, lam == 0
            XXt(IndXtCase3) = -((-XXt(IndXtCase3)+1).^(2-bt_tmp(IndXtCase3)).*(1-(2-bt_tmp(IndXtCase3)).*log(-XXt(IndXtCase3)+1))-1)./((2-bt_tmp(IndXtCase3)).^2);% X < 0, lam ~= 2
            XXt(IndXtCase4) = -0.5*log(-XXt(IndXtCase4)+1).^2;% X < 0, lam == 2
        end
    end
    XXa(:,NLTVariables,:) = Xt;
else	% No nonlinear transformation
    %     if EstimOpt.NumGrad == 0
    XXt = zeros(0,0,NP);
    %     end
end

%% Transformation from standard normal to normal(mu, sig2) draws
b0n = b0a + b0m*XXm;
b0n = reshape(b0n((1:size(b0n,1))'*ones(1,NRep),(1:size(b0n,2))'),[NVarA,NRep*NP]);  % NVarA x NRep*NP
b_mtx = b0n + VC*err;  % NVarA x NRep*NP;	transformation line
% draws for distributions 3, 4 and 5 remain standard normal in b_mtx
% and are adjusted below

%% Transformations for other distributions
if sum(Dist == 1) > 0 % Log-normal
    b_mtx(Dist == 1,:) = exp(b_mtx(Dist == 1,:));
end
if sum(Dist == 2) > 0 % Spike
    b_mtx(Dist == 2,:) = max(b_mtx(Dist == 2,:),0);
end
if sum(Dist == 3) > 0 % Triangular
    tmp = normcdf(b_mtx(Dist == 3,:));
    Triang = Triang(ones(NRep*NP,1),:)';
    Ftriang = (b0triang_c - Triang)./(b0triag_b - Triang);
    bmtx_triang = zeros(size(tmp));
    tmp2 = (b0triag_b - Triang).*(b0triang_c - Triang);
    bmtx_triang(tmp < Ftriang) = Triang(tmp < Ftriang)+ sqrt(tmp(tmp < Ftriang).*tmp2(tmp < Ftriang));
    tmp2 = (b0triag_b - Triang).*(b0triag_b-b0triang_c);
    %bmtx_triang(tmp >= Ftriang) = b0triag_b(tmp >= Ftriang)- sqrt((1-tmp(tmp >= Ftriang)).*tmp2(tmp >= Ftriang));
    bmtx_triang(tmp >= Ftriang) = b0triag_b- sqrt((1-tmp(tmp >= Ftriang)).*tmp2(tmp >= Ftriang));
    b_mtx(Dist == 3,:) = bmtx_triang;
end
if sum(Dist == 4) > 0 % Weibull
    tmpWeib = -log(1-normcdf(b_mtx(Dist == 4,:)));
    b_mtx(Dist == 4,:) = b0weibA.*(tmpWeib.^b0weibB);   % inverse CDF function
    
    if EstimOpt.NumGrad == 0
        tmpWeib = reshape(tmpWeib, [], NRep, NP);
    end
else
    tmpWeib = double.empty(0, 0, NP); % initialize variable for parfor loop (even thoough it is not used)
end
if sum(Dist >= 5) > 0 % Johnson
    if sum(Dist == 5) > 0 % Sinh-Arcsinh
        b_mtx(Dist == 5,:) = b0sinhA + b0sinhB.*asinh(b_mtx(Dist == 5,:));
        b_mtx(Dist == 5,:) = b0j(1:Johnson,ones(NRep*NP,1)) + exp(b0j(Johnson+1:end,ones(NRep*NP,1))).*sinh(b_mtx(Dist == 5,:));
        b_mtx(Dist == 5,:) = b0j(1:Johnson,:) + exp(b0j(Johnson+1:end,:)).*sinh(b_mtx(Dist == 5,:));
    end
    if sum(Dist == 6) > 0 % Johnson Sb
        tmp = exp(b_mtx(Dist == 6,:));
        b_mtx(Dist == 6,:) = tmp./(1+tmp);
        b_mtx(Dist == 6,:) = b0j(1:Johnson,:) + exp(b0j(Johnson+1:end,:)).*b_mtx(Dist == 6,:);
    end
    if sum(Dist == 7) > 0 % Johnson Su
        b_mtx(Dist == 7,:) = sinh(b_mtx(Dist == 7,:));
        b_mtx(Dist == 7,:) = b0j(1:Johnson,:) + exp(b0j(Johnson+1:end,:)).*b_mtx(Dist == 7,:);
    end
end

if WTP_space > 0
    b_mtx_grad = reshape(b_mtx,[NVarA,NRep,NP]); % needed for gradient calculation in WTP_space
    b_mtx(1:end-WTP_space,:) = b_mtx(1:end-WTP_space,:).*b_mtx(WTP_matrix,:);
else
    b_mtx_grad = zeros(0,0,NP);
end

if NVarS > 0
    cs = reshape(exp(Xs*b0s),[NAlt*NCT,1,NP]);
    XXa = XXa .* cs;
end

b_mtx = reshape(b_mtx,[NVarA,NRep,NP]);
p0 = zeros(NP,1); % simulated choice probabilities
% p0 = zeros([NP,1],'gpuArray');


%% main loop

if nargout == 1 % function value only
    
    if any(isnan(XXa(:))) == 0 % faster version for complete dataset
        % chosen alternative indicator
        YYy = (YY == 1);
        % calculation of the chosen alternative probability estimator for
        % each individual
        parfor n = 1:NP   % for each person
            % utility function
            U = reshape(XXa(:,:,n)*b_mtx(:,:,n),[NAlt,NCT,NRep]);
            U = exp(U - max(U,[],1)); % rescale utility to avoid exploding
            % denominator term of the conditional probability fraction
            % (denominator of the logit probability (sum over alternatives))
            U_sum = reshape(sum(U,1),[NCT,NRep]);
            YYy_n = YYy(:,n);
            % numerator term of the conditional probability fraction
            U_selected = reshape(U(YYy_n(:,ones(NRep,1))),[NCT,NRep]);
            % calculate probability estimator for the chosen alternative
            % prod for panel data, mean for simulated probability
            p0(n) = mean(prod(U_selected./U_sum,1),2);
        end
    else
        parfor n = 1:NP
            YnanInd = ~isnan(YY(:,n));
            XXa_n = XXa(:,:,n);
            NAltMissIndExp_n = NAltMissIndExp(:,n);
            NAltMissIndExp_n = NAltMissIndExp_n(YnanInd);
            if var(NAltMissIndExp_n(NAltMissIndExp_n > 0)) == 0 % if NAlt is constant per individual (but can vary between individuals)
                U = reshape(XXa_n(YnanInd,:)*b_mtx(:,:,n),[NAltMiss(n),NCTMiss(n),NRep]);
                U = exp(U - max(U,[],1));
                U_sum = reshape(sum(U,1),[NCTMiss(n),NRep]);
            else
                NAltMissInd_n = NAltMissInd(:,n);
                NAltMissInd_n = NAltMissInd_n(MissingCT(:,n) == 0);
                U = XXa_n(YnanInd,:)*b_mtx(:,:,n);
                Uniq = unique(NAltMissIndExp_n);
                U_sum = zeros(NCTMiss(n),NRep);
                if length(Uniq) == 2
                    U_tmp = U(NAltMissIndExp_n == Uniq(1),:);
                    U_tmp = reshape(U_tmp,[Uniq(1),size(U_tmp,1)/Uniq(1),NRep]);
                    U_tmp = exp(U_tmp - max(U_tmp));
                    U_sum(NAltMissInd_n == Uniq(1),:) = reshape(sum(U_tmp,1),[size(U_tmp,2),NRep]);
                    U(NAltMissIndExp_n == Uniq(1),:) = reshape(U_tmp,[size(U_tmp,2)*Uniq(1),NRep]);
                    U_tmp = U(NAltMissIndExp_n == Uniq(2),:);
                    U_tmp = reshape(U_tmp,[Uniq(2),size(U_tmp,1)/Uniq(2),NRep]);
                    U_tmp = exp(U_tmp - max(U_tmp));
                    U_sum(NAltMissInd_n == Uniq(2),:) = reshape(sum(U_tmp,1),[size(U_tmp,2),NRep]);
                    U(NAltMissIndExp_n == Uniq(2),:) = reshape(U_tmp,[size(U_tmp,2)*Uniq(2),NRep]);
                else
                    for i = 1:length(Uniq)
                        U_tmp = U(NAltMissIndExp_n == Uniq(i),:);
                        U_tmp = reshape(U_tmp,[Uniq(i),size(U_tmp,1)/Uniq(i),NRep]);
                        U_tmp = exp(U_tmp - max(U_tmp));
                        U_sum(NAltMissInd_n == Uniq(i),:) = reshape(sum(U_tmp,1),[size(U_tmp,2),NRep]);
                        U(NAltMissIndExp_n == Uniq(i),:) = reshape(U_tmp,[size(U_tmp,2)*Uniq(i),NRep]);
                    end
                end
            end
            YYy_n = YY(:,n)==1;
            U_selected = reshape(U(YYy_n(YnanInd,ones(NRep,1))),[NCTMiss(n),NRep]);
            p0(n) = mean(prod(U_selected./U_sum,1));
        end
    end
    
elseif nargout == 2 %% function value + gradient
   
    if NVarS > 0
        Xs_sliced = reshape(Xs,[NAlt*NCT,NP,NVarS]);
        Xs_sliced = permute(Xs_sliced(1:NAlt:end,:,:),[1,3,2]); % NCT x NVarS x NP
    else
        Xs_sliced = reshape(Xs,[NAlt*NCT,NP,0]);
        Xs_sliced = permute(Xs_sliced,[1,3,2]);
    end
    if FullCov == 0
%         g = zeros([NP,2*NVarA+NVarNLT+NVarS],'gpuArray');
        g = zeros(NP,2*NVarA+NVarNLT+NVarS);
        VC2 = reshape(err,[NVarA,NRep,NP,]);
%         VC2f = zeros([0,0,NP],'gpuArray');
        VC2f = zeros([0,0,NP]);
    else        
%         g = zeros([NP,2*NVarA+NVarA*(NVarA-1)/2+NVarNLT+NVarS],'gpuArray');
        g = zeros([NP,2*NVarA+NVarA*(NVarA-1)/2+NVarNLT+NVarS]);
%         VC2 = zeros(0,0,NP,'gpuArray');
        VC2 = zeros(0,0,NP);
        VC2f = reshape(err,[NVarA,NRep,NP]);
    end
    
    if any(isnan(XXa(:))) == 0 % faster version for complete dataset
        YYy = (YY == 1);              
        parfor n = 1:NP % for each person
            F3sum = [];
            b_mtx_n = b_mtx(:,:,n);
            XXa_n = XXa(:,:,n);
            U = reshape(XXa_n*b_mtx_n,[NAlt,NCT,NRep]);  % NAlt x NCT x NRep
            U = exp(U - max(U,[],1));  % rescale utility to avoid exploding
            % probability for each alternative
            U_prob = U./sum(U,1); % NAlt x NCT x NRep
            YYy_n = YYy(:,n);
            % probability for the chosen alternative (prod for panel data)
            U_prod = prod(reshape(U_prob(YYy_n(:,ones(NRep,1))),[NCT,NRep]),1);  % 1 x NRep
            % averaging for getting the chosen alternative probability estimator
            if RealMin == 1
                p0(n) = max(mean(U_prod,2),realmin);
            else
                p0(n) = mean(U_prod,2);
            end
            % calculations for gradient
            U_prob = reshape(U_prob,[NAlt*NCT,1,NRep]);  % NAlt*NCT x NVarA x NRep
            if WTP_space == 0
                % auxiliary variables for gradient calculation
                % summing over the first dimension -- alternatives
                X_hat = sum(reshape(U_prob.*XXa_n,[NAlt,NCT,NVarA,NRep]),1);
                F = XXa_n(YYy_n,:) - reshape(X_hat,[NCT,NVarA,NRep]);  %NCT x NVarA x NRep
                if NVarS > 0
                    FScale = sum(F.*reshape(b_mtx_n,[1,NVarA,NRep]),2);
                    % FScale = squeeze(sum(FScale,1));
                    % FScale = reshape(sum(FScale,1),[1,NRep]);
                    FScale = reshape(sum(FScale.*Xs_sliced(:,:,n),1),[NVarS,NRep]);
                end
                % sum over choice tasks/situations
                sumFsqueezed = reshape(sum(F,1),[NVarA,NRep]);  %NVarA x NRep
                if sum(Dist == 1) > 0
                    sumFsqueezed(Dist == 1,:) = sumFsqueezed(Dist == 1,:).*b_mtx_n(Dist == 1,:);
                end
                if NVarNLT == 1
                    XXt_n = XXt(:,:,n);
                    XXtt = XXt_n*b_mtx_n(NLTVariables,:); %NAlt*NCT x NRep
                    X_hat_lam = sum(reshape(reshape(U_prob(:,1,:),[NAlt*NCT,NRep]).*XXtt,[NAlt,NCT,NRep]),1);
                    F3 = XXtt(YY(:,n) == 1,:) - reshape(X_hat_lam,[NCT,NRep]); % CT x NRep
                    F3sum = sum(F3,1); % 1  x NRep
                elseif NVarNLT > 1
                    XXt_n = XXt(:,:,n);
                    XXtt = XXt_n.*permute(b_mtx_n(NLTVariables,:,ones(NCT*NAlt,1)),[3,1,2]); %NAlt*NCT x NVarNLT x NRep
                    X_hat_lam = sum(reshape(U_prob.*XXtt,[NAlt,NCT,NVarNLT,NRep]),1);
                    F3 = XXtt(YY(:,n) == 1,:,:) - reshape(X_hat_lam,[NCT,NVarNLT,NRep]); % CT x NVarNLT x NRep
                    F3sum = reshape(sum(F3,1),[NVarNLT,NRep]); % NVarNLT  x NRep                                        
                end
            else % WTP_space > 0
                Xalpha = XXa_n(:,1:end-WTP_space).*reshape(b_mtx_n(WTP_matrix,:),[1,NVarA-WTP_space,NRep]);
                % for non-cost variables
                X_hat1 = sum(reshape(U_prob.*Xalpha,[NAlt,NCT,NVarA-WTP_space,NRep]),1);
                F1 = Xalpha(YY(:,n) == 1,:,:) - reshape(X_hat1,[NCT,NVarA-WTP_space,NRep]);  %NCT x NVarA-WTP_space x NRep
                % for cost variables
                b_mtx_grad_n = b_mtx_grad(:,:,n);
                if WTP_space == 1 % without starting the loop
                    pX = XXa_n(:,NVarA) + XXa_n(:,1:end-WTP_space)*b_mtx_grad_n(1:end-WTP_space,:);
                    X_hat2 = sum(reshape(reshape(U_prob,[NAlt*NCT,NRep]).*pX,[NAlt,NCT,WTP_space,NRep]),1);
                    F2 = pX(YY(:,n) == 1,:) - reshape(X_hat2,[NCT,NRep]);  %NCT x WTP_space x NRep
                else
                    pX = zeros(NCT*NAlt,WTP_space,NRep);
                    for i = 1:WTP_space
                        pX(:,i,:) = XXa_n(:,NVarA-WTP_space+i) + XXa_n(:,WTP_matrix == NVarA-WTP_space+i)*b_mtx_grad_n(WTP_matrix == NVarA-WTP_space+i,:);
                    end
                    X_hat2 = sum(reshape(U_prob.*pX,[NAlt,NCT,WTP_space,NRep]),1);
                    F2 = pX(YY(:,n) == 1,:,:) - reshape(X_hat2,[NCT,WTP_space,NRep]);  %NCT x WTP_space x NRep
                end
                if NVarS > 0
                    FScale = sum(reshape(F2,[NCT,WTP_space,NRep]).*reshape(b_mtx_grad_n(NVarA-WTP_space+1:end,:),[1,WTP_space,NRep]),2);
                    Xs_tmp = Xs_sliced(:,:,n);
                    FScale = FScale.*Xs_tmp;
                    FScale = reshape(sum(FScale,1),[NVarS,NRep]);
                end
                sumFsqueezed = [reshape(sum(F1,1),[NVarA-WTP_space,NRep]);reshape(sum(F2,1),[WTP_space,NRep])];  %NVarA x NRep
                if sum(Dist == 1) > 0
                    sumFsqueezed(Dist == 1,:) = sumFsqueezed(Dist == 1,:).*b_mtx_grad_n(Dist == 1,:);
                end
                if NVarNLT == 1
                    XXt_n = XXt(:,:,n);
                    XXtt = XXt_n*b_mtx_n(NLTVariables,:,:); %NAlt*NCT x NRep
                    X_hat_lam = sum(reshape(squeeze(U_prob(:,1,:)).*XXtt,[NAlt,NCT,NRep]),1);
%                     F3 = XXtt(YY(:,n) == 1,:) - squeeze(X_hat_lam) ; % CT x NRep
                    F3 = XXtt(YY(:,n) == 1,:) - reshape(X_hat_lam,[NCT,NRep]) ; % CT x NRep                                        
                    F3sum = sum(F3,1); % 1  x NRep
                elseif NVarNLT > 1
                    XXt_n = XXt(:,:,n);
                    XXtt = XXt_n.*permute(b_mtx_n(NLTVariables,:,ones(NCT*NAlt,1)),[3 1 2]); %NAlt*NCT x NVarNLT x NRep
                    X_hat_lam = sum(reshape(U_prob.*XXtt,[NAlt,NCT,NVarNLT,NRep]),1);
%                     F3 = XXtt(YY(:,n) == 1,:,:) - squeeze(X_hat_lam); % CT x NVarNLT x NRep
                    F3 = XXtt(YY(:,n) == 1,:,:) - reshape(X_hat_lam,[NCT,NVarNLT,NRep]); % CT x NVarNLT x NRep                                        
                    F3sum = squeeze(sum(F3,1)); % NVarNLT  x NRep
                end
            end
            
            if FullCov == 0
                sumFtmp1 = ones(NVarA, NRep);
                sumVC2tmp2 = ones(NVarA, NRep);
                
                % check which attributes has normal/lognormal or weibull distribution
                normDist = logical(sum([Dist == 0; Dist == 1], 1));
                weibDist = (Dist == 4);
                % calculations for gradient columns for normal and lognormal distributions
                if any(normDist ~= 0)
                    sumFtmp1(normDist, :) = sumFsqueezed(normDist, :);
                    sumVC2tmp2(normDist, :) = ...
                        sumFsqueezed(normDist, :).*VC2(normDist,:,n);  % NVarA x NRep
                end
                % calculations for gradient columns for Weibull distr.
                if any(weibDist ~= 0)
                    tmpWeib_n = tmpWeib(:, :, n);
                    sumFtmp1(weibDist, :) = sumFsqueezed(weibDist, :).*tmpWeib_n.^b0weibB;
                    sumVC2tmp2(weibDist, :) = ...
                        sumFsqueezed(weibDist, :).*(-b_mtx_n(Dist == 4,:).*log(tmpWeib_n).*b0weibB.^2);
                end
                gtmp = -mean([sumFtmp1.*U_prod;sumVC2tmp2.*U_prod],2)./p0(n);
            else % FullCov = 1
                sumVC2tmp = sumFsqueezed(indx1,:).*VC2f(indx2,:,n);
                gtmp = -mean([sumFsqueezed.*U_prod;sumVC2tmp.*U_prod],2)./p0(n);
            end
            if NVarS > 0
                gtmp = [gtmp;-mean(FScale.*U_prod,2)./p0(n)];
            end
            if NVarNLT > 0
                gtmp = [gtmp;-mean(F3sum.*U_prod,2)./p0(n)];
            end
            g(n,:) = gtmp';
        end
        
    else % data has missing choices
        
% save tmp1
        parfor n = 1:NP
            YYy_n = YY(:,n) == 1;
            YnanInd = ~isnan(YY(:,n));
            b_mtx_n = b_mtx(:,:,n);
            XXa_n = XXa(:,:,n);
            NAltMissIndExp_n = NAltMissIndExp(:,n);
            NAltMissIndExp_n = NAltMissIndExp_n(YnanInd);
            if WTP_space > 0
                b_mtx_wtp = reshape(b_mtx_n,[1,NVarA,NRep]);
                Xalpha = XXa_n(YnanInd,1:end-WTP_space).*b_mtx_wtp(:,WTP_matrix,:);
                b_mtx_grad_n = b_mtx_grad(:,:,n);
            end
            if var(NAltMissIndExp_n(NAltMissIndExp_n > 0)) == 0 % if NAlt is constant per individual (but can vary between individuals)
                U = reshape(XXa_n(YnanInd,:)*b_mtx_n,[NAltMiss(n),NCTMiss(n),NRep]);
                U = exp(U - max(U,[],1));
                U_sum = reshape(sum(U,1),[1,NCTMiss(n),NRep]);
                U_prob = U./U_sum;  % NAlt x NCT x NRep
                U_prod = prod(reshape(U_prob(YYy_n(YnanInd,ones(NRep,1))),[NCTMiss(n),NRep]),1);  % 1 x NRep
                U_prob = reshape(U_prob,[NAltMiss(n)*NCTMiss(n),1,NRep]);  % NAlt*NCT x NVarA x NRep
                
                if WTP_space == 0
                    X_hat = sum(reshape(U_prob.*XXa_n(YnanInd,:),[NAltMiss(n),NCTMiss(n),NVarA,NRep]),1);
                    X_hat = reshape(X_hat,[NCTMiss(n),NVarA,NRep]);
                else
                    X_hat1 = sum(reshape(U_prob.*Xalpha,[NAltMiss(n),NCTMiss(n),NVarA-WTP_space,NRep]),1);
                    X_hat1 = reshape(X_hat1,[NCTMiss(n),NVarA-WTP_space,NRep]);
                    if WTP_space == 1
                        pX = XXa_n(YnanInd,NVarA) + XXa_n(YnanInd,1:end-WTP_space)*b_mtx_grad_n(1:end-WTP_space,:);
                        X_hat2 = sum(reshape(reshape(U_prob,[NAltMiss(n).*NCTMiss(n),NRep]).*pX,[NAltMiss(n),NCTMiss(n),WTP_space,NRep]),1);
                        X_hat2 = reshape(X_hat2,[NCTMiss(n),NRep]);
                    end
                end
            else
                NAltMissInd_n = NAltMissInd(:,n);
                NAltMissInd_n = NAltMissInd_n(MissingCT(:,n) == 0);
                U = XXa_n(YnanInd,:)*b_mtx_n;
                Uniq = unique(NAltMissIndExp_n);
                U_prob = zeros(size(U,1),1,NRep);
                XXa_tmp = XXa_n(YnanInd,:);
                if WTP_space == 0
                    X_hat = zeros(NCTMiss(n),NVarA,NRep);
                else
                    X_hat1 = zeros(NCTMiss(n),NVarA-WTP_space,NRep);
                    X_hat2 = zeros(NCTMiss(n),NRep);
                    pX = XXa_tmp(:,NVarA) + XXa_tmp(:,1:end-WTP_space)*b_mtx_grad_n(1:end-WTP_space,:);
                end
                if length(Uniq) == 2 % choice tasks with 2 different numbers of alternatives
                    U_tmp = U(NAltMissIndExp_n == Uniq(1),:);
                    U_tmp = reshape(U_tmp,[Uniq(1),size(U_tmp,1)/Uniq(1),NRep]);
                    U_max_tmp = max(U_tmp);
                    U_tmp = exp(U_tmp - U_max_tmp);
                    U_sum = reshape(sum(U_tmp,1),[1,size(U_tmp,2),NRep]);
                    U_tmp = U_tmp./U_sum;
                    U_prob(NAltMissIndExp_n == Uniq(1),:,:) = reshape(U_tmp,[size(U_tmp,2)*Uniq(1),1,NRep]);
                    if WTP_space == 0
                        X_hat_tmp = sum(reshape(U_prob(NAltMissIndExp_n == Uniq(1),:,:).*XXa_tmp(NAltMissIndExp_n == Uniq(1),:),[Uniq(1),size(U_tmp,2),NVarA,NRep]),1);
                        X_hat(NAltMissInd_n == Uniq(1),:,:) = reshape(X_hat_tmp,[size(U_tmp,2),NVarA,NRep]);
                    else
                        X_hat1_tmp = sum(reshape(U_prob(NAltMissIndExp_n == Uniq(1),:,:).*Xalpha(NAltMissIndExp_n == Uniq(1),:,:),[Uniq(1),size(U_tmp,2),NVarA-WTP_space,NRep]),1);
                        X_hat1(NAltMissInd_n == Uniq(1),:,:) = reshape(X_hat1_tmp,[size(U_tmp,2),NVarA-WTP_space,NRep]);
                        X_hat2_tmp = sum(reshape(reshape(U_prob(NAltMissIndExp_n == Uniq(1),:,:),[size(U_prob(NAltMissIndExp_n == Uniq(1),:,:),1),NRep]).*pX(NAltMissIndExp_n == Uniq(1),:),[Uniq(1),size(U_tmp,2),NRep]),1);
                        X_hat2(NAltMissInd_n == Uniq(1),:,:) = reshape(X_hat2_tmp,[size(U_tmp,2),NRep]);
                    end
                    U_tmp = U(NAltMissIndExp_n == Uniq(2),:);
                    U_tmp = reshape(U_tmp,[Uniq(2),size(U_tmp,1)/Uniq(2),NRep]);
                    U_tmp = exp(U_tmp - max(U_tmp));
                    U_sum = reshape(sum(U_tmp,1),[1,size(U_tmp,2),NRep]);
                    U_tmp = U_tmp./U_sum;
                    U_prob(NAltMissIndExp_n == Uniq(2),:,:)= reshape(U_tmp,[size(U_tmp,2)*Uniq(2),1,NRep]);
                    if WTP_space == 0
                        X_hat_tmp = sum(reshape(U_prob(NAltMissIndExp_n == Uniq(2),:,:).*XXa_tmp(NAltMissIndExp_n == Uniq(2),:),[Uniq(2),size(U_tmp,2),NVarA,NRep]),1);
                        X_hat(NAltMissInd_n == Uniq(2),:,:) = reshape(X_hat_tmp,[size(U_tmp,2),NVarA,NRep]);
                    else
                        X_hat1_tmp = sum(reshape(U_prob(NAltMissIndExp_n == Uniq(2),:,:).*Xalpha(NAltMissIndExp_n == Uniq(2),:,:),[Uniq(2),size(U_tmp,2),NVarA-WTP_space,NRep]),1);
                        X_hat1(NAltMissInd_n == Uniq(2),:,:) = reshape(X_hat1_tmp,[size(U_tmp,2),NVarA-WTP_space,NRep]);
                        X_hat2_tmp = sum(reshape(reshape(U_prob(NAltMissIndExp_n == Uniq(2),:,:),[size(U_prob(NAltMissIndExp_n == Uniq(2),:,:),1),NRep]).*pX(NAltMissIndExp_n == Uniq(2),:),[Uniq(2),size(U_tmp,2),NRep]),1);
                        X_hat2(NAltMissInd_n == Uniq(2),:,:) = reshape(X_hat2_tmp,[size(U_tmp,2),NRep]);
                    end
                else % choice tasks with 3+ different numbers of alternatives
                    for i = 1:length(Uniq)
                        U_tmp = U(NAltMissIndExp_n == Uniq(i),:);
                        U_tmp = reshape(U_tmp,[Uniq(i),size(U_tmp,1)/Uniq(i),NRep]);
                        U_tmp = exp(U_tmp - max(U_tmp));
                        U_sum = reshape(sum(U_tmp,1),[1,size(U_tmp,2),NRep]);
                        U_tmp = U_tmp./U_sum;
                        U_prob(NAltMissIndExp_n == Uniq(i),:,:) = reshape(U_tmp,[size(U_tmp,2)*Uniq(i),1,NRep]);
                        if WTP_space == 0
                            X_hat_tmp = sum(reshape(U_prob(NAltMissIndExp_n == Uniq(i),:,:).*XXa_tmp(NAltMissIndExp_n == Uniq(i),:),[Uniq(i),size(U_tmp,2),NVarA,NRep]),1);
                            X_hat(NAltMissInd_n == Uniq(i),:,:) = reshape(X_hat_tmp,[size(U_tmp,2),NVarA,NRep]);
                        else
                            X_hat1_tmp = sum(reshape(U_prob(NAltMissIndExp_n == Uniq(i),:,:).*Xalpha(NAltMissIndExp_n == Uniq(i),:,:),[Uniq(i),size(U_tmp,2),NVarA-WTP_space,NRep]),1);
                            X_hat1(NAltMissInd_n == Uniq(i),:,:) = reshape(X_hat1_tmp,[size(U_tmp,2),NVarA-WTP_space,NRep]);
                            X_hat2_tmp = sum(reshape(reshape(U_prob(NAltMissIndExp_n == Uniq(i),:,:),[size(U_prob(NAltMissIndExp_n == Uniq(i),:,:),1),NRep]).*pX(NAltMissIndExp_n == Uniq(i),:),[Uniq(i),size(U_tmp,2),NRep]),1);
                            X_hat2(NAltMissInd_n == Uniq(i),:,:) = reshape(X_hat2_tmp,[size(U_tmp,2),NRep]);
                        end
                    end
                end
                U_prod = prod(reshape(U_prob(YYy_n(YnanInd,ones(NRep,1))),[NCTMiss(n),NRep]),1);  % 1 x NRep
            end
            
            if RealMin == 1
                p0(n) = max(mean(U_prod),realmin);
            else
                p0(n) = mean(U_prod);
            end
                        
            % calculations for gradient
            if WTP_space == 0
                F = XXa_n(YYy_n,:) - X_hat;  %NCT x NVarA x NRep
                if NVarS > 0
                    FScale = sum(F.*reshape(b_mtx_n,[1,NVarA,NRep]),2);
                    Xs_tmp = Xs_sliced(:,:,n); %
                    Xs_tmp = Xs_tmp(MissingCT(:,n) == 0,:);
                    FScale = reshape(sum(FScale.*Xs_tmp,1),[NVarS,NRep]);
                end
                sumFsqueezed = reshape(sum(F,1),[NVarA,NRep]);  %NVarA x NRep
                if sum(Dist == 1) > 0
                    sumFsqueezed(Dist == 1,:) = sumFsqueezed(Dist == 1,:).*b_mtx_n(Dist == 1,:);
                end
            else % WTP_space > 0
                F1 = Xalpha(YYy_n(YnanInd) == 1,:,:) - X_hat1;  %NCT x NVarA-WTP_space x NRep
                % for cost variables
                if WTP_space == 1 % without starting the loop
                    F2 = pX(YYy_n(YnanInd) == 1,:,:) - X_hat2;
                else
                    pX = zeros(NCTMiss(n)*NAltMiss(n),WTP_space,NRep);
                    for i = 1:WTP_space
                        pX(:,i,:) = XXa_n(YnanInd,NVarA-WTP_space+i) + XXa_n(YnanInd,WTP_matrix == NVarA-WTP_space+i)*b_mtx_grad_n(WTP_matrix == NVarA-WTP_space+i,:);
                    end
                    X_hat2 = sum(reshape(U_prob.*pX,[NAltMiss(n),NCTMiss(n),WTP_space,NRep]),1);
                    F2 = pX(YYy_n(YnanInd) == 1,:,:) - reshape(X_hat2,[NCTMiss(n),WTP_space,NRep]);
                end
                sumFsqueezed = [reshape(sum(F1,1),[NVarA-WTP_space,NRep]);reshape(sum(F2,1),[WTP_space,NRep])];  %NVarA x NRep
                if NVarS > 0
                    FScale = sum(reshape(F2,[NCTMiss(n),WTP_space,NRep]).*reshape(b_mtx_grad_n(NVarA-WTP_space+1:end,:),[1,WTP_space,NRep]),2);
                    Xs_tmp = Xs_sliced(:,:,n); %
                    Xs_tmp = Xs_tmp(MissingCT(:,n) == 0,:);
                    FScale = FScale.*Xs_tmp;
                    FScale = reshape(sum(FScale,1),[NVarS,NRep]);
                end
                if sum(Dist == 1) > 0
                    sumFsqueezed(Dist == 1,:) = sumFsqueezed(Dist == 1,:).*b_mtx_grad_n(Dist == 1,:);
                end
            end
            
            if FullCov == 0
                sumFtmp1 = ones(NVarA, NRep);
                sumVC2tmp2 = ones(NVarA, NRep);
                
                % check which attributes has normal/lognormal or weibull distribution
                normDist = logical(sum([Dist == 0; Dist == 1], 1));
                weibDist = (Dist == 4);
                % calculations for gradient columns for normal and lognormal distributions
                if any(normDist ~= 0)
                    sumFtmp1(normDist, :) = sumFsqueezed(normDist, :);
                    sumVC2tmp2(normDist, :) = ...
                        sumFsqueezed(normDist, :).*VC2(normDist,:,n);  % NVarA x NRep
                end
                % calculations for gradient columns for Weibull distr.
                if any(weibDist ~= 0)
                    tmpWeib_n = tmpWeib(:, :, n);
                    sumFtmp1(weibDist, :) = sumFsqueezed(weibDist, :).*tmpWeib_n.^b0weibB;
                    sumVC2tmp2(weibDist, :) = ...
                        sumFsqueezed(weibDist, :).*(-b_mtx_n(Dist == 4,:).*log(tmpWeib_n).*b0weibB.^2);
                end
                gtmp = -mean([sumFtmp1.*U_prod;sumVC2tmp2.*U_prod],2)./p0(n);
            else % FullCov = 1
                sumVC2tmp = sumFsqueezed(indx1,:).*VC2f(indx2,:,n);
                gtmp = -mean([sumFsqueezed.*U_prod;sumVC2tmp.*U_prod],2)./p0(n);
            end
            if NVarS > 0
                gtmp = [gtmp;-mean(FScale.*U_prod,2)./p0(n)];
            end
            
            g(n,:) = gtmp';
            
        end
    end
    if NVarM > 0
        gm =  g(:,repmat(1:NVarA,[1,NVarM])).*(XXm(kron(1:NVarM,ones(1,NVarA)),:)');
        if FullCov == 0
            g = [g(:,1:2*NVarA),gm,g(:,2*NVarA+1:end)];
        else
            g = [g(:,1:NVarA*(NVarA/2+1.5)),gm,g(:,NVarA*(NVarA/2+1.5)+1:end)];
        end
    end  
    
elseif nargout == 3 % function value + gradient + hessian
    
    if FullCov == 0
        g = zeros(NP,2*NVarA);
        VC2 = reshape(2*diag(b0(NVarA+1:NVarA*2))*err,[NVarA,NRep,NP]);
        VC2f = zeros(0,0,NP);
        SIG = 4*b0(NVarA+1:NVarA*2)*b0(NVarA+1:NVarA*2)';
        %         SIG = SIG(:,:,ones(NRep,1));
        VC3 = reshape(2*err,[NVarA,NRep,NP]);
        NVarG = 2*NVarA; % general no. of parameters
        VCxx0 = diag(ones(NVarA,1));
        VCxx0 = VCxx0(:,:,ones(NRep,1));
    else
        g = zeros(NP,2*NVarA+NVarA*(NVarA-1)/2);
        VC2 = zeros(0,0,NP);
        VC2f = reshape(err,[NVarA,NRep,NP]);
        VC3 = zeros(0,0,NP); % for parfor
        SIG = []; % for parfor
        VCxx0 = []; % for parfor
        NVarG = 2*NVarA+NVarA*(NVarA-1)/2;
    end
    hx1 = zeros(NVarG,NVarG,NP);
    hx2 = zeros(NVarG,NVarG,NP);
    XXX = permute(mmx('square',permute(XXa,[2,4,1,3]),[]),[3,1,2,4]);
    
    if any(isnan(XXa(:))) == 0 % faster version for complete dataset
        YYy = (YY == 1);
        parfor n = 1:NP
            XXa_n = XXa(:,:,n);
            YYy_n = YYy(:,n);
            b_mtx_n = b_mtx(:,:,n);
            U = reshape(XXa_n*b_mtx_n,[NAlt,NCT,NRep]);  % NAlt x NCT x NRep
            %             U_max = max(U);
            %             U = exp(U - U_max(ones(NAlt,1),:,:));  % rescale utility to avoid exploding
            U = exp(U - max(U,[],1));  % rescale utility to avoid exploding
            %             U_sum = reshape(sum(U,1),1,NCT,NRep);
            %             U_prob = U./U_sum(ones(NAlt,1,1),:,:);  % NAlt x NCT x NRep
            U_prob = U./sum(U,1); % NAlt x NCT x NRep
            %             U_prod = prod(reshape(U_prob(YYy_n(:,ones(NRep,1))),NCT,NRep),1);  % 1 x NRep
            U_prod = prod(reshape(U_prob(YYy_n(:,ones(NRep,1))),[NCT,NRep]),1);  % 1 x NRep
            %             p0(n) = mean(U_prod);
            %             p0(n) = max(mean(U_prod),realmin);
            if RealMin == 1
                p0(n) = max(mean(U_prod),realmin);
            else
                p0(n) = mean(U_prod);
            end
            
            % calculations for gradient
            U_prob = reshape(U_prob,[NAlt*NCT,1,NRep]);  % NAlt*NCT x 1 x NRep
            X_hat = reshape(sum(reshape(U_prob.*XXa_n,[NAlt,NCT,NVarA,NRep]),1),[NCT,NVarA,NRep]); % NCT x NVarA x NRep
            F = XXa_n(YYy_n,:) - X_hat;  %NCT x NVarA x NRep
            sumFsqueezed = reshape(sum(F,1),[NVarA,NRep]);  %NVarA x NRep
            if FullCov == 0
                sumVC2tmp = sumFsqueezed.*VC2(:,:,n);  % NVarA x NRep
                g(n,:) = -mean([sumFsqueezed.*U_prod;sumVC2tmp.*U_prod],2)./p0(n);
            else % FullCov = 1
                sumVC2tmp = sumFsqueezed(indx1,:).*VC2f(indx2,:,n);
                g(n,:) = -mean([sumFsqueezed.*U_prod;sumVC2tmp.*U_prod],2)./p0(n);
            end

            % calculations for hessian
            % second term of hessian
            gtmp = [sumFsqueezed;sumVC2tmp];
            hx1(:,:,n) = ((gtmp.*U_prod)*gtmp')./(p0(n)*NRep); %NVarG x NVarG
            
            %third term of hessian
            U_prod = reshape(U_prod,[1,1,NRep]); % 1 x 1 x NRep
            U_prob = reshape(U_prob,[NAlt*NCT,1,1,NRep]); % NAlt*NCT x 1 x 1 x NRep
            X_hat_U = reshape(permute(X_hat.*U_prod,[1 3 2]),[NCT*NRep,NVarA]);
            X_hat = permute(X_hat,[1 3 2]); % NCT x NRep x NVarA
            
            hbbx = reshape(-sum(U_prob.*XXX(:,:,:,n),1),[NVarA,NVarA,NRep]);
            hbb = (mean(hbbx.*U_prod,3)+(X_hat_U')*reshape(X_hat,[NCT*NRep,NVarA])./NRep)./p0(n);
            
            % partial derivative of standard deviations and b
            
            if FullCov == 0
                VC2_n_tmp = reshape(permute(VC2(:,:,n),[2 1]),[1,NRep,NVarA]); % 1 x NRep x NVC
                X_hat_tmp_VC = reshape(X_hat.*VC2_n_tmp,[NCT*NRep,NVarA]);
                hsb = (mean(hbbx.*(U_prod.*reshape(VC2(:,:,n),[1,NVarA,NRep])),3) + X_hat_U'*X_hat_tmp_VC./NRep)./p0(n);
                VCxx = zeros(NVarA,NVarA,NRep);
                VCxx(VCxx0 == 1) = sumFsqueezed.*VC3(:,:,n);
                VCx_n = mmx('square',reshape(VC3(:,:,n)/2,[NVarA,1,NRep]),[]).*SIG;
                hss = (mean(U_prod.*(hbbx.*VCx_n + VCxx),3) + ...
                reshape(reshape(X_hat_U,[NCT,NRep,NVarA]).*VC2_n_tmp,[NCT*NRep,NVarA])'*X_hat_tmp_VC./NRep)./p0(n);
            else
                VC2_n_tmp = reshape(permute(VC2f(indx2,:,n),[2 1]),[1,NRep,NVarA*(NVarA-1)/2+NVarA]); % 1 x NRep x NVC
                X_hat_tmp_VC = reshape(X_hat(:,:,indx1).*VC2_n_tmp,[NCT*NRep,NVarA*(NVarA-1)/2+NVarA]);
                hsb = (mean((hbbx(:,indx1,:).*U_prod).*reshape(VC2f(indx2,:,n),[1,NVarA*(NVarA-1)/2+NVarA,NRep]),3) + X_hat_U'*X_hat_tmp_VC./NRep)./p0(n);
                hss = (mean((U_prod.*hbbx(indx1,indx1,:)).*mmx('square',reshape(VC2f(indx2,:,n),[NVarA*(NVarA-1)/2+NVarA,1,NRep]),[]),3) + ...
                    reshape(reshape(X_hat_U(:,indx1),[NCT,NRep,NVarA*(NVarA-1)/2+NVarA]).*VC2_n_tmp,[NCT*NRep,NVarA*(NVarA-1)/2+NVarA])'*X_hat_tmp_VC./NRep)./ p0(n);
            end
            hx2(:,:,n) = [hbb,hsb;hsb',hss];
        end
        
    else % this does not currently work
        %         parfor n = 1:NP
        %             U = reshape(XXa(~isnan(YY(:,n)),:,n)*b_mtx(:,:,n),NAlt,NCT-sum(isnan(YY(1:NAlt:end,n))),NRep);
        %             U_max = max(U);
        %             U = exp(U - U_max(ones(NAlt,1),:,:));  % rescale utility to avoid exploding
        %             U_sum = reshape(nansum(U,1),1,NCT-sum(isnan(YY(1:NAlt:end,n))),NRep);
        %             U_prob = U./U_sum(ones(NAlt,1,1),:,:);  % NAlt x NCT x NRep
        %             U_prod = prod(reshape(U_prob(YY(~isnan(YY(:,n)),n*ones(NRep,1))==1),NCT-sum(isnan(YY(1:NAlt:end,n))),NRep),1);  % 1 x NRep
        %             p0(n) = mean(U_prod);
        %             p0(n) = max(mean(U_prod),realmin);
        %
        %             % calculations for gradient
        %             U_prob = reshape(U_prob, NAlt*(NCT-sum(isnan(YY(1:NAlt:end,n)))),1, NRep);  % NAlt*NCT x NVarA x NRep
        %             X_hat = sum(reshape(U_prob(:,ones(1,NVarA,1),:).* XXa(:,:, n*ones(NRep,1)), NAlt, NCT-sum(isnan(YY(1:NAlt:end,n))), NVarA, NRep),1);
        %             F = XXa(YY(:,n) == 1,:,n*ones(NRep,1)) - squeeze(X_hat);  %NCT x NVarA x NRep
        %             sumFsqueezed = squeeze(sum(F,1));  %NVarA x NRep
        %             if FullCov == 0
        %                 sumVC2tmp = sumFsqueezed.*VC2(:,:,n);  % NVarA x NRep
        %                 g(n,:) = -mean([sumFsqueezed.*U_prod(ones(NVarA,1),:); sumVC2tmp.*U_prod(ones(NVarA,1),:)],2)./p0(n);
        %             else % FullCov = 1
        %                 sumVC2tmp = sumFsqueezed(indx1,:).*VC2(indx2,:,n);
        %                 g(n,:) = -mean([sumFsqueezed.*U_prod(ones(NVarA,1),:); sumVC2tmp.*U_prod(ones(NVarA*(NVarA-1)/2+NVarA,1),:)],2)./p0(n);
        %             end
        %         end;
    end
    h = g'*g - sum(hx1 + hx2,3);
end

% p0 = gather(p0);

if RealMin == 1
    f = -log(max(p0,realmin));
else
    f = -log(p0);
end
