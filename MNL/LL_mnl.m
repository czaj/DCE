function [f,g] = LL_mnl(y,Xa,Xs,EstimOpt,b0)

% save tmp_MNL_like
% return

B = b0(1:EstimOpt.NVarA + EstimOpt.NVarS);

if EstimOpt.NVarNLT > 0
    bt = b0(EstimOpt.NVarA + EstimOpt.NVarS + 1:end);
    % IndTransNon0 = (abs(bt) > 0.00001)';
    IndTransNon0 = (abs(bt) > eps)';

    Xt = Xa(:, EstimOpt.NLTVariables);

    
    
    
    if EstimOpt.NLTType == 1 % BC
        Xt(:, IndTransNon0) = -(Xt(:, IndTransNon0 == 1).^(bt(IndTransNon0 == 1,ones(size(Xa,1), 1))') - 1)./bt(IndTransNon0 == 1,ones(size(Xa,1), 1))';
        Xt(:, ~IndTransNon0) = -log(Xt(:, IndTransNon0 == 0));
    elseif EstimOpt.NLTType == 2 % YJ
        IndXtNon0 = (Xt >= 0);
        IndXtCase1 = IndXtNon0 & IndTransNon0(ones(size(IndXtNon0,1),1),:); % X >= 0, lam ~= 0
        IndXtCase2 = IndXtNon0 & ~IndTransNon0(ones(size(IndXtNon0,1),1),:); % X >= 0, lam = 0
    %     IndTransNon2 = (abs(bt - 2) < 0.00001)';
        IndTransNon2 = (abs(bt - 2) > eps)';
        IndXtCase3 = ~IndXtNon0 & IndTransNon2(ones(size(IndXtNon0,1),1),:);  % X < 0, lam ~= 2
        IndXtCase4 = ~IndXtNon0 & ~IndTransNon2(ones(size(IndXtNon0,1),1),:); % X < 0, lam = 2
        bt_tmp = bt(:,ones(size(Xa,1),1))';
        Xt(IndXtCase1) = ((Xt(IndXtCase1) + 1).^bt_tmp(IndXtCase1) - 1)./bt_tmp(IndXtCase1);
        Xt(IndXtCase2) = log(Xt(IndXtCase2) + 1);
        Xt(IndXtCase3) = -((-Xt(IndXtCase3) + 1).^(2 - bt_tmp(IndXtCase3)) - 1)./(2 - bt_tmp(IndXtCase3));
        Xt(IndXtCase4) = -log(-Xt(IndXtCase4) + 1);
    end
    
%     if EstimOpt.NumGrad == 0 % BC only
%        XXt =  Xa(:, EstimOpt.NLTVariables);
%     %    size(XXt(:, IndTransNon0 == 1))
%     %    size(bt(IndTransNon0 == 1,ones(size(Xa,1), 1))')
%        XXt(:, IndTransNon0 == 1) = (XXt(:, IndTransNon0 == 1).^(bt(IndTransNon0 == 1,ones(size(Xa,1), 1))').*(bt(IndTransNon0 == 1,ones(size(Xa,1), 1))'.*log(XXt(:, IndTransNon0 == 1))-1)+1)./(bt(IndTransNon0 == 1,ones(size(Xa,1), 1)).^2)';
%        XXt(:, IndTransNon0 == 0) = 0.5*log(XXt(:, IndTransNon0 == 0)).^2;
%        XXt = XXt.*b0(EstimOpt.NLTVariables, ones(size(Xa,1), 1))';
%     end
    if EstimOpt.NumGrad == 0 % 
       if EstimOpt.NLTType == 1 % BC
            XXt =  Xa(:, EstimOpt.NLTVariables);
            XXt(:, IndTransNon0 == 1) = -(XXt(:, IndTransNon0 == 1).^(bt(IndTransNon0 == 1,ones(size(Xa,1), 1))').*(bt(IndTransNon0 == 1,ones(size(Xa,1), 1))'.*log(XXt(:, IndTransNon0 == 1))-1)+1)./(bt(IndTransNon0 == 1,ones(size(Xa,1), 1)).^2)';
            XXt(:, IndTransNon0 == 0) = -0.5*log(XXt(:, IndTransNon0 == 0)).^2;
            
       elseif EstimOpt.NLTType == 2 % YJ
           XXt = Xa(:, EstimOpt.NLTVariables);
           XXt(IndXtCase1) = ((XXt(IndXtCase1)+1).^bt_tmp(IndXtCase1).*(bt_tmp(IndXtCase1).*log(XXt(IndXtCase1)+1)-1)+1)./(bt_tmp(IndXtCase1).^2);% X >= 0, lam ~= 0
           XXt(IndXtCase2) = 0.5*log(XXt(IndXtCase2)+1).^2;% X >= 0, lam == 0
           XXt(IndXtCase3) = -((-XXt(IndXtCase3)+1).^(2-bt_tmp(IndXtCase3)).*(1-(2-bt_tmp(IndXtCase3)).*log(-XXt(IndXtCase3)+1))-1)./((2-bt_tmp(IndXtCase3)).^2);% X < 0, lam ~= 2
           XXt(IndXtCase4) = -0.5*log(-XXt(IndXtCase4)+1).^2;% X < 0, lam == 2
       end
       XXt = XXt.*b0(EstimOpt.NLTVariables, ones(size(Xa,1), 1))';
    end
    
    Xa(:, EstimOpt.NLTVariables) = Xt;
end
    
N = nansum(y,1);

if EstimOpt.WTP_space > 0
    B(1:EstimOpt.NVarA-EstimOpt.WTP_space,:) = B(1:EstimOpt.NVarA-EstimOpt.WTP_space,:).*B(EstimOpt.WTP_matrix,:);
end

% computing utility levels

if EstimOpt.NVarS > 0 
    if EstimOpt.SCEXP==1
        cs = exp(Xs*B(end-EstimOpt.NVarS+1:end,1));
    else 
        cs = Xs(:,end-EstimOpt.NVarS+1:end)*B(end-EstimOpt.NVarS+1:end,1)+1;
    end
    betaX = Xa*B(1:end-EstimOpt.NVarS,1).*cs; %clear XX b_mtx
else
    betaX = Xa*B(1:EstimOpt.NVarA); %clear XX b_mtx
end

v = reshape(betaX,EstimOpt.NAlt,N);

maxv = max(v,[],1);
evdiff = exp(v- maxv(ones(EstimOpt.NAlt,1),:)); %clear v maxv
sum_evdiff = nansum(evdiff,1); %1 by N

P = evdiff./sum_evdiff(ones(EstimOpt.NAlt,1),:); % NAlt x N

probs = P(y == 1);

logprobs = log(max(probs,realmin)); %1 *N
% f = logprobs'.*EstimOpt.WT(1:end,:);
f = logprobs;

if nargout == 2 % function value + gradient

    if EstimOpt.WTP_space == 0
        XXa = reshape(Xa, EstimOpt.NAlt, N, EstimOpt.NVarA);
%         Xhat = squeeze(sum(P(:,:,ones(EstimOpt.NVarA,1)).*XXa,1))'; % N x NVarA
        Xhat = reshape(sum(P(:,:,ones(EstimOpt.NVarA,1)).*XXa,1),[N,EstimOpt.NVarA]); % N x NVarA
        g = Xa(y == 1, :) - Xhat;
        if EstimOpt.NVarNLT > 0
            XXtt = reshape(XXt, EstimOpt.NAlt, N, EstimOpt.NVarNLT);
            Xhatlam = squeeze(sum(P(:,:,ones(EstimOpt.NVarNLT,1)).*XXtt,1)); % N x NVarNLT
            if EstimOpt.NVarNLT == 1
                Xhatlam = Xhatlam';
            end
            g = [g, XXt(y == 1, :) - Xhatlam];
        end
    else
        % non - cost variables
        alphaX = Xa(:, 1:(EstimOpt.NVarA - EstimOpt.WTP_space)).*(B(EstimOpt.WTP_matrix,ones(EstimOpt.NAlt*N,1))');
        alphaXX = reshape(alphaX, EstimOpt.NAlt, N, EstimOpt.NVarA - EstimOpt.WTP_space);
        Xhat1 = squeeze(sum(P(:,:,ones(EstimOpt.NVarA- EstimOpt.WTP_space,1)).*alphaXX,1));
        g1 = alphaX(y == 1, :) - Xhat1;
        % cost variables

        if EstimOpt.WTP_space == 1
            pX = Xa(:, EstimOpt.NVarA) + Xa(:, 1:EstimOpt.NVarA-1)*b0(1:EstimOpt.NVarA-1);
            pXX = reshape(pX, EstimOpt.NAlt, N);
            Xhat2 = squeeze(sum(P.*pXX,1))'; % N x 1
            g2 = pX(y == 1, :) - Xhat2;
        else
            pX = zeros(EstimOpt.NAlt*N, EstimOpt.WTP_space);
            for i = 1:EstimOpt.WTP_space
                pX(:,i) = Xa(:, EstimOpt.NVarA - EstimOpt.WTP_space + i) + Xa(:, EstimOpt.WTP_matrix == EstimOpt.NVarA - EstimOpt.WTP_space + i)*b0(EstimOpt.WTP_matrix == EstimOpt.NVarA - EstimOpt.WTP_space + i);
            end
            pXX = reshape(pX, EstimOpt.NAlt, N, EstimOpt.WTP_space);
            Xhat2 = squeeze(sum(P(:,:, ones(EstimOpt.WTP_space,1)).*pXX,1)); % N x WTP_space
            g2 = pX(y == 1, :) - Xhat2;
        end
        g =[g1,g2];
        if EstimOpt.NVarNLT > 0
            XXtt = reshape(XXt, EstimOpt.NAlt, N, EstimOpt.NVarNLT);
            Xhatlam = squeeze(sum(P(:,:,ones(EstimOpt.NVarNLT,1)).*XXtt,1)); % N x NVarNLT
            if EstimOpt.NVarNLT == 1
                Xhatlam = Xhatlam';
            end
            g = [g,  XXt(y == 1, :) - Xhatlam];
        end
    end
end