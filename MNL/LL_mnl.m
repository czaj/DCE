function [f,g] = LL_mnl(y,Xa,Xm,Xs,EstimOpt,b0)

% save LL_mnl
% return

RealMin = EstimOpt.RealMin;
NVarA = EstimOpt.NVarA;
NVarM = EstimOpt.NVarM;
NVarS = EstimOpt.NVarS;
WTP_space = EstimOpt.WTP_space;
WTP_matrix = EstimOpt.WTP_matrix;
NP = EstimOpt.NP;
NAlt = EstimOpt.NAlt;
NCTMiss = EstimOpt.NCTMiss;
if isfield(EstimOpt,'XmIndx')
    XmIndx = EstimOpt.XmIndx;
end
if isfield(EstimOpt,'XmIndx2')
    XmIndx2 = EstimOpt.XmIndx2;
end
NVarNLT = EstimOpt.NVarNLT;
if NVarNLT > 0
    NLTVariables = EstimOpt.NLTVariables;
    NLTType = EstimOpt.NLTType;
end
NumGrad = EstimOpt.NumGrad;
ExpB = EstimOpt.ExpB;

if ~isempty(ExpB)
    b0(ExpB) = exp(b0(ExpB));
end

B = b0(1:NVarA*(1+NVarM)+NVarS);

if NVarNLT > 0
    bt = b0(NVarA*(1+NVarM)+NVarS+1:end);
    IndTransNon0 = (abs(bt) > eps)';
    Xt = Xa(:,NLTVariables);
    if NLTType == 1 % BC
        % Xt(:,IndTransNon0) = -(Xt(:,IndTransNon0 == 1).^(bt(IndTransNon0 == 1,ones(size(Xa,1),1))') - 1)./bt(IndTransNon0 == 1,ones(size(Xa,1),1))';
        Xt(:,IndTransNon0) = -(Xt(:,IndTransNon0 == 1).^(bt(IndTransNon0 == 1,:)') - 1)./bt(IndTransNon0 == 1,:)';
        Xt(:,~IndTransNon0) = -log(Xt(:,IndTransNon0 == 0));
    elseif NLTType == 2 % YJ
        IndXtNon0 = (Xt >= 0);
        IndXtCase1 = IndXtNon0 & IndTransNon0; % X >= 0, lam ~= 0
        IndXtCase2 = IndXtNon0 & ~IndTransNon0; % X >= 0, lam = 0
        IndTransNon2 = (abs(bt - 2) > eps)';
        IndXtCase3 = ~IndXtNon0 & IndTransNon2;  % X < 0, lam ~= 2
        IndXtCase4 = ~IndXtNon0 & ~IndTransNon2; % X < 0, lam = 2
        bt_tmp = bt(:,ones(size(Xa,1),1))';
        Xt(IndXtCase1) = ((Xt(IndXtCase1) + 1).^bt_tmp(IndXtCase1) - 1)./bt_tmp(IndXtCase1);
        Xt(IndXtCase2) = log(Xt(IndXtCase2) + 1);
        Xt(IndXtCase3) = -((-Xt(IndXtCase3) + 1).^(2 - bt_tmp(IndXtCase3)) - 1)./(2 - bt_tmp(IndXtCase3));
        Xt(IndXtCase4) = -log(-Xt(IndXtCase4) + 1);
    end
    if NumGrad == 0 %
        if NLTType == 1 % BC
            XXt = Xa(:,NLTVariables);
            %nie jestem pewien czy to jest dobrze
            %XXt2(:, IndTransNon0 == 1) = -(XXt(:, IndTransNon0 == 1).^(bt(IndTransNon0 == 1,ones(size(Xa,1), 1))').*(bt(IndTransNon0 == 1,ones(size(Xa,1), 1))'.*log(XXt(:, IndTransNon0 == 1))-1)+1)./(bt(IndTransNon0 == 1,ones(size(Xa,1), 1)).^2)';
            XXt(:,IndTransNon0 == 1) = -(XXt(:, IndTransNon0 == 1).^(bt(IndTransNon0 == 1,:)').*(bt(IndTransNon0 == 1,:)'.*log(XXt(:, IndTransNon0 == 1))-1)+1)./(bt(IndTransNon0 == 1,:).^2)';
            XXt(:,IndTransNon0 == 0) = -0.5*log(XXt(:, IndTransNon0 == 0)).^2;
        elseif NLTType == 2 % YJ
            XXt = Xa(:, NLTVariables);
            XXt(IndXtCase1) = ((XXt(IndXtCase1)+1).^bt_tmp(IndXtCase1).*(bt_tmp(IndXtCase1).*log(XXt(IndXtCase1)+1)-1)+1)./(bt_tmp(IndXtCase1).^2);% X >= 0, lam ~= 0
            XXt(IndXtCase2) = 0.5*log(XXt(IndXtCase2)+1).^2;% X >= 0, lam == 0
            XXt(IndXtCase3) = -((-XXt(IndXtCase3)+1).^(2-bt_tmp(IndXtCase3)).*(1-(2-bt_tmp(IndXtCase3)).*log(-XXt(IndXtCase3)+1))-1)./((2-bt_tmp(IndXtCase3)).^2);% X < 0, lam ~= 2
            XXt(IndXtCase4) = -0.5*log(-XXt(IndXtCase4)+1).^2;% X < 0, lam == 2
        end
        XXt = XXt.*b0(NLTVariables, ones(size(Xa,1), 1))';
    end
    Xa(:, NLTVariables) = Xt;
end

if NVarM > 0 && (WTP_space > 0 || NVarNLT > 0)
    ba = B(1:NVarA);
    b0m = reshape(B(NVarA+1:NVarA*(1+NVarM)),[NVarA,NVarM]);
    %ba = ba(:,ones(NP,1)) + b0m*Xm'; % NVarA x NP
    ba = ba + b0m*Xm';
    % ba = ba + b0m*Xm(XmIndx2,:)';
end

N = nansum(y,1);

if WTP_space > 0
    if NVarM == 0
        B(1:NVarA-WTP_space,:) = B(1:NVarA-WTP_space,:).*B(WTP_matrix,:);
    else
        ba_grad = ba;
        ba(1:NVarA-WTP_space,:) = ba(1:NVarA-WTP_space,:).*ba(WTP_matrix,:);
    end
end

% computing utility levels
if NVarM > 0 && (WTP_space > 0 || NVarNLT > 0)
    if EstimOpt.mCT == 0
        betaX = zeros(size(Xa,1),1);
        for i = 1:NP
            NCTno = sum(NCTMiss(1:i-1));
            betaX(NCTno*NAlt+1:(NCTno+NCTMiss(i))*NAlt) = Xa(NCTno*NAlt+1:(NCTno+NCTMiss(i))*NAlt,:)*ba(:,i);
        end
    else
        betaX = sum(Xa.*ba',2);
    end
else
    betaX = Xa*B(1:NVarA,1); %clear XX b_mtx
end

if NVarS > 0
    if EstimOpt.SCEXP == 1
        cs = exp(Xs*B(NVarA*(1+NVarM)+1:NVarA*(1+NVarM)+NVarS));
    else
        cs = Xs*B(NVarA*(1+NVarM)+1:NVarA*(1+NVarM)+NVarS)+1;
    end
    betaX = betaX.*cs;
end

% Computations of utility levels
v = reshape(betaX,[NAlt,N]);
maxv = max(v,[],1);

% Nominators of probabilities fractions
evdiff = exp(v - maxv); %clear v maxv

% Denominators of probabilities fractions
sum_evdiff = nansum(evdiff,1); %1 by N

% Probabilities for each alternative
P = evdiff./sum_evdiff; % NAlt x N

% Probabilities for chosen alternatives
probs = P(y == 1);

if RealMin == 1
    logprobs = log(max(probs,realmin));
else
    logprobs = log(probs);
    % f = logprobs'.*EstimOpt.WT(1:end,:);
end

f = logprobs;


if nargout == 2 % function value + gradient
    IsNaN = any(isnan(Xa(:)));
    if NVarS > 0
        Xs = reshape(Xs(1:NAlt:end,:),[N,NVarS]);
    end
    if WTP_space == 0
        if NVarS > 0
           Xa = Xa.*cs; 
        end
        XXa = reshape(Xa,[NAlt,N,NVarA]);
        if IsNaN == 0
            Xhat = reshape(sum(P.*XXa,1),[N,NVarA]); % N x NVarA
        else
            Xhat = reshape(nansum(P.*XXa,1),[N,NVarA]); % N x NVarA
        end
        g = Xa(y == 1,:) - Xhat;
        if NVarS > 0
            gScale = g*B(1:NVarA,1);
            gScale = gScale.*Xs;
            g = [g,gScale];
        end
        if NVarNLT > 0
            if NVarS > 0
               XXt = XXt.*Scale; 
            end
            XXtt = reshape(XXt,[NAlt,N,NVarNLT]);
            if IsNaN == 0
                Xhatlam = reshape(sum(P.*XXtt,1),[N,NVarNLT]); % N x NVarNLT
            else
                Xhatlam = reshape(nansum(P.*XXtt,1),[N,NVarNLT]); % N x NVarNLT
            end
            g = [g,XXt(y == 1,:) - Xhatlam];
            if NVarM > 0
                gm = g(:,repmat(1:NVarA,[1,NVarM])).*(Xm(XmIndx,kron(1:NVarM,ones(1,NVarA))));
                g = [g(:,1:NVarA),gm,g(:,NVarA+1:end)];
            end
        end               
    else
        if NVarS > 0
           Xa = Xa.*cs;
        end
        % non - cost variables
        if NVarM == 0
            %alphaX = Xa(:,1:(NVarA - WTP_space)).*(B(WTP_matrix,ones(NAlt*N,1))');
            alphaX = Xa(:,1:(NVarA - WTP_space)).*(B(WTP_matrix,:)');
        else
            alphaX = Xa(:,1:(NVarA - WTP_space)).*ba_grad(WTP_matrix,XmIndx2)';
        end
        alphaXX = reshape(alphaX,[NAlt,N,NVarA - WTP_space]);
        if IsNaN == 0
            %Xhat1 = squeeze(sum(P(:,:,ones(NVarA- WTP_space,1)).*alphaXX,1));
            Xhat1 = squeeze(sum(P.*alphaXX,1));
        else
            %Xhat1 = squeeze(nansum(P(:,:,ones(NVarA- WTP_space,1)).*alphaXX,1));
            Xhat1 = squeeze(nansum(P.*alphaXX,1));
        end
        g1 = alphaX(y == 1,:) - Xhat1;
        % cost variables
        if WTP_space == 1
            if NVarM == 0
                pX = Xa(:,NVarA) + Xa(:,1:NVarA-1)*b0(1:NVarA-1);
            else
                pX = Xa(:,NVarA) + sum(Xa(:,1:NVarA-1).*ba_grad(1:NVarA-1,XmIndx2)',2);
            end
            pXX = reshape(pX,[NAlt,N]);
            if IsNaN == 0
                Xhat2 = squeeze(sum(P.*pXX,1))'; % N x 1
            else
                Xhat2 = squeeze(nansum(P.*pXX,1))'; % N x 1
            end
            g2 = pX(y == 1,:) - Xhat2;
        else
            pX = zeros(NAlt*N,WTP_space);
            for i = 1:WTP_space
                if NVarM == 0
                    pX(:,i) = Xa(:,NVarA - WTP_space + i) + Xa(:,WTP_matrix == NVarA - WTP_space + i)*b0(WTP_matrix == NVarA - WTP_space + i);
                else
                    pX(:,i) = Xa(:,NVarA - WTP_space + i) + sum(Xa(:, WTP_matrix == NVarA - WTP_space + i).*ba_grad(WTP_matrix == NVarA - WTP_space + i,XmIndx2)',2);
                end
            end
            pXX = reshape(pX,[NAlt,N,WTP_space]);
            if IsNaN == 0
                %Xhat2 = squeeze(sum(P(:,:, ones(WTP_space,1)).*pXX,1)); % N x WTP_space
                Xhat2 = squeeze(sum(P.*pXX,1));
            else
                %Xhat2 = squeeze(nansum(P(:,:, ones(WTP_space,1)).*pXX,1)); % N x WTP_space
                Xhat2 = squeeze(nansum(P.*pXX,1)); % N x WTP_space
            end
            g2 = pX(y == 1,:) - Xhat2;
        end
        g = [g1,g2];
        if NVarS > 0
            if NVarM == 0
                gScale = g2*B(NVarA-WTP_space+1:NVarA,:);
            else
                gScale = sum(g2.*ba_grad(NVarA-WTP_space+1:end,XmIndx)',2);
            end
           gScale = gScale.*Xs;
           g = [g,gScale];
        end
        if NVarNLT > 0
            XXtt = reshape(XXt,[NAlt,N,NVarNLT]);
            if IsNaN == 0
                %Xhatlam = squeeze(sum(P(:,:,ones(NVarNLT,1)).*XXtt,1)); % N x NVarNLT
                Xhatlam = squeeze(sum(P.*XXtt,1)); % N x NVarNLT
            else
                %Xhatlam = squeeze(nansum(P(:,:,ones(NVarNLT,1)).*XXtt,1)); % N x NVarNLT
                Xhatlam = squeeze(nansum(P.*XXtt,1)); % N x NVarNLT
            end
            if NVarNLT == 1
                Xhatlam = Xhatlam';
            end
            g = [g,XXt(y == 1,:) - Xhatlam];
        end
        if NVarM > 0
            gm = g(:,repmat(1:NVarA,[1,NVarM])).*(Xm(XmIndx,kron(1:NVarM,ones(1,NVarA))));
            g = [g(:,1:NVarA),gm,g(:,NVarA+1:end)];
        end
    end
end
