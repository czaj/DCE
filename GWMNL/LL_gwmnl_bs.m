function [P,J, logprobs] = LL_gwmnl_bs(y,Xa,EstimOpt,B)

% save tmp_MNL_like
% return


N = nansum(y,1);

% computing utility levels
if EstimOpt.WTP_space > 0
    b0 = B;
    B(1:EstimOpt.NVarA-EstimOpt.WTP_space,:) = B(1:EstimOpt.NVarA-EstimOpt.WTP_space,:).*B(EstimOpt.WTP_matrix,:);
else
    b0 = [];
end

betaX = Xa*B;
v = reshape(betaX,EstimOpt.NAlt,N);

maxv = max(v,[],1);
evdiff = exp(v- maxv(ones(EstimOpt.NAlt,1),:)); %clear v maxv
sum_evdiff = nansum(evdiff,1); %1 by N

PX = evdiff./sum_evdiff(ones(EstimOpt.NAlt,1),:); % NAlt x N
P = reshape(PX, EstimOpt.NAlt*N,1);
probs = PX(y == 1);

logprobs = log(max(probs,realmin)); %1 *N
% f = logprobs'.*EstimOpt.WT(1:end,:);
%f = sum(reshape(logprobs, EstimOpt.NCT, EstimOpt.NP),1)';

if EstimOpt.WTP_space == 0
    XXa = reshape(Xa, EstimOpt.NAlt, N, EstimOpt.NVarA);
    Xhat = squeeze(sum(PX(:,:,ones(EstimOpt.NVarA,1)).*XXa,1)); % N x NVarA
    J = Xa(y == 1, :) - Xhat;% N x NVarA
else
    alphaX = Xa(:, 1:(EstimOpt.NVarA - EstimOpt.WTP_space)).*(B(EstimOpt.WTP_matrix,ones(EstimOpt.NAlt*N,1))');
    alphaXX = reshape(alphaX, EstimOpt.NAlt, N, EstimOpt.NVarA - EstimOpt.WTP_space);
    Xhat1 = squeeze(sum(PX(:,:,ones(EstimOpt.NVarA- EstimOpt.WTP_space,1)).*alphaXX,1));
    g1 = alphaX(y == 1, :) - Xhat1;
    % cost variables
    if EstimOpt.WTP_space == 1
        pX = Xa(:, EstimOpt.NVarA) + Xa(:, 1:EstimOpt.NVarA-1)*b0(1:EstimOpt.NVarA-1);
        pXX = reshape(pX, EstimOpt.NAlt, N);
        Xhat2 = squeeze(sum(PX.*pXX,1))'; % N x 1
        g2 = pX(y == 1, :) - Xhat2;
    else
        pX = zeros(EstimOpt.NAlt*N, EstimOpt.WTP_space);
        for i = 1:EstimOpt.WTP_space
            pX(:,i) = Xa(:, EstimOpt.NVarA - EstimOpt.WTP_space + i) + Xa(:, EstimOpt.WTP_matrix == EstimOpt.NVarA - EstimOpt.WTP_space + i)*b0(EstimOpt.WTP_matrix == EstimOpt.NVarA - EstimOpt.WTP_space + i);
        end
        pXX = reshape(pX, EstimOpt.NAlt, N, EstimOpt.WTP_space);
        Xhat2 = squeeze(sum(PX(:,:, ones(EstimOpt.WTP_space,1)).*pXX,1)); % N x WTP_space
        g2 = pX(y == 1, :) - Xhat2;
    end
    J =[g1,g2];
end
    

% g = sum(reshape(g, EstimOpt.NCT, EstimOpt.NP, EstimOpt.NVarA),1);
% g = squeeze(g);
