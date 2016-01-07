function [f,g] = LL_lc(YY,Xa,Xc,MissingInd,EstimOpt,B)

% save tmp_LL_lc
% return

if sum(EstimOpt.BActiveClass == 0,1) == 0
    Bclass = reshape(B(1:EstimOpt.NClass*EstimOpt.NVarA), EstimOpt.NVarA, EstimOpt.NClass);
else
   Bclass = B(1:EstimOpt.NVarA)*ones(1,EstimOpt.NClass);
   for i = 1:(EstimOpt.NClass - 1)
       Bclass(EstimOpt.BActiveClass == 1,i+1) = B(EstimOpt.NVarA + (i-1)*sum(EstimOpt.BActiveClass,1)+1:EstimOpt.NVarA + i*sum(EstimOpt.BActiveClass,1));
   end
end
if EstimOpt.WTP_space > 0
    Bclass_nonwtp = Bclass; % for gradient
    Bclass(1:end - EstimOpt.WTP_space,:) = Bclass(1:end - EstimOpt.WTP_space,:).*Bclass(EstimOpt.WTP_matrix,:);
end
U = exp(reshape(Xa * Bclass,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP*EstimOpt.NClass)); ...% NAlt x NCT*NP*NClass
U(MissingInd==1) = 0;... % do not include alternatives which were not available
U_sum = sum(U,1);... 
U_probs =  U ./ U_sum(ones(EstimOpt.NAlt,1),:); % NAlt x NCT*NP*NClass
P = reshape(sum(YY .*U_probs,1),EstimOpt.NCT,EstimOpt.NP*EstimOpt.NClass);... % NCT x NP*NClass
P(reshape(MissingInd(1:EstimOpt.NAlt:end),EstimOpt.NCT,EstimOpt.NP*EstimOpt.NClass)==1) = 1;... % do not include choice tasks which were not completed
probs = prod(P,1); ...
probs = reshape(probs,EstimOpt.NP,EstimOpt.NClass); ...
if sum(EstimOpt.BActiveClass == 0,1) == 0
    PClass = exp(Xc*reshape([B(EstimOpt.NClass*EstimOpt.NVarA+1:end);zeros(EstimOpt.NVarC,1)],EstimOpt.NVarC, EstimOpt.NClass));...
else
    PClass = exp(Xc*reshape([B((EstimOpt.NClass-1)*sum(EstimOpt.BActiveClass,1)+EstimOpt.NVarA+1:end);zeros(EstimOpt.NVarC,1)],EstimOpt.NVarC,EstimOpt.NClass));...
end
PClass_sum = sum(PClass,2);...
PClass = PClass./PClass_sum(:,ones(EstimOpt.NClass,1)); ...

probs_x = probs.*PClass;
f = max(sum(probs_x,2), realmin);
% f = sum(probs_x,2);
% f = -log(max(sum(probs.*PClass,2),realmin));

if nargout == 2 % function value + gradient

    U_probs = reshape(U_probs, EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP, 1, EstimOpt.NClass);
    probs_x = reshape(probs_x, EstimOpt.NP, 1, EstimOpt.NClass);
    if EstimOpt.WTP_space == 0
       
       U_probs = U_probs(:,:, ones(1,EstimOpt.NVarA),:);
       
       XXa = reshape(Xa(:,:, ones(EstimOpt.NClass,1)), EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP,EstimOpt.NVarA, EstimOpt.NClass);
       Xhat = reshape(sum(U_probs.*XXa,1), EstimOpt.NCT,EstimOpt.NP,EstimOpt.NVarA*EstimOpt.NClass) ; % NCT x NP x NClass*NVarA
       Xa_tmp = reshape(Xa(YY(:,1:EstimOpt.NCT*EstimOpt.NP) == 1, :, ones(EstimOpt.NClass,1)),EstimOpt.NCT,EstimOpt.NP,EstimOpt.NVarA*EstimOpt.NClass) ; % NCT*NP x NVarA
       F = squeeze(sum(Xa_tmp - Xhat,1)); % NP x NVarA*NClass
       
       probs_x = reshape(probs_x(:,ones(EstimOpt.NVarA,1),:),  EstimOpt.NP, EstimOpt.NVarA*EstimOpt.NClass);
       g1 = F.*probs_x; % gradient for choice parameters
       
    else
        % non cost variables
        alphaX = zeros(EstimOpt.NAlt*EstimOpt.NCT*EstimOpt.NP,EstimOpt.NVarA - EstimOpt.WTP_space, EstimOpt.NClass);
        for i = 1:EstimOpt.NClass %maybe this loop is unnecessary? 
            alphaX(:,:,i) = Xa(:, 1:(EstimOpt.NVarA - EstimOpt.WTP_space)).*(Bclass(EstimOpt.WTP_matrix,i*ones(EstimOpt.NAlt*EstimOpt.NCT*EstimOpt.NP,1))');
        end
     
        alphaXX = reshape(alphaX, EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP, EstimOpt.NVarA - EstimOpt.WTP_space, EstimOpt.NClass);
        Xhat1 = reshape(sum(U_probs(:,:, ones(1,EstimOpt.NVarA - EstimOpt.WTP_space),:).*alphaXX,1), EstimOpt.NCT,EstimOpt.NP,(EstimOpt.NVarA-EstimOpt.WTP_space)*EstimOpt.NClass) ; % NCT x NP x NClass*NVarA
        alphaX = reshape(alphaX(YY(:,1:EstimOpt.NCT*EstimOpt.NP) == 1, :,:), EstimOpt.NCT,EstimOpt.NP,(EstimOpt.NVarA-EstimOpt.WTP_space)*EstimOpt.NClass);
        F1 = squeeze(sum(alphaX - Xhat1,1)); % NP x (NVarA - WTP_space)*NClass
        g11 = F1.*reshape(probs_x(:,ones(1, EstimOpt.NVarA - EstimOpt.WTP_space),:),  EstimOpt.NP,(EstimOpt.NVarA-EstimOpt.WTP_space)*EstimOpt.NClass);
        
        % cost variables
        if EstimOpt.WTP_space == 1
            pX = Xa(:, EstimOpt.NVarA*ones(EstimOpt.NClass,1)) + Xa(:, 1:EstimOpt.NVarA-1)*Bclass_nonwtp(1:EstimOpt.NVarA-1,:);% NAlt*NCT*NP x NClass
            pXX = reshape(pX, EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP, EstimOpt.NClass );
            Xhat2 = squeeze(sum(squeeze(U_probs).*pXX,1)); % NCT*NP x NClass

            F2 = reshape(pX(YY(:,1:EstimOpt.NCT*EstimOpt.NP) == 1, :) - Xhat2, EstimOpt.NCT, EstimOpt.NP, EstimOpt.NClass);
            g12 = squeeze(sum(F2,1)).*squeeze(probs_x);
        else
            pX = zeros(EstimOpt.NAlt*EstimOpt.NCT*EstimOpt.NP, EstimOpt.WTP_space, EstimOpt.NClass);
            for i = 1:EstimOpt.WTP_space
                pX(:,i,:) = Xa(:, (EstimOpt.NVarA - EstimOpt.WTP_space + i)*ones(EstimOpt.NClass,1))...
                    + Xa(:, EstimOpt.WTP_matrix == EstimOpt.NVarA - EstimOpt.WTP_space + i)*Bclass_nonwtp(EstimOpt.WTP_matrix == EstimOpt.NVarA - EstimOpt.WTP_space + i,:);
            end
            pXX = reshape(pX, EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP, EstimOpt.WTP_space, EstimOpt.NClass);
            Xhat2 = squeeze(sum(U_probs(:,:, ones(EstimOpt.WTP_space,1),:).*pXX,1)); % N x WTP_space x NClass
            F2 = reshape(pX(YY(:,1:EstimOpt.NCT*EstimOpt.NP) == 1, :,:) - Xhat2, EstimOpt.NCT, EstimOpt.NP, EstimOpt.WTP_space*EstimOpt.NClass);
            probs_x = reshape(probs_x(:,ones(EstimOpt.WTP_space,1),:),  EstimOpt.NP, EstimOpt.WTP_space*EstimOpt.NClass);
            g12 = squeeze(sum(F2,1)).*probs_x;
        end
        g1 = zeros(EstimOpt.NP, EstimOpt.NVarA,EstimOpt.NClass);
        g1(:, 1:(EstimOpt.NVarA - EstimOpt.WTP_space),:) = reshape(g11,EstimOpt.NP,EstimOpt.NVarA - EstimOpt.WTP_space, EstimOpt.NClass)  ;
        g1(:, (EstimOpt.NVarA - EstimOpt.WTP_space+1):end,:) = reshape(g12,EstimOpt.NP,EstimOpt.WTP_space, EstimOpt.NClass)  ;
        g1 = reshape(g1, EstimOpt.NP, EstimOpt.NVarA*EstimOpt.NClass);
    end
   
   % class membership parameters
   PClass_tmp = reshape(PClass(:,1:EstimOpt.NClass-1), EstimOpt.NP, 1, EstimOpt.NClass-1);
   PClass_tmp = reshape(PClass_tmp(:, ones(EstimOpt.NClass-1,1),:),EstimOpt.NP,(EstimOpt.NClass-1)^2);
   PClass = reshape(PClass(:,:, ones(EstimOpt.NClass-1,1)), EstimOpt.NP, EstimOpt.NClass*(EstimOpt.NClass-1));
   indx1 = 1:EstimOpt.NClass-1;
   indx1 = indx1 + (indx1-1)*EstimOpt.NClass;
   PClass(:, indx1) = PClass(:, indx1).*(1-PClass(:, indx1));
   indx2 = 1:EstimOpt.NClass*(EstimOpt.NClass-1);
   indx2(indx1) = []; % 

   PClass(:,indx2) =-PClass(:,indx2).*PClass_tmp;
   PClass = reshape(PClass, EstimOpt.NP, EstimOpt.NClass, EstimOpt.NClass-1);
   g2 = sum(PClass.*probs(:,:, ones(EstimOpt.NClass-1,1)),2); % NP x 1 xNClass-1
   g2 = g2(:, ones(EstimOpt.NVarC,1),:).*Xc(:,:, ones(EstimOpt.NClass-1,1));
   g2 = reshape(g2, EstimOpt.NP, EstimOpt.NVarC*(EstimOpt.NClass-1)); % gradient for class probability parameters
   g = -[g1,g2];
   g = g./f(:, ones(length(B),1));
end

f = -log(f);

