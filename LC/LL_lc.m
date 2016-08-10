function [f,g] = LL_lc(YY,Xa,Xc,Xs,MissingInd,EstimOpt,B)

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
if EstimOpt.NVarS > 0
   bs = reshape(B(EstimOpt.NClass*EstimOpt.NVarA+1:EstimOpt.NClass*(EstimOpt.NVarA+EstimOpt.NVarS)), EstimOpt.NVarS, EstimOpt.NClass);
   Scale = exp(Xs*bs);
   U = exp(reshape((Xa * Bclass).*Scale,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP*EstimOpt.NClass)); ...% NAlt x NCT*NP*NClass
else
    U = exp(reshape(Xa * Bclass,EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP*EstimOpt.NClass)); ...% NAlt x NCT*NP*NClass
end

U(MissingInd==1) = 0;... % do not include alternatives which were not available
U_sum = sum(U,1);... 
U_probs =  U ./ U_sum(ones(EstimOpt.NAlt,1),:); % NAlt x NCT*NP*NClass
P = reshape(sum(YY .*U_probs,1),EstimOpt.NCT,EstimOpt.NP*EstimOpt.NClass);... % NCT x NP*NClass
P(reshape(MissingInd(1:EstimOpt.NAlt:end),EstimOpt.NCT,EstimOpt.NP*EstimOpt.NClass)==1) = 1;... % do not include choice tasks which were not completed
probs = prod(P,1); ...
probs = reshape(probs,EstimOpt.NP,EstimOpt.NClass); ...

if sum(EstimOpt.BActiveClass == 0,1) == 0
    PClass = exp(Xc*reshape([B(EstimOpt.NClass*(EstimOpt.NVarA+EstimOpt.NVarS)+1:end);zeros(EstimOpt.NVarC,1)],EstimOpt.NVarC, EstimOpt.NClass));...
else
    PClass = exp(Xc*reshape([B((EstimOpt.NClass-1)*sum(EstimOpt.BActiveClass,1)+EstimOpt.NVarA+EstimOpt.NClass*EstimOpt.NVarS+1:end);zeros(EstimOpt.NVarC,1)],EstimOpt.NVarC,EstimOpt.NClass));...
end
PClass_sum = sum(PClass,2);...
PClass = PClass./PClass_sum(:,ones(EstimOpt.NClass,1)); ...

probs_x = probs.*PClass;
% f = max(sum(probs_x,2), realmin);
f = sum(probs_x,2);
% f = -log(max(sum(probs.*PClass,2),realmin));

if nargout == 2 % function value + gradient

    U_probs = reshape(U_probs, EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP, 1, EstimOpt.NClass);
    probs_x = reshape(probs_x, EstimOpt.NP, 1, EstimOpt.NClass);
    if EstimOpt.NVarS > 0
       Scale = reshape(Scale, EstimOpt.NP*EstimOpt.NCT*EstimOpt.NAlt, 1, EstimOpt.NClass);
       
    end
    if EstimOpt.WTP_space == 0
       
       U_probs = U_probs(:,:, ones(1,EstimOpt.NVarA),:);
       if EstimOpt.NVarS > 0
           XXa = reshape(Xa(:,:, ones(EstimOpt.NClass,1)).*Scale(:, ones(1,EstimOpt.NVarA),:), EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP,EstimOpt.NVarA, EstimOpt.NClass);
           Xa_tmp = reshape(Xa(YY(:,1:EstimOpt.NCT*EstimOpt.NP) == 1, :, ones(EstimOpt.NClass,1)).*Scale(1:EstimOpt.NAlt:end, ones(1,EstimOpt.NVarA),:),EstimOpt.NCT,EstimOpt.NP,EstimOpt.NVarA*EstimOpt.NClass) ; % NCT*NP x NVarA
       else
            XXa = reshape(Xa(:,:, ones(EstimOpt.NClass,1)), EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP,EstimOpt.NVarA, EstimOpt.NClass);
            Xa_tmp = reshape(Xa(YY(:,1:EstimOpt.NCT*EstimOpt.NP) == 1, :, ones(EstimOpt.NClass,1)),EstimOpt.NCT,EstimOpt.NP,EstimOpt.NVarA*EstimOpt.NClass) ; % NCT*NP x NVarA
       end
       Xhat = reshape(sum(U_probs.*XXa,1), EstimOpt.NCT,EstimOpt.NP,EstimOpt.NVarA*EstimOpt.NClass) ; % NCT x NP x NClass*NVarA
       
       F = squeeze(sum(Xa_tmp - Xhat,1)); % NP x NVarA*NClass
       
       %probs_x = reshape(probs_x(:,ones(EstimOpt.NVarA,1),:),  EstimOpt.NP, EstimOpt.NVarA*EstimOpt.NClass);
       g1 = F.*reshape(probs_x(:,ones(EstimOpt.NVarA,1),:),  EstimOpt.NP, EstimOpt.NVarA*EstimOpt.NClass); % gradient for choice parameters
       if EstimOpt.NVarS > 0
           bss = reshape(Bclass,1,1, EstimOpt.NVarA*EstimOpt.NClass);
           Ftmp = Xa_tmp - Xhat;
           Ftmp = sum(reshape(Ftmp.*bss(ones(EstimOpt.NCT, 1), ones(EstimOpt.NP, 1),:), EstimOpt.NCT*EstimOpt.NP, EstimOpt.NVarA, EstimOpt.NClass),2);
           
           gs = reshape(Ftmp(:, ones(1, EstimOpt.NVarS),:).*Xs(1:EstimOpt.NAlt:end,:, ones(1,EstimOpt.NClass)), EstimOpt.NCT,EstimOpt.NP, EstimOpt.NVarS*EstimOpt.NClass);
           gs = squeeze(sum(gs,1)); % NP x NVarS*NClass
           gs = gs.*reshape(probs_x(:,ones(EstimOpt.NVarS,1),:),  EstimOpt.NP, EstimOpt.NVarS*EstimOpt.NClass); % gradient for choice parameters
           g1 = [g1, gs];
       end
    else
        % non cost variables
        alphaX = zeros(EstimOpt.NAlt*EstimOpt.NCT*EstimOpt.NP,EstimOpt.NVarA - EstimOpt.WTP_space, EstimOpt.NClass);
        for i = 1:EstimOpt.NClass %maybe this loop is unnecessary? 
            alphaX(:,:,i) = Xa(:, 1:(EstimOpt.NVarA - EstimOpt.WTP_space)).*(Bclass(EstimOpt.WTP_matrix,i*ones(EstimOpt.NAlt*EstimOpt.NCT*EstimOpt.NP,1))');
        end
        if EstimOpt.NVarS > 0
            alphaX = alphaX.*Scale(:, ones(EstimOpt.NVarA - EstimOpt.WTP_space,1),:);
        end
        alphaXX = reshape(alphaX, EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP, EstimOpt.NVarA - EstimOpt.WTP_space, EstimOpt.NClass);
        Xhat1 = reshape(sum(U_probs(:,:, ones(1,EstimOpt.NVarA - EstimOpt.WTP_space),:).*alphaXX,1), EstimOpt.NCT,EstimOpt.NP,(EstimOpt.NVarA-EstimOpt.WTP_space)*EstimOpt.NClass) ; % NCT x NP x NClass*NVarA
        alphaX = reshape(alphaX(YY(:,1:EstimOpt.NCT*EstimOpt.NP) == 1, :,:), EstimOpt.NCT,EstimOpt.NP,(EstimOpt.NVarA-EstimOpt.WTP_space)*EstimOpt.NClass);
        F1 = squeeze(sum(alphaX - Xhat1,1)); % NP x (NVarA - WTP_space)*NClass
        g11 = F1.*reshape(probs_x(:,ones(1, EstimOpt.NVarA - EstimOpt.WTP_space),:),  EstimOpt.NP,(EstimOpt.NVarA-EstimOpt.WTP_space)*EstimOpt.NClass);
        
        % cost variables
        if EstimOpt.WTP_space == 1
            pX = Xa(:, EstimOpt.NVarA*ones(EstimOpt.NClass,1)) + Xa(:, 1:EstimOpt.NVarA-1)*Bclass_nonwtp(1:EstimOpt.NVarA-1,:);% NAlt*NCT*NP x NClass
            if EstimOpt.NVarS > 0
               pX = pX.*squeeze(Scale); 
            end
            pXX = reshape(pX, EstimOpt.NAlt,EstimOpt.NCT*EstimOpt.NP, EstimOpt.NClass );
            Xhat2 = squeeze(sum(squeeze(U_probs).*pXX,1)); % NCT*NP x NClass

            F2 = reshape(pX(YY(:,1:EstimOpt.NCT*EstimOpt.NP) == 1, :) - Xhat2, EstimOpt.NCT, EstimOpt.NP, EstimOpt.NClass);
            if EstimOpt.NVarS > 0
                bss = reshape(Bclass(EstimOpt.NVarA,:),1,1, EstimOpt.NClass);
                gs = reshape(F2.*bss(ones(EstimOpt.NCT,1), ones(1, EstimOpt.NP), :),EstimOpt.NCT, EstimOpt.NP,1, EstimOpt.NClass);
                gs = gs(:,:, ones(1,EstimOpt.NVarS),:).*reshape(Xs(1:EstimOpt.NAlt:end,:, ones(EstimOpt.NClass,1)), EstimOpt.NCT, EstimOpt.NP, EstimOpt.NVarS,EstimOpt.NClass);
                gs = reshape(sum(gs,1),EstimOpt.NP, EstimOpt.NVarS, EstimOpt.NClass).*probs_x(:,ones(1, EstimOpt.NVarS),:);
            end
            g12 = squeeze(sum(F2,1)).*squeeze(probs_x);
        else
            pX = zeros(EstimOpt.NAlt*EstimOpt.NCT*EstimOpt.NP, EstimOpt.WTP_space, EstimOpt.NClass);
            for i = 1:EstimOpt.WTP_space
                pX(:,i,:) = Xa(:, (EstimOpt.NVarA - EstimOpt.WTP_space + i)*ones(EstimOpt.NClass,1))...
                    + Xa(:, EstimOpt.WTP_matrix == EstimOpt.NVarA - EstimOpt.WTP_space + i)*Bclass_nonwtp(EstimOpt.WTP_matrix == EstimOpt.NVarA - EstimOpt.WTP_space + i,:);
            end
            if EstimOpt.NVarS > 0
               pX = pX.*Scale(:, ones(EstimOpt.WTP_space,1),:); 
            end
            pXX = reshape(pX, EstimOpt.NAlt, EstimOpt.NCT*EstimOpt.NP, EstimOpt.WTP_space, EstimOpt.NClass);
            Xhat2 = squeeze(sum(U_probs(:,:, ones(EstimOpt.WTP_space,1),:).*pXX,1)); % N x WTP_space x NClass
            F2 = reshape(pX(YY(:,1:EstimOpt.NCT*EstimOpt.NP) == 1, :,:) - Xhat2, EstimOpt.NCT, EstimOpt.NP, EstimOpt.WTP_space*EstimOpt.NClass);
            probs_x = reshape(probs_x(:,ones(EstimOpt.WTP_space,1),:),  EstimOpt.NP, EstimOpt.WTP_space*EstimOpt.NClass);
            if EstimOpt.NVarS > 0
                bss = reshape(Bclass(EstimOpt.NVarA-EstimOpt.WTP_space+1:end,:),1,1,EstimOpt.WTP_space*EstimOpt.NClass);
                gs = sum(reshape(F2.*bss(ones(EstimOpt.NCT,1), ones(1, EstimOpt.NP), :),EstimOpt.NCT, EstimOpt.NP,EstimOpt.WTP_space, EstimOpt.NClass),3); % NCT x NP x 1 x NClass
                gs = gs(:,:, ones(1,EstimOpt.NVarS),:).*reshape(Xs(1:EstimOpt.NAlt:end,:, ones(EstimOpt.NClass,1)), EstimOpt.NCT, EstimOpt.NP, EstimOpt.NVarS,EstimOpt.NClass);
                gs = reshape(sum(gs,1),EstimOpt.NP, EstimOpt.NVarS, EstimOpt.NClass).*probs_x(:,ones(1, EstimOpt.NVarS),:);
            end
            
            g12 = squeeze(sum(F2,1)).*probs_x;
        end
        
        g1 = zeros(EstimOpt.NP, EstimOpt.NVarA,EstimOpt.NClass);
        g1(:, 1:(EstimOpt.NVarA - EstimOpt.WTP_space),:) = reshape(g11,EstimOpt.NP,EstimOpt.NVarA - EstimOpt.WTP_space, EstimOpt.NClass)  ;
        g1(:, (EstimOpt.NVarA - EstimOpt.WTP_space+1):end,:) = reshape(g12,EstimOpt.NP,EstimOpt.WTP_space, EstimOpt.NClass)  ;
        g1 = reshape(g1, EstimOpt.NP, EstimOpt.NVarA*EstimOpt.NClass);
        if EstimOpt.NVarS>0
           g1 = [g1, reshape(gs, EstimOpt.NP, EstimOpt.NVarS*EstimOpt.NClass)]; 
        end
    end
   
   % class membership parameters
   PClass_tmp = reshape(PClass(:,1:EstimOpt.NClass-1), EstimOpt.NP,1,EstimOpt.NClass-1);
   PClass_tmp = reshape(PClass_tmp(:,ones(EstimOpt.NClass-1,1),:),EstimOpt.NP,(EstimOpt.NClass-1)^2);
   PClass = reshape(PClass(:,:,ones(EstimOpt.NClass-1,1)),EstimOpt.NP,EstimOpt.NClass*(EstimOpt.NClass-1));
   indx1 = 1:EstimOpt.NClass-1;
   indx1 = indx1 + (indx1-1)*EstimOpt.NClass;
   PClass(:,indx1) = PClass(:,indx1).*(1-PClass(:,indx1));
   indx2 = 1:EstimOpt.NClass*(EstimOpt.NClass-1);
   indx2(indx1) = []; % 

   PClass(:,indx2) = -PClass(:,indx2).*PClass_tmp;
   PClass = reshape(PClass, EstimOpt.NP, EstimOpt.NClass, EstimOpt.NClass-1);
   g2 = sum(PClass.*probs(:,:, ones(EstimOpt.NClass-1,1)),2); % NP x 1 xNClass-1
   g2 = g2(:,ones(EstimOpt.NVarC,1),:).*Xc(:,:, ones(EstimOpt.NClass-1,1));
   g2 = reshape(g2,EstimOpt.NP, EstimOpt.NVarC*(EstimOpt.NClass-1)); % gradient for class probability parameters
   g = -[g1,g2];
   g = g./f(:, ones(length(B),1));
end

f = -log(f);

