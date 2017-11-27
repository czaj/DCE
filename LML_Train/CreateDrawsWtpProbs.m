
%Take draws from the parameter space. 
%These draws are held in BETAS which is size NPxNVxNDRAWS
%Then calculate the probability of each person's sequence of choices at
%each draw for that person. 
%These probabilities are held in PROBS which is size NPxNDRAWS

BETAS=randi(NGridPts,[NP,NV,NDRAWS]); %BETAS go from 1 to NGridPts
BETAS=(BETAS-1)./(NGridPts-1);   %Now BETAS go from zero to 1 inclusive
for r=1:NV
  BETAS(:,r,:)=COEF(r,1)+(COEF(r,2)-COEF(r,1)).*BETAS(:,r,:);  %Now BETAS go from lower limit to upper limit for each coefficient
end;

%Add up the non-price variables times their WTP
v=sum(bsxfun(@times,XMAT(:,IDV),BETAS(XMAT(:,1),2:end,:)),2);      %sum(NROWSx(NV-1)xNDRAWS,2) gives NROWSx1xNDRAWS
v=squeeze(v);                              %NROWSxNDRAWS
v=bsxfun(@minus,v,XMAT(:,IDPRICE)); %Subtract price
v=squeeze(BETAS(XMAT(:,1),1,:)).*v; %Multiply everything by price/scale coef
v=exp(v);
sparsematrix=bsxfun(@eq,sparse(1:NCS)',XMAT(:,2)'); %NCSxNROWS
denom=double(sparsematrix)*v;                           %NCSxNDRAWS
p=v(XMAT(:,3)==1,:)./denom;                     %NCSxNDRAWS
personid=XMAT(XMAT(:,3)==1,1);                  %NCSx1
sparsematrix=bsxfun(@eq,sparse(1:NP)',personid');    %NPxNCS
PROBS=double(sparsematrix)*log(p);                      %NPxNDRAWS
PROBS=exp(PROBS);                               %NPxNDRAWS

clear sparsematrix v denom p personid



