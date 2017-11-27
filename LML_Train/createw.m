function w=createw(alpha)

global Z NZ NDRAWS

% Z: NPxNDRAWSxNZ array of variables that explain distribution of coefficients
% NZ: scalar number of variables in Z
% NDRAWS: scalar number of draws in PROBS

% Input alpha is NZx1 vector of coefficients

w=sum(bsxfun(@times,Z,reshape(alpha,[1,1,NZ])),3);  %NPxNDRAWS
w(w>500)=500;  %As a precaution to prevent machine infty when exponentiate
w=squeeze(exp(w));
w=w./repmat(sum(w,2),1,NDRAWS); %NPxNDRAWS