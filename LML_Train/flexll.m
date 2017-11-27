function [ll,g]=flexll(alpha)

global PROBS Z NZ NDRAWS YesGPU

% PROBS: NPxNDRAWS matrix of probabilities for each person at each sampled coefficients
% Z:     NPxNDRAWSxNZ array of variables that explain distribution of coefficients
% NZ: scalar number of variables in Z
% NDRAWS: scalar number of draws in PROBS

% Input alpha is NZx1 vector of coefficients

w=sum(bsxfun(@times,Z,reshape(alpha,[1,1,NZ])),3);  %NPxNDRAWS
w(w<-500)=-500;  %As precaution against extreme parameters
w(w>500)=500;
w=squeeze(exp(w));
w=w./repmat(sum(w,2),1,NDRAWS); %NPxNDRAWS
logit_w=w.*PROBS; %NPxNDRAWS
mix_probs=sum(logit_w,2); %NPx1
ll=sum(log(mix_probs),1); %1x1
ll=-ll; %To minimize
if YesGPU==1
   ll=gather(ll); 
end;

if nargout>1;
  condw=bsxfun(@times,logit_w,1./mix_probs); %NPxNDRAWS
  g=bsxfun(@times,(condw-w),Z); %NPxNDRAWSxNZ
  g=sum(sum(g,1),2); %1x1xNZ
  g=reshape(g,NZ,1); %NZx1
  g=-g; %To minimize
  if YesGPU==1
     g=gather(g); 
  end
end;
  
  
