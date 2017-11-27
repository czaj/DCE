function [mn,stdv,cc,freqbins,midbins]=stats(alpha,NBins)

global NP NV Z NZ NDRAWS BETAS CrossCorr COEF

w=createw(alpha); %NPxNDRAWS
mn=zeros(NV,1);
v=zeros(NV,1);
cc=zeros(NV,NV);
freqbins=zeros(NV,NBins);
midbins=zeros(NV,NBins);
wtsum=sum(sum(w,1),2);
ww=reshape(w,NP*NDRAWS,1);
for r=1:NV
    thiscoef=squeeze(BETAS(:,r,:)); %NPxNDRAWS
    mn(r,1)=sum(sum(thiscoef.*w,1),2)./wtsum;
    v(r,1)=sum(sum(((thiscoef-mn(r,1)).^2).*w,1),2)./wtsum;
    cc(r,r)=v(r,1);
    if r>1 & CrossCorr==1;
      for thisc=2:NV;
         thiscoef2=squeeze(BETAS(:,thisc,:)); %NPxNDRAWS
         thismn2=sum(sum(thiscoef2.*w,1),2)./wtsum;
         cc(r,thisc)=sum(sum(((thiscoef-mn(r,1)).*(thiscoef2-thismn2)).*w,1),2)./wtsum;
         cc(thisc,r)=cc(r,thisc);
      end
    end
    [histw, midpt] = histwgt(reshape(thiscoef,NP*NDRAWS,1), ww, COEF(r,1), COEF(r,2),NBins); 
    freqbins(r,:)=histw./sum(histw);
    midbins(r,:)=midpt;
end
stdv=sqrt(v);