%Create the z variables that describe density over coefficient space
%The variables are given in an array that is NPxNDRAWSxNZ
% where NZ is the number of these variables

%ZTYPE determines the kind of variables to create
% 1=polynomial, 2= step function, 3=spline

NZ_Original=NZ;

if ZTYPE==1;  %Polynomial
    NZ=NV*PolyOrder;
    Z=zeros(NP,NDRAWS,NZ);
    count=0;
    for r=1:NV;
        hold=BETAS(:,r,:);
        hold=-1+(hold-COEF(r,1))./(COEF(r,2)-COEF(r,1));
        hold=reshape(hold,NP*NDRAWS,1);
        for k=1:PolyOrder;
            count=count+1;
            yy=legendre(k,hold');
            Z(:,:,count)=reshape(yy(1,:),[NP,NDRAWS,1]);
        end;
    end;
    if CrossCorr==1;
        NewNZ=NZ+((NV-1)*(NV-2)./2);
        count=0;
        for j=(PolyOrder+1):PolyOrder:(NV*PolyOrder);  %First coef is for price/scale coef, which is not correlated
            for i=(PolyOrder+1):PolyOrder:(j-1);
                cross=Z(:,:,j).*Z(:,:,i);
                count=count+1;
                Z(:,:,NZ+count)=cross;
            end;
        end
    NZ=size(Z,3);
    if NZ ~= NewNZ; disp('Error in NZ'); return; end;
    if NZ ~= NZ_Original; disp('Error in NZ'); return; end;
    end

elseif ZTYPE==2;  %Step functions
    NZ=(NLevels-1)*NV;  %One parameter for each level minus one by normalization
    Z=zeros(NP,NDRAWS,NZ);
    CUTS=zeros(NV,NLevels+1);
    CUTS(:,1)=COEF(:,1);
    cutsize=(COEF(:,2)-COEF(:,1))./NLevels;
    for k=1:NLevels;
     CUTS(:,k+1)=CUTS(:,1)+k.*cutsize;
    end
    count=0;
    for r=1:NV;
      for k=1:NLevels-1;
        count=count+1;
        Z(:,:,count)=(BETAS(:,r,:) >= CUTS(r,k)) & (BETAS(:,r,:) < CUTS(r,k+1));
      end;
    end;
    if CrossCorr==1;
        NewNZ=NZ+2.*(NV-1)+((NV-2)*(NV-1)./2);
        for r=1:NP
            for s=1:NDRAWS
                newbeta=reshape(BETAS(r,2:end,s),NV-1,1);  %Do not include price/scale coef
                Z(r,s,(NZ+1):(NZ+2*(NV-1)))=[newbeta ; newbeta.^2];
            end
        end
        NNZZ=size(Z,3);
        count=0;
        for j=1:(NV-1);  %First coef is for price/scale coef, which is not correlated
            for i=1:(j-1);
                cross=Z(:,:,NZ+j).*Z(:,:,NZ+i);
                count=count+1;
                Z(:,:,NNZZ+count)=cross;
            end;
        end;
        NZ=size(Z,3);
        if NZ ~= NewNZ; disp('Error in NZ'); return; end;
        if NZ ~= NZ_Original; disp('Error in NZ'); return; end;
    end
    
elseif ZTYPE==3; %Linear spline
    NZ=(NKnots+1)*NV;  %Number of parameters/z-variables is 1 startpoint+ 1 endpoint+ NKnots - 1 for normalization
    %Height at endpoint is normalized to zero.
    Z=zeros(NP,NDRAWS,NZ);
    CUTS=zeros(NV,NKnots+2); %Knots and endpoints
    CUTS(:,1)=COEF(:,1);
    cutsize=(COEF(:,2)-COEF(:,1))./(NKnots+1);
    for k=1:(NKnots+1);
     CUTS(:,k+1)=CUTS(:,1)+k.*cutsize;
    end
    count=0;
    for r=1:NV;
      for k=1:(NKnots+1);
        count=count+1;
        inseg=(BETAS(:,r,:) >= CUTS(r,k)) & (BETAS(:,r,:) < CUTS(r,k+1)); %NPx1xNDRAWS
        m=(BETAS(:,r,:)-CUTS(r,k))./(CUTS(r,k+1)-CUTS(r,k)); %NPx1xNDRAWS
        Z(:,:,count)=Z(:,:,count)+squeeze((1-m).*inseg);
        if k<(NKnots+1)
          Z(:,:,count+1)=Z(:,:,count+1)+ squeeze(m.*inseg);
        end;
      end;
    end;
    
    if CrossCorr==1;
        NewNZ=NZ+2.*(NV-1)+((NV-2)*(NV-1)./2);
        for r=1:NP
            for s=1:NDRAWS
                newbeta=reshape(BETAS(r,2:end,s),NV-1,1);  %Do not include price/scale coef
                Z(r,s,(NZ+1):(NZ+2*(NV-1)))=[newbeta ; newbeta.^2];
            end
        end
        NNZZ=size(Z,3);
        count=0;
        for j=1:(NV-1);  %First coef is for price/scale coef, which is not correlated
            for i=1:(j-1);
                cross=Z(:,:,NZ+j).*Z(:,:,NZ+i);
                count=count+1;
                Z(:,:,NNZZ+count)=cross;
            end;
        end;
        NZ=size(Z,3);
        if NZ ~= NewNZ; disp('Error in NZ'); return; end;
        if NZ ~= NZ_Original; disp('Error in NZ'); return; end;
    end
    
end;

    