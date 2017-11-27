ParamSe=zeros(NZ,1); %Standard error of estimated coefficients of Z variables
MeanSE=zeros(NV,1); %Standard error of Mean
StdSE=zeros(NV,1);  %  " of Std Devs
CorrSE=zeros(NV,NV);  % " of Correlation matrix
FreqSE=zeros(NV,NBins); % " of Freq in each bin
MidSE=zeros(NV,NBins); % " of Midpoint for each bin, (std err should be zero)

paramhold=zeros(NZ,NReps);
mnhold=zeros(NV,NReps);
stdhold=zeros(NV,NReps);
corrhold=zeros(NV-1,NV-1,NReps);
freqhold=zeros(NV,NBins,NReps);
midhold=zeros(NV,NBins,NReps);

XMAT_Original=XMAT;
PID_Original=XMAT(:,1);
CSID_Original=XMAT(:,2);
clear i r
for xi=1:NReps
    btsample=randi(NP,[NP,1]);
    XMAT=[];
        for r=1:NP;
              thisp=btsample(r,1);
              thisx=XMAT_Original(PID_Original==thisp,:);
              thisx(:,1)=r.*ones(size(thisx,1),1);
              XMAT=[XMAT ; thisx];   
        end
        XMAT(:,2)=CSID_Original;
    CreateDrawsWtpProbs
    CreateZ
    if YesGPU==1;
        PROBS=gpuArray(PROBS);
        BETAS=gpuArray(BETAS);
        Z=gpuArray(Z);
    end
    param=StartB;
    options=optimset('LargeScale','off','Display','iter','GradObj','on',...
    'MaxFunEvals',10000,'MaxIter',MAXITERS,'TolX',PARAMTOL,'TolFun',LLTOL,'DerivativeCheck','off');
    [paramhat,fval,exitflag]=fminunc(@flexll,param,options);
    if YesGPU==1;
        PROBS=gather(PROBS);
        BETAS=gather(BETAS);
        Z=gather(Z);
    end;
    paramhold(:,xi)=paramhat;
    [mnhold(:,xi),stdhold(:,xi),cc,freqhold(:,:,xi),midhold(:,:,xi)]=stats(paramhat,NBins);
    corrhold(:,:,xi)=corrcov(cc(2:end,2:end));
end
ParamSE=std(paramhold,0,2);
MeanSE=std(mnhold,0,2);
StdSE=std(stdhold,0,2);
CorrSE=std(corrhold,0,3); 
FreqSE=std(freqhold,0,3); 
MidSE=std(midhold,0,3);
XMAT=XMAT_Original;
clear XMAT_Original;

disp('Utility Parameters');
disp(' ');
disp('                    Mean                   Std Dev');
disp('              ------------------   -----------------------');
disp('                 Est     SE            Est         SE');
for r=1:length(NAMES);
    fprintf('%-10s %10.4f %10.4f %10.4f %10.4f\n', NAMES{r,1}, [MeanEst(r,1) MeanSE(r,1) StdEst(r,1) StdSE(r,1)]);
end

disp(' ');
disp('Correlations of WTPs');
disp(corrcov(CovMatEst(2:end,2:end)));
disp('T-stats on Correlations');
disp('Diagonal has Inf because diagonals of correlation matrix are always 1.');
disp(corrcov(CovMatEst(2:end,2:end))./CorrSE);
disp(' ');
disp('To obtain histogram for utility parameter k, type: bar(MidEst(k,:),FreqEst(k,:))');
