
if YesGPU==1
    PROBS=gpuArray(PROBS);
    BETAS=gpuArray(BETAS);
    Z=gpuArray(Z);
end

param=StartB;

disp('Start estimation');
disp('The negative of the log-likelihood is minimized,');
disp('which is the same as maximizing the log-likelihood.');
tic;
options=optimset('LargeScale','off','Display','iter','GradObj','on',...
    'MaxFunEvals',10000,'MaxIter',MAXITERS,'TolX',PARAMTOL,'TolFun',LLTOL,'DerivativeCheck','off');
if WantHessian==1;
    [paramhat,fval,exitflag,output,grad,hessian]=fminunc(@flexll,param,options);
else
    [paramhat,fval,exitflag]=fminunc(@flexll,param,options);
end
disp(' ');
disp(['Estimation took ' num2str(toc./60) ' minutes.']);
disp(' ');
disp('Calculating summary statistics for random utility parameters.');

if YesGPU==1;
    PROBS=gather(PROBS);
    BETAS=gather(BETAS);
    Z=gather(Z);
end;
[MeanEst,StdEst,CovMatEst,FreqEst,MidEst]=stats(paramhat,NBins);
disp('Estimated coefficients of Z variables are held in paramhat.');
if WantHessian==1
    disp('Hessian at convergence is held in hessian.');
end;
disp('Means, Standard Deviations, and Covariance Matrix of Utility Parameters');
disp('are held as MeanEst, StdEst, and CovMatEst.');
disp('Share of density in each bin and midpoint for each bin are held');
disp('in FreqEst and MidEst; each is size NV x NBins');
disp('Means and StdDevs of Random Utility Paramaters');
disp(' ');
disp('                Mean     Std Dev');
disp('              -------------------');
for r=1:length(NAMES);
    fprintf('%-10s %10.4f %10.4f\n', NAMES{r,1}, [MeanEst(r,1) StdEst(r,1)]);
end
disp(' ');
disp('Correlation Matrix for WTPs');
disp(corrcov(CovMatEst(2:end,2:end)));
disp('To obtain histogram for utility parameter k, type: bar(MidEst(k,:),FreqEst(k,:))');


