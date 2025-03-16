function Bhat = EM_LC(Y,Xa,Xc,Xs, MissingInd, EstimOpt,OptimOpt,B)

global B_backup;

NP = EstimOpt.NP;
NAlt = EstimOpt.NAlt;
NCT = EstimOpt.NCT;
NClass = EstimOpt.NClass;
NVarA = EstimOpt.NVarA;
% NVarC = EstimOpt.NVarC;
NVarS = EstimOpt.NVarS;
% WTP_space = EstimOpt.WTP_space;
% WTP_matrix = EstimOpt.WTP_matrix;

YY = reshape(Y,[NAlt,NCT*NP]);
YY = YY(:,(1:size(YY,2))'*ones(1,NClass)); %NAlt x NCT*NP*NClass
YY(isnan(YY)) = 0;
Xm = zeros(size(Y,1),0);
EstimOpt.NVarM = 0;
EstimOpt.NVarNLT = 0;
EstimOpt.SCEXP = 1;
EstimOpt.ExpB = [];
OptimOpt.Display = 'off';
OptimOpt.OutputFcn = [];

b_diff = 1;
LL_diff = 1; 
TolFun = OptimOpt.TolFun;
TolX = OptimOpt.TolX;


Bclass = reshape(B(1:NClass*NVarA),[NVarA,NClass]);
Bs = reshape(B(NClass*NVarA+1:NClass*(NVarA+NVarS)),[NVarS,NClass]);
Bmnl = [Bclass; Bs];
Theta = B(NClass*(NVarA+NVarS)+1:end);

LL = LL_lc(YY,Xa,Xc,Xs,MissingInd,EstimOpt,B);
LL = -sum(LL);
iter = 1;
while (b_diff > TolX) && (LL_diff > TolFun)
    
    % Calculating posterior probability
    Eta = BayesProbs(YY,Xa,Xc, Xs, MissingInd,EstimOpt,B); % NP x NClass

    % Estimating MNLs for Betas
    Bmnl_new = zeros(NVarA + NVarS,NClass);
    Eta_mnl = reshape(permute(Eta(:,:, ones(1, NCT)), [3, 1 ,2]), [NCT*NP, NClass]);
    for i = 1:NClass
        LLfun = @(Beta) LL_mnl_MATlike(Y,Xa,Xm,Xs,Eta_mnl(:,i),EstimOpt,OptimOpt,Beta);
        [bi] = fminunc(LLfun,Bmnl(:,i),OptimOpt);
        Bmnl_new(:,i) = bi;
    end
    % Estimating FMNLs for Thetas
    LLfun = @(Beta) LL_fmnl_MATlike(Eta, Xc, EstimOpt, OptimOpt, Beta);
    Theta_new = fminunc(LLfun,Theta,OptimOpt);
    
    % Updating coefficients and Calculating LL
    b_diff = sum(sum(abs(Bmnl_new - Bmnl),1),2) + sum(abs(Theta_new -Theta),1);
    Bmnl = Bmnl_new;
    Theta = Theta_new;

    Bclass = Bmnl(1:NVarA,:);
    Bs = Bmnl(NVarA+1:end,:);
    B = [Bclass(:); Bs(:); Theta(:)];
    B_backup = B;
    LLnew = LL_lc(YY,Xa,Xc,Xs,MissingInd,EstimOpt,B);
    LLnew = -sum(LLnew);
    LL_diff = abs(LLnew - LL);
    LL = LLnew;

    disp(num2str([iter, b_diff, LL_diff, LL],'Iteration: %1.0f | Parameters change: %1.2f | LL change: %1.2f | LL: %1.3f '))
    iter = iter+1;
end
disp(' ')
Bclass = Bmnl(1:NVarA,:);
Bs = Bmnl(NVarA+1:end,:);
Bhat = [Bclass(:); Bs(:); Theta(:)];

end

function LL = LL_fmnl(Eta, Xc, EstimOpt, B)

    NClass = EstimOpt.NClass;
%     NVarA = EstimOpt.NVarA;
    NVarC = EstimOpt.NVarC;
%     NVarS = EstimOpt.NVarS;

    PClass = exp(Xc*reshape([B;zeros(NVarC,1)],[NVarC,NClass]));
    PClass = PClass./sum(PClass,2); % NP x NClass
    LL = sum(log(PClass).*Eta,2);
end

function [f,g,h] = LL_fmnl_MATlike(Eta, Xc, EstimOpt, OptimOpt, b0)
    
    LLfun = @(B) LL_fmnl(Eta, Xc, EstimOpt, B);
    if isequal(OptimOpt.GradObj,'on')
        LL = LLfun(b0);
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            j = -numdiff(LLfun,LL,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
            g = sum(j,1)';
            h = j'*j;
        else
            j = -numdiff(LLfun,LL,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
            g = sum(j,1)';
        end
        f = -sum(LL);
    else
        EstimOpt.NumGrad = 1;
        LL = LLfun(b0);
        f = -sum(LL);   
    end
end
