function [OUTPUT]=MATLAB_Ordered_Probit_MLE(y,X,Optimiser_Settings);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%   Author: Ian Gregory
%
%   Functionality:  To obtain maximum likelihood estimates for basic ordered probit.
%
%   INPUTS:
%       y                                     :       categorical responses
%       X                                     :       Independent variables in columns.
%       start                                 :       starting vector of parameters (theta0 in our notes)
%       Convergence.Start_Criteria            :       Starting value for convergence criterion
%       Convergence.End_Criteria              :       End value for convergence criterion
%       Step_Size                             :       Step size
%       Max_Iterations                        :       Maximum number of allowed iterations
%       h                                     :       Perturbation for first derivative
%       dh                                    :       Perturbation for 2nd derivative
%       sw                                    :       Number of iterations beyond which we want to switch from BHHH to Hessian
%                                                     for direction vector to the Hessian
%
%    OUTPUTS:
%       Beta                :       the MLE solution vector
%       LLV                 :       Loglikelihood value at optimisation termination 
%       L                   :       full n by 1 vector of individual log-lh's (li's) at convergence 
%       Cut_Points          :       
%       First_Derivative    :       n by k matrix of (transposed) individual gradients (gi's)
%       Hessian             :       Hessian matrix at optimisation termination
%
%   DEPENDENCIES:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    %Data size:
    n=size(X,1);
    k=length(Optimiser_Settings.Initial_Coefficient_Guesses);
    k1=size(X,2);
    
    %Optimizer settings:
    Initial_Coefficient_Guesses=Optimiser_Settings.Initial_Coefficient_Guesses;    
		Start_Criteria=Optimiser_Settings.Convergence.Start_Criteria;  
    End_Criteria=Optimiser_Settings.Convergence.End_Criteria;  
    Step_Size=Optimiser_Settings.Step_Size;
    Max_Iterations=Optimiser_Settings.Max_Iterations;
    
    h=Optimiser_Settings.h;
    dh=Optimiser_Settings.dh;
    dh2=dh^2;    
    der=zeros(n,k);
    sw=Optimiser_Settings.sw;
    

	b=Initial_Coefficient_Guesses';                   % Beta coefficients    
    jj=0;
    while (Start_Criteria>End_Criteria)&(jj<Max_Iterations);
        jj=jj+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L = MATLAB_Ordered_Probit_Likelihood(b,y,X,k1);      % call your problem-specific function that 
                                        % computes the n by 1 vector of individual log-likelihood components
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        v0=L; 
        b0=b;
        %Following comutes numerical first derivatives of log likelihood
        for i = 1:k
            b(i)=b(i)-h;
            L = MATLAB_Ordered_Probit_Likelihood(b,y,X,k1);                                       % Ordered Probit likelihood function.
            vn = L; 
            b=b0;
            b(i)=b(i)+h; 
            L = MATLAB_Ordered_Probit_Likelihood(b,y,X,k1);                                       % Ordered Probit likelihood function.
            vp=L; 
            b=b0;
            der(:,i)=.5*(vp-vn)/h;
        end;
        % Following computes numerical second derivatives of log likelihood
        hess=zeros(k,k);   
        i=1; 
        j=1;  
        while j<=k;
            while i<=j;
                if i==j;
                    b(i)=b(i)+dh; 
                    L = MATLAB_Ordered_Probit_Likelihood(b,y,X,k1);                                       % Ordered Probit likelihood function.
                    vp=L; 
                    b=b0;
                    b(i)=b(i)-dh; 
                    L = MATLAB_Ordered_Probit_Likelihood(b,y,X,k1);                                       % Ordered Probit likelihood function.
                    vn=L; 
                    b=b0; 
                    hess(i,j)=sum((vp+vn-2*v0)/dh2);
                else; 
                    c=[i j];
                    b(c)=b(c)+dh;  
                    L = MATLAB_Ordered_Probit_Likelihood(b,y,X,k1);                                       % Ordered Probit likelihood function.
                    vpp=L; 
                    b=b0; 
                    b(c)=b(c)-dh; 
                    L = MATLAB_Ordered_Probit_Likelihood(b,y,X,k1);                                       % Ordered Probit likelihood function.
                    vnn=L; 
                    b=b0;
                    b(c)=b(c)+[dh;-dh]; 
                    L = MATLAB_Ordered_Probit_Likelihood(b,y,X,k1);                                       % Ordered Probit likelihood function.
                    vpn=L; 
                    b=b0; 
                    b(c)=b(c)+[-dh;dh]; 
                    L = MATLAB_Ordered_Probit_Likelihood(b,y,X,k1);                                       % Ordered Probit likelihood function.
                    vnp=L; 
                    b=b0;
                    hess(i,j)=.25*sum(vpp+vnn-vpn-vnp)/dh2;
                end
                i=i+1;
            end
            i=1; 
            j=j+1;
        end
        hess=eye(k).*hess+(1-eye(k)).*(hess+hess');
        Start_Criteria=sum(abs(sum(der)));

        % For the direction vector.  It is better to use BHHH algorithm for the beginning iterations as it is generally believed the Hessian reaction is not stable.  
        dbBHHH=inv(der'*der)*(sum(der)');
        dbhess=inv(-hess)*(sum(der)');

        %At "sw" iterations stopping BHHH and using Hessian updates.
        db=(1-min(1,max(0,jj-sw)))*dbBHHH + min(1,max(0,jj-sw))*dbhess;

        % Update coefficients by step-size db
        b=b+Step_Size*db; 
        
        iter=[jj sum(L) Start_Criteria];

        if jj==Max_Iterations
            disp('maximum number of optimisation iterations reached');
        end
    end

    % Return a Boolean for convergence.
    if jj<Max_Iterations
        Results.Convergence=1;
    else
        Results.Convergence=0;
    end
    
    Results.Beta=b(1:k1);                 % Independent variable coefficients.
    Results.Cut_Points=b(k1+1:end);       % Cut-off coefficients.
    Results.Likelihood.LLV=sum(L);        % Log likelihood value at termination.
    Results.Likelihood.LLVs=L;            % Individual likelihood values.
    Results.First_Derivative=der;         % Derivatives at termination.
    Results.Hessian=hess;                 % Hessian at termination.

    OUTPUT=Results;
    
end
