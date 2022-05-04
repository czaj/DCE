%-------------------------------------------------------------------------------------------------------------------------------------
%   Author              :   Ian Gregory
%   Date                :   25th January 2009
%   Functionality       :   To estimate ordered probit of the form y=Xb+e   OR  y=C+Xb+e.
%   
%   NOTE:                   There is an accompanying .m file: "Test_MATLAB_Ordered_Probit_Estimate" for testing and demonstrating
%                           how the function "MATLAB_Ordered_Probit_Estimate1" works.
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%           *Usage:*
%                               RESULTS=MATLAB_Ordered_Probit_Estimate1(INPUTS)
%                           
%           Where INPUTS are of the following form:    
%           
%   
%                                            =================        INPUTS:           =================
%            ---------------         *REQUIRED:*           ---------------        
%
%            INPUTS.DATA            :       Is a comma separated raw data file excluding column headings for estimating an ordered probit.
%                                   :       Column   1:  The y-response values (dependent)
%                                   :       Columns >1:  Are the raw explanatory variables  (X).  (independent variables)
%                                           Eg.
%                                           2,    0.8644,  0.57501
%                                           1,  0.094203,  -0.86613
%                                           1  ,-0.85191,  -2.1165
%
%            ---------------         *OPTIONAL:*           ---------------        
%           INPUTS.Display_Output_Switch  : Used to indicate if should output the results to MATLABs standard output (to the screen).
%                                           Take values 0 or 1.  Default is 1.
%                                           Eg.  
%                                           INPUTS.Display_Output_Switch=1;         % Display results output.
%
%           INPUTS.Confidence_Interval    : Confidence interval for the coefficients b.  Default is 97.5%.
%                                           Eg.
%                                           INPUTS.Confidence_Interval=0.975;
% 
%           INPUTS.INDC_SE                : Indicator variable to choose which method to calculate the standard error.  
%                                           Values are 0,1,2.  
%                                           Eg.
%                                           INDC_SE=0;          % Variance using inverted Hessian (default)
%                                           INDC_SE=1;          % Variance using the Outer Product of Gradient 
%                                           INDC_SE=2;          % Robust Variance % Covariance (VCOV) 
% 
% 
%           % Optimisation Settings:
%           INPUTS.Optimiser_Settings.Convergence.End_Criteria      :    This is the change in covergence tolerance at termination of optimisation.
%                                                                        Generally a value between 0.01 and 0.0001.  Default is 0.0001.
%                                                                        Eg. 
%                                                                        INPUTS.Optimiser_Settings.Convergence.End_Criteria=0.000001;                       
% 
%           INPUTS.Optimiser_Settings.Step_Size                     :   This is the step size in perturbing the direction of the derivative in the optimisation.
%                                                                       Generally a value between .1 and 1, default is 0.5
%                                                                       Eg. 
%                                                                       INPUTS.Optimiser_Settings.Step_Size=0.5;                                           
% 
%           INPUTS.Optimiser_Settings.Max_Iterations                :   This is the number of iterations to calculate in the optimisier.
%                                                                       Default is 20,000
%                                                                       Eg.
%                                                                       INPUTS.Optimiser_Settings.Max_Iterations=20000;
% 
%           INPUTS.Optimiser_Settings.h                             :   Perturbation for first derivative
%                                                                       Default is .000002.
%                                                                       Eg.
%                                                                       INPUTS.Optimiser_Settings.h=.000002;
% 
%           INPUTS.Optimiser_Settings.dh                            :   This is the perturbation for 2nd derivative
%                                                                       Default value is .0002;    
%                                                                       Eg.
%                                                                       INPUTS.Optimiser_Settings.dh=.0002;
% 
%           INPUTS.Optimiser_Settings.sw                                This is a switch to change the optimisation from BHHH to Hessian
%                                                                       for direction vector to the Hessian
%                                                                       Default is 15.    
%                                                                       Eg.
%                                                                       INPUTS.Optimiser_Settings.sw=15;       
%  
%                                            =================        OUTPUTS:           =================
%           RESULTS.Convergence             :  0 or 1, indicating successful convergence of the optimisation routine.
%           RESULTS.Beta                    :  Vector of estimated coefficients
%           RESULTS.Cut_Points              :  Vector of cut-points
%           RESULTS.Likelihood              :  The likelihood value returned from the optimisation routine.
%           RESULTS.First_Derivative        :  The value of the first derivatives returned from the optimisation routine.
%           RESULTS.Hessian                 :  The value of the second derivatives returned from the optimisation routine.
%           RESULTS.y                       :  The data used for the response varaible.  
%           RESULTS.X                       :  The data used for the explanatory varaibles.  
%           RESULTS.t_value                 :  The t statistic for the estimated values.
%           RESULTS.Standard_Error          :  The standrd error of the estimated values.
%           RESULTS.Optimiser_Settings      :  Settings used in the optimisation.
%------------------------------------------------------------------------------------

function RESULTS=MATLAB_Ordered_Probit_Estimate1(INPUTS)

%------------------------------------------------------------------------------------------
% BEGIN OF CHECKING INPUTS:
    switch nargin
        case 1
        % The amount of inputs is correct.
        otherwise
            error('MATLAB_Ordered_Probit_Estimate1:TooManyInputs' , ' Too many inputs specified.');
    end
    

    % Declaring switches to 0 for user specifying values.   These will be set to 1 if the user specifies the value.    
    User_Specified.INPUTS.Display_Output_Switch=0;    
    User_Specified.INPUTS.Confidence_Interval=0;
    User_Specified.INPUTS.INDC_SE=0;
    User_Specified.INPUTS.DATA=0;    
    User_Specified.INPUTS.Optimiser_Settings.Convergence.End_Criteria=0;
    User_Specified.INPUTS.Optimiser_Settings.Step_Size=0;
    User_Specified.INPUTS.Optimiser_Settings.Max_Iterations=0;
    User_Specified.INPUTS.Optimiser_Settings.h=0;
    User_Specified.INPUTS.Optimiser_Settings.dh=0;
    User_Specified.INPUTS.Optimiser_Settings.sw=0;

    
    Field_Names=fieldnames(INPUTS);
    for i=1:length(Field_Names)
        
        if strcmp(Field_Names{i},'DATA')
    		if isnumeric(INPUTS.DATA) && size(INPUTS.DATA,2)>1 && size(INPUTS.DATA,1)>5
            % Data is of numeric format.  ok to continue.
                User_Specified.INPUTS.DATA=1;
            else
                error('MATLAB_Ordered_Probit_Estimate1:UnspecifiedInputData' , ' Input DATA must be of numeric format and atleast two columns and greater than five rows of data.');
            end
        end
            
        if strcmp(Field_Names{i},'Display_Output_Switch')
            if (INPUTS.Display_Output_Switch==1 || INPUTS.Display_Output_Switch==0)
            % Display output switch is of the correct format.
                User_Specified.Display_Output_Switch=1;
            else
            	error('MATLAB_Ordered_Probit_Estimate1:UnspecifiedInputDisplay_Output_Switch' , ' Input ''Display_Output_Switch'' must take a value of 0 or 1.');            
            end
        end
        
        if strcmp(Field_Names{i},'Confidence_Interval')                
            if INPUTS.Confidence_Interval>=0 || INPUTS.Confidence_Interval<=0
            % confidence interval for the reported standard errors is of the correct format.     
                User_Specified.INPUTS.Confidence_Interval=1;
            else
            	error('MATLAB_Ordered_Probit_Estimate1:UnspecifiedInputConfidence_Interval' , ' Input ''Confidence_Interval'' must take a value between 0 and 1.');            
            end        
        end
        
        if strcmp(Field_Names{i},'INDC_SE')                        
            if INPUTS.INDC_SE==0 || INPUTS.INDC_SE==1 || INPUTS.INDC_SE==2
                 User_Specified.INPUTS.INDC_SE=1;
            else
            	error('MATLAB_Ordered_Probit_Estimate1:UnspecifiedInputINDC_SE' , ' Input ''INDC_SE'' must take a value of 0 or 1 or 2.');                                    
            end
        end
        
        % Handling the optimiser settings.
        if strcmp(Field_Names{i},'Optimiser_Settings')
            Optimiser_Field_Names=fieldnames(eval(['INPUTS.',Field_Names{i}]));

            for j=1:length(Optimiser_Field_Names)
                if strcmp(Field_Names{i},'End_Criteria')                                       
                    if INPUTS.Optimiser_Settings.End_Criteria>0 && INPUTS.Optimiser_Settings.End_Criteria<1
                    % Value ok.
                        User_Specified.INPUTS.Optimiser_Settings.Convergence.End_Criteria=1;
                    else
                    	error('MATLAB_Ordered_Probit_Estimate1:UnspecifiedInputEnd_Criteria' , ' The ''End_Criteria'' must take a value bigger than 0 and smaller than 1.');
                    end
                end
                if strcmp(Field_Names{i},'Step_Size')                                        
                    if INPUTS.Optimiser_Settings.Step_Size>0 && INPUTS.Optimiser_Settings.Step_Size<1
                        User_Specified.INPUTS.Optimiser_Settings.Step_Size=1;              
                    else
                        error('MATLAB_Ordered_Probit_Estimate1:UnspecifiedInputStep_Size' , ' The ''Step_Size'' must take a value bigger than 0 and smaller than 1.');
                    end
                end
                if strcmp(Field_Names{i},'Max_Iterations')                                        
                    if INPUTS.Optimiser_Settings.Max_Iterations>0
                        User_Specified.INPUTS.Optimiser_Settings.Max_Iterations=1;
                    else
                        error('MATLAB_Ordered_Probit_Estimate1:UnspecifiedInputMax_Iterations' , ' The ''Max_Iterations'' must be bigger than 0.');        
                    end
                end
                if strcmp(Field_Names{i},'h')                                        
                    if INPUTS.Optimiser_Settings.Optimiser_Settings.h>0 && INPUTS.Optimiser_Settings.Optimiser_Settings.h<1
                        User_Specified.INPUTS.Optimiser_Settings.h=1;      
                    else
                        error('MATLAB_Ordered_Probit_Estimate1:UnspecifiedInputMax_Iterations' , ' The ''h'' must be bigger than 0 and less than 1.');                
                    end
                end
                if strcmp(Field_Names{i},'dh')                                        
                    if INPUTS.Optimiser_Settings.Optimiser_Settings.dh>0
                        User_Specified.INPUTS.Optimiser_Settings.dh=1;
                    else
                        error('MATLAB_Ordered_Probit_Estimate1:UnspecifiedInputMax_Iterations' , ' The ''dh'' must be bigger than 0.');                
                    end
                end                    
                if strcmp(Field_Names{i},'sw')
                    if INPUTS.Optimiser_Settings.Optimiser_Settings.sw>0
                        User_Specified.INPUTS.Optimiser_Settings.sw=1;
                    else
                        error('MATLAB_Ordered_Probit_Estimate1:UnspecifiedInputMax_Iterations' , ' The ''sw'' must be bigger than 0.');        
                    end
                end                    
            end
        end
        
    end           
    
    if User_Specified.INPUTS.DATA==0
    % User did not specify any data.  Need to throw an error.
        error('MATLAB_Ordered_Probit_Estimate1:UnspecifiedInputData' , ' Need to specify some input DATA.');
    end
    
%----------------------------------------------------------------------------------------------------------------         
% END OF CHECKING INPUTS.
%----------------------------------------------------------------------------------------------------------------         



%----------------------------------------------------------------------------------------------------------------         
% BEGIN OF HANDLE INPUTS:
%----------------------------------------------------------------------------------------------------------------         
    % Extract the data.
    y=INPUTS.DATA(:,1);
    X=INPUTS.DATA(:,[2:end]);


    % Starting values.
    b=inv(X'*X)*X'*y;  % OLS solution.

    % Initial guesses for the cut-off.  Guesses are a random number for the number of categories-1):
    alpha=sort(randn(length(unique(y))-1,1));

    %----------------------------------------------------------------------------------------------------------------
    % Optimisation Settings:
    Optimiser_Settings.Initial_Coefficient_Guesses=[b',alpha'];                 % Initial guesses for estimates.  b are the dependent coefficients, alpha are the cut-off's.
    Optimiser_Settings.Convergence.Start_Criteria=1;
    

    if User_Specified.INPUTS.Optimiser_Settings.Convergence.End_Criteria==0;
        Optimiser_Settings.Convergence.End_Criteria=0.000001;                       % Generally a value between 0.01 and 0.0001        
    end    
    if User_Specified.INPUTS.Optimiser_Settings.Step_Size==0;
        Optimiser_Settings.Step_Size=0.5;                                           % Generally a value between .1 and 1        
    end    
    if User_Specified.INPUTS.Optimiser_Settings.Max_Iterations==0;
        Optimiser_Settings.Max_Iterations=20000;
    end
    if User_Specified.INPUTS.Optimiser_Settings.h==0;
        Optimiser_Settings.h=.000002;                                               % Perturbation for first derivative        
    end        
    if User_Specified.INPUTS.Optimiser_Settings.dh==0;
        Optimiser_Settings.dh=.0002;                                                % Perturbation for 2nd derivative        
    end        
    if User_Specified.INPUTS.Optimiser_Settings.sw==0;
        Optimiser_Settings.sw=15;       
    end        
%----------------------------------------------------------------------------------------------------------------


%----------------------------------------------------------------------------------------------------------------
    % Output variable settings:
    if User_Specified.INPUTS.Confidence_Interval==1
    % User has specifed the value
    else
        INPUTS.Confidence_Interval=0.975;       % The confidence bounds to specify around the estimate.
                                                % 0.975 ie. 1.96 (default).
    end
    
    if User_Specified.INPUTS.INDC_SE==1
        
    else
    % Standard error calculation:
        INPUTS.INDC_SE=0;% Variance using inverted Hessian (default)
%     INDC_SE=1;% Variance using the Outer Product of Gradient 
%     INDC_SE=2;% Robust Variance % Covariance (VCOV) 
    end

    % Set the default values for variables the user has not specified:
     if User_Specified.INPUTS.Display_Output_Switch==0
         INPUTS.Display_Output_Switch=1;
     end    
%----------------------------------------------------------------------------------------------------------------         
% END OF HANDLE INPUTS.
%----------------------------------------------------------------------------------------------------------------         



%----------------------------------------------------------------------------------------------------------------         
% BEGIN OF PERFORM ESTIMATION:
%----------------------------------------------------------------------------------------------------------------         
    % Obtain maximum likelihood estimates
    RESULTS=MATLAB_Ordered_Probit_MLE(y,X,Optimiser_Settings);
    
    % Run the estimation a 2nd time to 'improve' on the result.
    Optimiser_Settings.Initial_Coefficient_Guesses=[RESULTS.Beta',RESULTS.Cut_Points'];
    RESULTS=MATLAB_Ordered_Probit_MLE(y,X,Optimiser_Settings);    
    
    RESULTS.y=y;
    RESULTS.X=X;
%----------------------------------------------------------------------------------------------------------------         
% END OF PERFORM ESTIMATION.
%----------------------------------------------------------------------------------------------------------------         





%----------------------------------------------------------------------------------------------------------------         
% BEGIN OF OUTPUT RESULTS:
%----------------------------------------------------------------------------------------------------------------         


    % Standard error calculation:
    switch INPUTS.INDC_SE        
        case 0
        %%%%%%%%%%%%%%%%%%%%%%%%%% Variance using inverted Hessian %%%%%%%%%%%%%%%%%%%%%%%%%% 
            se = sqrt(diag(inv(-1*RESULTS.Hessian)));
            Heading_Output='                               MLE Output using inverted Hessian';
        case 1
             H1=inv(RESULTS.First_Derivative'*RESULTS.First_Derivative); 
             d1 = diag(H1);
             se = sqrt(d1);
            Heading_Output='                               MLE Output using Outer Product Gradient';

        case 2
            H1=inv(-RESULTS.Hessian)*(RESULTS.First_Derivative'*RESULTS.First_Derivative)*inv(-RESULTS.Hessian); 
            d1 = diag(H1);
            se = sqrt(d1);                    
            Heading_Output='                               MLE Output using robust Variance Covariance (VCOV)';
        otherwise
            disp('Unknown display type for the standard error')
    end    

    z_statistic = [RESULTS.Beta;RESULTS.Cut_Points]./se;
%     t_statistic = [RESULTS.Beta;RESULTS.Cut_Points]./se;

    Estimated_Coefficient=[RESULTS.Beta;RESULTS.Cut_Points];

    Confidence_Interval_Lower=Estimated_Coefficient-norminv(INPUTS.Confidence_Interval)*se;
    Confidence_Interval_Upper=Estimated_Coefficient+norminv(INPUTS.Confidence_Interval)*se;
    
    if INPUTS.Display_Output_Switch==1
    % User would like to see some output.
    
        %%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT RESULTS: %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Column heading for reporting results.
    	fprintf('                                  Ordered Probit Results Output:\n')
        fprintf('LLH value at termination:\t%6.4f\n',RESULTS.Likelihood.LLV);
        Row_Name_Count=1;
        for i=1:length(RESULTS.Beta)
            Row_Heading{Row_Name_Count}=['independent_variable_',num2str(i-1)];  % Stating the 'true' value in the column heading 
            Row_Name_Count=Row_Name_Count+1;
        end
        for i=1:length(RESULTS.Cut_Points)
            Row_Heading{Row_Name_Count}=['Cut_Point_',num2str(i-1)];  % Stating the 'true' value in the column heading    
            Row_Name_Count=Row_Name_Count+1;    
        end
    
        fprintf('----------------------------------------------------------------------------------------------------------------------')
        fprintf('\n');
        fprintf(Heading_Output); 
        % NOTE: These error calculations are the same in STATA 9.2 command 'oprobit' uses for its standard output for after estimation.
        fprintf('\n');
        fprintf('%23s','variable','coeff','s.e.','t-value');
        fprintf('\n');
        for i=1:Row_Name_Count-1    
            fprintf('%23s',Row_Heading{i}); %print ith variable name    
            fprintf('%23.4f',[Estimated_Coefficient(i,:) se(i,:) z_statistic(i,:)     ]);%print results for variable i        
            fprintf('\n');
        end
        fprintf('----------------------------------------------------------------------------------------------------------------------')
        fprintf('\n');
    end    
    RESULTS.t_value=z_statistic;                           % Adding the t-beta values to the results structure.
    RESULTS.Standard_Error=se;                             % Adding the beta standard errors to the results structure
    RESULTS.Optimiser_Settings=Optimiser_Settings;         % Adding the settings used for the optimiser to the results structure
%----------------------------------------------------------------------------------------------------------------         
% END OF OUTPUT RESULTS.
%----------------------------------------------------------------------------------------------------------------         

    % Clear all redundant variables.
    clear('Confidence_Interval','Confidence_Interval_Lower','Confidence_Interval_Upper','DATA','Estimated_Coefficient', ...
            'Heading_Output','INDC_SE','Optimiser_Settings','Row_Heading','Row_Name_Count','X','alpha','b','i','se','y','z_statistic')
end