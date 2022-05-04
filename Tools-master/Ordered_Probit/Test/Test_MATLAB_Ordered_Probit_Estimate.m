function Test_MATLAB_Ordered_Probit_Estimate
%------------------------------------------------------------------------------------
%   Author: Ian Gregory
%   Date:   25th January 2009
%   Functionality:  To be a test script for the function "MATLAB_Ordered_Probit_Estimate".
%   
%   REQUIRED FILES:
%   1.) MATLAB_Ordered_Probit_Estimate.m 
%   2.) Test_MATLAB_Ordered_Probit_Estimate.csv
%   3.) MATLAB_Ordered_Probit_Likelihood.m
%   4.) MATLAB_Ordered_Probit_MLE.m
% 
%   EXPECTED OUTPUT:
%
%
%                                   Start of testing function: "MATLAB_Ordered_Probit_Estimate"
%
%                                                              Ordered Probit Results Output:
%                            LLH value at termination:	-4.6059
%                            ----------------------------------------------------------------------------------------------------------------------
%                                                                   MLE Output using inverted Hessian
%                                         variable                  coeff                   s.e.                t-value
%                           independent_variable_0                 2.1410                 1.3812                 1.5501
%                           independent_variable_1                 2.9463                 1.4006                 2.1036
%                                      Cut_Point_0                 0.1580                 0.5582                 0.2832
%                           ----------------------------------------------------------------------------------------------------------------------
% 
%                            ans = 
% 
%                                       Convergence: 1
%                                              Beta: [2x1 double]
%                                        Cut_Points: 0.1580
%                                        Likelihood: [1x1 struct]
%                                  First_Derivative: [20x3 double]
%                                           Hessian: [3x3 double]
%                                                 y: [20x1 double]
%                                                 X: [20x2 double]
% 
%
%                                End of testing function: "MATLAB_Ordered_Probit_Estimate"
%------------------------------------------------------------------------------------

    clc                                                                         % Clear the MATLAB output screen.
    disp('Start of testing function: "MATLAB_Ordered_Probit_Estimate"')

    % REQUIRED INPUT:
    % Open DATA file 
    DATA=csvread('Test_MATLAB_Ordered_Probit_Estimate.csv',1,0);       % First row is column headings.   csvread uses base 0.
    
    % OPTIONAL INPUTS:
    INPUTS.DATA=DATA;
    INPUTS.Display_Output_Switch=1;
    INPUTS.Optimiser_Settings.Step_Size=0.1;     
    INPUTS.Confidence_Interval=0.975;
    INPUTS.INDC_SE=0;
    INPUTS.Optimiser_Settings.Convergence.End_Criteria=0.0001;
    INPUTS.Optimiser_Settings.Step_Size=0.5;
    INPUTS.Optimiser_Settings.Max_Iterations=20000;
    INPUTS.Optimiser_Settings.h=.000002;
    INPUTS.Optimiser_Settings.dh=.0002;
    INPUTS.Optimiser_Settings.sw=15;

    % Call the function.
    MATLAB_Ordered_Probit_Estimate(INPUTS);
    
    
    disp('End of testing function: "MATLAB_Ordered_Probit_Estimate"')
end
















