function Test_MATLAB_Ordered_Probit_Simulate1
%------------------------------------------------------------------------------------
%   Author: Ian Gregory
%   Date:   25th January 2009
%   Functionality:  To be a test script for the function "MATLAB_Ordered_Probit_Simulate1".
%   
%   REQUIRED FILES:
%   1.) MATLAB_Ordered_Probit_Simulate.m 
% 
%   EXPECTED OUTPUT:
%
%               SCREEN OUTPUT:
%                                   Start of testing function: "MATLAB_Ordered_Probit_Simulate"
%
%                                               Finished simulating ordered probit data.
% 
%                                   End of testing function: "MATLAB_Ordered_Probit_Estimate"  
% 
%               FILE OUTPUT:
% 
%------------------------------------------------------------------------------------

    clc                                                                         % Clear the MATLAB output screen.
    disp('Start of testing function: "MATLAB_Ordered_Probit_Simulate"')

    % REQUIRED INPUT:

   INPUTS.Beta=[0.5,2,3];
   INPUTS.No_Data_Points=1000;
   INPUTS.Cut_Points=[-0.5,0.5];
   
    % OPTIONAL INPUTS
   INPUTS.Constant_Value=3;
   INPUTS.Include_Const_Switch=1;
   INPUTS.fNameOutput='My_Sim_OProbit_DATA.csv';
   
   % Call the function.   
    MATLAB_Ordered_Probit_Simulate(INPUTS)

   disp('End of testing function: "MATLAB_Ordered_Probit_Simulate"')
end
















