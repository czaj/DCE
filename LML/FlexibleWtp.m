% Matlab code to estimate a flexible mixed logit model in wtp space.
% Written by Kenneth Train, first version Oct 21, 2015, 

clear all

% Do not change the next line. It sets the global variables.
global NP PROBS Z NZ NDRAWS NV BETAS ZTYPE CrossCorr YesGPU COEF

% OUTPUT FILE
% Put the name you want for your output file (including full path if not the current 
% working directory) after words "delete" and "diary".
% The 'diary off' and 'delete filename' commands close and delete the previous version 
% of the file created during your current matlab session or in any previous sessions. 
% If you want to append the new output to the old output, then 
% put % in front of the 'diary off' and 'delete filename' commands (or erase them).

diary off
delete MyResults.out
diary MyResults.out

% TITLE
% Put a title for the run in the quotes below, to be printed at the top of the output file.
disp 'Flexible mixed logit.'

% DATA

% Number of people (decision-makers) in dataset 
NP=100; 


% Number of choice situations in dataset. This is the number faced by all the people combined.
% There are 11 choice situations per person in the example dataset.
% You do not need to have the same number of choice situations for each person. 
NCS=NP*11;   

% Total number of alternatives faced by all people in all choice situations combined.
% This is the number of rows of data in XMAT below.
% In the example data, there are 5 alternatives in each choice situation.
% You do not need to have the same number of alternatives in each choice situation.
NROWS=NCS*5;

% Load and/or create XMAT, a matrix that contains the data.
%
% XMAT must contain one row of data for each alternative in each choice situation for each person.
% The rows are grouped by person, and by choice situations faced by each person.
% The number of rows in XMAT must be NROWS, specified above.
% The columns in XMAT are variable that describe the alternative.
% 
% The *first* column of XMAT identifies the person who faced this alternative. 
% The people must be numbered sequentially from 1 to NP, in ascending order.
% All alternatives for a given person must be grouped together.
% The *second* column of XMAT identifies the choice situation. The choice
% situations must be numbered sequentially from 1 to NCS.
% All alternatives for a given choice situation must be grouped together.
% The *third* column of XMAT identifies the chosen alternatives (1 for
% chosen, 0 for not). One and only one alternative must be chosen for each
% choice situation.
% The remaining columns of XMAT can be any variables.

load('videodata100.mat');  %This loads XMAT with the variables below

% To help you keep up with the variables, list the variables in XMAT here.
% NOTES for XMAT for example run:
% This dataset is for people's choice among video streaming services in stated-preference
% experiments. Each person faced 11 experiments, and each
% experiment contained 5 alternatives representing four different services
% whose price and other attributes were described, plus a fifth option
% of not obtaining any service. The person stated which
% of the options he/she would choose if facing this choice in the real world.
% The variables in XMAT are:
% 1. Person number (1-NP)            MUST BE THIS. DO NOT CHANGE.
% 2. Choice situation number (1-NCS) MUST BE THIS. DO NOT CHANGE.
% 3. Chosen alternative (1/0)        MUST BE THIS. DO NOT CHANGE.
% 4. Price in dollars per month 
% 5. share.usage
% 6. share.all
% 7. commercials.yes
% 8. content.fast
% 9. more.TV.fewer.movies
% 10. more.content
% 11. no service

% MODEL SPECIFICATION 

% RANDOM COEFFICIENTS
%Identify the variable in XMAT that is the price variable for the model in wtp space
IDPRICE=4;

% Identify the NON-price variables in XMAT that enter the model.
% Put semicolons between the numbers (so that IDV is a column vector).
IDV=[7;8;9;10;5;6;11];

%Do not change the following line. It gives the number of random parameters.
NV=1+size(IDV,1); 

% Give a name to each of the variables in the model, with the price variable first. They can 
% have up to ten characters including spaces. Put the names in single quotes and separate 
% the quotes with semicolons. If IDV=[], then set NAMES=[];
NAMES={'price';'commer';'fast';'more TV';'more both';'share use'; 'share all';'no service'}; 

% Utility in the model is U=-Alpha*Price+(Alpha*Wtp)*NonPriceVariable
% Give the range of values for Alpha, which is the price/scale coefficient.
% Give the minimum value first, and then the maximum value, 
% *without* a semicolon between.
% Since Alpha*Priceis subtracted from utility and price is positive,
% Alpha must be positive. That is: both values of P_Range must be >= 0.
P_Range=[0 2];

% Give the range for the WTP for each non-price variable.
% WTP_Range is a matrix with one row for each variable in IDV, with 2 columns. 
% The first column gives the lower limit and the second column gives the upper limit.
% Limits are inclusive, in that value of the limit is included.
% Put colons between range for each variable.
WTP_Range=[-9.44    6.32;
           -3.31   11.21;
          -10.42    9.02;
           -2.08    8.00;
           -5.62    4.38;
          -16.20   10.80;
          -66.10   11.58];  %I set these as 2 std dev above and below mean WTP,
                            %with mean and std dev estimated by a mixed
                            %logit in wtp space with all normals for WTPs.

%Do not change the next line
COEF=[P_Range;WTP_Range];


%Give the number of points in each dimension to define the grid for the
%parameter space. The grid includes the two endpoints given in P_Range and WTP_Range. 
%So if you want a grid with eg 1000 intervals between these endpoints, then setNGridPts=1001.
%The total number of points in the grid is NGridPts^NV
NGridPts=1000;

%Number of random draws from the parameter space to use in simulation for each person
NDRAWS=2000;

%Specify the seed to use in the random number generator for simulation
ThisSeed=1234;

%No not change the next line
rng(ThisSeed);

% Specify the type of variables that you want to describe the marginal densities of
% the random parameters. The options are:
% ZTYPE=1 for polynomials.
%      =2 for step functions
%      =3 for splines
ZTYPE=1;

%If ZTYPE=1, specify the order of the polynomial.
%The number of Z variables that will be created is PolyOrder*NV where 
%where NV is the number of random utility parameters.
PolyOrder=12;

%If ZTYPE=2, specify the number of levels for the step functions.
%The height of each level is estimated, with the height of the final step 
%normalized to zero. The number of Z variables that will be 
%created is (NLevels-1)*NV 
NLevels=6;  

%If ZTYPE=3, specify the number of knots for each spline.
%The spline starts at the low end of the range given above,
%changes slopes at each knot, and stops at the high end of the range.
%The number of parameters for a spline with K knots is K+1, i.e.,
%the height at the low end and each knot, with the height
%at the high end normalized to zero. The number of Z variables that 
%will be created is (NKnots+1)*NV 
NKnots=5;

%Do you want to include cross-products for correlations among WTPs?
%Set CrossCorr=1 for yes, =0 for no.
%Note that the Price/scale coef is not allowed to be correlated with WTPs
%   If ZTYPE=1, then the number of extra Z variables that are created with
%CrossCorr=1 is (NV-1)*(NV-2)/2. (There are NV-1 WTPs, and so
%(NV-1)*(NV-2)/2 pairs of WTPs.)
%   If ZTYPE=2 or 3, then the number of extra Z variables is
% 2*(NV-1)+[(NV-1)*(NV-2)/2], i.e, a second order polynomial in each WTP by
% itself (for 2*(NV-1) extra Z variables) plus cross-product terms for each pair
% of WTPs.
CrossCorr=1;

%Do not change the following lines. They calculate the number of Z variables
if ZTYPE==1; 
    NZ=PolyOrder*NV;
    if CrossCorr==1;
        NZ=NZ+(NV-1)*(NV-2)/2;
    end
elseif ZTYPE==2
    NZ=(NLevels-1)*NV;
    if CrossCorr==1;
        NZ=NZ+2*(NV-1)+(NV-1)*(NV-2)/2;
    end
elseif ZTYPE==3
    NZ=(NKnots+1)*NV;
    if CrossCorr==1;
        NZ=NZ+2*(NV-1)+(NV-1)*(NV-2)/2;
    end
end

%Set the starting values for the coefficients of the Z variables
StartB=zeros(NZ,1);  

%The code will estimate the coefficients of the Z variables and create
% a histogram for each random utility parameter based on the estimated
% coefficients of the Z variables.
%Specify the number of bins you want in the histogram.
%More bins gives more detailed shapes but results in more simulation noise in each bin.
NBins=12;

%Do you want the Hessian at convergence? WantHessian=1 for yes, =0 for no
%Setting WantHessian=0 saves computation time
WantHessian=0;

%Do you want to bootstrap the standard errors for estimated coefficients
%and summary statistics?
%Set WantBoot=1 for yes, =0 for no.
%Bootstrapping takes much longer than original estimation, and so
%you might want to bootstrap only after you are pretty sure of your model
WantBoot=1;
 
%If WantBoot=1, specify the number of resamples:
NReps=4;

%Do you have the ability to use Matlab's gpu processor functions?
%YesGPU=1 for yes, and =0 for no.
%To use the gpu processor, you must have Matlab's parallel processing
%toolbox installed on your machine. It helps to have a fast gpu.
%The gpu processor speeds up estimation and bootstrapping considerably.
YesGPU=1;

% ITERATIONS 
% Specify the maximum number of iterations for estimation.
% The code will abort after ITERMAX iterations, even if convergence has
% not been achieved. The default is 1000, which is used when MAXITERS=[];
MAXITERS=2000;

%CONVERENCE TOLERANCE.
% Specify the convergence tolerance for change in parameters as PARAMTOL
% and the change in log-likelihood as LLTOL. Set =[] to use matlab
% defaults.
PARAMTOL=[];
LLTOL=[];

%Do not change the next lines. It runs the model and prints out the results.

tic;
Check;  
if check_ok==1
    CreateDrawsWtpProbs;
    CreateZ;
    disp(' ');
    disp(['Data setup took ' num2str(toc./60) ' minutes.']);
    disp(' ');

    DoEstimation
    if WantBoot==1
       Boot
    end
end
diary off
