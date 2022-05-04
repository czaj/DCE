%-------------------------------------------------------------------------------------------------------------------------------------
%   Author: Ian Gregory
%   Date:   25th January 2009
%   Functionality:  To simulate the basic ordered probit equation.
%                   (Default Model)  - Assume:  y*=X*b+e          , where e approx N(0,1).             
%                   (User specified) - Assume:  y*=Const+X*b+e    , where e approx N(0,1).             
%
%   Dependencies:   None.
%                                            =================        INPUTS:           =================
%            ---------------         *REQUIRED:*           ---------------        
%       INPUTS.Beta                                   :       The independent variables to simulate.  Eg. INPUTS.Independent_Variables=[0.5,2,3];.  ie. y*=0.5X_1 + 2X_2 + 3X_3
%       INPUTS.No_Data_Points                         :       number of data points.  ie. n=1000.
%       INPUTS.Cut_Points                             :       Ordered list of cut-points to split the data into. Eg. INPUTS.Cut_Points=[-0.5,0.5].
%
%            ---------------         *OPTIONAL:*           ---------------        
%       INPUTS.Constant_Value                         :    Value of a constant if wishing to include one.  ie. INPUTS.Constant_Value=3;      
%       INPUTS.Include_Const_Switch                   :    Switch to indicate, want to simulate with a constant.  Take value of 1.  Eg.  INPUTS.Constant.Include_Const_Switch=1; 
%
%       INPUTS.fNameOutput                            :     User specified output filename.  Default filename: "DATA_Ordered_Probit_Simulated.csv". Eg. INPUTS.fNameOutput='My_Sim_OProbit_DATA.csv';
%                                                               OUTPUT FORMAT:
%                                 eg. dependent_variable,independent_variable_1 (-0.43256),independent_variable_2 (-1.6656),independent_variable_3 (0.12533),independent_variable_4 (0.28768)                                
%                                     2.504,0.8644,-1.7486,-0.27884,0.0014914
%
%                                            =================        OUTPUTS:           =================
% 
%       The function outputs a .csv flat file with a default filename: "DATA_Ordered_Probit_Simulated.csv" or a user specified filename.
%  
%                                            =================        Example Usage:           ================= 
%    INPUTS.Beta=[0.5,2,3];
%    INPUTS.No_Data_Points=1000;            % Number of data points.
%    INPUTS.Cut_Points=[-0.5,0.5];          % Cut-points..  ie. 2 cut-points gives 3 bins.
%    INPUTS.Constant_Value=3;
%    INPUTS.Include_Const_Switch=1;
%    INPUTS.fNameOutput='My_Sim_OProbit_DATA.csv';
%    MATLAB_Ordered_Probit_Simulate(INPUTS)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

function MATLAB_Ordered_Probit_Simulate1(varargin)
    if nargin<1
        disp('Error: No inputs were specified.')
        disp('The function "MATLAB_Ordered_Probit_Simulate" requires syntax as follows')
        help MATLAB_Ordered_Probit_Simulate
        return
    end
    if ~isstruct(varargin{1})
        disp('Error: Input must be a structure.');
    end    
    if length(varargin)>1
       disp('Error: Expecting only one input (which is as struct).') 
    end


    %
    % Create an empty shell for the INPUT specification structure.
    %
    options = struct('No_Data_Points'    , []   ,   'Cut_Points'              , [] , ...
                     'Beta'              , []   ,   'Include_Const_Switch'    , [] , ...
                     'Constant_Value'    , []   ,   'fNameOutput'             , []       );


    Required_Struct_Elements=strvcat('No_Data_Points','Cut_Points','Beta');
    Optional_Struct_Elements=strvcat('options.Include_Const_Switch','options.Constant_Value','options.fNameOutput');

    % NOTE: The following constraint conditions need to line up with the names specified in the struct 'options'.    
    Value_Lower_Bound=[0,-inf,-inf,0,-inf];
    Value_Uppper_Bound=[inf,inf,inf,2,inf];
    
    Struct_Elements=fieldnames(options);
    for i=1:length(Struct_Elements)
    % For all the struct elements
        if isfield(varargin{1},Struct_Elements{i})
        % The user specified the coefficients container for the structure.
            Values=eval((['varargin{1}.',eval('Struct_Elements{i}')]));
            if ~isempty(Values)
                % Checking the value the user has supplied.
                if isnumeric(Values) && sum(isnan(Values))==0
                % Values are numeric and don't contain any NaN's.
                	if sum(Values>Value_Lower_Bound(i))==length(Values) && sum(Values<Value_Uppper_Bound(i))==length(Values)
                    % Assigning the user inputs here.                   
                    	options=setfield(options,Struct_Elements{i},eval(['varargin{1}.',eval('Struct_Elements{i}')]));
                    else
                    	disp(['Error: There is a constraint on  ',['varargin{1}.',eval('Struct_Elements{i}')]])
                        disp(['Lower value must be greater than: ',num2str(Value_Lower_Bound(i)),' and smaller than: ',num2str(Value_Uppper_Bound(i))])
                    	return
                    end
                else
                	if strcmp(eval('Struct_Elements{i}'),'fNameOutput')
                    	if ischar(Values)
                        % File name is a character array (string).  Ok to set the value.
                            options=setfield(options,Struct_Elements{i},eval(['varargin{1}.',eval('Struct_Elements{i}')]));
                        else                                
                        	disp('Error: Output filename must be a valid character array (string).') 
                            return                                
                        end
                    else
                    	disp('Error: All Beta values must be numeric type and not contain any NaN.') 
                        return
                    end
                end
            else
                disp('Error: You did not specify any Beta to simulate.') 
                return
            end
        else
        % Check if this field is required or optional.
            if ~isempty(strmatch(Struct_Elements{i},Required_Struct_Elements))
                disp(['Error: The required field ',Struct_Elements{i},' was not assigned by user.'])
                return
            else
            % Assign the default value (if applicable).

            end
        end
    end

    Beta=options.Beta;
    if size(Beta,1)==1 || size(Beta,2)==1
    % Betas are a vector, ok to continue.
    else
    % Betas are a matrix.  Betas, mispecified.  Notifying user of the error and quitting.
        disp('Error: Coefficients must be a vector not a matrix.')

    end
    
    
    Cut_Points=sort(options.Cut_Points);            % Ordering the cut-offs from smallest to largest.      
    n=length(options.Beta);
    No_Independent_Variables=n;
    m=options.No_Data_Points;

    
    %-------------------------------------------------------------   
    % Begin of: Globl variables.
    % Declare the state to simulate from.
    rand('state',37); % set arbitrary seed for uniform draws
    randn('state',37); % set arbitrary seed for normal draws    
    
    if isempty(options.fNameOutput)
    % Use default output filename.
        fNameOutput='DATA_Ordered_Probit_Simulated.csv';
    else
    % User has specified a valid output file name.  Using this one.
       fNameOutput=options.fNameOutput;
    end

    % End of: Globl variables.    
    %-------------------------------------------------------------          
    
    %-------------------------------------------------------------   
    % Begin of: Create the DATA.   
    randn('state',0);          % State of data generation for the independent coefficients.        
        
    if exist('Beta')~=1
    % User did not specify the coefficients, will generate some for them.
        b=randn(n,1);              % declare simulation coefficients values.
    else
        if size(Beta,2)~=No_Independent_Variables
            disp('Incorrect number of independent variables specified.')
            disp(['You specified number of independent variables as: ',num2str(No_Independent_Variables)])
            disp(strcat('You specified the coefficients "Betas" as: ',num2str(Beta')))
            disp('Exited function with error.')
            return
        else
            if isempty(Beta)
                disp(['You did not specify the beta coefficients.'])
                disp('Exited function with error.')            
                
            else
                if Beta==0
                    disp(['You specified the beta coefficient as zero.  This does not make sense.'])
                    disp('Exited function with error.')                            
                else                
                   b=Beta;
                end
            end
        end
    end
    randn('state',1);       % State of data generation for the independent variables.       
    % Create the independent component of the equation.
    X=randn(m,n);   

    randn('state',2);          % State of data generation for error.        
    e=randn(m,1);      
    
    % Observed y:
    if isempty(options.Include_Const_Switch)==1 || options.Include_Const_Switch==0
        if size(Beta,1)>1 || size(Beta,2)>1
            if size(Beta,1)>1
                y_star=X*b+e;    
            else
                y_star=X*b'+e;
            end
        else
            y_star=X*b'+e;                
        end
    else
%         Constant_Vector=repmat(options.Constant,size(X,1),1);
        if size(b,2)==size(X,2)
            y_star=options.Constant_Value+X*b'+e;
        else
            y_star=options.Constant_Value+X*b+e;    
        end
    end
    
    % Creating the unobserved y:
    y(1:m,1)=0;
	for i=1:length(Cut_Points)
        if length(Cut_Points)==1
        % If only one cut-point, splitting into two bins.
            eval([['f',num2str(i)],'=find(y_star<=Cut_Points(i));'])
            eval([[['y(f',num2str(i),',1)=',num2str(i)]],';']);
                
            eval([['f',num2str(i+1)],'=find(y_star>Cut_Points(i));'])
            eval([[['y(f',num2str(i+1),',1)=',num2str(i+1)]],';']);
        else
        % Multiple cut-points generate multiple bins.
        % The following code performs similar to the following.
%             f0=find(y_star<=cut_point_1);
%             f1=find(y_star>cut_point_1 & y_star<=cut_point_2);
%             f2=find(y_star>cut_point_2 & y_star<=cut_point_3);
%             f3=find(y_star>cut_point_3);                        
%
% And then assigns an integer to the data in the cut-off points.  Similar to the following:
%             y(Position,1)=0;
%             y(Position,1)=1;
%             y(Position,1)=2;
%             y(Position,1)=3;  
    
            if i==1
                eval([['f',num2str(i)],'=find(y_star<=Cut_Points(i));'])
                eval([[['y(f',num2str(i),',1)=',num2str(i)]],';']);
            else                                
                if i==length(Cut_Points)
                    
                    eval([['f',num2str(i)],'=find(y_star>Cut_Points(i-1) & y_star<=Cut_Points(i));'])
                    eval([[['y(f',num2str(i),',1)=',num2str(i)]],';']);
                    
                    eval([['f',num2str(i+1)],'=find(y_star>Cut_Points(i));'])
                    eval([[['y(f',num2str(i+1),',1)=',num2str(i+1)]],';']);
                else
                    eval([['f',num2str(i)],'=find(y_star>Cut_Points(i-1) & y_star<=Cut_Points(i));'])
                    eval([[['y(f',num2str(i),',1)=',num2str(i)]],';']);
                end
            end            
        end
    end
    
    % End of: Create the DATA.      
    %-------------------------------------------------------------       
        
    %-------------------------------------------------------------   
    % Begin of: Output the column headings and data:
	fid_output = fopen(fNameOutput,'wt');                        
    
    if fid_output==-1;
       warning('File could not be opened.  Probably because file already exsits and already opened.');
       return;
    end
    
    for j=1:(m+1)
    % For the length of all the data plus a column heading.
        for i=1:(n+1)
        % For all the independent variables plus the dependent variable.
            if i==1 
                if j==1
                    output_str=['dependent_variable'];  
                else
                    output_str=num2str(y(j-1,i));
                end
            else
                if j==1
%                     output_str=[output_str,',',['independent_variable_',num2str(i-1)],' (',num2str(b(i-1)),')'];  % Stating the 'true' value in the column heading
                    output_str=[output_str,',',['independent_variable_',num2str(i-1)]];  % Stating the 'true' value in the column heading
                else
                    output_str=[output_str,',',num2str(X(j-1,i-1))];
                end
            end
        end   
        fprintf(fid_output,output_str);
    	fprintf(fid_output,'\n');
    end   
    status = fclose(fid_output);
    % End of: Output the column headings and data.
    %-------------------------------------------------------------               

    disp('Finished simulating ordered probit data.')
end
    
    function CloseFiles
    % On an error.  This function closes all open files.   
        FIDs = fopen('all');
        for i=1:length(FIDs)
            fclose(FIDs(i));
        end
    end