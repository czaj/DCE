function LLV = MATLAB_Ordered_Probit_Likelihood(b,y,X,k)
%------------------------------------------------------------------------------------
% Author        :       Ian Gregory
% Date:         :       14th November 2009
% Functionality :       To produce a vector of the ordered probit log likelihood
%
%   INPUTS:   
%   b           :       beta coefficients
%   y           :       categories
%   X           :       Independent variables
%   k           :       Number of independent variables
%
%   OUTPUTS:
%   LLV         :       Log likelihood values.  Vector size is 1 by size(X,1).
%------------------------------------------------------------------------------------
    Xb=X*b([1:k],:);         % y_star
    alpha=b([k+1:end],:);

    LLV=zeros(size(X,1),1);   % Declare likelihood vector

    Unique_Dependent_Variables=unique(y);           % ie.  Assuming the categories are conditioned as integers and sorting them.  ie. 0,1,2,3...               
        
	% 1st section
	d_start=find(y==Unique_Dependent_Variables(1));     % Distribution start
 
    LLV(d_start)=log(normcdf(alpha(1)-Xb(d_start)));

    % Last section
    d_end=find(y==Unique_Dependent_Variables(end));     % Distribution end
    LLV(d_end)=log(1-normcdf(alpha(end)-Xb(d_end)));

	% Middle sections (if applicable)
    if length(Unique_Dependent_Variables)>1
    % If there are atleast two cut-off points.
        for i=2:length(Unique_Dependent_Variables)-1      

            d_middle=find(y==Unique_Dependent_Variables(i));     % Distribution middle
            LLV(d_middle)=log(normcdf(alpha(i)-Xb(d_middle)) - normcdf(alpha(i-1)-Xb(d_middle))); 
        end 
    end    
end