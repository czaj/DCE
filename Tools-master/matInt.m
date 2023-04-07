function MAT = matInt(X,I,varargin) 
% XfirstI==1 or missing - first element of X multiplied by each element of I
% XfirstI==0 - each element of X multiplied by first element of I

if nargin == 2
    XfirstI = 1;
elseif nargin == 3
    XfirstI = varargin{1};   
end

% save tmp1

if iscell(X) && iscell(I) % this is not thorougly tested, especially for size(X,1) > 1
    if XfirstI == 1
        for i = 1:size(X,2)
            for j = 1:size(I,2)
                MAT{:,(i-1)*size(I,2)+j} = [X{:,i},I{:,j}];
            end
        end
    elseif XfirstI == 0
        for j = 1:size(I,2)
            for i = 1:size(X,2)
                MAT{:,(j-1)*size(X,2)+i} = [X{:,i},I{:,j}];
            end
        end
    end
else
    if XfirstI == 1
        MAT = X(:,ceil((1:size(I,2)*size(X,2))/size(I,2))) .* I(:,[1:size(I,2)]' * ones(1,size(X,2)));
    elseif XfirstI == 0
        MAT = X(:,[1:size(X,2)]' * ones(1,size(I,2))) .* I(:,ceil((1:size(X,2)*size(I,2))/size(X,2)));
    end
end