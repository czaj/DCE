function Y = ReplicateRows(X,no,varargin)

% Y = ReplicateRows(X,no,varargin)
% only tested for 2D matrices 

if (size(X,2) == 1) % X is a column vector
    Y = repmat(X,1,no)';
    Y = Y(:);
elseif size(X,2) > 1
%     Y = repmat(X,[1,1,no]);
%     Y = permute(X
    Y = repelem(X,no,1);
end
    