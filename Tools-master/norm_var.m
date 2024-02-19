function X = norm_var(INPUT,varargin) % option 1 for 0-1 normalization

if nargin == 1
    norm_version = 0;
else
    norm_version = varargin{1};
end

if norm_version == 0 % normalize for mean = 0, s.d. = 1
    X = INPUT - nanmean(INPUT);
    X = X ./ nanstd(X);
elseif norm_version == 1
    X = INPUT - min(INPUT);
    X = X ./ max(X);
end
