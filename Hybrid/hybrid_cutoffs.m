function [threshold,se,p] = hybrid_cutoffs(raw,V,maxExp)
% Ordered thresholds: first raw cutoff + cumulative positive increments.
if nargin < 3, maxExp = inf; end
raw = raw(:); n = numel(raw);
assert(n >= 1 && isequal(size(V),[n n]),'DCE:Hybrid:Cutoffs','Cutoff covariance dimensions do not match.');
increments = exp(min(max(raw(2:end),-maxExp),maxExp));
threshold = cumsum([raw(1);increments]);
derivative = [1;increments.*(abs(raw(2:end)) < maxExp)];
J = tril(ones(n)).*derivative';
variance = sum((J*V).*J,2);
se = nan(n,1); ok = isfinite(variance) & variance > 0;
se(ok) = sqrt(variance(ok));
p = erfc(abs(threshold./se)/sqrt(2));
end
