function X = norm_var(INPUT)

X = INPUT - nanmean(INPUT);
X = X ./ nanstd(X);
