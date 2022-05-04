function Y = VarNorm(X)
X = X - nanmean(X);
X = X ./ nanstd(X);
Y = X;

