function [W, p, dF] = waldtestx(hfun, b0, hess)

% save tmp1

h = feval(@(b) hfun(b), b0);
dF = length(h);
H = jacobianest(@(b) hfun(b), b0);
I = inv(hess);
X = inv(H*I*H');
W = h'*X*h;
p =  1 - chi2cdf(W,dF);
