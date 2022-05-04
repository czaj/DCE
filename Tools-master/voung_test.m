function T = voung_test(LL1, LL2, p,q)
disp('Voung test for strictly non-nested models') 
if nargin < 3
   disp('Assuming equal no. of parameters')
   p = 0; q = 0;
end
if length(LL1) ~= length(LL2)
   error('Likelihoods vectors should have the same length')
end

N = length(LL1);

LR = sum(LL1 - LL2,1) - 0.5*(p-q)*log(N);

Voung = LR./std(LL1 - LL2);
Voung = Voung/sqrt(N);

pval = 2*(1 - normcdf(abs(Voung),0,1));

disp(num2str(Voung,'Voung test statistics: %8.4f'))
disp(num2str(pval,'pvalue: %1.4f'))
disp(' ')
T = [Voung, pval];