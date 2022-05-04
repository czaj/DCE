function [h,p,chi2stat,df] = prop_test2(X1,X2,varargin)

% save tmp1

% check no. of inputs
if nargin < 2
    error('Too few input arguments - need to provide two vectors to compare')
elseif nargin == 2
    correct = false;
elseif nargin == 3
    correct = varargin{1}; % should probably test if this is ok (true/false)
end

% Test if X1 and X2 are vectors (should also check if matrix?)
if ~isvector(X1)
    X1 = X1(:);
end
if ~isvector(X2)
    X2 = X2(:);
end

% Check for missing values (also Infs? test... )
if any(isnan(X1)) || any(isnan(X2))
    cprintf(rgb('DarkOrange'), 'WARNING: Input vectors include NaN values - they will be skipped \n')
    % NaNs could be treated as if they were a yet another category - add option, for later)
    X1 = X1(~isnan(X1));
    X2 = X2(~isnan(X2));
end

% Search for unique values:
uX1 = unique(X1);
uX2 = unique(X2);
u = unique([uX1;uX2]);

% Count number of successes
nY1 = zeros(size(u)-1);
nY2 = zeros(size(u)-1);
for i = 1:size(u,1)-1
    nY1(i) = sum(X1(:) == u(i));
    nY2(i) = sum(X2(:) == u(i));
end

% Count number of cases
N1 = size(X1,1);
N2 = size(X2,1);

% Caucluate expected proportions
p0 = (nY1 + nY2) / (N1 + N2);

% Calculate expected counts under H0
n10 = N1 * p0;
n20 = N2 * p0;
observed = [nY1 N1-nY1 nY2 N2-nY2];
expected = [n10 N1-n10 n20 N2-n20];

% degrees of freedom
df = size(u,1)-1;

if correct == false
    % Standard Chi-square test
    chi2stat = sum((observed-expected).^2 ./ expected);
    p = 1 - chi2cdf(chi2stat,df);
else
    % Yates continuity correction
    chi2stat = sum((abs(observed - expected) - 0.5).^2 ./ expected);
    p = 1 - chi2cdf(chi2stat,df);
end

h=0;
if p<0.05
    h=1;
end

% n1 = 19; N1 = 135;
% n2 = 16; N2 = 119;
% % Pooled estimate of proportion
% p0 = (n1+n2) / (N1+N2)
% % Expected counts under H0 (null hypothesis)
% n10 = N1 * p0;
% n20 = N2 * p0;
% % Chi-square test, by hand
% observed = [n1 N1-n1 n2 N2-n2];
% expected = [n10 N1-n10 n20 N2-n20];
% chi2stat = sum((observed-expected).^2 ./ expected)
% p = 1 - chi2cdf(chi2stat,1)

% [h,p, chi2stat,df] = prop_test(X , N, correct)
%
% A simple Chi-square test to compare two proportions
% It is a 2 sided test with alpha=0.05
%
% Input:
% X = vector with number of success for each sample (e.g. [20 22])
% N = vector of total counts for each sample (e.g. [48 29])
% correct = true/false : Yates continuity correction for small samples?
%
% Output:
% h = hypothesis (H1/H0)
% p = p value
% chi2stat= Chi-square value
% df = degrees of freedom (always equal to 1: 2 samples)
%
% Needs chi2cdf from the Statistics toolbox
% Inspired by prop.test() in "R" but much more basic
%
% Example: [h,p,chi]=prop_test([20 22],[48 29], true)
% The above example tests if 20/48 differs from 22/29 using Yate's correction

% if (length(X)~= 2)||(length(X)~=length(N))
%     disp('Error: bad vector length')
% elseif (X(1)>N(1))|| (X(2)>N(2))
%     disp('Error: bad counts (X>N)')
% else
%     df=1; % 2 samples
%     
%     % Observed data
%     n1 = X(1);
%     n2 = X(2);
%     N1 = N(1);
%     N2 = N(2);
%     
%     % Pooled estimate of proportion
%     p0 = (n1+n2) / (N1+N2);
%     
%     % Expected counts under H0 (null hypothesis)
%     n10 = N1 * p0;
%     n20 = N2 * p0;
%     observed = [n1 N1-n1 n2 N2-n2];
%     expected = [n10 N1-n10 n20 N2-n20];
%     
%     if correct == false
%         % Standard Chi-square test
%         chi2stat = sum((observed-expected).^2 ./ expected);
%         p = 1 - chi2cdf(chi2stat,1);
%     else
%         % Yates continuity correction
%         chi2stat = sum((abs(observed - expected) - 0.5).^2 ./ expected);
%         p = 1 - chi2cdf(chi2stat,1);
%     end
%     
%     h=0;
%     if p<0.05
%         h=1;
%     end
% end
% end