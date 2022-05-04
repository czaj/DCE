function f = RatioDeltaSE(B,VARB,varargin)

% Calculate s.e. of the ratio using the delta method in the following specification_cases:
%   0 - coefficients enter linearly
%   1 - cost (last coefficient in B) enters exponentially
%   2 - all coefficients enter exponentially
%   3 - non-cost is log-normally distributed (parameters of the underlying normal provided), cost (last coefficient in B) enters linearly
%   4 - non-cost is log-normally distributed (parameters of the underlying normal provided), cost (last coefficient in B) enters exponentially

if nargin == 2
    specification_case = 0;
else
    specification_case = varargin{1};
end

if specification_case < 3
    
    m1 = B(1:end-1);
    m2 = B(end);
    s1 = VARB(1,1)^0.5;
    s2 = VARB(2,2)^0.5;
    cr = VARB(1,2);
    
else
    
    m1m = B(1);
    m1s = ((B(2))^2)^0.5;
    m2 = B(3);
    s1m = VARB(1,1)^0.5;
    s1s = VARB(2,2)^0.5;
    s2 = VARB(3,3)^0.5;
    c1 = VARB(1,2);
    c2 =  VARB(1,3);
    c3 =  VARB(2,3);
    
end

switch specification_case
    case 0
        f = (m2.^(-4).*((-2).*cr.*m1.*m2+m2.^2.*s1.^2+m1.^2.*s2.^2))^0.5;
    case 1
        f = (exp(1).^((-2).*m2).*((-2).*cr.*m1+s1.^2+m1.^2.*s2.^2))^0.5;
    case 2
        f = (exp(1).^(2.*m1+(-2).*m2).*((-2).*cr+s1.^2+s2.^2))^0.5;
    case 3
        f = (exp(1).^(2.*m1m+m1s.^2).*m2.^(-4).*((-2).*c2.*m2+(-2).*c3.*m1s.* ...
            m2+2.*c1.*m1s.*m2.^2+m2.^2.*s1m.^2+m1s.^2.*m2.^2.*s1s.^2+s2.^2))^0.5; % delta se for mean b / c
    case 4
        f = (exp(1).^(2.*m1m+m1s.^2+(-2).*m2).*((-2).*c2+2.*c1.*m1s+(-2).*c3.* ...
            m1s+s1m.^2+m1s.^2.*s1s.^2+s2.^2))^0.5;
end

