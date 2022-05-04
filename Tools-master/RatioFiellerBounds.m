function [a1,a2,a3,a4] = RatioFiellerBounds(B,VARB,varargin)

% Calculate 95% Fieller bounds of a ratio of X1/X2 where (X1,X2)~BVN(m1,m2;s1,s2;r)
% in the following specification_cases:
%   0 - coefficients enter linearly
%   1 - cost (last coefficient in B) enters exponentially
%   2 - all coefficients enter exponentially
%   3 - non-cost is log-normally distributed (parameters of the underlying normal provided), cost (last coefficient in B) enters linearly
%   4 - non-cost is log-normally distributed (parameters of the underlying normal provided), cost (last coefficient in B) enters exponentially


if nargin > 3 % case 0-2
    
    m1 = B;
    m2 = varb;
    s1 = varargin{1};
    s2 = varargin{2};
    r = varargin{3};
    if nargin == 5
        specification_case = 0;
    else
        specification_case = varargin{4};
    end    
    VARB = [s1^2 r*s1*s2; r*s1*s2 s2^2];
    cr = VARB(1,2);

else % case 3-4
    
    m1m = B(1);
    m1s = ((B(2))^2)^0.5;
    m2 = B(3);
    s1m = VARB(1,1)^0.5;
    s1s = VARB(2,2)^0.5;
    s2 = VARB(3,3)^0.5;
    c1 = VARB(1,2);
    c2 =  VARB(1,3);
    c3 =  VARB(2,3);
    if nargin == 2
        specification_case = 3;
    else
        specification_case = varargin{1};
    end
    
end

z0 = norminv(0.975,0,1); % critical value for the 95% c.i. (change for other c.i.)

switch specification_case
    
    case 0
        
        q025f_tmp = (-1).*(m2.^2+(-1).*s2.^2.*z0.^2).^(-1).*((-1).*m1.*m2+cr.*z0.^2+( ...
            z0.^2.*((-2).*cr.*m1.*m2+m2.^2.*s1.^2+cr.^2.*z0.^2+s2.^2.*(m1.^2+( ...
            -1).*s1.^2.*z0.^2))).^(1/2));
        q975f_tmp = (m2.^2+(-1).*s2.^2.*z0.^2).^(-1).*(m1.* ...
            m2+(-1).*cr.*z0.^2+(z0.^2.*((-2).*cr.*m1.*m2+m2.^2.*s1.^2+cr.^2.* ...
            z0.^2+s2.^2.*(m1.^2+(-1).*s1.^2.*z0.^2))).^(1/2));
        Fposit = s2.^2.*z0.^2 - m2.^2; % g^2 coefficient, if < 0 - interval is bounded
        % x_tmp = -(2.*m1.*m2+(-2).*cr.*z0.^2) / (2*((-1).*m2.^2+s2.^2.*z0.^2));
        % Dposit = (-1).*m1.^2+s1.^2.*z0.^2+x_tmp.*(2.*m1.*m2+(-2).*cr.*z0.^2)+x_tmp.^2.*((-1).*m2.^2+s2.^2.*z0.^2);
        
    case 1
        
        q025f_tmp = [exp(1).^((-2).*m2).*((-1)+s2.^2.*z0.^2).^(-1).*(exp(1).^m2.*((-1) ...
            .*m1+cr.*z0.^2)+(-1).*(exp(1).^(2.*m2).*z0.^2.*((-2).*cr.*m1+s1.^2+ ...
            m1.^2.*s2.^2+cr.^2.*z0.^2+(-1).*s1.^2.*s2.^2.*z0.^2)).^(1/2))];
        q975f_tmp = [exp(1).^((-2).*m2).*((-1)+s2.^2.*z0.^2).^(-1).*(exp(1).^m2.*((-1).*m1+ ...
            cr.*z0.^2)+(exp(1).^(2.*m2).*z0.^2.*((-2).*cr.*m1+s1.^2+m1.^2.* ...
            s2.^2+cr.^2.*z0.^2+(-1).*s1.^2.*s2.^2.*z0.^2)).^(1/2))];
        Fposit = s2.^2.*z0.^2 - 1; % g^2 coefficient, if < 0 - interval is bounded
        
    case 2
        
        q025f_tmp = [exp(1).^((-2).*m2).*((-1)+s2.^2.*z0.^2).^(-1).*(exp(1).^(m1+m2).*( ...
            (-1)+cr.*z0.^2)+(-1).*(exp(1).^(2.*(m1+m2)).*z0.^2.*((-2).*cr+s1.^2+ ...
            s2.^2+cr.^2.*z0.^2+(-1).*s1.^2.*s2.^2.*z0.^2)).^(1/2))];
        q975f_tmp = [exp(1).^((-2).*m2).*((-1)+s2.^2.*z0.^2).^(-1).*(exp(1).^(m1+m2).*((-1)+cr.*z0.^2) ...
            +(exp(1).^(2.*(m1+m2)).*z0.^2.*((-2).*cr+s1.^2+s2.^2+cr.^2.*z0.^2+( ...
            -1).*s1.^2.*s2.^2.*z0.^2)).^(1/2))];
        Fposit = s2.^2.*z0.^2 - 1; % g^2 coefficient, if < 0 - interval is bounded
        
    case 3
        
        q025f_tmp = (-1).*((-1).*m2.^2+s2.^2.*z0.^2).^(-1).*(exp(1).^(m1m+(1/2).* ...
            m1s.^2).*(m2+(-1).*(c2+c3.*m1s).*z0.^2)+(exp(1).^(2.*m1m+m1s.^2).*( ...
            (m2+(-1).*(c2+c3.*m1s).*z0.^2).^2+(-1).*((-1)+2.*c1.*m1s.*z0.^2+ ...
            s1m.^2.*z0.^2+m1s.^2.*s1s.^2.*z0.^2).*((-1).*m2.^2+s2.^2.*z0.^2))).^(1/2));
        q975f_tmp = (1/2).*((-1).*m2.^2+s2.^2.*z0.^2).^(-1).*((-2).*exp(1).^(m1m+ ...
            (1/2).*m1s.^2).*(m2+(-1).*(c2+c3.*m1s).*z0.^2)+2.*(exp(1).^(2.*m1m+ ...
            m1s.^2).*((m2+(-1).*(c2+c3.*m1s).*z0.^2).^2+(-1).*((-1)+2.*c1.* ...
            m1s.*z0.^2+s1m.^2.*z0.^2+m1s.^2.*s1s.^2.*z0.^2).*((-1).*m2.^2+s2.^2.* ...
            z0.^2))).^(1/2));
        Fposit = (-1)+s2.^2.*z0.^2;
        
    case 4
        
        q025f_tmp = (1/2).*exp(1).^((-2).*m2).*((-1)+s2.^2.*z0.^2).^(-1).*(2.*exp(1) ...
            .^(m1m+(1/2).*m1s.^2+m2).*((-1)+c2.*z0.^2+c3.*m1s.*z0.^2)+(-2).*( ...
            exp(1).^(2.*m1m+m1s.^2+2.*m2).*(((-1)+c2.*z0.^2+c3.*m1s.*z0.^2).^2+( ...
            -1).*((-1)+2.*c1.*m1s.*z0.^2+s1m.^2.*z0.^2+m1s.^2.*s1s.^2.*z0.^2).*(( ...
            -1)+s2.^2.*z0.^2))).^(1/2));
        q975f_tmp = exp(1).^((-2).*m2).*((-1)+s2.^2.*z0.^2) ...
            .^(-1).*(exp(1).^(m1m+(1/2).*m1s.^2+m2).*((-1)+c2.*z0.^2+c3.*m1s.* ...
            z0.^2)+(exp(1).^(2.*m1m+m1s.^2+2.*m2).*(((-1)+c2.*z0.^2+c3.*m1s.* ...
            z0.^2).^2+(-1).*((-1)+2.*c1.*m1s.*z0.^2+s1m.^2.*z0.^2+m1s.^2.* ...
            s1s.^2.*z0.^2).*((-1)+s2.^2.*z0.^2))).^(1/2));
        Fposit = (-1)+s2.^2.*z0.^2;
        
end



% q025f_tmp = q025f_tmp*isreal(q025f_tmp);
% q975f_tmp = q975f_tmp*isreal(q975f_tmp);
% q025f = min(q025f_tmp, q975f_tmp);
% q975f = max(q025f_tmp, q975f_tmp);
if isreal(q025f_tmp)
    q025f = min(q025f_tmp,q975f_tmp);
else
    q025f = NaN;
end
if isreal(q975f_tmp)
    q975f = max(q025f_tmp,q975f_tmp);
else
    q975f = NaN;
end


a2 = q025f;
a3 = q975f;

if Fposit < 0 % confidence set is bounded
    a1 = NaN;
    a4 = NaN;
else % confidence set is unbounded
    a1 = - Inf;
    a4 = Inf;
end
