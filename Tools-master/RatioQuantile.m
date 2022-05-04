function Y = RatioQuantile(m1,m2,s1,s2,r,cen)

% Calculate percentile 'cen' of a ratio distribution X1/X2 where
% (X1,X2)~BVN(m1,m2;s1,s2;r).

% I could not integrate expression using symbolic math toolbox - 
% CDF is calculated numerically at this point

% my solution using erf function
% f0 = @(w) sqrt(1-r^2)*s1*s2/(2*pi*K2)*exp(-(m1-w*m2)^2/(2*K2))*...
%     (sqrt(2*pi)/(s1*s2)*sqrt(K1)*erf(sqrt(K1)/(sqrt(2)*s1*s2))+...
%     2*exp(-K1/(2*s1^2*s2^2)));
%     K1 = ((m1*s2*(w*s2-r*s1)+m2*s1*(s1-w*r*s2))^2/((1-r^2)*K2));
%     K2 = (w^2*s2^2-2*w*r*s2*s1+s1^2);
% 
% f0 = @(w) sqrt(1-r.^2).*s1.*s2./(2.*pi.*(w.^2.*s2.^2-2.*w.*r.*s2.*s1+s1.^2)).*exp(-(m1-w.*m2).^2./(2.*(w.^2.*s2.^2-2.*w.*r.*s2.*s1+s1.^2))).*...
%     (sqrt(2.*pi)./(s1.*s2).*sqrt(((m1.*s2.*(w.*s2-r.*s1)+m2.*s1.*(s1-w.*r.*s2)).^2./((1-r.^2).*(w.^2.*s2.^2-2.*w.*r.*s2.*s1+s1.^2)))).*erf(sqrt(((m1.*s2.*(w.*s2-r.*s1)+m2.*s1.*(s1-w.*r.*s2)).^2./((1-r.^2).*(w.^2.*s2.^2-2.*w.*r.*s2.*s1+s1.^2))))./(sqrt(2).*s1.*s2))+...
%     2.*exp(-((m1.*s2.*(w.*s2-r.*s1)+m2.*s1.*(s1-w.*r.*s2)).^2./((1-r.^2).*(w.^2.*s2.^2-2.*w.*r.*s2.*s1+s1.^2)))./(2.*s1.^2.*s2.^2))); 
    
% Pham-Gia et al. (2006) - very slow
% f0 = @(w) (1./(2.*pi.*s1.*s2.*(1-r.^2).^0.5) .* exp (-(s2.^2.*m1.^2 - 2.*r.*s1.*s2.*m1.*m2+m2.^2.*s1.^2)./(2.*(1-r.^2).*s1.^2.*s2.^2))) .* (2.*(1-r.^2).*s1.^2.*s2.^2) ./ (s2.^2.*w.^2-2.*r.*s1.*s2.*w+s1.^2) .* hypergeom(1,0.5,((-s2.^2.*m1.*w+r.*s1.*s2.*(m2.*w+m1)-m2.*s1.^2).^2 ./ (2.*s1.^2.*s2.^2.*(1-r.^2).*(s2.^2.*w.^2-2.*r.*s1.*s2.*w+s1.^2))));

% Shanmugalingam(1982)
% f0 = @(w) ...
% K3*K4/(sqrt(2*pi)*s1*s2*K1^3) *...
% (normcdf(K3/(sqrt(1-r^2)*K1),0,1) - normcdf(-K3/(sqrt(1-r^2)*K1),0,1)) + ...
% sqrt(1-r^2)/(pi*s1*s2*K1^2)*...
% exp(-K2/(2*(1-r^2)));

% K1 = (sqrt(w.^2./s1.^2-2.*r.*w./(s1.*s2)+1./(s2).^2));
% K2 = (m1.^2./s1.^2-2.*r.*m1.*m2./(s1.*s2)+m2.^2./s2.^2);   
% K3 = ((m1./s1-r.*m2./s2).*w./s1 + (m2./s2-r.*m1./s1).*1./s2);
% K4 = (exp(1./(2.*(1-r.^2)).*(K3.^2./K1.^2-K2)));
% 
% f0 = @(w) ...
% ((m1./s1-r.*m2./s2).*w./s1 + (m2./s2-r.*m1./s1).*1./s2).*(exp(1./(2.*(1-r.^2)).*(((m1./s1-r.*m2./s2).*w./s1 + (m2./s2-r.*m1./s1).*1./s2).^2./(sqrt(w.^2./s1.^2-2.*r.*w./(s1.*s2)+1./(s2).^2)).^2-(m1.^2./s1.^2-2.*r.*m1.*m2./(s1.*s2)+m2.^2./s2.^2))))./(sqrt(2.*pi).*s1.*s2.*(sqrt(w.^2./s1.^2-2.*r.*w./(s1.*s2)+1./(s2).^2)).^3) .*...
% (normcdf(((m1./s1-r.*m2./s2).*w./s1 + (m2./s2-r.*m1./s1).*1./s2)./(sqrt(1-r.^2).*(sqrt(w.^2./s1.^2-2.*r.*w./(s1.*s2)+1./(s2).^2))),0,1) - normcdf(-((m1./s1-r.*m2./s2).*w./s1 + (m2./s2-r.*m1./s1).*1./s2)./(sqrt(1-r.^2).*(sqrt(w.^2./s1.^2-2.*r.*w./(s1.*s2)+1./(s2).^2))),0,1)) + ...
% sqrt(1-r.^2)./(pi.*s1.*s2.*(sqrt(w.^2./s1.^2-2.*r.*w./(s1.*s2)+1./(s2).^2)).^2).*...
% exp(-(m1.^2./s1.^2 - 2.*r.*m1.*m2./(s1.*s2)+ m2.^2./s2.^2)./(2.*(1-r.^2)));

% Hinkley (1969) - same as Shanmugalingam(1982), differences in notation
% f0 = @(w) ...
%      K3*K4 / (sqrt(2*pi)*s1*s2*K1^3)*...
%      (normcdf(K3/(sqrt(1-r^2)*K1),0,1) - normcdf(-K3/(sqrt(1-r^2)*K1),0,1))+...
%      sqrt(1-r^2)/(pi*s1*s2*K1^2)*...
%      exp(-K2/(2*(1-r^2)));
% 
% K1 = (sqrt(w^2/s1^2 - 2*r*w/(s1*s2) + 1/s2^2)); % dawniej a, tak samo
% K2 = (m1^2/s1^2 - 2*r*m1*m2/(s1*s2) + m2^2/s2^2); % dawniej c, tak samo
% K3 = (m1*w/s1^2 - r*(m1+m2*w)/(s1*s2) + m2/s2^2); % dawniej b, tak samo
% K4 = (exp((K3^2 - K2*K1^2) / (2*(1-r^2)*K1^2))); % dawniej d

f0 = @(w) ...
    (m1.*w./s1.^2 - r.*(m1+m2.*w)./(s1.*s2) + m2./s2.^2).*(exp(((m1.*w./s1.^2 - r.*(m1+m2.*w)./(s1.*s2) + m2./s2.^2).^2 - (m1.^2./s1.^2 - 2.*r.*m1.*m2./(s1.*s2) + m2.^2./s2.^2).*(sqrt(w.^2./s1.^2 - 2.*r.*w./(s1.*s2) + 1./s2.^2)).^2) ./ (2.*(1-r.^2).*(sqrt(w.^2./s1.^2 - 2.*r.*w./(s1.*s2) + 1./s2.^2)).^2))) ./ (sqrt(2.*pi).*s1.*s2.*(sqrt(w.^2./s1.^2 - 2.*r.*w./(s1.*s2) + 1./s2.^2)).^3).*...
    (normcdf((m1.*w./s1.^2 - r.*(m1+m2.*w)./(s1.*s2) + m2./s2.^2)./(sqrt(1-r.^2).*(sqrt(w.^2./s1.^2 - 2.*r.*w./(s1.*s2) + 1./s2.^2))),0,1) - normcdf(-(m1.*w./s1.^2 - r.*(m1+m2.*w)./(s1.*s2) + m2./s2.^2)./(sqrt(1-r.^2).*(sqrt(w.^2./s1.^2 - 2.*r.*w./(s1.*s2) + 1./s2.^2))),0,1))+...
    sqrt(1-r.^2)./(pi.*s1.*s2.*(sqrt(w.^2./s1.^2 - 2.*r.*w./(s1.*s2) + 1./s2.^2)).^2).*...
    exp(-(m1.^2./s1.^2 - 2.*r.*m1.*m2./(s1.*s2) + m2.^2./s2.^2)./(2.*(1-r.^2)));


t0 = -1e6;
f = @(x) quad(@(w) f0(w),t0,x) - cen;
options = optimset('display','off','TolFun',1e-12,'TolX',1e-12);
Y = fsolve(f,m1/m2,options);

% Y = fsolve(f,1); for 0/0

% If optimization toolbox not available - use the following (much slower):
% bsolve = fzero(f,1); 





