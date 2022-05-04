function result = pv(input1,varargin)

% save tmp1

if nargin < 1 % check no. of inputs
    error('Too few input arguments')
elseif nargin > 1
    %     m = m(:,1); %data must be in columns
    m = input1;
    s = varargin{1};
elseif nargin == 1
    if size(input1,2) > 2
        cprintf(rgb('DarkOrange'),'WARNING: using the first 2 columns of input only, assuming they are means and standard errors, respectively  \n');
        s = input1(:,2);
        m = input1(:,1);
    elseif size(m,2) == 2
        s = input1(:,2);
        m = input1(:,1);
    else
        error('Input must include means and standard errors (as two separate inputs, or one 2-column matrix)')
    end
end

s(logical(imag(s)) | s < 0) = NaN;

result = (1-normcdf(abs(m)./s,0,1))*2;

% result = NaN(size(m));
% result(~logical(imag(s)) & s>=0) = (1-normcdf(abs(m(~logical(imag(s)) & s>=0))./real(s(~logical(imag(s)) & s>=0)),0,1))*2;    

% if any(s < 0) || ~isreal(s)
%     result = zeros(size(m));    
%     for i = 1:length(m)
%         if s(i) < 0 || ~isreal(s(i))
%             result(i) = NaN;
%         else
%             result(i) = (1-normcdf(abs(m(i))./real(s(i)),0,1))*2;
%         end
%     end
% else
%     result = (1-normcdf(abs(m)./s,0,1))*2;
% end
