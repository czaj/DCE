function der = numdiff(FUN, f, b0, central, BActive)

% save res_numdiff
% return

b0 = b0(:);

DiffMaxChange = 1e-1;
DiffMinChange = 1e-8;
TypicalX = ones(size(b0,1),1);

if nargin < 5 || isempty(BActive) % no fixed parameters
   BActive = ones(size(b0,1), 1);
   if nargin < 4 %assuming forward differencial
       central = 0;
   end
end

if central ~= 1 %forward differencial
    signX = (b0 >= 0) - (b0 < 0);
%     signX = sign(b0);
    deltaX = sqrt(eps)*signX.*max(abs(b0),abs(TypicalX));
    deltaX = signX.*min(max(abs(deltaX),DiffMinChange),DiffMaxChange); % checking if disturbance is within the bounds
    b0_mtx = b0(:, ones(1,size(b0,1))) + diag(deltaX);
    deltaX = diag(b0_mtx) - b0;
    der = zeros(size(f,1),size(b0,1));
    for j = 1:length(b0)
        if  BActive(j) == 1
            der(:,j) = (feval(FUN,b0_mtx(:,j)) - f)/deltaX(j);
        end
    end
else % central differencial
    deltaX = (eps^(1/3))*max(abs(b0),abs(TypicalX));
    deltaX = min(max(abs(deltaX),DiffMinChange),DiffMaxChange); % checking if disturbance is within the bounds
    b0_mtx1 = b0(:, ones(1,size(b0,1))) + diag(deltaX);
    b0_mtx2 = b0(:, ones(1,size(b0,1))) - diag(deltaX);
    deltaX = diag(b0_mtx1) - diag(b0_mtx2);
    der = zeros(size(f,1),size(b0,1));
    for j = 1:length(b0)
        if  BActive(j) == 1
            der(:,j) = (feval(FUN,b0_mtx1(:,j)) - feval(FUN,b0_mtx2(:,j)))/deltaX(j);
        end
    end
end 