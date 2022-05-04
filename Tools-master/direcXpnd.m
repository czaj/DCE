function new_direc = direcXpnd(direc, bactive)

% direc is K-no.of constranted parameters by 1 parameters
% bactive is 1 by K vector
bactive = bactive(:)';

if size(direc,2) ~= (size(bactive,2)-sum(bactive==0))
    error('Number of columns of direc and bactive are not consistent.')
end

new_direc = zeros(size(direc,1),size(bactive,2));

for i = 1:size(bactive,2)
    if bactive(i)==1
        tmp = sum(bactive(1:i)==0);
        new_direc(:,i) = direc(:,i-tmp);
    
    elseif bactive(i) == 0
        new_direc(:,i)= 0;
    end
end




