function MAT = matInt(X,I)

if iscell(X) && iscell(I) % this is not thorougly tested, especially for size(X,1) > 1
    for i = 1:size(X,2)
        for j = 1:size(I,2)
            MAT{:,(i-1)*size(I,2)+j} = [X{:,i},I{:,j}];
        end
    end    
else
    MAT = X(:,ceil((1:size(I,2)*size(X,2))/size(I,2))) .* I(:,[1:size(I,2)]' * ones(1,size(X,2)));
end