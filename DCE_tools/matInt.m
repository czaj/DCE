function MAT = matInt(X,I);

MAT = X(:,ceil((1:size(I,2)*size(X,2))/size(I,2))) .* I(:,[1:size(I,2)]' * ones(1,size(X,2)));