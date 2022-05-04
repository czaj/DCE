function MAT = reprc(X,n,k)

MAT = X(ceil((1:n*size(X,1))/n),:);...
MAT = MAT(:,ceil((1:k*size(MAT,2))/k));
