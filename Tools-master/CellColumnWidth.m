function [w,mw] = CellColumnWidth(X)

CellNumTrue = cellfun(@(x) isnumeric(x), X);
X(CellNumTrue) = cellfun(@(x) num2str(floor(x)),X(CellNumTrue),'UniformOutput',0);
CellNaNTrue = cellfun(@(x) strcmp(x,'NaN'), X);
X(CellNaNTrue) = {' '};
w = cellfun(@(x) numel(x), X);
mw = max(w,[],1);
