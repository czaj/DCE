function desctiptive_stats = stats(input,mode)

if nargin < 2 || isempty(mode) || (mode ~=  0 && mode ~= 1) % no mode specified
    mode = 0;
end

input = double(input);

% inputname(1)
% inputname(2)
% desctiptive_stats2 = @(x) inputname(1)

% [la,idx] = ismember(id_target,id_source); ...
% data_target = data_source(idx(la),:);


for k = 1 : size(input,2)
    desctiptive_stats(1:10,k*2-1) = {'mean';'median';'std';'min';'max';'quantile 0.025';'quantile 0.975'; 'sum NaN'; 'sum Inf'; 'cases'};
    desctiptive_stats(1:10,k*2) = num2cell([nanmean(input(:,k));nanmedian(input(:,k));nanstd(input(:,k));min(input(:,k));max(input(:,k));quantile(input(:,k),0.025);quantile(input(:,k),0.975);sum(isnan(input(:,k)));sum(input(:,k) == -Inf  | input(:,k) == Inf);size(input(:,k),1)]);
    if mode == 1
        desctiptive_stats(11,k*2-1:k*2) = {'levels', 'shares'};
        levels = [unique(input(~isnan(input(:,k)),k)); NaN];
        for i = 1:numel(levels)-1
            share(i,1) = sum(input(:,k)==levels(i),1)./numel(input(:,k));
        end
        share(numel(levels),1) = sum(isnan(input(:,k)))./numel(input(:,k));
        desctiptive_stats(12:11 + numel(levels),k*2-1:k*2) = num2cell([levels, share]);
        clear share levels
    end
end

