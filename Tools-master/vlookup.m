function data_target = vlookup(id_target, id_source, data_source)
% INPUT: id_target, id_source, data_source

% save tmp1

[~,idx] = ismember(id_target,id_source);

data_target = NaN(size(idx,1),size(data_source,2));
data_target(idx~=0,:) = data_source(idx(idx~=0),:);



% this was supposed to work even in the case of missing id in source / target, but does not... 
% [~,idx2] = ismember(id_source,id_target);
% data_target = data_source(idx(idx~=0),:);
% data_target = zeros(size(id_target,1),size(data_source,2));
% data_target(idx2,:) = data_source(idx(idx~=0),:);

% data_target = data_source(idx(idx~=0),:);



% [la,idx] = ismember(id_target,id_source); ...
% data_target = data_source(idx(la),:);