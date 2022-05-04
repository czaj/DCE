function output_matrix = ExtendData(input_matrix,resp_id_idx)

% Restructure dataset so that each respondent has equal number of rows (choices * alternatives), as required by the DCE package

if resp_id_idx > size(input_matrix,2)
    error('resp_id_idx must indicate unique respondent ID column number in input_matrix')
end

id_unique = unique(input_matrix(:,resp_id_idx));
id_cases = histc(input_matrix(:,resp_id_idx),id_unique);
output_matrix = NaN(max(id_cases),size(input_matrix,2),length(id_unique));
for i = 1:length(id_unique)
    output_matrix(1:id_cases(i),:,i) = input_matrix(input_matrix(:,resp_id_idx) == id_unique(i),:);
end
output_matrix = permute(output_matrix,[1,3,2]);
output_matrix = reshape(output_matrix,[max(id_cases).*length(id_unique),size(input_matrix,2)]);

