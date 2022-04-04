function subcombos = subCombos(combos)

% Arranges combos into sub arrays depending on their set_output index.
% 
% INPUTS
%
% combos: cell array containing prev advisory, output set, current advisory,
%         state set.
%
% OUTPUTS
%
% subcombos: combos divided depending on output set number. Length equal 
%            to output set length. 

last_ro_id = combos{end,2}; % number of total output_sets
subcombos = cell(1,last_ro_id);
for j = 1:size(combos,1)
    ro_id = combos{j,2};
    subcombos{ro_id} = [subcombos{ro_id};combos(j,:)];
end

end