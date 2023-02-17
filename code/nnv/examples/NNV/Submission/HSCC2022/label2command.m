function star_com = label2command(label,listofcommands)

% Selects which physical command send to the plant depending on the label 
% obtained on the current step. Equivalent to PostProcessing function but
% with another "structure treatment" for the labels.
% 
% INPUTS
%
% label: vector or scalar containing the label.
% lisofcommands: list of possible commands the controller might take.
%
% OUTPUTS
%
% star_com: star containing all of the selected commands which are grouped
% in array form. length = numagents.

array_com = [];
for i=1:length(label)  
    array_com = [array_com;listofcommands{label(i)}]; 
end
star_com = Star(array_com,array_com);

end