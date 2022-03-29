function new_cmbs = preSetCombos(combos)

% Adapts the values of combos when agent > 1. For each new agent, Ro is
% "initialized" so that the splits from previous agents do not interfere
% with the possible new splits. 
%
% INPUTS
%
% combos: vector containing prev advisory, output set, current advisory,
%         state set.
%
% OUTPUTS
%
% new_cmbs: update combos.

new_cmbs = {};
for j = 1:size(combos,1)
    % Get init_set of departure for this branch
    set_id = combos{j,4};
    
    % Updates variables
    ccc = {combos{j,1}, set_id, combos{j,3}, set_id};
    new_cmbs = [new_cmbs;ccc];
    
end

end