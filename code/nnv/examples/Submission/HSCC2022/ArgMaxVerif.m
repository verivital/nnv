function [idxout,new_cmbs] = ArgMaxVerif(yNN,combos,agent)

% Get max index for networks (advisory command)
% 
% INPUTS
%
% yNN: cell array with the weighted vector of all possible branches. 
% combos: vector containing prev advisory, output set, current advisory,
%         state set.
% agent: id of the agent
%
% OUTPUTS
%
% idxs: cell with the indexes of the max value in the weigth vector for
%       each init_set.
% new_cmbs: update combos.

new_cmbs = {};
idxout = {};

% Divide combos into sub combos depending on their set_id
subcmbs = subCombos(combos);

for ro_id = 1:size(combos,1)
    xSs = yNN{ro_id}; % takes the array from the current branch
    if length(yNN) > 1
        Xi = Star.merge_stars(xSs,1,'single');
    else
        Xi = xSs;
    end

    [mm,MM] = Xi.getRanges;
    [~, idxmax] = max(MM); 
    new_idx = [];
    % position of the maximum value among all intervals
    % checks if there are overlapses with the interval that contains the
    % maximum value. Fixes the minimum value of the interval containing the
    % maximum value and checks if any local maximum is higher than such
    % value.
    for idxvar = 1:length(MM)
        if MM(idxvar) >= mm(idxmax)
            new_idx = [new_idx idxvar]; 
        end
    end
    
    idxout = [idxout new_idx];
    
    for i = 1:length(new_idx)
        for j = 1:size(subcmbs{ro_id},1)
            subcmbs{ro_id}{j,3}(agent) = new_idx(i);
        end
        new_cmbs = [new_cmbs; subcmbs{ro_id}];
    end  
end
end