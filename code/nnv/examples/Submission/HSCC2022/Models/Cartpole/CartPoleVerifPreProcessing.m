function [Routf,new_cmbs] = CartPoleVerifPreProcessing(init_set,outputMat,combos,agent)

% Adapts init_set to NN entry format.
%
% INPUTS
%
% init_set: state of the physical component of the system in Star form.
% outputMat: observer matrix.
% combos: vector containing prev advisory, output set, current advisory,
%         state set.
% agent: id of the agent.
%
% OUTPUTS
%
% Routf: array with the state in every possible branch.
% new_cmbs: update combos.

Rout = [];
new_cmbs = {};
Routf = [];

visited_ro_id = [];

for j = 1:size(combos,1)
    
    % Get init_set of departure for this branch
    set_id = combos{j,4};
    
    % Get Ro of departure for this branch
    ro_id = combos{j,2};
    
    if ~ismember(ro_id,visited_ro_id)
        
        visited_ro_id = [visited_ro_id ro_id];
    
        % Gets variables of interest from init_set
        Rout = init_set(set_id).affineMap(outputMat,[]);

        % Gets lb, ub of x, th, v and omega
        [mx,Mx] = Rout.getRange(1);
        [mth,Mth] = Rout.getRange(2);
        [mv,Mv] = Rout.getRange(3);
        [momega,Momega] = Rout.getRange(4);

        % Reconstruct Rout to make it suitable for NN
        Routf = [Routf Star([mx;mth;mv;momega],[Mx;Mth;Mv;Momega])];

        % Updates auxiliar variables
        ccc = {combos{j,1}, ro_id, combos{j,3}, set_id};
        new_cmbs = [new_cmbs;ccc];
        
    end
    
end
end