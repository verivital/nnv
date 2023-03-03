function [Routf,new_cmbs] = VCASVerifPreProcessing(init_set,outputMat,combos,agent)

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

        % Gets lb, ub of h, h'own, h'int and tau (variables of interest)
        [mh,Mh] = Rout.getRange(1);
        [mdhown,Mdhown] = Rout.getRange(2);
        [mdhint,Mdhint] = Rout.getRange(3);
        [mtau,Mtau] = Rout.getRange(4);

        % Reconstruct Rout to make it suitable for NN
        if agent == 1
            Routf = [Routf Star([mh;mdhown;mdhint;mtau],[Mh;Mdhown;Mdhint;Mtau])];
        elseif agent == 2
            Routf = [Routf Star([-Mh;mdhint;mdhown;mtau],[-mh;Mdhint;Mdhown;Mtau])];
        end

        % Updates auxiliar variables
        ccc = {combos{j,1}, ro_id, combos{j,3}, set_id}; % Combos(2) does not change because we don't split at this preproc function
        new_cmbs = [new_cmbs;ccc];
        
    end
         

    
end
end