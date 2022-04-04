function [Routf,new_cmbs] = ACASXuVerifPreProcessing(init_set,outputMat,combos,agent)

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
ro_id = 1; % new subranches
last_set_id = combos{end,4}; % number of total init_sets

if agent == 1
    
    for j = 1:size(combos,1)
        
        % Get init_set of departure for this branch
        set_id = combos{j,4};
        
        % Gets variables of interest from init_set: rho, theta, psi, vown, vint
        Rout = init_set(set_id).affineMap(outputMat,[]);

        % Gets lower and upper bounds
        [x_own_lb,x_own_ub] = Rout.getRange(1);
        [y_own_lb,y_own_ub] = Rout.getRange(2);
        [psi_own_lb,psi_own_ub] = Rout.getRange(3);
        [x_int_lb,x_int_ub] = Rout.getRange(4);
        [y_int_lb,y_int_ub] = Rout.getRange(5);
        [psi_int_lb,psi_int_ub] = Rout.getRange(6);
        [v_own_lb,v_own_ub] = Rout.getRange(7);
        [v_int_lb,v_int_ub] = Rout.getRange(8);

        x_own = [x_own_lb,x_own_ub];
        y_own = [y_own_lb,y_own_ub];
        psi_own =[psi_own_lb,psi_own_ub];
        x_int = [x_int_lb,x_int_ub];
        y_int = [y_int_lb,y_int_ub];
        psi_int = [psi_int_lb,psi_int_ub];
        v_own = [v_own_lb,v_own_ub];
        v_int = [v_int_lb,v_int_ub];

        % calculate rho
        delta_x = [x_int(1) - x_own(2), x_int(2) - x_own(1)];
        delta_y = [y_int(1) - y_own(2), y_int(2) - y_own(1)];
        delta_x_squared = [max(0.0, min(delta_x(1) * delta_x(1), min(delta_x(2) * delta_x(2), delta_x(1) * delta_x(2)))), ...
            max(delta_x(1) * delta_x(1), delta_x(2) * delta_x(2))];
        delta_y_squared = [max(0.0, min(delta_y(1) * delta_y(1), min(delta_y(2) * delta_y(2), delta_y(1) * delta_y(2)))), ...
            max(delta_y(1) * delta_y(1), delta_y(2) * delta_y(2))];
        sum_delta_squared = [delta_x_squared(1) + delta_y_squared(1), delta_x_squared(2) + delta_y_squared(2)];
        rho = [sqrt(sum_delta_squared(1)),sqrt(sum_delta_squared(2))];

        % calculate theta_o (from ownship POV)
        if delta_x(1) < 0.0 && delta_x(2) > 0.0 && delta_y(1) < 0.0 && delta_y(2) > 0.0
            angle = [-pi,pi];
        elseif delta_x(1) < 0.0 && delta_x(2) > 0.0 && delta_y(1) < 0.0
            angle_lb = min(-atan2(delta_x(1),delta_y(1)), min(-atan2(delta_x(1),delta_y(2)),...
                min(-atan2(delta_x(2),delta_y(1)), -atan2(delta_x(2),delta_y(2)))));
            angle_ub = max(-atan2(delta_x(1),delta_y(1)), max(-atan2(delta_x(1),delta_y(2)),...
                max(-atan2(delta_x(2),delta_y(1)), -atan2(delta_x(2),delta_y(2)))));
            angle = [angle_ub, angle_lb + 2*pi];
        else
            angle_lb = min(-atan2(delta_x(1),delta_y(1)), min(-atan2(delta_x(1),delta_y(2)),...
                min(-atan2(delta_x(2),delta_y(1)), -atan2(delta_x(2),delta_y(2)))));
            angle_ub = max(-atan2(delta_x(1),delta_y(1)), max(-atan2(delta_x(1),delta_y(2)),...
                max(-atan2(delta_x(2),delta_y(1)), -atan2(delta_x(2),delta_y(2)))));
            angle = [angle_lb, angle_ub];
        end
        theta_o = [angle(1) + pi/2 - psi_own(2), angle(2) + pi/2 - psi_own(1)];

        % calculate psi_o (from ownship's POV)
        psi_o = [psi_int(1) - psi_own(2), psi_int(2) - psi_own(1)];

        % keep angles in [-pi,pi] and split if necessary
        % theta_o
        theta_o_to_pi = [set_angleRange(theta_o(1)), set_angleRange(theta_o(2))];
        if theta_o_to_pi(1) > theta_o_to_pi(2)
            list_theta_o = {[-pi, theta_o_to_pi(2)],[theta_o_to_pi(1), pi]};
        else
            list_theta_o = {[theta_o_to_pi(1),theta_o_to_pi(2)]};
        end
        % psi_o
        psi_o_to_pi = [set_angleRange(psi_o(1)), set_angleRange(psi_o(2))];
        if psi_o_to_pi(1) > psi_o_to_pi(2)
            list_psi_o = {[-pi, psi_o_to_pi(2)],[psi_o_to_pi(1), pi]};
        else
            list_psi_o = {[psi_o_to_pi(1),psi_o_to_pi(2)]};
        end

        % Gets lb, ub of theta branch by branch
        for k_theta_o = 1:length(list_theta_o)    
            theta_o = list_theta_o{k_theta_o};
            % Gets lb, ub of psi branch by branch
            for k_psi_o = 1:length(list_psi_o)
                psi_o = list_psi_o{k_psi_o};

                % Reconstruct Rout to make it suitable for NN
                Routf = [Routf Star([rho(1);theta_o(1);psi_o(1);v_own(1);v_int(1)],...
                    [rho(2);theta_o(2);psi_o(2);v_own(2);v_int(2)])];

                % Updates auxiliar variables
                ccc = {combos{j,1}, ro_id, combos{j,3}, set_id};
                ro_id = ro_id + 1;
                new_cmbs = [new_cmbs;ccc];

            end
        end
        
    end
    
elseif agent == 2
    
    % Divide combos into sub combos depending on their set_id
    subcmbs = subCombos(combos);
    
    % Compute Routf
    for set_id = 1:last_set_id
        
        
        % Gets variables of interest from init_set: rho, theta, psi, vown, vint
        Rout = init_set(set_id).affineMap(outputMat,[]);
        
        % Gets lower and upper bounds
        [x_own_lb,x_own_ub] = Rout.getRange(1);
        [y_own_lb,y_own_ub] = Rout.getRange(2);
        [psi_own_lb,psi_own_ub] = Rout.getRange(3);
        [x_int_lb,x_int_ub] = Rout.getRange(4);
        [y_int_lb,y_int_ub] = Rout.getRange(5);
        [psi_int_lb,psi_int_ub] = Rout.getRange(6);
        [v_own_lb,v_own_ub] = Rout.getRange(7);
        [v_int_lb,v_int_ub] = Rout.getRange(8);
        
        x_own = [x_own_lb,x_own_ub];
        y_own = [y_own_lb,y_own_ub];
        psi_own =[psi_own_lb,psi_own_ub];
        x_int = [x_int_lb,x_int_ub];
        y_int = [y_int_lb,y_int_ub];
        psi_int = [psi_int_lb,psi_int_ub];
        v_own = [v_own_lb,v_own_ub];
        v_int = [v_int_lb,v_int_ub];
        
        % calculate rho
        delta_x = [x_int(1) - x_own(2), x_int(2) - x_own(1)];
        delta_y = [y_int(1) - y_own(2), y_int(2) - y_own(1)];
        delta_x_squared = [max(0.0, min(delta_x(1) * delta_x(1), min(delta_x(2) * delta_x(2), delta_x(1) * delta_x(2)))), ...
            max(delta_x(1) * delta_x(1), delta_x(2) * delta_x(2))];
        delta_y_squared = [max(0.0, min(delta_y(1) * delta_y(1), min(delta_y(2) * delta_y(2), delta_y(1) * delta_y(2)))), ...
            max(delta_y(1) * delta_y(1), delta_y(2) * delta_y(2))];
        sum_delta_squared = [delta_x_squared(1) + delta_y_squared(1), delta_x_squared(2) + delta_y_squared(2)];
        rho = [sqrt(sum_delta_squared(1)),sqrt(sum_delta_squared(2))];

        % calculate theta_i (from intruder's POV)
        if delta_x(1) < 0.0 && delta_x(2) > 0.0 && delta_y(1) < 0.0 && delta_y(2) > 0.0
            angle = [-pi,pi];
        elseif delta_x(1) < 0.0 && delta_x(2) > 0.0 && delta_y(1) < 0.0
            angle_lb = min(-atan2(delta_x(1),delta_y(1)), min(-atan2(delta_x(1),delta_y(2)),...
                min(-atan2(delta_x(2),delta_y(1)), -atan2(delta_x(2),delta_y(2)))));
            angle_ub = max(-atan2(delta_x(1),delta_y(1)), max(-atan2(delta_x(1),delta_y(2)),...
                max(-atan2(delta_x(2),delta_y(1)), -atan2(delta_x(2),delta_y(2)))));
            angle = [angle_ub, angle_lb + 2*pi];
        else
            angle_lb = min(-atan2(delta_x(1),delta_y(1)), min(-atan2(delta_x(1),delta_y(2)),...
                min(-atan2(delta_x(2),delta_y(1)), -atan2(delta_x(2),delta_y(2)))));
            angle_ub = max(-atan2(delta_x(1),delta_y(1)), max(-atan2(delta_x(1),delta_y(2)),...
                max(-atan2(delta_x(2),delta_y(1)), -atan2(delta_x(2),delta_y(2)))));
            angle = [angle_lb, angle_ub];
        end
        theta_i = [angle(1) + pi + pi/2 - psi_int(2), angle(2) + pi + pi/2 - psi_int(1)];
        
        % calculate psi_o (from ownship POV)
        psi_i = [psi_own(1) - psi_int(2), psi_own(2) - psi_int(1)];
        
        % keep angles in [-pi,pi] and split if necessary
        % theta_o
        theta_i_to_pi = [set_angleRange(theta_i(1)), set_angleRange(theta_i(2))];
        if theta_i_to_pi(1) > theta_i_to_pi(2)
            list_theta_i = {[-pi, theta_i_to_pi(2)],[theta_i_to_pi(1), pi]};
        else
            list_theta_i = {[theta_i_to_pi(1),theta_i_to_pi(2)]};
        end
        % psi_o
        psi_i_to_pi = [set_angleRange(psi_i(1)), set_angleRange(psi_i(2))];
        if psi_i_to_pi(1) > psi_i_to_pi(2)
            list_psi_i = {[-pi, psi_i_to_pi(2)],[psi_i_to_pi(1), pi]};
        else
            list_psi_i = {[psi_i_to_pi(1),psi_i_to_pi(2)]};
        end
        
        % Gets lb, ub of theta branch by branch
        for k_theta_i = 1:length(list_theta_i)    
            theta_i = list_theta_i{k_theta_i};
            % Gets lb, ub of psi branch by branch
            for k_psi_i = 1:length(list_psi_i)
                psi_i = list_psi_i{k_psi_i};
                
                % Reconstruct Rout to make it suitable for NN
                Routf = [Routf Star([rho(1);theta_i(1);psi_i(1);v_own(1);v_int(1)],...
                    [rho(2);theta_i(2);psi_i(2);v_own(2);v_int(2)])];
                
                % Updates ro_id
                for j = 1:size(subcmbs{set_id},1)
                    subcmbs{set_id}{j,2} = ro_id;
                end
                
                % Updates auxiliar variables
                ro_id = ro_id + 1;
                new_cmbs = [new_cmbs;subcmbs{set_id}];
                
            end
        end
        
    end
end
end
