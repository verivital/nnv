function results = verify_voltage_spec(GS_out, model_data, v_min, v_max)
% verify_voltage_spec - Verify voltage magnitude bounds on GNN output
%
% Syntax:
%   results = verify_voltage_spec(GS_out, model_data, v_min, v_max)
%
% Inputs:
%   GS_out     - Output GraphStar from GNN reachability analysis
%   model_data - Loaded model struct with normalization stats
%   v_min      - Physical voltage lower bound (e.g., 0.95 p.u.)
%   v_max      - Physical voltage upper bound (e.g., 1.05 p.u.)
%
% Outputs:
%   results - Per-node array:
%              1 = verified safe (bounds within spec)
%              0 = violated (bounds completely outside spec)
%              2 = unknown (bounds cross spec boundary)
%             -1 = N/A (not a voltage-output node)
%
% Author: Anne Tumlin
% Date: 01/13/2026

    voltage_idx = 3;   % Index of voltage magnitude in output features
    bus_type_idx = 4;  % Index of bus_type in input features

    % Normalize physical voltage bounds to model space
    v_min_norm = (v_min - model_data.global_mean_labels(voltage_idx)) / ...
                 model_data.global_std_labels(voltage_idx);
    v_max_norm = (v_max - model_data.global_mean_labels(voltage_idx)) / ...
                 model_data.global_std_labels(voltage_idx);

    % Identify voltage-output nodes (bus_type == 1)
    % Denormalize input features to check bus type
    X_physical = model_data.X_test_g{1} .* model_data.global_std + model_data.global_mean;
    voltage_mask = (X_physical(:, bus_type_idx) == 1);

    numNodes = size(GS_out.V, 1);
    results = zeros(numNodes, 1);

    % Get output bounds from GraphStar
    [lb, ub] = GS_out.getRanges();

    for i = 1:numNodes
        if ~voltage_mask(i)
            results(i) = -1;  % Not applicable (not a voltage-output node)
            continue;
        end

        v_lb = lb(i, voltage_idx);
        v_ub = ub(i, voltage_idx);

        % Check if output bounds satisfy specification
        if v_lb >= v_min_norm && v_ub <= v_max_norm
            results(i) = 1;  % Verified safe
        elseif v_ub < v_min_norm || v_lb > v_max_norm
            results(i) = 0;  % Definitely violated
        else
            results(i) = 2;  % Unknown (bounds cross spec boundary)
        end
    end
end
