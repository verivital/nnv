function results = verify_voltage_spec(GS_out, model_data, v_min, v_max)
% verify_voltage_spec - Verify voltage magnitude bounds on GNN output
%
% Uses LP-based verification (verify_specification) for precise results,
% with fallback to interval bounds when LP returns unknown.
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

    % Print normalization parameters for debugging
    fprintf('\n--- Normalization Debug ---\n');
    fprintf('global_mean_labels: %s\n', mat2str(model_data.global_mean_labels', 4));
    fprintf('global_std_labels:  %s\n', mat2str(model_data.global_std_labels', 4));
    fprintf('Voltage (idx=%d): mean=%.4f, std=%.4f\n', voltage_idx, ...
            model_data.global_mean_labels(voltage_idx), model_data.global_std_labels(voltage_idx));

    % Normalize physical voltage bounds to model space
    v_min_norm = (v_min - model_data.global_mean_labels(voltage_idx)) / ...
                 model_data.global_std_labels(voltage_idx);
    v_max_norm = (v_max - model_data.global_mean_labels(voltage_idx)) / ...
                 model_data.global_std_labels(voltage_idx);

    fprintf('Physical spec: [%.2f, %.2f] p.u.\n', v_min, v_max);
    fprintf('Normalized spec: [%.4f, %.4f]\n', v_min_norm, v_max_norm);
    fprintf('Spec width (normalized): %.4f\n', v_max_norm - v_min_norm);

    % Identify voltage-output nodes (bus_type == 1)
    % Denormalize input features to check bus type
    X_physical = model_data.X_test_g{1} .* model_data.global_std + model_data.global_mean;
    voltage_mask = (X_physical(:, bus_type_idx) == 1);

    fprintf('Voltage nodes (bus_type==1): %d of %d\n', sum(voltage_mask), length(voltage_mask));

    numNodes = size(GS_out.V, 1);
    numFeatures = size(GS_out.V, 2);
    results = zeros(numNodes, 1);

    % Convert GraphStar to Star for verification
    Y_star = GS_out.toStar();

    % Print per-node voltage bounds comparison
    fprintf('\n--- Per-Node Voltage Verification (first 10 voltage nodes) ---\n');
    fprintf('%-6s %-12s %-12s %-10s\n', 'Node', 'LB', 'UB', 'Status');
    voltage_node_count = 0;

    for i = 1:numNodes
        if ~voltage_mask(i)
            results(i) = -1;  % Not applicable (not a voltage-output node)
            continue;
        end

        % Extract the i-th node's features as a Star set
        % The Star has dimension [numNodes * numFeatures]
        % We need to select node i's voltage feature
        matIdx = zeros(1, numNodes * numFeatures);

        % Index for node i, feature voltage_idx in flattened Star
        % GraphStar.toStar() uses reshape() which is column-major in MATLAB
        % So flattening is: node1_feat1, node2_feat1, ..., nodeN_feat1, node1_feat2, ...
        % flat_idx = (feature - 1) * numNodes + node
        flat_idx = (voltage_idx - 1) * numNodes + i;
        matIdx(flat_idx) = 1;

        % Extract 1D Star for this node's voltage feature
        Y_node = Y_star.affineMap(matIdx, []);

        % Create halfspace constraints for voltage bounds
        % v_min <= v <= v_max is equivalent to:
        %   v <= v_max  (upper bound)
        %   -v <= -v_min (lower bound)
        G = [1; -1];
        g = [v_max_norm; -v_min_norm];
        Hs = [HalfSpace(G(1,:), g(1)); HalfSpace(G(2,:), g(2))];

        % LP-based verification
        res = verify_specification(Y_node, Hs);

        % Fallback to interval bounds if LP returns unknown
        if res == 2
            [lb, ub] = Y_node.getRanges;
            if lb(1) >= v_min_norm && ub(1) <= v_max_norm
                res = 1;  % Verified safe
            elseif ub(1) < v_min_norm || lb(1) > v_max_norm
                res = 0;  % Violated
            end
            % Otherwise keep res = 2 (unknown)
        end

        results(i) = res;

        % Determine status string and get bounds for display
        [lb, ub] = Y_node.getRanges;
        if res == 1
            status = 'SAFE';
        elseif res == 0
            status = 'VIOLATED';
        else
            status = 'UNKNOWN';
        end

        % Print details for first 10 voltage nodes
        voltage_node_count = voltage_node_count + 1;
        if voltage_node_count <= 10
            fprintf('%-6d [%-10.4f, %-10.4f] %-10s\n', i, lb(1), ub(1), status);
        end
    end
    fprintf('Spec range: [%.4f, %.4f], width: %.4f\n', v_min_norm, v_max_norm, v_max_norm - v_min_norm);
end
