function results = verify_voltage_spec(GS_out, model_data, v_min, v_max, timeout_flag)
% verify_voltage_spec - Verify voltage magnitude bounds on GNN output
%
% Uses LP-based verification (verify_specification) for precise results,
% with fallback to interval bounds when LP returns unknown.
%
% Syntax:
%   results = verify_voltage_spec(GS_out, model_data, v_min, v_max)
%   results = verify_voltage_spec(GS_out, model_data, v_min, v_max, timeout_flag)
%
% Inputs:
%   GS_out       - Output GraphStar from GNN reachability analysis
%   model_data   - Loaded model struct with normalization stats
%   v_min        - Physical voltage lower bound (e.g., 0.95 p.u.)
%   v_max        - Physical voltage upper bound (e.g., 1.05 p.u.)
%   timeout_flag - (optional) Name of base workspace variable to check for timeout
%                  If provided, checks evalin('base', timeout_flag) between nodes
%
% Outputs:
%   results - Per-node array:
%              1 = verified safe (bounds within spec)
%              0 = violated (bounds completely outside spec)
%              2 = unknown (LP returned unknown / bounds cross spec boundary)
%              3 = unknown (timeout before verification)
%             -1 = N/A (not a voltage-output node)
%
% Author: Anne Tumlin
% Date: 01/13/2026

    if nargin < 5
        timeout_flag = '';
    end

    voltage_idx = 3;   % Index of voltage magnitude in output features
    bus_type_idx = 4;  % Index of bus_type in input features

    % Normalize physical voltage bounds to model space
    v_min_norm = (v_min - model_data.global_mean_labels(voltage_idx)) / ...
                 model_data.global_std_labels(voltage_idx);
    v_max_norm = (v_max - model_data.global_mean_labels(voltage_idx)) / ...
                 model_data.global_std_labels(voltage_idx);

    % Identify voltage-output nodes (bus_type == 1)
    X_physical = model_data.X_test_g{1} .* model_data.global_std + model_data.global_mean;
    voltage_mask = (X_physical(:, bus_type_idx) == 1);

    numNodes = size(GS_out.V, 1);
    numFeatures = size(GS_out.V, 2);
    results = zeros(numNodes, 1);

    % Convert GraphStar to Star for verification
    Y_star = GS_out.toStar();

    for i = 1:numNodes
        % Check for timeout BEFORE verifying this node
        if ~isempty(timeout_flag) && evalin('base', timeout_flag)
            % Mark remaining voltage nodes as timeout-unknown
            for j = i:numNodes
                if voltage_mask(j)
                    results(j) = 3;  % Timeout before verification
                else
                    results(j) = -1;
                end
            end
            break;
        end

        if ~voltage_mask(i)
            results(i) = -1;  % Not applicable (not a voltage-output node)
            continue;
        end

        % Extract the i-th node's features as a Star set
        matIdx = zeros(1, numNodes * numFeatures);
        flat_idx = (voltage_idx - 1) * numNodes + i;
        matIdx(flat_idx) = 1;

        % Extract 1D Star for this node's voltage feature
        Y_node = Y_star.affineMap(matIdx, []);

        % Create halfspace constraints for voltage bounds
        G = [1; -1];
        g = [v_max_norm; -v_min_norm];
        Hs = [HalfSpace(G(1,:), g(1)); HalfSpace(G(2,:), g(2))];

        % LP-based verification
        res = verify_specification(Y_node, Hs);

        % Fallback to interval bounds if LP returns unknown
        if res == 2
            [lb, ub] = Y_node.getRanges;
            if lb(1) >= v_min_norm && ub(1) <= v_max_norm
                res = 1;  % Verified safe - bounds fully within spec
            elseif ub(1) < v_min_norm || lb(1) > v_max_norm
                res = 0;  % Violated - bounds fully outside spec
            else
                res = 2;  % Unknown - bounds cross spec boundary
            end
        end

        results(i) = res;
    end
end
