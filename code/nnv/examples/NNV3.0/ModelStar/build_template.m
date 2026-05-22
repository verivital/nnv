function data = build_template()
% build_template — programmatically construct the empty ModelStar
% experiment configuration. Replaces the legacy YAML template
% (MNIST_MLP_empty_results_file.yaml) so the example has zero external
% MATLAB toolbox dependencies.
%
% Returns a struct matching the YAML's structure: title/network/n_images
% scalar fields plus `layers`, a cell array of per-layer specs. Each
% layer carries a `fracs` cell array of perturbation magnitudes and the
% NaN-initialized `percents`/`times` slots that conv_expt_any_layer.m
% populates per fraction.

    data.title    = '';
    data.network  = 'MNIST_MLP';
    data.n_images = 100;
    data.position = [128, 79, 1038, 439];

    % Per-layer configuration: {name, n_weights, n_neurons, [fracs]}
    layer_specs = { ...
        'fc_6', 2560,    10,   [0.005, 0.01,  0.015, 0.02 ]; ...
        'fc_5', 65536,   256,  [0.001, 0.002, 0.003, 0.004]; ...
        'fc_4', 65536,   256,  [0.001, 0.002, 0.003, 0.004]; ...
        'fc_3', 131072,  256,  [5e-4,  0.001, 0.0015, 0.002]; ...
        'fc_2', 524288,  512,  [4e-4,  6e-4,  8e-4,  0.001 ]; ...
        'fc_1', 802816,  1024, [2e-4,  3e-4,  4e-4,  5e-4  ]; ...
    };

    data.layers = cell(size(layer_specs, 1), 1);
    for i = 1:size(layer_specs, 1)
        L = struct();
        L.name      = layer_specs{i, 1};
        L.n_weights = layer_specs{i, 2};
        L.n_neurons = layer_specs{i, 3};
        fracs       = layer_specs{i, 4};
        L.fracs     = cell(length(fracs), 1);
        for j = 1:length(fracs)
            % conv_expt_any_layer writes back num2cell([..., ..., ...])
            % into these slots, so initialize as 1×3 cell of NaN.
            L.fracs{j} = struct( ...
                'frac',     fracs(j), ...
                'percents', {{NaN, NaN, NaN}}, ...
                'times',    {{NaN, NaN, NaN}});
        end
        data.layers{i} = L;
    end
end
