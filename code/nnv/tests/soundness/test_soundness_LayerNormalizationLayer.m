% test_soundness_LayerNormalizationLayer
% Tests for Layer Normalization layer
% To run: results = runtests('test_soundness_LayerNormalizationLayer')

%% Test 1: LayerNormalizationLayer constructor
L = LayerNormalizationLayer('Name', 'ln_test', ...
    'NumChannels', 3, ...
    'Epsilon', 1e-5, ...
    'Scale', [1.0, 1.0, 1.0], ...
    'Offset', [0.0, 0.0, 0.0]);

assert(strcmp(L.Name, 'ln_test'), 'Name should be set correctly');
assert(L.NumChannels == 3, 'NumChannels should be 3');
assert(L.Epsilon == 1e-5, 'Epsilon should match');

%% Test 2: LayerNormalizationLayer with default values
L2 = LayerNormalizationLayer('Name', 'ln_default', 'NumChannels', 4);

assert(~isempty(L2.Scale), 'Scale should be initialized');
assert(~isempty(L2.Offset), 'Offset should be initialized');
assert(all(L2.Scale == 1), 'Default scale should be ones');
assert(all(L2.Offset == 0), 'Default offset should be zeros');

%% Test 3: LayerNormalizationLayer properties
L3 = LayerNormalizationLayer('Name', 'ln_props', ...
    'NumChannels', 8, ...
    'Scale', ones(1, 8) * 2.0, ...
    'Offset', ones(1, 8) * 0.5);

assert(L3.NumChannels == 8, 'NumChannels should be 8');
assert(all(L3.Scale == 2.0), 'Scale should be 2.0');
assert(all(L3.Offset == 0.5), 'Offset should be 0.5');

