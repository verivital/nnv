% test_load_vnnlib
% Unit tests for VNN-LIB specification parsing
% Tests load_vnnlib function with ACAS Xu properties
% To run: results = runtests('test_load_vnnlib')

%% Test 1: Load prop_1.vnnlib
rng(42);

% Get path to VNN-LIB files
acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
vnnlib_path = [acas_path, 'vnnlib', filesep];

% Load property 1
vnnlib_file = [vnnlib_path, 'prop_1.vnnlib'];
property = load_vnnlib(vnnlib_file);

% Verify structure has required fields
assert(isfield(property, 'lb'), 'Property should have lower bounds');
assert(isfield(property, 'ub'), 'Property should have upper bounds');
assert(isfield(property, 'prop'), 'Property should have output constraints');

%% Test 2: Verify input bounds dimensions
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
vnnlib_path = [acas_path, 'vnnlib', filesep];

vnnlib_file = [vnnlib_path, 'prop_1.vnnlib'];
property = load_vnnlib(vnnlib_file);

% ACAS Xu has 5 inputs
assert(length(property.lb) == 5, 'Should have 5 input lower bounds');
assert(length(property.ub) == 5, 'Should have 5 input upper bounds');

%% Test 3: Verify bounds ordering
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
vnnlib_path = [acas_path, 'vnnlib', filesep];

vnnlib_file = [vnnlib_path, 'prop_1.vnnlib'];
property = load_vnnlib(vnnlib_file);

% Lower bounds should be <= upper bounds
assert(all(property.lb <= property.ub), 'Lower bounds should be <= upper bounds');

%% Test 4: Verify HalfSpace output constraint
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
vnnlib_path = [acas_path, 'vnnlib', filesep];

vnnlib_file = [vnnlib_path, 'prop_1.vnnlib'];
property = load_vnnlib(vnnlib_file);

% Should have output constraints
assert(~isempty(property.prop), 'Property should have output constraints');
assert(iscell(property.prop), 'Output constraints should be a cell array');

% First constraint should have HalfSpace
assert(isfield(property.prop{1}, 'Hg'), 'Property should have HalfSpace');
hs = property.prop{1}.Hg;
assert(isa(hs, 'HalfSpace'), 'Constraint should be a HalfSpace');

%% Test 5: Load prop_2.vnnlib
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
vnnlib_path = [acas_path, 'vnnlib', filesep];

vnnlib_file = [vnnlib_path, 'prop_2.vnnlib'];
property = load_vnnlib(vnnlib_file);

assert(isfield(property, 'lb'), 'Property 2 should have lower bounds');
assert(isfield(property, 'ub'), 'Property 2 should have upper bounds');
assert(length(property.lb) == 5, 'Property 2 should have 5 inputs');

%% Test 6: Load prop_3.vnnlib
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
vnnlib_path = [acas_path, 'vnnlib', filesep];

vnnlib_file = [vnnlib_path, 'prop_3.vnnlib'];
property = load_vnnlib(vnnlib_file);

assert(isfield(property, 'lb'), 'Property 3 should have lower bounds');
assert(isfield(property, 'ub'), 'Property 3 should have upper bounds');
assert(~isempty(property.prop), 'Property 3 should have output constraints');

%% Test 7: Bounds are finite
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
vnnlib_path = [acas_path, 'vnnlib', filesep];

vnnlib_file = [vnnlib_path, 'prop_1.vnnlib'];
property = load_vnnlib(vnnlib_file);

assert(all(isfinite(property.lb)), 'Lower bounds should be finite');
assert(all(isfinite(property.ub)), 'Upper bounds should be finite');

%% Test 8: Multiple properties loaded correctly
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
vnnlib_path = [acas_path, 'vnnlib', filesep];

% Load multiple properties
props = dir([vnnlib_path, '*.vnnlib']);
assert(~isempty(props), 'Should find VNN-LIB files');

% Load first 3 properties
for i = 1:min(3, length(props))
    vnnlib_file = fullfile(props(i).folder, props(i).name);
    property = load_vnnlib(vnnlib_file);

    assert(isfield(property, 'lb'), sprintf('Property %d should have lb', i));
    assert(isfield(property, 'ub'), sprintf('Property %d should have ub', i));
    assert(length(property.lb) == 5, sprintf('Property %d should have 5 inputs', i));
end

%% Test: a constraint after an (or ...) is CONJOINED into each disjunct
% Regression for the load_vnnlib bug where a standalone assertion that follows an
% (or ...) was APPENDED as a new OR disjunct instead of AND-ed into every disjunct.
% That dropped the conjunction, so a witness satisfying one disjunct but violating the
% conjunct was wrongly accepted as SAT (a false counterexample -> -150). Real-world
% case: lsnc_relu / quadrotor2d_state, where the unsafe region is (OR of Y conditions)
% AND (Y_1 in a narrow box).
rng(42);
tf = [tempname '.vnnlib'];
fid = fopen(tf, 'w');
fprintf(fid, '(declare-const X_0 Real)\n(declare-const Y_0 Real)\n(declare-const Y_1 Real)\n');
fprintf(fid, '(assert (<= X_0 1.0))\n(assert (>= X_0 -1.0))\n');
fprintf(fid, '(assert (or (and (>= Y_0 1.0)) (and (<= Y_1 -1.0))))\n');
fprintf(fid, '(assert (<= Y_0 5.0))\n');
fclose(fid);
property = load_vnnlib(tf);
delete(tf);
Hg = property.prop{1}.Hg;
assert(numel(Hg) == 2, 'two OR disjuncts should yield two HalfSpaces');
assert(size(Hg(1).G, 1) == 2, 'disjunct 1 must conjoin the trailing Y_0<=5 constraint (2 rows)');
assert(size(Hg(2).G, 1) == 2, 'disjunct 2 must conjoin the trailing Y_0<=5 constraint (2 rows)');
assert(any(Hg(1).G(:, 1) ~= 0), 'the Y_0 conjunct must be present in disjunct 1');
assert(any(Hg(2).G(:, 1) ~= 0), 'the Y_0 conjunct must be present in disjunct 2');

%% Summary
disp('test_load_vnnlib: all sections passed');

