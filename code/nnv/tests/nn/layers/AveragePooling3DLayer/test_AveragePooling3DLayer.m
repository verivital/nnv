% Test AveragePooling3DLayer
% To run: results = runtests('test_AveragePooling3DLayer')

%% Test 1: AveragePooling3DLayer Constructor (4 args)
% Test basic constructor
name = 'test_avgpool3d';
poolSize = [2 2 2];
stride = [2 2 2];
paddingSize = [0 0 0; 0 0 0];

L = AveragePooling3DLayer(name, poolSize, stride, paddingSize);

assert(strcmp(L.Name, name));
assert(isequal(L.PoolSize, poolSize));
assert(isequal(L.Stride, stride));
assert(isequal(L.PaddingSize, paddingSize));

%% Test 2: AveragePooling3DLayer Constructor (8 args)
% Test full constructor with all parameters
name = 'test_avgpool3d_full';
poolSize = [2 3 2];
stride = [1 1 1];
paddingSize = [1 1 1; 1 1 1];
numInputs = 1;
inputNames = {'in'};
numOutputs = 1;
outputNames = {'out'};

L = AveragePooling3DLayer(name, poolSize, stride, paddingSize, ...
    numInputs, inputNames, numOutputs, outputNames);

assert(strcmp(L.Name, name));
assert(isequal(L.PoolSize, poolSize));
assert(isequal(L.Stride, stride));
assert(isequal(L.PaddingSize, paddingSize));
assert(L.NumInputs == numInputs);
assert(L.NumOutputs == numOutputs);

%% Test 3: AveragePooling3DLayer Reach with VolumeStar
% Test reachability with small 3D volume
name = 'avgpool3d_reach';
poolSize = [2 2 2];
stride = [2 2 2];
paddingSize = [0 0 0; 0 0 0];

L = AveragePooling3DLayer(name, poolSize, stride, paddingSize);

% Create a small VolumeStar input (4x4x4 volume, 1 channel)
vol_size = [4, 4, 4, 1];
vol_center = rand(vol_size);
disturbance = 0.1 * ones(vol_size);

inputVolume = VolumeStar(vol_center, -disturbance, disturbance);

% Perform reachability analysis
outputVolume = L.reach_star_single_input(inputVolume);

% Check output is VolumeStar
assert(isa(outputVolume, 'VolumeStar'));

%% Test 4: AveragePooling3DLayer Different Pool Sizes
% Test with non-square pool size
name = 'avgpool3d_nonsquare';
poolSize = [2 3 2];
stride = [2 3 2];
paddingSize = [0 0 0; 0 0 0];

L = AveragePooling3DLayer(name, poolSize, stride, paddingSize);

vol_size = [6, 6, 6, 1];
vol_center = rand(vol_size);
disturbance = 0.05 * ones(vol_size);

inputVolume = VolumeStar(vol_center, -disturbance, disturbance);
outputVolume = L.reach_star_single_input(inputVolume);

assert(isa(outputVolume, 'VolumeStar'));

%% Test 5: AveragePooling3DLayer with Padding
% Test with non-zero padding
name = 'avgpool3d_padding';
poolSize = [2 2 2];
stride = [1 1 1];
paddingSize = [1 1 1; 1 1 1];

L = AveragePooling3DLayer(name, poolSize, stride, paddingSize);

vol_size = [3, 3, 3, 1];
vol_center = rand(vol_size);
disturbance = 0.1 * ones(vol_size);

inputVolume = VolumeStar(vol_center, -disturbance, disturbance);
outputVolume = L.reach_star_single_input(inputVolume);

assert(isa(outputVolume, 'VolumeStar'));
