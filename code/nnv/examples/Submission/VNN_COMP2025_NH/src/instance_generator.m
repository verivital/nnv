function case_instances = instance_generator(directory)

% List ONNX files
onnx_dir = fullfile(directory, 'onnx');
onnx_files_struct = dir(fullfile(onnx_dir, '*.onnx'));
onnx_files = {onnx_files_struct.name};

% List VNNLIB files
vnnlib_dir = fullfile(directory, 'vnnlib');
vnnlib_files_struct = dir(fullfile(vnnlib_dir, '*.vnnlib'));
vnnlib_files = {vnnlib_files_struct.name};

% Initialize cell array to hold pairs
instance_pairs = {};

% Generate all combinations
for i = 1:numel(onnx_files)
    for j = 1:numel(vnnlib_files)
        % Only relative paths
        onnx_path = "onnx/" + onnx_files{i};
        vnnlib_path = "vnnlib/" + vnnlib_files{j};
        instance_pairs(end+1, :) = {onnx_path, vnnlib_path};
    end
end

% Convert to string array
case_instances = string(instance_pairs);

end