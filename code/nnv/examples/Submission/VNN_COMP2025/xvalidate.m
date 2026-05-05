function diff = xvalidate(onnx_path, mat_path, n_samples)
%XVALIDATE  Compare onnxruntime vs NNV on N random inputs.
% Returns the max absolute difference across all samples and outputs.
% Returns NaN if either side fails.
%
% Calls Python's onnxruntime via system() since MATLAB doesn't have a
% built-in ONNX runner that bypasses importNetworkFromONNX (the very thing
% we're trying to bypass).

if nargin < 3, n_samples = 5; end

try
    net = load_nnv_from_mat(mat_path);
catch ME
    fprintf('  load FAIL: %s\n', ME.message);
    diff = NaN;
    return;
end

% Determine input size
in_layer = net.Layers{1};
if isa(in_layer, 'FeatureInputLayer')
    in_shape = in_layer.InputSize;
    if isscalar(in_shape), in_shape = [in_shape 1]; end
elseif isa(in_layer, 'ImageInputLayer')
    in_shape = in_layer.InputSize;
else
    error('Unknown input layer: %s', class(in_layer));
end

n = prod(in_shape);
rng(0);

% Use a single Python invocation that processes all samples via an inputs.npy
% file to avoid Windows command-line-length limits on large input dims.
inputs_file = [tempname '.npy'];
output_file = [tempname '.npy'];

% Generate samples already in MATLAB-native shape, then write a 4D array
% that Python can permute into ONNX BCHW order:
%   MATLAB ImageInputLayer expects [H, W, C, B]
%   ONNX expects [B, C, H, W]
% For each sample we keep the MATLAB order on disk and transpose in Python.
is_image = isa(net.Layers{1}, 'ImageInputLayer');
if is_image
    sz = in_layer.InputSize;   % [H W C]
    samples_mat = single(randn([sz n_samples], 'single'));   % [H W C N]
else
    samples_mat = single(randn([n n_samples], 'single'));    % [n N]
end

% NNV side
y_nnv_all = cell(1, n_samples);
for k = 1:n_samples
    if is_image
        x = samples_mat(:,:,:,k);
    else
        x = samples_mat(:, k);
    end
    try
        y = double(net.evaluate(x));
        % If output is multi-dim (image), normalize to a canonical
        % ordering that matches Python's `y.flatten()` (row-major BCHW).
        % NNV's image output is HWC; permute to CHW then flatten in
        % MATLAB column-major (which equals C-order on a CHW tensor
        % only when we reshape correctly).
        if ~ismatrix(y)
            % Numpy row-major flatten of BCHW [1,C,H,W] iterates
            % W-fastest, then H, then C. MATLAB column-major flatten
            % of [W,H,C] also iterates W-fastest, then H, then C, so
            % these two flat vectors match elementwise.
            y_perm = permute(y, [2 1 3]);     % HWC -> WHC
            y_nnv_all{k} = y_perm(:);
        else
            y_nnv_all{k} = y(:);
        end
    catch ME
        fprintf('  nnv evaluate FAIL: %s\n', ME.message);
        diff = NaN;
        return;
    end
end

% Reshape samples_mat to a 2D matrix [N, prod(per-sample-dims)] for npy output
if is_image
    samples = reshape(permute(samples_mat, [4 1 2 3]), n_samples, []);  % N x (H*W*C)
else
    samples = samples_mat.';   % N x n
end

% Save inputs as flat-row matrix (N x flat_features). Each row's flat layout
% MATLAB-native HWC if image, plain feature vector if FC.
np_save_npy(inputs_file, samples);

% Build a Python driver that reshapes the flat row appropriately.
% For image inputs, samples row k holds H*W*C floats in HWC layout (matlab
% column-major flatten); reshape via numpy to [H, W, C] (using order='F'
% which matches MATLAB column-major), then transpose to ONNX [C, H, W] and
% prepend batch axis.
if is_image
    % Peek at the ONNX input shape to decide BCHW vs BHWC layout. If the
    % model declares [B, H, W, C] (TF-export style), feed it BHWC instead
    % of BCHW.
    py_code = sprintf(['import sys, numpy as np, onnxruntime as ort\n', ...
        'sess = ort.InferenceSession(r"%s")\n', ...
        'inp = sess.get_inputs()[0]\n', ...
        'samples = np.load(r"%s").astype(np.float32)\n', ...
        'H, W, C = %d, %d, %d\n', ...
        'outs = []\n', ...
        'sh = inp.shape\n', ...
        'is_bhwc = (len(sh) == 4 and isinstance(sh[3], int) and sh[3] == C and isinstance(sh[1], int) and sh[1] == H)\n', ...
        'for s in samples:\n', ...
        '    hwc = s.reshape((H, W, C), order="F")\n', ...
        '    if is_bhwc:\n', ...
        '        feed = hwc[None, ...]\n', ...
        '    else:\n', ...
        '        feed = hwc.transpose(2, 0, 1)[None, ...]\n', ...
        '    y = sess.run(None, {inp.name: feed.astype(np.float32)})[0]\n', ...
        '    outs.append(np.asarray(y).flatten())\n', ...
        'np.save(r"%s", np.stack(outs).astype(np.float32))\n'], ...
        onnx_path, inputs_file, sz(1), sz(2), sz(3), output_file);
else
    py_code = sprintf(['import sys, numpy as np, onnxruntime as ort\n', ...
        'sess = ort.InferenceSession(r"%s")\n', ...
        'inp = sess.get_inputs()[0]; sh = inp.shape\n', ...
        'sh = [d if isinstance(d,int) and d>0 else 1 for d in sh]\n', ...
        'samples = np.load(r"%s").astype(np.float32)\n', ...
        'outs = []\n', ...
        'for s in samples:\n', ...
        '    x = s.reshape(sh)\n', ...
        '    y = sess.run(None, {inp.name: x})[0]\n', ...
        '    outs.append(np.asarray(y).flatten())\n', ...
        'np.save(r"%s", np.stack(outs).astype(np.float32))\n'], ...
        onnx_path, inputs_file, output_file);
end
script_file = [tempname '.py'];
fidp = fopen(script_file, 'w'); fprintf(fidp, '%s', py_code); fclose(fidp);

[s, o] = system(sprintf('python "%s"', script_file));
if s ~= 0
    fprintf('  onnxruntime FAIL: %s\n', strtrim(o));
    diff = NaN; cleanup(inputs_file, output_file, script_file); return;
end
if ~isfile(output_file)
    diff = NaN; cleanup(inputs_file, output_file, script_file); return;
end
y_ort_all = np_load_npy(output_file);

diff = 0;
for k = 1:n_samples
    y_ort = double(y_ort_all(k, :)).';
    y_nnv = y_nnv_all{k};
    if numel(y_ort) ~= numel(y_nnv)
        fprintf('  shape mismatch: nnv=%d, ort=%d\n', numel(y_nnv), numel(y_ort));
        diff = NaN;
        cleanup(inputs_file, output_file, script_file); return;
    end
    d = max(abs(y_ort - y_nnv));
    if d > diff, diff = d; end
end
cleanup(inputs_file, output_file, script_file);
end


function cleanup(varargin)
    for i = 1:nargin
        if isfile(varargin{i}), delete(varargin{i}); end
    end
end

function np_save_npy(path, arr)
% Minimal .npy writer for 2D float32 arrays
arr = single(arr);
shape = size(arr);
header = sprintf('{''descr'': ''<f4'', ''fortran_order'': False, ''shape'': (%d, %d), }', shape(1), shape(2));
% pad to 16-byte boundary including 10-byte preamble + 1 newline
total = 10 + numel(header) + 1;
pad_len = 16 - mod(total, 16); if pad_len == 16, pad_len = 0; end
header = [header, repmat(' ', 1, pad_len)];
header_bytes = uint8(header);
header_len_field = uint16(numel(header_bytes) + 1);
fid = fopen(path, 'w');
fwrite(fid, [uint8(147) uint8('NUMPY')], 'uint8');     % magic
fwrite(fid, uint8([1 0]), 'uint8');                     % version
fwrite(fid, header_len_field, 'uint16');                % header_len (little endian)
fwrite(fid, header_bytes, 'uint8');
fwrite(fid, uint8(10), 'uint8');                        % \n
% Numpy default order is C-style (row-major); MATLAB writes column-major.
% Write transposed so the on-disk layout matches NumPy's expected row-major.
fwrite(fid, single(arr.'), 'single');
fclose(fid);
end

function arr = np_load_npy(path)
% Minimal .npy reader for 2D float32 arrays
fid = fopen(path, 'r');
magic = fread(fid, 6, 'uint8=>uint8');
if ~(magic(1)==147 && all(char(magic(2:6)).' == 'NUMPY'))
    fclose(fid); error('Not a NPY file');
end
ver = fread(fid, 2, 'uint8=>uint8'); %#ok<NASGU>
hlen = fread(fid, 1, 'uint16');
header = char(fread(fid, hlen, 'uint8=>uint8').');
% parse shape (we know format)
tok = regexp(header, '''shape'':\s*\((\d+),\s*(\d+)', 'tokens', 'once');
nr = str2double(tok{1}); nc = str2double(tok{2});
data = fread(fid, nr*nc, 'single=>single');
fclose(fid);
arr = reshape(data, nc, nr).';   % numpy row-major -> MATLAB
end
