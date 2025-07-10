function Reachableset = test_Prob_class()

model_name = 'm2nist_dilated_72iou_24layer.mat';
dir0 = pwd;


parentDir = fileparts(dir0);
parentDir = fileparts(parentDir);
nnv_dir = fileparts(parentDir);
addpath(genpath(nnv_dir))
files_dir = fullfile(nnv_dir, 'engine', 'nn', 'Prob_reach', 'Temp_files_mid_run');
pe = pyenv;
py_dir = pe.Executable;
load(model_name)
layers = net.Layers;
newLayers = layers(1:end-2);
Net = SeriesNetwork(newLayers);

load("m2nist_6484_test_images.mat")
image_number = 1;
img = im_data(:,:,image_number);



delta_rgb = 1;
delta = 0.99;%0.9999;
confidence = 0.9;%0.9995;
height = 64;
width = 84;
n_channel = 1;
n_class = 11;
epochs = 50;

[N_dir , N , Ns] = CP_specification(delta, confidence, height*width*n_class , 'gpu', 'single');


im_lb = img - delta_rgb;
im_ub = img + delta_rgb;

ct = 0;
indices = zeros(height*width*n_channel , 2);
for i=1:height
    for j=1:width
        ct =ct+1;
        indices(ct , :) = [i,j];
    end
end

output_dim = [height, width, n_class];
mode = 'Linear';
model = Net;


params = struct;
params.epochs = epochs;
params.trn_batch = floor(N_dir/3);
params.dims = [-1 -1];
params.N_dir = N_dir;

params.Nt = N;
params.Ns = Ns;
params.files_dir = files_dir;

params.threshold_normal = 1e-5;
params.guarantee = delta;
params.py_dir = py_dir;


IS = ImageStar(im_lb , im_ub);

obj = ProbReach_ImageStar(model,IS,indices,output_dim,mode, params);
Reachableset = obj.ProbReach();

end