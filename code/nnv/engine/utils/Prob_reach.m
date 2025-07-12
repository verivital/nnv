function Out_ImS = Prob_reach(Net, In_ImS, reachOptions)

pe = pyenv;
py_dir = pe.Executable;


if isfield(reachOptions, 'coverage')
    coverage = reachOptions.coverage;
else
    coverage = 0.99;
end
if isfield(reachOptions, 'confidence')
    confidence = reachOptions.confidence;
else
    confidence = 0.99;
end
if isempty(In_ImS.im_lb)
    error('We assume the input ImageStar is a box, and also contains the feature, im_lb. In case your input is not a box but contains im_lb, then your reachset will be more conservative as we assume a box with lower bound im_lb and upper bound im_ub.')
end

if isfield(reachOptions, 'train_device')
    train_device = reachOptions.train_device;
else
    train_device = 'gpu';
end
if isfield(reachOptions, 'train_epochs')
    train_epochs = reachOptions.train_epochs;
else
    train_epochs = 50;
end

if isfield(reachOptions, 'train_lr')
    train_lr = reachOptions.train_lr;
else
    train_lr = 0.01;
end

[N_dir , N , Ns] = CP_specification(coverage, confidence, numel(In_ImS.im_lb) , train_device, 'single');

SizeIn = size(In_ImS.im_lb);
SizeOut = size(forward(Net, In_ImS.im_lb));
height = SizeIn(1);
width = SizeIn(2);
if length(SizeOut) == 2
    SizeOut = [ SizeOut , 1];
end


if isfield(reachOptions, 'indices')
    indices = reachOptions.indices;
else
    [J,I] = ndgrid(1:width,1:height);
    indices = [I(:), J(:)];
end

if isfield(reachOptions, 'train_mode')
    train_mode = reachOptions.train_mode;
else
    train_mode = 'Linear';
end

if isfield(reachOptions, 'surrogate_dim')
    surrogate_dim = reachOptions.surrogate_dim;
else
    surrogate_dim = [-1, -1];
end

if isfield(reachOptions, 'threshold_normal')
    threshold_normal = reachOptions.threshold_normal;
else
    threshold_normal = 1e-5;
end

params = struct;
params.epochs = train_epochs;
params.lr = train_lr;
params.trn_batch = floor(N_dir/3);
params.dims = surrogate_dim;
params.N_dir = N_dir;
params.Nt = N;
params.Ns = Ns;
params.threshold_normal = threshold_normal;
params.guarantee = coverage;
params.py_dir = py_dir;


obj = ProbReach_ImageStar(Net,In_ImS,indices,SizeOut,train_mode, params);
Out_ImS = obj.ProbReach();

end





function out = forward(model, x)
    
    model_source = class(model);
    
    switch model_source
    
        case 'SeriesNetwork'
            out = model.predict(x);
    
        case 'DAGNetwork'
            out = model.predict(x);
    
        case 'dlnetwork'
            dlX = dlarray(X, 'SSC'); % 'SSC' = Spatial, Spatial, Channel
            out = model.predict(dlX);
    
        case 'NN'
            out = model.evaluate(x);
    
        otherwise
            error("Unknown model source: " + model_source + ". We only cover NN, SeriesNetwork, dlnetwork and DAGNetwork.");
    end
end