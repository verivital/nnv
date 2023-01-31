function [networks, names2idxs] = load_rl_NNs()
    %% 1) Load networks
    benchmarkFolder = "onnx/";
    listNN = dir(benchmarkFolder);
    networks = {}; % create a cell array of neural networks
    names = [];
    idxs = [];
    count = 1;
    t = tic;
    for h = 1:length(listNN) % generlize NN loading options for all benchmarks
        if endsWith(listNN(h).name, ".onnx")
            net = importONNXNetwork(benchmarkFolder+string(listNN(h).name), InputDataFormats="BC");
            % transform for matlab
            if ~contains(listNN(h).name, "dubins")
                Layers = net.Layers([1,4:end-1]);
                net = dlnetwork(Layers);
            else
                net = dlnetwork(net.Layers(1:end-1));
            end
            nn = matlab2nnv(net);
            % store networks
            NNs.matlab = net;
            NNs.nnv = nn;
            networks{count} = NNs;
            names{count} = listNN(h).name;
            idxs{count} = count;
            count = count + 1;
        end
    end
    t = toc(t);
    names2idxs = containers.Map(names,idxs);
    disp("All networks are loaded in " + string(t) + " seconds");
    % Remove extra files
    try 
        rmdir +cartpole s
    end
    try
        rmdir +lunarlander s
    end
end
