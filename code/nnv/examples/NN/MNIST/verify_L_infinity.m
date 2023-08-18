function verify_L_infinity()

% This mnist network was trained with a single output in a regression
% manner. The output of the network is continuous, and the idx
% (label/classification) is determined by comparing idx to bounds
% i.e. if output is in between 0.5 and 1.5, corresponding label is 1
% i.e. if output is in between 1.5 and 2.5, corresponding label is 2
    
    % Load network info
    NNinfo = load('fc_6_141_mnist.mat');
    NNinfo = NNinfo.mnist_model;
    
    % Create NN
    Layers = {};
    k = 1; % layer idx in NN
    for i=1:length(NNinfo.W)
        Layers{k} = FullyConnectedLayer(NNinfo.W{i}, NNinfo.b{i}');
        Layers{k+1} = ReluLayer;
        k = k+2;
    end
    net = NN(Layers);
    
    % Load mnist data
    data = load('digits/ones.mat');
    
    % Create input
    input_vec = data.ones(5, :)';
    disturbance = 1/255;
    lb = input_vec-disturbance;
    ub = input_vec+disturbance;
    I = Star(lb, ub);
    
    % Define unsafe/not robust region
    G1 = 1;
    g1 = 0.5; 
    G2 = -1;
    g2 = -1.5;
    U1 = HalfSpace(G1,g1);
    U2 = HalfSpace(G2,g2);
    U = [U1, U2]; % unrobust region is y < 0.5 or y > 1.5
    
    % Evaluate input image
    y = net.evaluate(input_vec); % (expected: y = ?)
    
    % Verification
    reachOptions.reachMethod = 'approx-star';
    res_approx = net.verify_robustness(I, reachOptions, U);
    reachOptions.reachMethod = 'exact-star';
    res_exact = net.verify_robustness(I, reachOptions, U);
    
    fprintf('\n--------- RESULTS -----------\n');
    disp("Disturbance = "+string(disturbance));
    disp("Exact -> "+string(res_exact));
    disp("Approx -> "+string(res_approx));
    
end


