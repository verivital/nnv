function falsify_mnist_L_infinity()

% This mnist network was trained with a single output in a regression
% manner. The output of the network is continuous, and the idx
% (label/classification) is determined by comparing idx to bounds
% i.e. if output is in between 0.5 and 1.5, corresponding label is 1
% i.e. if output is in between 1.5 and 2.5, corresponding label is 2
    
    t = tic;
    
    % Set random seed for reproducibility purposes
    rng(0); 

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
    input_vec = data.ones(10, :)';
    disturbance = 1/255;
    lb = input_vec-disturbance;
    ub = input_vec+disturbance;
%     I = Star(lb, ub);
    I = Box(lb,ub); % much faster to gnerate samples than with Star
    n_samples = 1000; % number of samples to attempt falsification
    
    % Define unsafe/not robust region
    G1 = 1;
    g1 = 0.5; 
    G2 = -1;
    g2 = -1.5;
    U1 = HalfSpace(G1,g1);
    U2 = HalfSpace(G2,g2);
    U = [U1, U2]; % unrobust region is y < 0.5 or y > 1.5
    
    % Evaluate input image
    y = net.evaluate(input_vec); % (expected: y = 0.8657)
    
    % Falsification
    counter_inputs = net.falsify(I, U, n_samples); % can prove robustness of NN with approx-star in 4-5 seconds
    
    fprintf('\n--------- RESULTS -----------\n');
    if isempty(counter_inputs)
        disp("No counter examples were found.")
    else
        disp("Property is falsified, with "+string(length(counter_inputs))+ " counter examples found");
    end
    toc(t);

% -------------  Try again with larger disturbance ------------------

    % Create input
    disturbance = 20/255; 
    lb = input_vec-disturbance;
    ub = input_vec+disturbance;
    %     I = Star(lb, ub);
    I = Box(lb,ub); % much faster to gnerate samples than with Star
    n_samples = 1000; % number of samples to attempt falsification
    
    % Falsification
    counter_inputs = net.falsify(I, U, n_samples);
    
    if isempty(counter_inputs)
        disp("No counter examples were found.")
    else
        disp("Property is falsified. A total of "+string(width(counter_inputs))+ " counter examples found");
    end
    
    toc(t);
    
end

