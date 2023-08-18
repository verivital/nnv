function get_robustness_bound_mnist()

% This mnist network was trained with a single output in a regression
% manner. The output of the network is continuous, and the idx
% (label/classification) is determined by comparing idx to bounds
% i.e. if output is in between 0.5 and 1.5, corresponding label is 1
% i.e. if output is in between 1.5 and 2.5, corresponding label is 2

    t = tic;

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
    input_vec = data.ones(1, :)';
    n = length(input_vec);
    
    % Define parameters to find robustness bounds
    init_dis_bound = 1/255; % initial bound of disturbance
    tol = 1/255; % accuracy in finding the maximum bound
    % ------- why do we choose 1/255 ? ----------------------
    % mnist images are defined to be uint8 with range 0 to 255, so this
    % 1/255 value means we are "perturbing" the image the minimum possible
    % change that we can to a pixel
    %--------------------------------------------------------
    max_steps = 10; % maximum number of searching step
    lb_allowable = zeros(n, 1); % input value must be within 0 and 1
    ub_allowable = ones(n, 1);
    
    % Define unsafe/not robust region
    G1 = 1;
    g1 = 0.5; 
    G2 = -1;
    g2 = -1.5;
    U1 = HalfSpace(G1,g1);
    U2 = HalfSpace(G2,g2);
    un_robust_reg = [U1, U2]; % unrobust region is y < 0.5 or y > 1.5
    
    % Evaluate input image
    y = net.evaluate(input_vec); % (expected: y = 0.9394)
    
    % Define reachability parameters
    reachOptions.reachMethod = 'approx-star';
    
    % Search for robustness bound based on the approx-star method
    robustness_bound = net.get_robustness_bound(input_vec, init_dis_bound, tol, max_steps,...
        lb_allowable, ub_allowable, un_robust_reg, reachOptions);
    
    disp('--------- APPROX-STAR RESULTS -----------');
    if isempty(robustness_bound)
        disp("Robustness bound not found.")
    else
        disp("Robustness bound:")
        disp(robustness_bound)
    end
    
    
    % Define reachability parameters
    reachOptions.reachMethod = 'approx-zono';
    
    % Search for robustness bound based on the approx-star method
    robustness_bound = net.get_robustness_bound(input_vec, init_dis_bound, tol, max_steps,...
        lb_allowable, ub_allowable, un_robust_reg, reachOptions);
    
    disp('--------- APPROX-ZONO RESULTS -----------');
    if isempty(robustness_bound)
        disp("Robustness bound not found.")
    else
        disp("Robustness bound:")
        disp(robustness_bound)
    end
    
    toc(t);

end
