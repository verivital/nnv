function verify_robustness_classification_rnn()

    %% Construct the network
    dense = load('dense.mat');
    W = dense.W;
    b = dense.b;

    simple_rnn = load('simple_rnn.mat');
    rnn.bh = double(simple_rnn.bias);
    rnn.Wi = double(simple_rnn.kernel);
    rnn.Wh = double(simple_rnn.recurrent_kernel);
    rnn.fh = 'poslin';
    rnn.Wo = eye(2); % outputs equal to hidden states
    rnn.bo = zeros(2,1);
    rnn.fo = 'purelin';
    
    L1 = RecurrentLayer(rnn); % recurrent layer
    L2 = LayerS(double(W{1}),double(b{1}), 'poslin'); % feedfoward
    L3 = LayerS(double(W{2}),double(b{2}), 'poslin'); % feedfoward
    L4 = LayerS(double(W{3}),double(b{3}), 'poslin'); % feedfoward
    L5 = LayerS(double(W{4}),double(b{4}), 'poslin'); % feedfoward
    L6 = LayerS(double(W{5}),double(b{5}), 'poslin'); % feedfoward
    L7 = LayerS(double(W{6}),double(b{6}), 'purelin'); % feedfoward
    
    L = {L1, L2, L3, L4, L5, L6, L7}; % all layers of the networks
    
    net = NN(L);
    
    %% Create the input points & Verify the network
    data = load('points.mat');
    M = 5; % number of tested input points
    x = data.pickle_data(1:M,:); % load first M datapoints
    x = x';
    
    eps = 0.01; % adversarial disturbance bound: |xi' - xi| <= eps
    Tmax = [5 10 15 20];
    N = length(Tmax);
    result = cell(M,N);
    reachOptions.reachMethod = 'approx-star';
    net.OutputSize = 20; % set for reachability (checkRobust)
    
    % Using Approximate Reachability
    for k=1:M
        for i=1:N
            input_points = [];
            for j=1:Tmax(i)
                input_points = [input_points x(:, k)];
            end
            y = net.evaluate(input_points);
            [~, target] = max(y);
            result{k,i} = net.verify_sequence_robustness(input_points, eps, target, reachOptions);
        end
    end

    disp('--------  RESULTS  --------');
    for k=1:M
        for i=1:N
            disp("Robustness of sequence #"+string(k)+ " with Tmax of "+string(Tmax(i)) + " = "+string(result{k,i}(end)));
        end
    end

end