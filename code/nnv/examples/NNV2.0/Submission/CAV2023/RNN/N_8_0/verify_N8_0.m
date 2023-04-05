function verify_N8_0()

    %% Construct the network
    dense = load("dense.mat");
    W = dense.W;
    b = dense.b;

    simple_rnn = load("simple_rnn_3.mat");
    
    rnn.bh = double(simple_rnn.bias);
    rnn.Wi = double(simple_rnn.kernel);
    rnn.Wh = double(simple_rnn.recurrent_kernel);
    rnn.fh = 'poslin';
    
    rnn.Wo = eye(8); % outputs equal to hidden states
    rnn.bo = zeros(8,1);
    rnn.fo = 'purelin';
    
    L1 = RecurrentLayer(rnn); % recurrent layer
    L2 = FullyConnectedLayer(double(W{1}),double(b{1}));
    L2a = ReluLayer();
    L3 = FullyConnectedLayer(double(W{2}),double(b{2}));
    L3a = ReluLayer();
    L4 = FullyConnectedLayer(double(W{3}),double(b{3}));
    L4a = ReluLayer();
    L5 = FullyConnectedLayer(double(W{4}),double(b{4}));
    L5a = ReluLayer();
    L6 = FullyConnectedLayer(double(W{5}),double(b{5}));
    L6a = ReluLayer();
    L7 = FullyConnectedLayer(double(W{6}),double(b{6}));
    
    L = {L1, L2, L2a, L3, L3a, L4, L4a, L5, L5a, L6, L6a, L7}; % all layers of the networks
    net = NN(L);
    
    
    %% Create the input points & Verify the network
    data = load("points.mat");
    M = 5; % number of tested input points
    x = data.pickle_data(1:M,:); % load first M datapoints
    x = x';
    
    eps = 0.01; % adversarial disturbance bound: |xi' - xi| <= eps
    Tmax = [5 10 15 20];
    N = length(Tmax);
    rb1 = cell(M,N);
    vt1 = Inf(M,N);
    
    % Using Approximate Reachability
    reachOptions.reachMethod = 'approx-star';
    for k=1:M
        for i=1:N
            % Create input
            input_points = [];
            for j=1:Tmax(i)
                input_points = [input_points x(:, k)];
            end

            % verify reach sets compare with groundtruth, i.e., non-attacked signal
            y = net.evaluate(input_points);
            [~,max_id] = max(y); % find the classified output

            % reach + robustness verification
            t = tic;
            rb1{k, i} = net.verify_sequence_robustness(input_points, eps, max_id, reachOptions);
            vt1(k, i) = toc(t); % time verification
        end
    end
    save('N8_0_results.mat','rb1',"vt1");

end
