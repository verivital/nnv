function verify_N4_4()

    %% Construct the network
    dense = load("dense.mat");
    W = dense.W;
    b = dense.b;

    simple_rnn = load("simple_rnn_8.mat");
    rnn.bh = double(simple_rnn.bias);
    rnn.Wi = double(simple_rnn.kernel);
    rnn.Wh = double(simple_rnn.recurrent_kernel);
    rnn.fh = 'poslin';
    rnn.Wo = eye(4); % outputs equal to hidden states
    rnn.bo = zeros(4,1);
    rnn.fo = 'purelin';
    
    L1 = RecurrentLayer(rnn); % recurrent layer
    
    simple_rnn = load("simple_rnn_9.mat");
    
    rnn.bh = double(simple_rnn.bias);
    rnn.Wi = double(simple_rnn.kernel);
    rnn.Wh = double(simple_rnn.recurrent_kernel);
    rnn.fh = 'poslin';
    
    rnn.Wo = eye(4); % outputs equal to hidden states
    rnn.bo = zeros(4,1);
    rnn.fo = 'purelin';
    
    L2 = RecurrentLayer(rnn); % recurrent layer
    
    L3 = LayerS(double(W{1}),double(b{1}), 'poslin'); % feedfoward
    L4 = LayerS(double(W{2}),double(b{2}), 'poslin'); % feedfoward
    L5 = LayerS(double(W{3}),double(b{3}), 'poslin'); % feedfoward
    L6 = LayerS(double(W{4}),double(b{4}), 'poslin'); % feedfoward
    L7 = LayerS(double(W{5}),double(b{5}), 'poslin'); % feedfoward
    L8 = LayerS(double(W{6}),double(b{6}), 'purelin'); % feedfoward
    
    L = {L1, L2, L3, L4, L5, L6, L7, L8}; % all layers of the networks
    
    net = VanillaRNN(L, 'N_4_4');
    
    
    %% Create the input points & Verify the network
    data = load('points.mat');
    M = 5; % number of tested input points
    x = data.pickle_data(1:M,:); % load first M datapoints
    x = x';
    
    eps = 0.01; % adversarial disturbance bound: |xi' - xi| <= eps
    Tmax = [5 10 15 20];
    N = length(Tmax);
    rb1 = cell(M,N);
    vt1 = Inf(M,N);
    
    % Using Approximate Reachability
    for k=1:M
        for i=1:N
            input_points = [];
            for j=1:Tmax(i)
                input_points = [input_points x(:, k)];
            end
            [rb1{k, i}, vt1(k, i)] = net.verifyRBN(input_points, eps);
        end
    end
    save('N4_4_results.mat','rb1',"vt1");

end