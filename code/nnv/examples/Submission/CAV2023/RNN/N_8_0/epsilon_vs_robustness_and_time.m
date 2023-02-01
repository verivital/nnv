%% Construct the network
load dense.mat;
load simple_rnn_3.mat;

rnn.bh = double(bias);
rnn.Wi = double(kernel);
rnn.Wh = double(recurrent_kernel);
rnn.fh = 'poslin';

rnn.Wo = eye(8); % outputs equal to hidden states
rnn.bo = zeros(8,1);
rnn.fo = 'purelin';

L1 = RecurrentLayer(rnn); % recurrent layer
L2 = LayerS(double(W{1}),double(b{1}), 'poslin'); % feedfoward
L3 = LayerS(double(W{2}),double(b{2}), 'poslin'); % feedfoward
L4 = LayerS(double(W{3}),double(b{3}), 'poslin'); % feedfoward
L5 = LayerS(double(W{4}),double(b{4}), 'poslin'); % feedfoward
L6 = LayerS(double(W{5}),double(b{5}), 'poslin'); % feedfoward
L7 = LayerS(double(W{6}),double(b{6}), 'purelin'); % feedfoward

L = {L1, L2, L3, L4, L5, L6, L7}; % all layers of the networks

net = VanillaRNN(L, 'N_8_0');


%% Create the input points & Verify the network
load points.mat;
M = 25; % number of tested input points
x = pickle_data(1:M,:); % load first M datapoints
x = x';


RF = [0 0.25 0.5 0.75 1];
eps = [0.01 0.02 0.03 0.04 0.05]; % adversarial disturbance bound: |xi' - xi| <= eps
Tmax = 5;
N = length(eps);
Q = length(RF);
rb = cell(M,N,Q);
vt = Inf(M,N,Q);
t = tic;

for k=1:M
    input_points = [];
    for j=1:Tmax
        input_points = [input_points x(:, k)];
    end
    for i=1:N
        for l=1:Q
            if l==5
                [rb{k, i, l}, vt(k, i, l)] = net.verifyRBN(input_points, eps(i), 'fast', 1, RF(l), 'relax-star-area');
            else
                [rb{k, i, l}, vt(k, i, l)] = net.verifyRBN(input_points, eps(i), 'normal', 1, RF(l), 'relax-star-area');
            end
        end
    end

end

vt1 = zeros(N,Q);
for i=1:N
    for l=1:Q
        vt1(i, l) = sum(vt(:,i, l))/M;
    end
end

rb1 = zeros(N,Q);

for k=1:M
    for i=1:N
        for j=1:Q
            rb0 = rb{k, i, j};
            if rb0(Tmax) ==1
                rb1(i,j) = rb1(i,j) + 1;
            end
        end
    end
end


%% Plot figures

vt_RF_0 = vt1(:,1);
vt_RF_025 = vt1(:,2);
vt_RF_05 = vt1(:,3);
vt_RF_075 = vt1(:,4);
vt_RF_1 = vt1(:,5); 

rb_RF_0 = rb1(:,1);
rb_RF_025 = rb1(:,2);
rb_RF_05 = rb1(:,3);
rb_RF_075 = rb1(:,4);
rb_RF_1 = rb1(:,5);

load RnnVerify_eps_result.mat;

fig = figure;
subplot(1,2,1);
plot(eps, rb_RF_0, '-x');
hold on;
plot(eps, rb_RF_025, '-*');
hold on;
plot(eps, rb_RF_05, '-o');
hold on;
plot(eps, rb_RF_075, '-s');
hold on;
plot(eps, rb_RF_1, '-v');
hold on; 
plot(eps, RnnVerify_eps_rb, '-d');
legend('RF = 0', 'RF = 0.25', 'RF = 0.5', 'RF = 0.75', 'RF = 1', 'RnnVerify');
xlabel('$\epsilon$', 'Interpreter', 'latex');
ylabel('Robustness');


subplot(1,2,2);
plot(eps, vt_RF_0, '-x');
hold on;
plot(eps, vt_RF_025, '-*');
hold on;
plot(eps, vt_RF_05, '-o');
hold on;
plot(eps, vt_RF_075, '-s');
hold on;
plot(eps, vt_RF_1, '-v');
hold on; 
plot(eps, RnnVerify_eps_rb, '-d');
legend('RF = 0', 'RF = 0.25', 'RF = 0.5', 'RF = 0.75', 'RF = 1', 'RnnVerify_eps_vt');
xlabel('$\epsilon$', 'Interpreter', 'latex');
ylabel('Verification Time (s)');
saveas(fig, 'eps_rf.fig');

save eps_rf.mat vt1 rb1;
total_time = toc(t);

