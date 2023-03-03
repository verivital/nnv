%% Construct the network
load dense.mat;
load simple_rnn_1.mat;

rnn.bh = double(bias);
rnn.Wi = double(kernel);
rnn.Wh = double(recurrent_kernel);
rnn.fh = 'poslin';

rnn.Wo = eye(4); % outputs equal to hidden states
rnn.bo = zeros(4,1);
rnn.fo = 'purelin';

L1 = RecurrentLayer(rnn); % recurrent layer
L2 = LayerS(double(W{1}),double(b{1}), 'poslin'); % feedfoward
L3 = LayerS(double(W{2}),double(b{2}), 'poslin'); % feedfoward
L4 = LayerS(double(W{3}),double(b{3}), 'poslin'); % feedfoward
L5 = LayerS(double(W{4}),double(b{4}), 'poslin'); % feedfoward
L6 = LayerS(double(W{5}),double(b{5}), 'poslin'); % feedfoward
L7 = LayerS(double(W{6}),double(b{6}), 'purelin'); % feedfoward

L = {L1, L2, L3, L4, L5, L6, L7}; % all layers of the networks

net = VanillaRNN(L, 'N_4_0');

%% Create the input points & Verify the network
load points.mat;
M = 25; % number of tested input points
x = pickle_data(1:M,:); % load first M datapoints
x = x';

eps = 0.01; % adversarial disturbance bound: |xi' - xi| <= eps
Tmax = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
N = length(Tmax);
rb1 = cell(M,N);
vt1 = Inf(M,N);
t = tic;

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

vt1 = sum(vt1);
vt1 = vt1/25; % average verification time in seconds
vt1 = vt1';

rb11 = zeros(N,1);
non_rb_id_1 = cell(N,1);
fprintf('============ Verification using Approximate Reachability ==============\n');
for i=1:N
    nrb_id = [];
    for k=1:M
        rb0 = rb1{k, i};
        if ~isempty(rb0) && rb0(Tmax(i)) == 1
            rb11(i) = rb11(i) + 1;
        else
            fprintf('No robustness for input point %d \n', k);
            nrb_id = [nrb_id k];
        end
    end
    non_rb_id_{i} = nrb_id;
end


% Using Relax Reachability RF = 0.5
rb2 = cell(M,N);
vt2 = Inf(M,N);
RF = 0.5;
for k=1:M
    for i=1:N
        input_points = [];
        for j=1:Tmax(i)
            input_points = [input_points x(:, k)];
        end
        
        [rb2{k, i}, vt2(k, i)] = net.verifyRBN(input_points, eps, 'normal', 1, RF, 'relax-star-area');
    
    end
end

vt2 = sum(vt2);
vt2 = vt2/25; % average reachability time in seconds
vt2 = vt2';

rb22 = zeros(N,1);
non_rb_id_2 = cell(N,1);
fprintf('============ Verification using Relaxed Reachability with RF = 0.5 ==============\n');
for i=1:N
    nrb_id = [];
    for k=1:M
        rb0 = rb2{k, i};
        if ~isempty(rb0) && rb0(Tmax(i)) == 1
            rb22(i) = rb22(i) + 1;
        else
            fprintf('No robustness for input point %d \n', k);
            nrb_id = [nrb_id k];
        end
    end
    non_rb_id_2{i} = nrb_id;
end

% Using Relax Reachability RF = 1
rb3 = cell(M,N);
vt3 = Inf(M,N);
RF = 1;
for k=1:M
    for i=1:N
        input_points = [];
        for j=1:Tmax(i)
            input_points = [input_points x(:, k)];
        end
        [rb3{k, i}, vt3(k, i)] = net.verifyRBN(input_points, eps, 'fast', 1, RF, 'relax-star-area');
    end
end

vt3 = sum(vt3);
vt3 = vt3/25; % average reachability time in seconds
vt3 = vt3';

rb33 = zeros(N,1);
non_rb_id_3 = cell(N,1);
fprintf('============ Verification using Relaxed Reachability with RF = 1 ==============');
for i=1:N
    nrb_id = [];
    for k=1:M
        rb0 = rb3{k, i};
        if rb0(Tmax(i)) == 1
            rb33(i) = rb33(i) + 1;
        else
            fprintf('No robustness for input point %d \n', k);
            nrb_id = [nrb_id k];
        end
    end
    non_rb_id_3{i} = nrb_id;
end


%% Print table
% load RnnVerify Result
load RnnVerify_result.mat;

RnnVerify_VT = T.rnn4_5fc32_avg_time;
RnnVerify_robust = T.rnn4_5fc32_result;
RnnV_rb = RnnVerify_robust;
RnnV_vt = RnnVerify_VT;
N_4_0 = table;
N_4_0.Tmax = Tmax';
N_4_0.RnnVerify_Robust = RnnV_rb;
N_4_0.RnnVerify_VT = RnnV_vt;
N_4_0.NNV_Robust = rb11; 
N_4_0.rc = rb11 - RnnV_rb; % conservativeness improvement
N_4_0.NNV_VT = vt1;
N_4_0.rt = RnnV_vt./vt1;% time improvement
N_4_0.NNV_RF_05_Robust = rb22; % relaxed reachability with RF = 0.5
N_4_0.NNV_RF_05_rc = rb22 - RnnV_rb; % conservativeness improvement
N_4_0.NNV_RF_05_VT = vt2; % relaxed reachability with RF = 0.5
N_4_0.NNV_RF_05_rt = RnnV_vt./vt2; % conservativeness improvement

N_4_0.NNV_RF_1_Robust = rb33; % relaxed reachability with RF = 0.5
N_4_0.NNV_RF_1_rc = rb33 - RnnV_rb; % conservativeness improvement
N_4_0.NNV_RF_1_VT = vt3; % relaxed reachability with RF = 0.5
N_4_0.NNV_RF_1_rt = RnnV_vt./vt3; % conservativeness improvement


N_4_0

%% print latex table
fileID = fopen('N_4_0_full_tab.tex','w');
formatSpec1 = '\\multirow{4}{*}{$\\mathcal{N}_{4,0}$} & $%d$ & $%d$ & $%1.2f$ & $%d$ & $%d$ & $%1.2f$ & $%1.1f\\times$ & $%d$ & $%d$ & $%1.2f$ & $%1.1f\\times$ & $%d$ & $%d$ & $%1.2f$ & $%1.1f\\times$ \\\\ \n';
formatSpec2 = ' & $%d$ & $%d$ & $%1.2f$ & $%d$ & $%d$ & $%1.2f$ & $%1.1f\\times$ & $%d$ & $%d$ & $%1.2f$ & $%1.1f\\times$ & $%d$ & $%d$ & $%1.2f$ & $%1.1f\\times$ \\\\ \n';
for i=1:N
    if i==1
        fprintf(fileID, formatSpec1, Tmax(i), RnnV_rb(i), RnnV_vt(i), rb11(i), N_4_0.rc(i), vt1(i), N_4_0.rt(i), rb22(i), N_4_0.NNV_RF_05_rc(i), vt2(i), N_4_0.NNV_RF_05_rt(i), rb33(i), N_4_0.NNV_RF_1_rc(i), vt3(i), N_4_0.NNV_RF_1_rt(i));
    else
        fprintf(fileID, formatSpec2, Tmax(i), RnnV_rb(i), RnnV_vt(i), rb11(i), N_4_0.rc(i), vt1(i), N_4_0.rt(i), rb22(i), N_4_0.NNV_RF_05_rc(i), vt2(i), N_4_0.NNV_RF_05_rt(i), rb33(i), N_4_0.NNV_RF_1_rc(i), vt3(i), N_4_0.NNV_RF_1_rt(i));
    end
    
end
fclose(fileID);

total_time = toc(t);

save N_4_0_result_full.mat N_4_0;