
%% Construct the network
load dense.mat;
load simple_rnn.mat;

rnn.bh = double(bias);
rnn.Wi = double(kernel);
rnn.Wh = double(recurrent_kernel);
rnn.fh = 'poslin';

rnn.Wo = double(W); % outputs 
rnn.bo = double(b);
rnn.fo = 'poslin';

L1 = RecurrentLayer(rnn); % recurrent layer
L = {L1}; % all layers of the networks

net = VanillaRNN(L, 'Stanford_net_1');

%% Construct the input sets and input points
load points.mat;
M = 25; % number of tested input points
x = pickle_data(1:M,:); % load first M datapoints
x =  x';
eps = 0.01; % adversarial disturbance bound: |xi' - xi| <= eps


% input sets corresponding to each random input points
X = [];
for i=1:M
    X1 = Star(x(:,i) - eps, x(:,i) + eps);
    X = [X X1];
end

% compute reachable set
%Tmax = [2]; % for testing only
Tmax = [2 5 10 15 20];
%Tmax = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20, 21, 22, 23, 24, 25];
N = length(Tmax);
inputSet = cell(M, N); % input sets for reachability
inputPoints = cell(M, N); % input points for ground-truth evaluation
reachSet = cell(M, N); % reachability set
rt = zeros(M, N); 
vt = zeros(M, N);

for i=1:M
    for j=1:N
        inputSet1 = [];
        inputPoints1 = [];
        for k=1:Tmax(j)
            inputSet1 = [inputSet1 X(i)];
            inputPoints1 = [inputPoints1 x(:, i)];
        end
        inputSet{i, j} = inputSet1;
        inputPoints{i,j} = inputPoints1; 
    end
end



%% Perform Verification with Exact Reachability
% Reachability Analysis
numCores = 1;
outputSet = cell(M, N);

for i=1:M
    for j=1:N
        t = tic;
        R = L1.reach(inputSet{i,j}, 'exact-star');
        outputSet{i,j} = R{Tmax(j)};
        rt(i, j) = toc(t);
    end
end

% checking robustness
ct = zeros(M, N); % checking time
ex_rb = zeros(M,N); % 1 = robust, 0 = not robust
output_CE = cell(M,N); % counter examples
NumCE = zeros(M, N); % number of counter examples;
NStars = zeros(M, N); % number of stars in the output set
ground_truth_ids = zeros(M, N); % ground truth label index
for i=1:M
    for j=1:N
        t = tic;
        O = outputSet{i,j};
        K = length(O);
        NStars(i,j) = K;
        y = net.evaluate(inputPoints{i,j});
        l = 0;
        [~,max_id] = max(y(:, Tmax(j))); % find the classified output
        ground_truth_ids(i,j) = max_id; % groundtruth label
        CE1 = [];
        for k=1:K
            max_cands = O(k).get_max_point_candidates;
            has_CE = 0;
                for l=1:length(max_cands)
                    if (max_cands(l) ~= max_id) && O(k).is_p1_larger_than_p2(max_cands(l), max_id)
                        fprintf("\n xi = x%d, Tmax = %d: Counterexample is found, label = % d, ground-truth label = %d", i, Tmax(j), max_cands(l), max_id);
                        has_CE = 1;
                    end
                end
                if has_CE
                    CE1 = [CE1 O(k)];
                end
        end
        ex_rb(i,j) = isempty(CE1); % if no counter example set is found, the network is robust
        output_CE{i,j} = CE1;
        NumCE(i,j) = length(CE1);
        ct(i,j) = toc(t);
    end
end

vt1 = ct + rt; % total verification time using exact analysis;


%% Perform Verification with Approximate Reachability

% verify with approx-star
approx_rb = zeros(M, N);
vt2 = zeros(M, N); % verification time using approximate analysis
t = tic;
for k=1:M
    for i=1:N
        input_points = [];
        for j=1:Tmax(i)
            input_points = [input_points x(:, k)];
        end
        [rb, vt2(k, i)] = net.verifyRBN(input_points, eps);
        approx_rb(k, i) = rb(Tmax(i));
    end
end
total_time = toc(t);



%% Print Tables of Results

x = {'x1'; 'x2'; 'x3'; 'x4'; 'x5'; 'x6'; 'x7'; 'x8'; 'x9'; 'x10'; 'x11'; 'x12'; 'x13'; 'x14'; 'x15'; 'x16'; 'x17'; 'x18'; 'x19'; 'x20'; 'x21'; 'x22'; 'x23'; 'x24'; 'x25'};

%Tmax = 2
Robust = ex_rb(:,1);
VT = vt1(:,1);
N = NStars(:,1);
CE = NumCE(:,1);

Tmax_2 = table;
Tmax_2.x = x;
Tmax_2.Robust = Robust;
Tmax_2.VT = VT;
Tmax_2.N = N;
Tmax_2.CE = CE;

%Tmax = 5
Robust = ex_rb(:,2);
VT = vt1(:,2);
N = NStars(:,2);
CE = NumCE(:,2);

Tmax_5 = table;
Tmax_5.x = x;
Tmax_5.Robust = Robust;
Tmax_5.VT = VT;
Tmax_5.N = N;
Tmax_5.CE = CE;

%Tmax = 10
Robust = ex_rb(:,3);
VT = vt1(:,3);
N = NStars(:,3);
CE = NumCE(:,3);

Tmax_10 = table;
%Tmax_10.x = x;
Tmax_10.Robust = Robust;
Tmax_10.VT = VT;
Tmax_10.N = N;
Tmax_10.CE = CE;

%Tmax = 15
Robust = ex_rb(:,4);
VT = vt1(:,4);
N = NStars(:,4);
CE = NumCE(:,4);

Tmax_15 = table;
Tmax_15.x = x;
Tmax_15.Robust = Robust;
Tmax_15.VT = VT;
Tmax_15.N = N;
Tmax_15.CE = CE;

%Tmax = 20
Robust = ex_rb(:,5);
VT = vt1(:,5);
N = NStars(:,5);
CE = NumCE(:,5);

Tmax_20 = table;
Tmax_20.x = x;
Tmax_20.Robust = Robust;
Tmax_20.VT = VT;
Tmax_20.N = N;
Tmax_20.CE = CE;

fprintf("=======================EXACT VERIFICATION RESULTS WITH DIFFERENT Tmax ===============================\n");
Tmax_2
Tmax_5
Tmax_10
Tmax_15
Tmax_20

N = length(Tmax);
for i=1:M
    for j=1:N
        NStars(i,j) = length(outputSet{i,j});
        NumCE(i,j) = length(output_CE{i,j});
    end
end

%% Print latex table
fileID = fopen('exact_verification_tab.tex','w');
formatSpec = '$x_{%d}$ & $%d$ & $%2.2f$ & $%d$ & $%d$ & $%d$ & $%2.2f$ & $%d$ & $%d$ & $%d$ & $%2.2f$ & $%d$ & $%d$ & $%d$ & $%2.2f$ & $%d$ & $%d$ \\\\ \n';
for i=1:25
    rb5 = ex_rb(i, 2); vt5 = vt1(i, 2); N5 = NStars(i, 2); CE5 = NumCE(i, 2);
    rb10 = ex_rb(i, 3); vt10 = vt1(i, 3); N10 = NStars(i, 3); CE10 = NumCE(i, 3);
    rb15 = ex_rb(i, 4); vt15 = vt1(i, 4); N15 = NStars(i, 4); CE15 = NumCE(i, 4);
    rb20 = ex_rb(i, 5); vt20 = vt1(i, 5); N20 = NStars(i, 5); CE20 = NumCE(i, 5);
    fprintf(fileID, formatSpec, i, rb5, vt5, N5, CE5, rb10, vt10, N10, CE10, rb15, vt15, N15, CE15, rb20, vt20, N20, CE20);
end
fclose(fileID);

%% Plot ranges of counter example output sets

CE5 = output_CE{1,5};
max_id = ground_truth_ids(1,5);
CE = CE5(1); 
[lb, ub] = CE.getRanges; 
center = (lb + ub)/2;
err = (ub-lb)/2;
x = 1:1:20;
y = center; 

y1 = lb(max_id)*ones(20,1); % lower bound of ground truth output label
y2 = ub(10)*ones(20, 1); % upper bound of counter output label

fig = figure;
e = errorbar(x,y,err);
e.LineStyle = 'none';
e.LineWidth = 1;
e.Color = 'red';
xlabel('Output', 'FontSize', 11);
ylabel('Ranges', 'FontSize', 11);
xlim([0 21]);
xticks([1 10 20]);
xticklabels({'1', '10', '20'});
set(gca, 'FontSize', 10);

hold on; 
plot(x, y1, '--', 'color', 'blue');
hold on;
plot(x, y2, '--', 'color', 'blue');

p1 = nsidedpoly(1000, 'Center', [max_id center(max_id)], 'Radius', 0.5);
hold on;
plot(p1, 'FaceColor', 'blue');

p2 = nsidedpoly(1000, 'Center', [10 center(10)], 'Radius', 0.3);
hold on;
plot(p2);
saveas(fig, 'counterExampleRange.fig');

%% Plot verification time performance of different approaches

load model_marabou_rnsverify_compare_26_20210915-123159.mat

RnnVerify = [our(1) our(4) our(9) our(14) our(19)]; % Tmax = 2, Tmax = 5, Tmax = 10, Tmax = 15, Tmax = 20
RNSVerify = [rns(1) rns(4) rns(9) rns(14) rns(19)]; % Tmax = 2, Tmax = 5, Tmax = 10, Tmax = 15, Tmax = 20

ex_vt = sum(vt1)/25;
approx_vt = sum(vt2)/25;

Tmax = [2 5 10 15 20];
fig = figure;
plot(Tmax, ex_vt, '--*'); % our average exact verification time
hold on;
plot(Tmax, RNSVerify, '--x') % RNSVerify time
hold on; 
plot(Tmax, approx_vt, '-s'); % our approximate verification time
hold on;
plot(Tmax, RnnVerify, '-o'); % RNNVerify Time
legend('exact-star', 'RNSVerify', 'approx-star', 'RnnVerify');
xlabel('$T_{max}$', 'FontSize', 13, 'interpreter', 'latex');
ylabel('Verification Time (s)', 'FontSize', 13);
saveas(fig, 'verificationTime.fig');
