load 2L_pollution_nets.mat;

W1 = nnetwork.W{1, 1};
b1 = nnetwork.b{1, 1};

W2 = nnetwork.W{1, 2};
b2 = nnetwork.b{1, 2};

W3 = nnetwork.W{1, 3};
b3 = nnetwork.b{1, 3};

L1 = Layer(W1, b1, 'ReLU');
L2 = Layer(W2, b2, 'ReLU');
L3 = Layer(W3, b3, 'ReLU');


F = FFNN([L1, L2, L3]);


lb = nnetwork.min';
ub = nnetwork.max';

% Input set
% x[i] = 0, i = 1:22
% lb(i) <= x[i] <= ub(i), i=23, 24

n = length(lb);

Ae = eye(n);
Ae(n-1:n, :) = [];
be = zeros(n-2, 1);

A1 = zeros(2, n);
A1(1, n) = 1;
A1(2, n) = -1;

A2 = zeros(2, n);
A2(1, n-1) = 1;
A2(2, n-1) = -1;

A = vertcat(A1, A2);
b = vertcat(ub(n), -lb(n), ub(n-1), -lb(n-1));

% Input Set
I = Polyhedron('A', A, 'b', b, 'Ae', Ae, 'be', be);

% exact range analysis
[R1, t1] = F.reach(I, 'exact', 4, []); % exact scheme
R11 = Reduction.hypercubeHull(R1);
range1 = [R11.lb, R11.ub];


% lazy-approximate + input partition method for range analysis

n = 4;
CSV3 = cell(1,n);
N_partitions = [];
T3 = [];
for k=1:n
    I1 = Partition.partition_box(I,k);
    N_partitions = [N_partitions length(I1)];
    [R3, t3] = F.reach(I1, 'approx', 4, []); % lazy-approximate scheme
    R31 = Reduction.hypercubeHull(R3);
    range3 = [R31.lb R31.ub];
    CSV = zeros(3,1);
    for i=1:3
        CSV(i, 1) = max(abs(range3(i,1) - range1(i,1)), abs(range3(i,2) - range1(i,2))) / (range1(i,2) - range1(i,1));
    end
    CSV3{1, k} = 100*CSV;
    T3 = [T3 t3];
end


y1 = [];
y2 = [];
y3 = [];
for i=1:n
    y1 = [y1 CSV3{1, i}(1)];
    y2 = [y2 CSV3{1, i}(2)];
    y3 = [y3 CSV3{1, i}(3)];
end

x = N_partitions;



% mixing scheme for output range analysis
N_max = [4 6 8 10];
CSV4 = cell(1,5);
T4 = [];
for k=1:4
    [R4, t4] = F.reach(I, 'mix', 4, N_max(k));
    R41 = Reduction.hypercubeHull(R4);
    range4 = [R41.lb, R41.ub];
    CSV = zeros(3,1);
    for i=1:3
        CSV(i, 1) = max(abs(range4(i,1) - range1(i,1)), abs(range4(i,2) - range1(i,2))) / (range1(i,2) - range1(i,1));
    end
    CSV4{1, k} = 100*CSV;
    T4 = [T4 t4];    
end


y11 = [];
y21 = [];
y31 = [];

for i=1:4
    y11 = [y11 CSV4{1, i}(1)];
    y21 = [y21 CSV4{1, i}(2)];
    y31 = [y31 CSV4{1, i}(3)];
end

x1 = N_max;



fig = figure;
subplot(2,1,1);
plot(x, y1, '-*');
hold on;
plot(x, y2, '-x');
hold on;
plot(x, y3, '-o');
title({'Lazy-approximate + input partition scheme'; '(Conservativeness vs. number of partitions)'}, 'FontSize', 30);
xlabel('Number of partitions', 'FontSize', 25);
ylabel('Conservativeness (%)', 'FontSize', 25);
legend('y_1', 'y_2', 'y_3');
set(gca,'FontSize',25);

subplot(2,1,2);
plot(x1, y11, '-*');
hold on;
plot(x1, y21, '-x');
hold on;
plot(x1, y31, '-o');
title({'Mixing scheme'; '(Conservativeness vs. maximum, allowable number of polyhedra)'}, 'FontSize', 30);
xlabel('Maximum, allowable number of polyhedra', 'FontSize', 25);
ylabel('Conservativeness (%)', 'FontSize', 25);
legend('y_1', 'y_2', 'y_3');
set(gca,'FontSize',25);
