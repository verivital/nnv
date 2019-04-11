
fprintf('\nLoading data...');
load counterInputSet.mat;
load ACASXU_run2a_2_8_batch_2000.mat;
fprintf('\nLoading data done!, start ploting...');

% normalized counter input set
% input = [x1 = \rho; x2 = \Theta; x3 = \psi; x4 = \v_own; x5 = \v_in]
maps1 = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0]; % plot a projection on \rho, \theta, \psi
maps2 = [1 0 0 0 0; 0 0 0 1 0; 0 0 0 0 1]; % plot a projection on \rho, \v_own, \v_in
R = counterInputSet;
n = length(R);
R1 = [];
R2 = [];
parfor i=1:n
    fprintf('\nProjecting %d^th counter input set on three-dimensonal subspace', i)
    X1 = R(i).affineMap(maps1, []);
    X2 = R(i).affineMap(maps2, []); 
    R1 = [R1 X1]; % normalized counter input set on \rho, \theta, \psi
    R2 = [R2 X2]; % normalized counter input set on \rho, \v_own, \v_in
end

P1 = [];
P2 = [];
parfor i=1:n
    fprintf('\nConverting %d^th star set to polyhedron', i);
    X1 = R1(i).toPolyhedron;
    X2 = R2(i).toPolyhedron;
    P1 = [P1 X1];
    P2 = [P2 X2];
end

% normalized input set
lb = [1500; -0.06; 3.1; 1000; 700];
ub = [1800; 0.06; 3.14; 1200; 800];
for i=1:5
    lb(i) = (lb(i) - means_for_scaling(i))/range_for_scaling(i);
    ub(i) = (ub(i) - means_for_scaling(i))/range_for_scaling(i);   
end

I = Box(lb, ub);
I = I.toStar;

fprintf('\nPlotting the counter input set...');
% plot reachable set
fig = figure;
subplot(1, 2, 1);
Star.plotBoxes_3D(I, 1, 2, 3, 'red'); % plot input set on \rho, \theta and \psi
hold on;
P1.plot; % plot counter input set on \rho, \theta, and \psi
xlabel('\rho', 'Fontsize', 16);
ylabel('\theta', 'Fontsize', 16);
zlabel('\psi', 'Fontsize', 16);
set(gca, 'Fontsize', 16);

subplot(1, 2, 2);
Star.plotBoxes_3D(I, 1, 4, 5, 'red'); % plot input set on \rho, \v_own and \v_in
hold on;
P2.plot; % plot counter input set on \rho, \v_own, and \v_in
xlabel('\rho', 'Fontsize', 16);
ylabel('v_{own}', 'Fontsize', 16);
zlabel('v_{int}', 'Fontsize', 16);
set(gca, 'Fontsize', 16);

set(gca, 'Fontsize', 16);
saveas(gcf, 'counterInputSet_on_N2_8.pdf');
