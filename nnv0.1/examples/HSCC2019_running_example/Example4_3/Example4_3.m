
load NeuralNetwork7_3.mat;
R = reachSet;
N1 = 16;
N2 = 128;

P1 = Reduction.merge_box(R, N1, 'single');
P2 = Reduction.merge_box(R, N2, 'single');

fig = figure;
subplot(1, 3, 1);
R.plot;
title('Exact reachable set, N = 1250', 'FontSize', 20);
xlabel('x_1', 'FontSize', 16);
ylabel('x_2', 'FontSize', 16);

subplot(1,3, 2);
Box.plots(P1);
title('Merged reachable set, N = 16', 'FontSize', 20);
xlabel('x_1', 'FontSize', 16);
ylabel('x_2', 'FontSize', 16);

subplot(1,3, 3);
Box.plots(P2);
title('Merged reachable set, N = 128', 'FontSize', 20);
xlabel('x_1', 'FontSize', 16);
ylabel('x_2', 'FontSize', 16);
