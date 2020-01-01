load 'verify_nonlinear_ACC.mat'

cI = counterExamples{1};
cI = cell2mat(cI);
d_rel = [1 0 0 -1 0 0]*cI;
d_safe = [0 0 0 1.4 0 0]*cI + 10;

figure; 
T = 0:1:50;
plot(T, d_rel, 'blue');
hold on;
plot(T, d_safe, 'red');

xlabel('Control Time Steps', 'FontSize', 13);
ylabel('Distance', 'FontSize', 13);
xticks([0:5:50]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',13)
title('Actual Distance (blue) vs. Safe Distance (red)');
