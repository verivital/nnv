%% this file generates fig. 3 captioned
% "Percentage Robustness and Runtime plots w.r.t increasing noise"
figure;

subplot(1,3,1);
plot(100*[0 perturbations],100*[1, PR_n1{3,1},PR_n1{3,2},PR_n1{3,3},PR_n1{3,4}]);
hold on
plot(100*[0 perturbations],100*[1, PR_n2{3,1},PR_n2{3,2},PR_n2{3,3},PR_n2{3,4}]);
hold on
plot(100*[0 perturbations],100*[1, PR_n3{3,1},PR_n3{3,2},PR_n3{3,3},PR_n3{3,4}]);
hold on
plot(100*[0 perturbations],100*[1, PR_n4{3,1},PR_n4{3,2},PR_n4{3,3},PR_n4{3,4}]);
xlabel('Pecentage Noise (%) ');
ylabel('Percentage Robustness by Samples (%)')
set(gca, 'FontSize', 22);


subplot(1,3,2);
plot(100*[0 perturbations],100*[1, POR_n1{3,1},POR_n1{3,2},POR_n1{3,3},POR_n1{3,4}]);
hold on
plot(100*[0 perturbations],100*[1, POR_n2{3,1},POR_n2{3,2},POR_n2{3,3},POR_n2{3,4}]);
hold on
plot(100*[0 perturbations],100*[1, POR_n3{3,1},POR_n3{3,2},POR_n3{3,3},POR_n3{3,4}]);
hold on
plot(100*[0 perturbations],100*[1, POR_n4{3,1},POR_n4{3,2},POR_n4{3,3},POR_n4{3,4}]);
xlabel('Pecentage Noise (%) ');
ylabel('Percentage Overlap Robustness (%)')
set(gca, 'FontSize', 22);

subplot(1,3,3);
plot(100*[0 perturbations],100*[0,T_sum_n1{3,1},T_sum_n1{3,2},T_sum_n1{3,3},T_sum_n1{3,4}]);
hold on
plot(100*[0 perturbations],100*[0,T_sum_n2{3,1},T_sum_n2{3,2},T_sum_n2{3,3},T_sum_n2{3,4}]);
hold on
plot(100*[0 perturbations],100*[0,T_sum_n3{3,1},T_sum_n3{3,2},T_sum_n3{3,3},T_sum_n3{3,4}]);
hold on
plot(100*[0 perturbations],100*[0,T_sum_n4{3,1},T_sum_n4{3,2},T_sum_n4{3,3},T_sum_n4{3,4}]);
xlabel('Pecentage Noise (%) ');
ylabel('Total Runtime (sec)')
set(gca, 'FontSize', 22);
legend('SFSI','SFAI','MFSI');
