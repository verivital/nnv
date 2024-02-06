%% this file generates fig. 5 captioned
% "Percentage Robustness and Runtime plots w.r.t increasing noise"

load japanesevowel_cnnlstmClassifier_results2.mat
figure;
perturbations = percent;
subplot(1,2,1);
plot(100*[0 perturbations],100*[1, PR_MFAI]);
hold on
plot(100*[0 perturbations],100*[1, PR_MFSI]);
hold on
plot(100*[0 perturbations],100*[1, PR_SFAI]);
hold on
plot(100*[0 perturbations],100*[1, PR_SFSI]);
xlabel('Pecentage Noise (%) ');
ylabel('Percentage Robustness (%)')
set(gca, 'FontSize', 22);
legend('MFAI','MFSI','SFAI','SFSI');

subplot(1,2,2);
plot(100*[0 perturbations],100*[0,T_sum_MFAI]);
hold on
plot(100*[0 perturbations],100*[0,T_sum_MFSI]);
hold on
plot(100*[0 perturbations],100*[0,T_sum_SFAI]);
hold on
plot(100*[0 perturbations],100*[0,T_sum_SFSI]);
xlabel('Pecentage Noise (%) ');
ylabel('Total Runtime (sec)')
set(gca, 'FontSize', 22);
legend('MFAI','MFSI','SFAI','SFSI');
