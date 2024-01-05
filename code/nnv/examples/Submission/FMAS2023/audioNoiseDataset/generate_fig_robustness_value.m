%% this file generates fig. 5 captioned
% "Percentage Robustness and Runtime plots w.r.t increasing noise"

load audioNoiseClassifier_results.mat
figure;
perturbations = percent;
subplot(1,2,1);
plot(100*[0 perturbations],100*[1, GPR_MFAI]);
hold on
plot(100*[0 perturbations],100*[1, GPR_MFSI]);
hold on
plot(100*[0 perturbations],100*[1, GPR_SFAI]);
hold on
plot(100*[0 perturbations],100*[1, GPR_SFSI]);
xlabel('Pecentage Noise (%) ');
ylabel('Percentage Robustness (%)')
set(gca, 'FontSize', 22);
legend('MFAI','MFSI','SFAI','SFSI');

subplot(1,2,2);
plot(100*[0 perturbations],100*[0,GTsum_MFAI]);
hold on
plot(100*[0 perturbations],100*[0,GTsum_MFSI]);
hold on
plot(100*[0 perturbations],100*[0,GTsum_SFAI]);
hold on
plot(100*[0 perturbations],100*[0,GTsum_SFSI]);
xlabel('Pecentage Noise (%) ');
ylabel('Total Runtime (sec)')
set(gca, 'FontSize', 22);
legend('MFAI','MFSI','SFAI','SFSI');
