figure
subplot(1,3,1);
plot(audionoise.layerGraph);
title("Audio Noise LSTM Classifier")
set(gca, 'FontSize', 22);

subplot(1,3,2);
plot(audionoise.layerGraph);
title("Japanese Vowel LSTM Classifier")
set(gca, 'FontSize', 22);

subplot(1,3,3);
plot(cnnlstm.layerGraph);
title("Japanese Vowel CNN+LSTM Classifier")
set(gca, 'FontSize', 22);