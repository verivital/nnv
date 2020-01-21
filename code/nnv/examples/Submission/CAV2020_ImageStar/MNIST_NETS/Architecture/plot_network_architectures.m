
load Medium_ConvNet.mat;
load Large_ConvNet.mat;

figure;

subplot(1,3,1);
load Small_ConvNet.mat;
layers = net.Layers;
lgraph = layerGraph(layers);
plot(lgraph);
title('MNIST\_Small', 'FontSize', 11);
xlabel('98% accuracy');
set(gca, 'FontSize', 11);


subplot(1,3,2);
load Medium_ConvNet.mat;
layers = net.Layers;
lgraph = layerGraph(layers);
plot(lgraph);
title('MNIST\_Medium', 'FontSize', 11);
xlabel('99.7% accuracy');
set(gca, 'FontSize', 11);

subplot(1,3,3);
load Large_ConvNet.mat;
layers = net.Layers;
lgraph = layerGraph(layers);
plot(lgraph);
title('MNIST\_Large', 'FontSize', 11);
xlabel('99.9% accuracy');
set(gca, 'FontSize', 11);
