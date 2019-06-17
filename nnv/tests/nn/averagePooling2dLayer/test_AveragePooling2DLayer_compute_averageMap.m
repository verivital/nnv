
% this example is from (see max pooling layer)
%https://en.wikipedia.org/wiki/Convolutional_neural_network

I = [1 0 2 3; 4 6 6 8; 3 1 1 0; 1 2 2 4]; % input
L = AveragePooling2DLayer([2 2], [2 2], [0 0 0 0]);
averageMap = L.compute_averageMap(I);
display(averageMap);
display(I);


