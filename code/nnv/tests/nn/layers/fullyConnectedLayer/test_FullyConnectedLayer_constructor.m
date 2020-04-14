% construct FullyConnectedLayer objects
W = rand(4, 48);
b = rand(4,1);
fc = FullyConnectedLayer(W,b);
fc1 = FullyConnectedLayer('fc1', W, b);