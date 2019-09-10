% construct a FullyConnectedLayer object
rl = ReluLayer();
% image input set
IM(:,:,1) = [-1 1 0 -1; 0 0 1 -1; 1 0 -1 0; 1 -1 -1 1]; % center image channel 1
IM(:,:,2) = [0 1 0 0; 1 0 0 -1; 0 1 -1 0; 0 1 0 -1]; % center image channel 2
IM(:,:,3) = [1 -1 1 1; 1 -1 0 1; 0 1 -1 0; 1 0 -1 0]; % center image channel 3

output = rl.evaluate(IM);
