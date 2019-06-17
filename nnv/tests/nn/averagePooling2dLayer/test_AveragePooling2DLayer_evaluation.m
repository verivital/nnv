% original input volume: color image with 3 channels
inputVol(:, :, 1) = [0 0 2 0 0; 1 2 0 2 0; 0 0 2 2 0; 0 2 2 2 2; 2 2 2 1 1]; % channel 1 input matrix
inputVol(:, :, 2) = [1 2 2 1 2; 2 1 2 0 2; 2 2 2 0 1; 1 1 1 0 0; 1 0 2 2 1]; % channel 2 input matrix
inputVol(:, :, 3) = [0 0 2 2 1; 0 2 1 1 2; 0 2 0 0 1; 0 2 1 0 1; 1 2 1 0 0]; % channel 3 input matrix

L = AveragePooling2DLayer([3 3], [2 2], [0 0 0 0]);

y = L.evaluate(inputVol);

for i=1:3
    display(y(:,:,i));
    display(inputVol(:,:,i));
end





