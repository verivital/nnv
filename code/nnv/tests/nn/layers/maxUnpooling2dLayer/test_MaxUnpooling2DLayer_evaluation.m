% original input volume: color image with 3 channels
inputVol(:, :, 1) = [0 0 2 0; 1 2 0 2 ; 0 0 2 2; 0 2 2 2]; % channel 1 input matrix
inputVol(:, :, 2) = [1 2 2 1; 2 1 2 0; 2 2 2 0; 1 1 1 0]; % channel 2 input matrix
inputVol(:, :, 3) = [0 0 2 2; 0 2 1 1; 0 2 0 0; 0 2 1 0]; % channel 3 input matrix

L = MaxPooling2DLayer([2 2], [2 2], [0 0 0 0]);

y = L.evaluate(inputVol);

L1 = MaxUnpooling2DLayer();


y1 = L1.evaluate(y, L.MaxIndx, L.InputSize);




