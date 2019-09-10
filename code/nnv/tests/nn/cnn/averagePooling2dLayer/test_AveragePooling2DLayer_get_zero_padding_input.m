
% original input volume: color image with 3 channels
inputVol(:, :, 1) = [2 0 1 2 1; 1 0 2 2 2; 1 2 2 0 2; 1 2 0 0 1; 1 0 1 1 2]; % channel 1 input matrix
inputVol(:, :, 2) = [0 0 1 0 1; 0 0 2 1 1; 1 1 0 1 1; 1 1 0 2 2; 2 1 2 0 0]; % channel 2 input matrix
inputVol(:, :, 3) = [1 2 2 1 0; 2 0 0 2 0; 0 0 1 0 1; 1 2 0 2 0; 1 0 2 1 0]; % channel 3 input matrix


% construct input with padding operation
paddingSize = [1 1 1 1];

L = AveragePooling2DLayer();
L.set_padding(paddingSize);

I = L.get_zero_padding_input(inputVol);

