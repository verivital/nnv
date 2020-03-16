% /* Example of manually create nnv FFNN object */
W1 = [1 -1; 0.5 2; -1 1]; % first layer weight matrix
b1 = [-1; 0.5; 0];        % first layer bias vector

W2 = [-2 1 1; 0.5 1 1];   % second layer weight matrix
b2 = [-0.5; -0.5];        % second layer bias vector

L1 = LayerS(W1, b1, 'poslin'); % construct first layer object
L2 = LayerS(W2, b2, 'purelin');   % construct second layer object

F = FFNNS([L1 L2]); % construct Feedforward neural network
