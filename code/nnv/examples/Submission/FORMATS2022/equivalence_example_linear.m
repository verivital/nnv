% Show equivalence transformation with one example
% Weight matrices
w1 = [1 3;2 4]; 
w2 = [5 6;7 8]; 
w3 = [1 10;1 0];
% Bias matrices
b1 = [1;10];
b2 = [5;2];
b3 = [3;1];

% Initial state
x0 = [1;2]; % initial state

% State-space representation
A = w3*w2*w1; 
B = 0;
c = b3+w3*b2+w3*w2*b1;
xF = A*x0+c; 

% Neural network (concatenate functions)
x1 = w1*x0+b1; % layer 1
x2 = w2*x1+b2; % layer 2
x3 = w3*x2+b3; % layer 3

all(x3 == xF)