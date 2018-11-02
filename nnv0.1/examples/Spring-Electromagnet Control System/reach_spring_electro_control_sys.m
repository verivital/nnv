load controller.mat;
load sysd.mat;

L1 = Layer(W{1,1}, b{1,1}, 'ReLU');
L2 = Layer(W{1,2}, b{1,2}, 'Linear');

NN_Controller = FFNN([L1 L2]); % feedforward neural network controller
Plant = DLinearODE(sysd.A, sysd.B, sysd.C, sysd.D, 0.1);
feedbackMap = [1;2]; % feedback map

ncs = NNCS(NN_Controller, Plant, feedbackMap); % the neural network control system

B = Box([0.8; -1],[1; -0.8]);
% initial set of state of the plant 
I = B.toPolyhedron();


N = 2; % number of step


Rx = Plant.stepReachPolyhedron(I, []); % the reachable set of state
Ry = Rx.affineMap(Plant.C);
Y = Ry;
X = Rx;


% construct input set for neural net controller in each step


%next_I = NNCS.nextInputSet(feedbackMap, I, Ry);

