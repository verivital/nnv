%/* An example of automatically verifying an FFNN */
%/* construct an NNV network
W1 = [1 -1; 0.5 2; -1 1]; 
b1 = [-1; 0.5; 0];        
W2 = [-2 1 1; 0.5 1 1];   
b2 = [-0.5; -0.5];        
L1 = LayerS(W1, b1, 'poslin'); 
L2 = LayerS(W2, b2, 'purelin');   
F = FFNNS([L1 L2]); % construct an NNV FFNN
%/* construct input set
lb = [-1; -2]; % lower bound vector
ub = [1; 0]; % upper bound vector
I = Star(lb, ub); % star input set
B = Box(lb, ub); % a box input set
I_Zono = B.toZono; % convert to a zonotope
%/* Properties   
P = HalfSpace([-1 0], -1.5); % P: y1 >= 1.5
% P = HalfSpace([-1 0], -0.4); % P: y1 >= 0.4
%/* verify the network
nC = 1; % number of cores
nS = 100; % number of samples
[safe1, t1, cE1] = F.verify(I, P, 'exact-star', nC, nS);
[safe2, t2, cE2] = F.verify(I, P, 'approx-star', nC, nS);
[safe3, t3, cE3] = F.verify(I_Zono, P, 'approx-zono', nC, nS);
[safe4, t4, cE4] = F.verify(I, P, 'abs-dom', nC, nS);

