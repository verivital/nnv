%/* An example of computing reachable sets of an FFNN */
%/* construct an NNV network
W1 = [1 -1; 0.5 2; -1 1]; 
b1 = [-1; 0.5; 0];        
W2 = [-2 1 1; 0.5 1 1];   
b2 = [-0.5; -0.5];        
L1 = LayerS(W1, b1, 'poslin'); 
L2 = LayerS(W2, b2, 'purelin');   
F = FFNNS([L1 L2]); % construct an NNV FFNN

%/* choose the number of cores
numCores = 2; 

%/* construct input set
lb = [-1; -2]; % lower bound vector
ub = [1; 0]; % upper bound vector
I = Star(lb, ub); % star input set
I_Poly = Polyhedron('lb', lb, 'ub', ub); % polyhedron input set
B = Box(lb, ub); % a box input set
I_Zono = B.toZono; % convert to a zonotope

%/* compute the reachable sets with a selected method
[R1, t1] = F.reach(I, 'exact-star', numCores); 
[R2, t2] = F.reach(I_Poly, 'exact-polyhedron', numCores); 
[R3, t3] = F.reach(I, 'approx-star');
[R4, t4] = F.reach(I_Zono, 'approx-zono');
[R5, t5] = F.reach(I, 'abs-dom');
