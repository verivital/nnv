% /* An example of visualizing verification results */

%/* An example of visualizing the verification results of an FFNN */
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
% P = HalfSpace([-1 0], -1.5); % P: y1 >= 1.5
P = HalfSpace([-1 0], -0.4); % P: y1 >= 0.4
%/* verify the network
nC = 1; % number of cores
nS = 0; % number of samples

map_mat = eye(2); % mapping matrix
map_vec = []; % mapping vector
P_poly = Polyhedron('A', P.G, 'b', P.g); % polyhedron obj

fig1 = figure;
subplot(2, 2, 1);
[safe1, t1, cE1] = F.verify(I, P, 'exact-star', nC, nS);
F.visualize(map_mat, map_vec); % plot y1 y2
hold on;
plot(P_poly); % plot unsafe region
title('exact-star', 'FontSize', 13);

subplot(2,2,2);
[safe2, t2, cE2] = F.verify(I, P, 'approx-star', nC, nS);
F.visualize(map_mat, map_vec); % plot y1 y2
hold on;
plot(P_poly); % plot unsafe region
title('approx-star', 'FontSize', 13);

subplot(2,2,3);
[safe3, t3, cE3] = F.verify(I_Zono, P, 'approx-zono', nC, nS);
F.visualize(map_mat, map_vec); % plot y1 y2
hold on;
plot(P_poly); % plot unsafe region
title('approx-zono', 'FontSize', 13);

subplot(2, 2, 4);
[safe4, t4, cE4] = F.verify(I, P, 'abs-dom', nC, nS);
F.visualize(map_mat, map_vec); % plot y1 y2
hold on;
plot(P_poly); % plot unsafe region
title('abs-dom', 'FontSize', 13);
