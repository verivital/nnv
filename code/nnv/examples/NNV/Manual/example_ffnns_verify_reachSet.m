%/* An example of manually verifying an FFNN */
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

% Properties   
P = HalfSpace([-1 0], -1.5); % P: y1 >= 1.5
P_Poly = Polyhedron('A', P.G, 'b', P.g); 

rs1 = cell(1, length(R1));
rs2 = cell(1, length(R2));

% verifying R1
for i=1:length(R1)
    M = R1(i).intersectHalfSpace(P.G, P.g);
    if isempty(M)
        rs1{i} = 'UNSAT';
    else
        rs1{i} = 'SAT';
    end
end

% verifying R2
for i=1:length(R2)
    M = R2(i).intersect(P_Poly);
    if M.isEmptySet
        rs2{i} = 'UNSAT';
    else
        rs2{i} = 'SAT';
    end
end

% verify R3
M = R3.intersectHalfSpace(P.G, P.g);
if isempty(M)
    rs3 = 'UNSAT';
else
    rs3 = 'SAT';
end

% verify R4
M = R4.intersectHalfSpace(P.G, P.g);
if isempty(M)
    rs4 = 'UNSAT';
else
    rs4 = 'SAT';
end

% verify R5
M = R5.intersectHalfSpace(P.G, P.g);
if isempty(M)
    rs5 = 'UNSAT';
else
    rs5 = 'SAT';
end
