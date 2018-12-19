Center = [1 0 0; 1 1 0; 1 0 1]; % center matrix
Basis = [1 0 0; 0 0 0; 0 0 0]; % basic matrix

V = cell(1, 2);
V{1} = Center;
V{2} = Basis;

% constraint: -1<= a <= 1
Constr_mat = [1; -1];
Constr_vec = [1; 1]; 

S = Star2D(V, Constr_mat, Constr_vec);
