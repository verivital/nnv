% /* An example of constructing a polyhedron input set */

lb = [-1; 1]; % lower bound vector
ub = [1; 2];  % upper bound vector

% polyhedron from input ranges
I1 = Polyhedron('lb', lb, 'ub', ub); 

A = [2 1; 1 0; -1 0; 0 1; 0 -1; 1 1]; % inequality matrix
b = [2; 1; 1; 0; 1; 1]; % inequality vector

I2 = Polyhedron('A', A, 'b', b);% polyhedron without equalities

Ae = [2 3]; % equality matrix
be = 1.5;   % equality vector

% polyhedron with one equality
I3 = Polyhedron('A', A, 'b', b, 'Ae', Ae, 'be', be); 

