
% /* An example of constructing a star input set */

c1 = [1; -1]; % center vector
v1 = [1; 0];  % basis vector 1
v2 = [0; 0.5]; % basis vector 2

V = [c1 v1 v2]; % basis matrix

% predicate constraint: P = C*[a] <= d
% -1<= a1 <= 1
% 0 <= a2 <= 1
% a1 + a2 <= 1
C = [1 0; -1 0; 0 1; 0 -1; 1 1]; % constraint matrix 
d = [1; 1; 1; 0; 1];             % constraint vector

I1 = Star(V, C, d); % star input set

% -2 <= x1 <= 2
% 0 <= x2 <= 1
lb = [-2; 0]; % lower bound vector
ub = [2; 1]; % upper bound vector

I2 = Star(lb, ub); % star input set
