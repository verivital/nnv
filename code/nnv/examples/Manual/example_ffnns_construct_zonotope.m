% /* An example of constructing a zonotope input set */

c = [-1; 0]; % center vector
v1 = [2; 1]; % 1st basis vector
v2 = [1; 1]; % 2nd basis vector
v3 = [-1; 1]; % 3rd basis vector
V = [v1 v2 v3]; % basis matrix

% zonotope from center vector and basis matrix
I1 = Zono(c, V);  

lb = [-1; 1]; % lower bound vector
ub = [1; 2];  % upper bound vector

% zonotope from input ranges
I2 = Box(lb, ub); % a box object
I2 = I2.toZono;   % transform to zonotope