% Test Sign functionality
% To run: results = runtests('test_Sign')

%% Test 1: Sign evaluate with polar_zero_to_pos_one mode
x = [-2; -1; 0; 1; 2];
mode = "polar_zero_to_pos_one";
y = Sign.evaluate(x, mode);

expected = [-1; -1; 1; 1; 1];  % zero maps to 1
assert(isequal(y, expected), 'Sign evaluate with polar_zero_to_pos_one failed');

%% Test 2: Sign evaluate with nonnegative_zero_to_pos_one mode
x = [-2; -1; 0; 1; 2];
mode = "nonnegative_zero_to_pos_one";
y = Sign.evaluate(x, mode);

expected = [0; 0; 1; 1; 1];  % negative maps to 0, zero maps to 1
assert(isequal(y, expected), 'Sign evaluate with nonnegative_zero_to_pos_one failed');

%% Test 3: Sign get_sign_lb
sign_lb = Sign.get_sign_lb('polar_zero_to_pos_one');
assert(sign_lb == -1, 'get_sign_lb for polar_zero_to_pos_one should be -1');

sign_lb = Sign.get_sign_lb('nonnegative_zero_to_pos_one');
assert(sign_lb == 0, 'get_sign_lb for nonnegative_zero_to_pos_one should be 0');

%% Test 4: Sign get_sign_ub
sign_ub = Sign.get_sign_ub('polar_zero_to_pos_one');
assert(sign_ub == 1, 'get_sign_ub for polar_zero_to_pos_one should be 1');

sign_ub = Sign.get_sign_ub('nonnegative_zero_to_pos_one');
assert(sign_ub == 1, 'get_sign_ub for nonnegative_zero_to_pos_one should be 1');

%% Test 5: Sign stepReach with constant positive value
lb = [1];
ub = [1];
I = Star(lb, ub);

index = 1;
mode = 'polar_zero_to_pos_one';
S = Sign.stepReach(I, index, 'linprog', mode);

% Should return star with sign value of 1
assert(isa(S, 'Star'), 'stepReach should return Star');
assert(S.V(index, 1) == 1, 'stepReach with constant +1 should set output to 1');

%% Test 6: Sign stepReach with constant negative value  
lb = [-1];
ub = [-1];
I = Star(lb, ub);

index = 1;
mode = 'polar_zero_to_pos_one';
S = Sign.stepReach(I, index, 'linprog', mode);

% Should return star with sign value of -1
assert(isa(S, 'Star'), 'stepReach should return Star');
assert(S.V(index, 1) == -1, 'stepReach with constant -1 should set output to -1');

%% Test 7: Sign stepReach with positive range
V = [0.5 1; 0.5 1];
C = [1; -1];
d = [1; 1];
lb = [0];
ub = [2];
I = Star(V, C, d, lb, ub);

index = 1;
mode = 'polar_zero_to_pos_one';
S = Sign.stepReach(I, index, 'linprog', mode);

% Should return star with sign value of 1
assert(isa(S, 'Star'), 'stepReach should return Star');
assert(S.V(index, 1) == 1, 'stepReach with positive range should set output to 1');

%% Test 8: Sign reach with Star - exact method
lb = [-1; -1];
ub = [1; 1];
% B = Box(lb, ub);
I_star = Star(lb,ub);

% mode = 'polar_zero_to_pos_one';
S = Sign.reach(I_star, 'exact-star');
assert(~isempty(S), 'Sign reach exact-star should return result');

%% Test 9: Sign reach with Star - approx method
lb = [-1; -1];
ub = [1; 1];
%B = Box(lb, ub);
I_star = Star(lb,ub);

%mode = 'polar_zero_to_pos_one';
S = Sign.reach(I_star, 'approx-star');
assert(~isempty(S), 'Sign reach approx-star should return result');
