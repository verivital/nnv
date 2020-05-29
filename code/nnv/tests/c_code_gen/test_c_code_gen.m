%to run this as a test, use results_c_code_gen=runtests('test_c_code_gen')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables





%___________________________________________________________________________________________________
%tests below originally taken from test.m

%% test 1: basic vertice test
lb = [1; 1];
ub = [2; 1];
test_box_getVertices_c_gen([1;1], [1.5; 1.5]);


%___________________________________________________________________________________________________
%tests below originally taken from test_box_getVertices.m

%% test 2: box get vertices



lb = [-0.1; -0.2; -1];
ub = [1; 0.5; 1];
B = Box(lb, ub);

V = B.getVertices();











