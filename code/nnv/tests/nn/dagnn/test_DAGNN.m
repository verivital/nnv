%to run this as a test, use results_DAGNN=runtests('test_DAGNN')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables



%___________________________________________________________________________________________________
%tests below originally taken from test_DAGNN_parse.m

%% test 1: DAGNN
load SegNet.mat;
nnvNet = DAGNN.parse(net, 'SetNet');

