% Load network with computed reach sets
nn = load("SthToTest.mat");
nn = nn.nn;

% Get all the reach sets
set1 = nn.reachSet{1};
set2 = nn.reachSet{2};
set3 = nn.reachSet{3};
set4 = nn.reachSet{4};
set5 = nn.reachSet{5};
set6 = nn.reachSet{6};
set7 = nn.reachSet{7};
set8 = nn.reachSet{8};
set9 = nn.reachSet{9};

% DepthConcatenation concatenates two inputs along the channel dimension
% (3rd dimension in ImageStars.V)
% new_V = 

concatSet = ImageStar(new_V, set9.C, set9.d, set9.pred_lb, set9.pred_ub);  