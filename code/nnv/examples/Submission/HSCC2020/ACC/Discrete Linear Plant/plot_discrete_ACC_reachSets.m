% Plot reachable sets for Discrete Linear ACC model
% Dung Tran: 10/22/2019



%% Load verification results

load discrete_ACC_verification_results.mat;


%% Plot output reach sets: actual distance vs. safe distance

% plot reachable set of the distance between two cars d = x1 - x4 vs. ego
% car's velocity
ncs = dis_ACC_approx{1};
figure; 
map_mat = [1 0 0 -1 0 0 0];
map_vec = [];
ncs.plotOutputReachSets('blue', map_mat, map_vec);

hold on;

% plot safe distance between two cars: d_safe = D_default + t_gap * v_ego;
% D_default = 10; t_gap = 1.4 
% d_safe = 10 + 1.4 * x5; 

map_mat = [0 0 0 0 1.4 0 0];
map_vec = [10];

ncs.plotOutputReachSets('red', map_mat, map_vec);
title('Actual Distance versus. Safe Distance');

%% END OF SCRIPT