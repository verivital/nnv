% Plot reachable sets for Discrete Linear ACC model
% Dung Tran: 10/22/2019



%% Plot output reach sets: actual distance vs. safe distance

% plot reachable set of the distance between two cars d = x1 - x4 vs. ego
% car's velocity

figure;
h1 = subplot(2,1,1);
load dis_ACC_ncs_1.mat;
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
title(h1,'Actual Distance (blue) vs. Safe Distance (red)');
xlabel(h1, 'Control Time Steps');
ylabel(h1, 'Distance');
xticks(h1, [0:5:50])

h2 = subplot(2,1,2);
load dis_ACC_ncs_6.mat;
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
title(h2,'Actual Distance (blue) vs. Safe Distance (red)');
xlabel(h2, 'Control Time Steps');
ylabel(h2, 'Distance');
xticks(h2, [0:5:50])


%% END OF SCRIPT