load Engine_Toy_Tansig_net.mat;

nnvNet = FFNNS.parse(net); % parse the network trained from matlab
lb = [-1; -1];
ub = [1; 1];
I = Star(lb, ub);

fprintf('\nReachable set of the toy network:\n');
reachSet = nnvNet.reach(I, 'approx-star') % perform reachability analyis using approx-star

fprintf('\nRange of the output of the toy network:\n');
B = reachSet.getBox % get the output range 