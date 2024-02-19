%% Let's take a look into this

load('results/verification_multipleAttacks_nodulemnist3d_debug.mat');

% Only get verification results
results = results(1,:,:,:);
results = squeeze(results);
results = results';
outputSets = outputSets';

% first discrepancy (5)
R1a = outputSets{1};
R1 = R1a{5}{end};
[lb1, ub1] = R1.getRanges;
R2a = outputSets{2};
R2 = R2a{5}{end};
[lb2, ub2] = R2.getRanges;
R3a = outputSets{3};
R3 = R3a{5}{end};
[lb3, ub3] = R3.getRanges;

% The ranges are different, but the 3 bounds should be containing the other
% 2... We are probably doing somethig wrong, either when defining the input
% set, or something with the reachability (most likely the first...)

% Let's visualize the output sets (exact)
figure;
Star.plot(R1,'r');
figure;
Star.plot(R2,'r');
figure;
Star.plot(R3,'r');

