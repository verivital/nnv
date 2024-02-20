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


% Let's see what the differences in the layer-by-layer sets are...

% Layer 1
r1 = R1a{5}{1};
r2 = R2a{5}{1};
r3 = R3a{5}{1};

rs1 = r1.toStar;
rs2 = r2.toStar;
rs3 = r3.toStar;

% The next 3 statements should return true
% bool1 = rs1.isSubSet(rs2);
% bool2 = rs1.isSubSet(rs3);
% bool3 = rs2.isSubSet(rs3);
% this takes too long as stars are high dimensional (polyhedron not very
% efficient)

% get ranges (overapprox) and see if the smaller bounds are within te
% larger ones (rs1 < rs2 < rs3)
[rs1L, rs1U] = rs1.estimateRanges;
[rs2L, rs2U] = rs2.estimateRanges;
[rs3L, rs3U] = rs3.estimateRanges;

b1L = all(rs1L >= rs2L);
b1U = all(rs1U <= rs2U);
b2L = all(rs1L >= rs3L);
b2U = all(rs1U <= rs3U);
b3L = all(rs2L >= rs3L);
b3U = all(rs2U <= rs3U);

% Lower bounds seem fine, but ot the upper bounds...

% where do these assertions fail?
u1 = find(rs1U > rs2U);
u2 = find(rs1U > rs3U);
u3 = find(rs2U > rs3U);

% these fail on all the perturbed pixels...
uu11 = rs1L(u1);
uu12 = rs2L(u1);

% we are not defining the constraints on the input set correctly, some
% pixel values are outside the valid ranges...
