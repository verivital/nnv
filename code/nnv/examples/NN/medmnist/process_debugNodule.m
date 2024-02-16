%% Let's take a look into this

load('results/verification_multipleAttacks_nodulemnist3d_debug.mat');

% Only get verification results
results = results(1,:,:,:);
results = squeeze(results);
results = results';
outputSets = outputSets';

% first discrepancy (5)
R1 = outputSets(1,5); % 
[lb1, ub1] = R1.getRanges;
R2 = outputSets(2,5);
[lb2, ub2] = R1.getRanges;
R3 = outputSets(3,5);
[lb3, ub3] = R1.getRanges;

% So the ranges are the same, same results should be getting...
% What is happening?

% Nodulemnist is a binary classification problem
% The bounds are clearly apart, so no unkown should be happening...
% Are we representing the robustness halfspace incorrectly?

% Verification property
G = [1,-1]; g = 0;
Hs = HalfSpace(G,g);

res1 = verify_specification(R1, Hs);
res2 = verify_specification(R2, Hs);
res3 = verify_specification(R3, Hs);

% bounds are the same, but the verification results are not...
% why is there no intersections? Does the intersection method have a bug on
% it?

% Let's visualize the output sets (exact)
% figure;
% Star.plot(R1,'r');
% figure;
% Star.plot(R2,'r');
% figure;
% Star.plot(R3,'r');

% Are the predicate bounds incorrect? They show 254 and 255... 
% This seems just wrong... Can we simply remove the predicate bounds in 
% all cases? I don't think so...
% Then we need to fix either the how the predicate bounds get computed, or
% do not use them for certain operations such as plotting or computing
% intersections/empty sets.

% Remove predicate bounds
R1.predicate_lb = [];
R1.predicate_ub = [];
R2.predicate_lb = [];
R2.predicate_ub = [];
R3.predicate_lb = [];
R3.predicate_ub = [];

% Check verification results now
res11 = verify_specification(R1, Hs);
res22 = verify_specification(R2, Hs);
res33 = verify_specification(R3, Hs);

% We get the same results, so I guess the predicate bounds do not affect
% this verification result (intersection of Star sets)

% The difference between R3 and the other 2 is that there is an extra
% constraint. Is this messing things up?
% Let's manually visualize the stars...

% get centers
c1 = R1.V(:,1);
c2 = R2.V(:,1);
c3 = R3.V(:,1);

% get basis vectors
V1 = R1.V(:,2:end);
V2 = R2.V(:,2:end);
V3 = R3.V(:,2:end);

% get constraints
C1 = R1.C;
d1 = R1.d;
C2 = R2.C;
d2 = R2.d;
C3 = R3.C;
d3 = R3.d;

% can we get the bounds of the stars from here?
R4 = Star(lb1,ub1);
Star.plot(R4, 'r');
