%% Example of how to verify a property using Multiple Star sets and HalfSpaces

% The property is the following
% y2 > y1
% However, we represent the properties as counterexamples, so the Halfspace
% will be represented as 
% y2 - y1 <= 0
% G = [-1 1]
% g = 0
G = [-1 1];
g = 0;
Hs = HalfSpace(G,g);

% Set1 (Star)
XLower = [0; 10];
XUpper = [2; 11];
Set1 = Star(XLower, XUpper);

% Set1 (Star)
XLower = [-10; 1];
XUpper = [0 ; 10];
Set2 = Star(XLower, XUpper);

% Sets to verify
Sets = [Set1 Set2];

% Verification
res = verify_specification(Sets, Hs); % res = 1, which means property is UNSAT

% How we are verifying properties returns the property is UNSAT...
% But let's visualize it... 
figure;
Star.plots(Sets);
hold on;
% the point [1.5; 1.2] is within the boundaries of Sets, which will violate the property...
plot(1.5,1.2,'*r');
% However, we can see the point is not contained in the sets

% If we were to compute the hypercube hull of the sets (overapprox), the
% set would violate the property
Sets_approx = Star.get_hypercube_hull(Sets);
figure;
Box.plot(Sets_approx)
hold on;
Star.plots(Sets);
hold on;
plot(1.5, 1.2, '*b');
