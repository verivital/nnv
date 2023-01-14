%% Example of how to verify a property using Stars and HalfSpaces

% Set (Star)
X = [0;0;0];
perturbation = 0.1;
XLower = X - perturbation;
XUpper = X + perturbation;
Set = Star(XLower, XUpper);

% Properties
% Property 1 (HalfSpace)
H1 = [1 0 0; 0 1 0; 0 0 1];
b1 = [0; 0; 0];
% Hs1 = HalfSpace(H1,b1);
% Property 2 (HalfSpace)
H2 = [1 0 0; 0 1 0; 0 0 1];
b2 = [-2; -2; -2];
% Hs2 = HalfSpace(H2,b2);
% Property 3 (HalfSpace)
H3 = [1 0 0; 0 1 0; 0 0 1];
b3 = [2; 2; 2];
% Hs3 = HalfSpace(H3,b3);

% Verification
verifyProperty(Set, H1, b1); % Property 1
verifyProperty(Set, H2, b2); % Property 2
verifyProperty(Set, H3, b3); % Property 3

% Comparing indexes
% Set (Star)
X = [0;1;-1];
perturbation = 0.1;
XLower = X - perturbation;
XUpper = X + perturbation;
Set = Star(XLower, XUpper);

% Properties
% Property 1 (HalfSpace)
H1 = [1 -1 0]; % verified (x1 - x2 <= 0)
b1 = [0];
% Hs1 = HalfSpace(H1,b1);
% Property 2 (HalfSpace)
H2 = [-1 1 0]; % falsified (-x1 + x2 <= 0)
b2 = [0];
% Hs2 = HalfSpace(H2,b2);
% Property 3 (HalfSpace)
H3 = [-1 0 1; 1 -1 0]; % verified [(-x1 + x3 <= 0) && (x1 - x2 <= 0)]
b3 = [0; 0];
% Hs3 = HalfSpace(H3,b3);

% Verification
verifyProperty(Set, H1, b1); % Property 1
verifyProperty(Set, H2, b2); % Property 2
verifyProperty(Set, H3, b3); % Property 3


%% Helper Functions

function verifyProperty(Set, H, b)
    S = Set.intersectHalfSpace(H,b);
    if isempty(S)
        disp("Violated")
    elseif S.isSubSet(Set) && Set.isSubSet(S)
        disp("Verified")
    else
        disp("Unknown if method is approx, violated otherwise")
    end
end
