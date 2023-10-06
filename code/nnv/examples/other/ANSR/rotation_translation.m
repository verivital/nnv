% Create translation and rotation examples

rng(3);

% Create random set
I = ExamplePoly.randVrep;
I.outerApprox;
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b, I.Internal.lb, I.Internal.ub); % input star

% Visualize set
figure;
Star.plot(I,'r')
xlim([-2 2])
ylim([-2 2])


%% Translation
% Translate all possible directions 3 units

% 1) pos x axis 
xpos = I.affineMap(eye(2), [3;0]);
figure;
Star.plot(xpos,'b'); hold on

% 2) neg x axis
xneg = I.affineMap(eye(2), [-3;0]);
figure;
Star.plot(xneg,'b'); hold on

% 3) pos y axis
ypos = I.affineMap(eye(2), [0;3]);
figure;
Star.plot(ypos,'b'); hold on

% 4) neg y axis
yneg = I.affineMap(eye(2), [0;-3]);
figure;
Star.plot(yneg,'b')


%% Rotation
% rotate a few differet angles, show that 360 degrees is equal to original

% 1) rot +30 degrees
R1 = StarRotate(I, 30);
figure;
Star.plot(R1,'b')
xlim([-2 2])
ylim([-2 2])

% 2) rot -45 degrees
R2 = StarRotate(I,-45);
figure;
Star.plot(R2,'b')
xlim([-2 2])
ylim([-2 2])

% 3) rot 90 degrees
R3 = StarRotate(I, 90);
figure;
Star.plot(R3,'b')
xlim([-2 2])
ylim([-2 2])

% 4) rot 180 degrees
R4 = StarRotate(I, 180);
figure;
Star.plot(R4,'b')
xlim([-2 2])
ylim([-2 2])

% 5) rot 360 degrees
R5 = StarRotate(R4, 180);
figure;
Star.plot(R5,'b')
xlim([-2 2])
ylim([-2 2])


%% Helper functions

function Y = StarRotate(X, angle)
    % rotation function for 2D sets
    % X: set to rotate
    % angle: angle in degrees to rotate set
    % Y: output (rotated) set

    % rotate angle degrees
    angle = deg2rad(angle);

    % rotation matrix
    rot_mat = [cos(angle), -sin(angle);
        sin(angle), cos(angle)];
    
    % Compute rotated set
    Y = X;
    Y.V(:,2:end) = rot_mat * Y.V(:,2:end);
end