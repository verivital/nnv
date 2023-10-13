% Create translation and rotation examples

rng(3);

% Create random set
I = ExamplePoly.randVrep;
I.outerApprox;
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b, I.Internal.lb, I.Internal.ub); % input star


%% diagram with steps for the different transformations

% 0) original

R0 = I;

% 1) Translate
R1 = R0.affineMap(eye(2), [0;3]);

% 2) Translate and rotate
R2 = R1.affineMap(eye(2), [3;0]);
R2 = StarRotate(R2, 30);

% 3) Translate and scale
R3 = R2.affineMap(eye(2), [0;-3]);
R3.V(:,2:end) = R3.V(:,2:end) * 0.5;

% 4) Translate and rotate
R4 = R3.affineMap(eye(2), [0;-3]);
R4 = StarRotate(R4, 60);

% 5) Translate, rotate and scale
R5 = R4.affineMap(eye(2), [-3;0]);
R5 = StarRotate(R5, 90);
R5.V(:,2:end) = R5.V(:,2:end) * 2;

% 6) Translate and scale
R6 = R5.affineMap(eye(2), [-6;0]);
R6.V(:,2:end) = R6.V(:,2:end) * 2;

% 7) Translate and rotate
R7 = R6.affineMap(eye(2), [0; 3]);
R7 = StarRotate(R7, 90);

% 8) Translate, and scale
R8 = R7.affineMap(eye(2), [3;3]);
R8.V(:,2:end) = R8.V(:,2:end) * 0.5;


%% Visualize transformations

figure;
Star.plot(I, 'k');
hold on;
Star.plot(R1,'b');
hold on;
Star.plot(R2,'c');
hold on;
Star.plot(R3,'r');
hold on;
Star.plot(R4,'c');
hold on;
Star.plot(R5,'m');
hold on;
Star.plot(R6,'r');
hold on;
Star.plot(R7,'c');
hold on;
Star.plot(R8,'r');
% saveas(gcf, 'affineTransformations.png');


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

