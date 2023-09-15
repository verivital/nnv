% Create set representations
rng(0);

%% Star and ImageStar

% ImageStar
lower_bound = -rand([32 32 3]);
upper_bound = rand([32 32 3]);
IMS = ImageStar(lower_bound, upper_bound);

% Star
x_lb = -rand(15, 1);
x_ub = rand(15, 1);
X = Star(x_lb, x_ub);


%% Zonotope / ImageZono

% Zono
center = rand(3,1);
generators = [1 0 0; 0 1 0; 1 0 -1];
Z = Zono(center, generators);
figure; Zono.plot(Z, 'm');
ZS = Z.toStar;
figure; Star.plot(ZS);

% ImageZono
IMZ = ImageZono(lower_bound, upper_bound);

IMS1 = IMZ.toImageStar; % convert to ImageStar


%% Polyhedron
I = ExamplePoly.randVrep; % create random set
figure; I.plot('color', 'r');

I1 = Conversion.toStar(I); % convert to Star
figure; Star.plot(I1, 'c');

% Can also define the Star set using the center, basis vectors and constraints directly
I.outerApprox;
V = [0 0; 1 0; 0 1];
I1b = Star(V', I.A, I.b, I.Internal.lb, I.Internal.ub); % input star
figure; Star.plot(I1b, 'c');

I2 = I1.toPolyhedron; % convert back to polyhedron
figure; I2.plot('color', 'b')


%% From a vnnlib file
vnnlibfile = "../../../../../data/ACASXu/vnnlib/prop_2.vnnlib";
property = load_vnnlib(vnnlibfile);

% vnnlib files define lower and upper bounds for every variable 
lb = property.lb
ub = property.ub

S = Star(lb, ub);

% We can only plot sets with <= 3 dimensions
% For higher dimensions, there are some options:

% 1) project set into 2d or 3D and then plot

S1 = S.affineMap([1 0 0 0 0; 0 1 0 0 0], []); % 2D, dims 1 and 2
figure;
Star.plot(S1);

S3 = S.affineMap([1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0], []); % 3D, dims 1,2 and 3
figure;
Star.plot(S3);

% 2) Plot an overapproximation (box around Star) of the set and plot

% 2D
figure;
Star.plotBoxes_2D(S,1,2,'c');
figure;
Star.plotBoxes_2D_noFill(S,3,5,'b');

% 3D
figure;
Star.plotBoxes_3D(S,1,2,4,'g');

