
%shared variables
% image input set
IM(:,:,1) = [1 1 0 1; 0 0 1 1; 1 0 1 0; 1 1 1 1]; % center image channel 1
IM(:,:,2) = [0 1 0 0; 1 0 0 1; 0 1 1 0; 0 0 0 1]; % center image channel 2
IM(:,:,3) = [1 1 1 1; 1 1 0 1; 0 1 1 0; 1 0 1 0]; % center image channel 3

LB(:,:,1) = [-0.1 -0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % attack on pixel (1,1) and (1,2)
LB(:,:,2) = [-0.1 -0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; 
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image_star = ImageStar(IM, LB, UB);
image_zono = ImageZono(LB, UB);

% image input set
IMs(:,:,1) = single([1 1 0 1; 0 0 1 1; 1 0 1 0; 1 1 1 1]); % center image channel 1
IMs(:,:,2) = single([0 1 0 0; 1 0 0 1; 0 1 1 0; 0 0 0 1]); % center image channel 2
IMs(:,:,3) = single([1 1 1 1; 1 1 0 1; 0 1 1 0; 1 0 1 0]); % center image channel 3

LBs(:,:,1) = single([-0.1 -0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]); % attack on pixel (1,1) and (1,2)
LBs(:,:,2) = single([-0.1 -0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]); 
LBs(:,:,3) = LB(:,:,2);

UBs(:,:,1) = single([0.1 0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]);
UBs(:,:,2) = single([0.1 0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]);
UBs(:,:,3) = UB(:,:,2);

image_star_s = ImageStar(IMs, LBs, UBs);


%% test 1: FullyConnected Constructor

% construct FullyConnectedLayer objects
W = rand(4, 48);
b = rand(4,1);
fc = FullyConnectedLayer(W,b);
fc1 = FullyConnectedLayer('fc1', W, b);


%% test 2: FullyConnected reach star exact

% construct a FullyConnectedLayer object
W = rand(4, 48);
b = rand(4,1);
fc = FullyConnectedLayer(W,b);

fc.reach_star_single_input(image_star);
fc.reach_zono(image_zono);


%% test 3: FullyConnected reach zono exact

% construct a FullyConnectedLayer object
W = rand(4, 48);
b = rand(4,1);
fc = FullyConnectedLayer(W,b);

output = fc.reach_zono(image_zono);

%% test 4: Fullyconnected precision check (parameters)

rng(0);
W = rand(4, 48);
b = rand(4, 1);
fc = FullyConnectedLayer(W,b);
Rd = fc.reach_star_single_input(image_star);

W = single(W);
b = single(b);
fc = FullyConnectedLayer(W,b);
Rs = fc.reach_star_single_input(image_star);

V_error = abs(Rs.V - Rd.V);
max(V_error, [], 'all')

[lb_s, ub_s] = Rs.getRanges;
[lb_d, ub_d] = Rd.getRanges;

error_tolerance = 1e-5;
all(abs(lb_s-lb_d) <= error_tolerance)
all(abs(ub_s-ub_d) <= error_tolerance)

%% test 5: Fullyconnected precision check (set)

rng(0);
W = rand(4, 48);
b = rand(4, 1);
fc = FullyConnectedLayer(W,b);

Rd = fc.reach_star_single_input(image_star);
IS = image_star.changeVarsPrecision('single');
Rs = fc.reach_star_single_input(image_star);

V_error = abs(Rs.V - Rd.V);
max(V_error, [], 'all')

[lb_s, ub_s] = Rs.getRanges;
[lb_d, ub_d] = Rd.getRanges;

error_tolerance = 1e-5;
all(abs(lb_s-lb_d) <= error_tolerance)
all(abs(ub_s-ub_d) <= error_tolerance)


%% test 6: Fullyconnected evaluate precision check

rng(0);
W = rand(4, 48);
b = rand(4, 1);
fc = FullyConnectedLayer(W,b);
ydd = fc.evaluate(IM);
yds = fc.evaluate(single(IM));
fc = fc.changeParamsPrecision('single');
ysd = fc.evaluate(IM);
yss = fc.evaluate(single(IM));

% somehow these 2 are the same
assert(all(yss == yds)) 
assert(all(ysd == yss))

% But they are not if both are double precision
diff_sd = abs(ydd - yds);
max(diff_sd, [], 'all')


%% test 7: contain (inference in reach set) - Double

rng(0);
W = rand(4, 48);
b = rand(4, 1);
fc = FullyConnectedLayer(W,b);
% inference
y = fc.evaluate(IM);
% reach set
R = fc.reach_star_single_input(image_star);
R = R.toStar;

assert(R.contains(y));

[lb,ub] = R.estimateRanges;

assert(all(y > lb));
assert(all(y < ub));


%% test 8: contain (inference in reach set) - Single

rng(0);
W = rand(4, 48, 'single');
b = rand(4, 1, 'single');
fc = FullyConnectedLayer(W,b);
% inference
y = fc.evaluate(IMs);
y = double(y);
% reach set
R2 = fc.reach_star_single_input(image_star_s);
R2 = R2.toStar;
% R = R.changeVarsPrecision('double');

[lb,ub] = R2.estimateRanges;

assert(all(y > lb));
assert(all(y < ub));

