% test constructor for ImageStar class
IM(:,:,1) = [1 1 0 1; 0 0 1 1; 1 0 1 0; 1 1 1 1]; % center image channel 1
IM(:,:,2) = [0 1 0 0; 1 0 0 1; 0 1 1 0; 0 0 0 1]; % center image channel 2
IM(:,:,3) = [1 1 1 1; 1 1 0 1; 0 1 1 0; 1 0 1 0]; % center image channel 3

LB(:,:,1) = [-0.1 -0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % attack on pixel (1,1) and (1,2)
LB(:,:,2) = [-0.1 -0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; 
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

%% 1) Get Range
image = ImageStar(IM, LB, UB);
[xmin, xmax] = image.getRange(1,1,1);
display(xmin);
display(xmax);

%% 2) Get Ranges
image = ImageStar(IM, LB, UB);
[xmin, xmax] = image.getRanges;

%% 3) Estimate Ranges
image = ImageStar(IM, LB, UB);
[xmin, xmax] = image.estimateRanges;

%% 4) Compare ranges (estimate >= get)
image1 = ImageStar(IM, LB, UB);
t = tic;
[xmin1, xmax1] = image1.estimateRanges;
t1 = toc(t);

% estimate ranges should be faster, but an overapproximation of getRanges
image2 = ImageStar(IM, LB, UB);
t = tic;
[xmin2, xmax2] = image2.getRanges;
t2 = toc(t);

disp("Estimate took " + string(t1) +" seconds vs getRanges, that run on " +string(t2) +"seconds");

assert(all(xmin1 - xmin2 <= eps, 'all'));
assert(all(xmax1 - xmax2 >= -eps, 'all'));

%% 5) Test from issues

V(1,1,1,1) = 0;
V(1,1,1,2) = -1;
V(1,1,1,3) = 1;
C = [0 0];
d = 0;
ub = [1; 1];
lb = -ub;

I1 = ImageStar(V, C, d, lb, ub);
I2 = ImageStar(V, C, d, lb, ub);

[xmin1, xmax1] = I1.estimateRanges;
[xmin2, xmax2] = I2.getRanges;


assert(all(xmin1 - xmin2 <= eps, 'all'));
assert(all(xmax1 - xmax2 >= -eps, 'all'));