

%% test 1: constructor
L1 = GlobalAveragePooling2DLayer();


%% test 2: inference
L1 = GlobalAveragePooling2DLayer();
x = load('one_image.mat');
x = x.one_image;
L1.evaluate(x);


%% test 3: equivalence (inference)
L1 = GlobalAveragePooling2DLayer();
x = load('one_image.mat');
x = x.one_image;
y = L1.evaluate(x);

dlX = dlarray(x, 'SSBC');
dlY = avgpool(dlX,'global');

assert(all(dlY == y, 'all'));

%% test 4: inference, higher dimension

miniBatchSize = 10;
inputSize = [5 5];
numChannels = 3;
X = rand(inputSize(1),inputSize(2),numChannels,miniBatchSize);

L1 = GlobalAveragePooling2DLayer();
Y = L1.evaluate(X);

dlX = dlarray(X,'SSCB');
dlY = avgpool(dlX,'global');
dlY = extractdata(dlY);

assert(all(dlY == Y, 'all'));

%% test 5: reachability

x = load('one_image.mat');
X = x.one_image;

lb = X - 0.1;
ub = X + 0.1;
IS = ImageStar(lb,ub);

L1 = GlobalAveragePooling2DLayer();
Y = L1.evaluate(X);
Yset = L1.reach(IS,'approx-star');

[LB,UB] = Yset.estimateRanges;

assert(all(LB <= Y,'all'))
assert(all(UB >= Y,'all'))

%% test 6: reach (sound)

N = 100; % random samples

x = load('one_image.mat');
X = x.one_image;

lb = X - 0.1;
ub = X + 0.1;
IS = ImageStar(lb,ub);
x_samples = IS.sample(N);

L1 = GlobalAveragePooling2DLayer();
Yset = L1.reach(IS,'approx-star');
[LB,UB] = Yset.estimateRanges;

for i=1:N
    xi = x_samples{i};
    Yi = L1.evaluate(xi);
    assert(all(LB <= Yi,'all'))
    assert(all(UB >= Yi,'all'))
end


%% test 7: reachability

miniBatchSize = 1;
inputSize = [5 5];
numChannels = 3;
X = rand(inputSize(1),inputSize(2),numChannels,miniBatchSize);

lb = X - 0.1;
ub = X + 0.1;
IS = ImageStar(lb,ub);

L1 = GlobalAveragePooling2DLayer();
Y = L1.evaluate(X);
Yset = L1.reach(IS,'approx-star');

[LB,UB] = Yset.estimateRanges;

assert(all(LB <= Y,'all'))
assert(all(UB >= Y,'all'))

%% test 8: reach (sound)

N = 200; % random samples

miniBatchSize = 1;
inputSize = [5 5];
numChannels = 3;
X = rand(inputSize(1),inputSize(2),numChannels,miniBatchSize);

lb = X - 0.1;
ub = X + 0.1;
IS = ImageStar(lb,ub);

x_samples = IS.sample(N);

L1 = GlobalAveragePooling2DLayer();
Yset = L1.reach(IS,'approx-star');
[LB,UB] = Yset.estimateRanges;

for i=1:N
    xi = x_samples{i};
    Yi = L1.evaluate(xi);
    assert(all(LB <= Yi,'all'))
    assert(all(UB >= Yi,'all'))
end

