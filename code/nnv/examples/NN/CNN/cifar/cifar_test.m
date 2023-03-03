% download and extract cifar data set, then load one of them
% eg
% https://www.cs.toronto.edu/~kriz/cifar.html
% cifar-10: https://www.cs.toronto.edu/~kriz/cifar-10-matlab.tar.gz
% cifar-100: https://www.cs.toronto.edu/~kriz/cifar-100-matlab.tar.gz
% extract
% load data_batch_1.mat

% todo: cache and check if already downloaded more generally for the files that will be loaded, this only checks if folder exists
if ~isfolder('cifar-10-batches-mat')
  websave('cifar10.tar.gz', 'https://www.cs.toronto.edu/~kriz/cifar-10-matlab.tar.gz');
  untar('cifar10.tar.gz');
end

% variable with images is named data
load('cifar-10-batches-mat/data_batch_1.mat');

% each row of data is 3072 long: RGB channels 3072/3 = 1024, the 32x32
% pixel color values

num_show = 12;
sp_w = ceil(sqrt(num_show));
sp_h = floor(sqrt(num_show));

for idx_n = 1 : num_show

idx = randi([1 length(data)]);

width = 32;
height = 32;
pixels = width * height;

rgb = [data(idx,1:pixels); data(idx,1+pixels:2*pixels); data(idx,1+2*pixels:3*pixels)];

red = [];
for i = 1 : height
    red(i,:) = rgb(1,(i-1) * width + 1 : i * width);
    blue(i,:) = rgb(2,(i-1) * width + 1 : i * width);
    green(i,:) = rgb(3,(i-1) * width + 1 : i * width);
end

% image is a width x height x color_channel matrix,
% so 32x32x3
% concatenate each color in the 3rd matrix dimension
rgb_rg = cat(3,red,green);
rgb_rgb = cat(3,rgb_rg,blue);

  subplot(sp_w,sp_h,idx_n)
  image(rgb_rgb)
  axis equal
end
