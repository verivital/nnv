% load from: data\examples\CNN\VGG19\DEEPFOOL_Attack
path_base = [nnvroot(), filesep, 'data', filesep, 'examples', filesep, 'CNN', filesep, 'VGG19', filesep, 'DEEPFOOL_Attack', filesep];

%% Construct input sets
dif_images = load([path_base, 'pepper_dif_images.mat']);
ori_images = load([path_base, 'pepper_ori_images.mat']);

dif_images = struct2cell(dif_images);
ori_images = struct2cell(ori_images);

N = 20; % choose 20 cases

l = [0.96 0.97 0.98];
delta = [0.0000001 0.0000002];

%l = 0.96;
%delta = 0.0000002;

P = length(l);
M = length(delta);

inputSetStar = cell(P, M);
correct_labels = cell(P,M);
for i=1:P
    for j=1:M
        lb = l(i);
        ub = l(i) + delta(j);
        IS(N) = ImageStar;
        for k=1:N
            V(:,:,:,1) = double(ori_images{k});
            V(:,:,:,2) = double(dif_images{k});
            IS(k) = ImageStar(V, [1;-1], [ub; -lb], lb, ub);
        end
        inputSetStar{i, j} = IS;
        correct_labels{i, j} = 946 * ones(1, N);
    end
end

%% parse VGG19
fprintf('\n\n=============================LOAD VGG19 ======================\n');

% Load the trained model 
if is_codeocean()
    % the vgg support packages require a gui to install
    % so, load alternatively from the data repository 
    % when running on codeocean / other headless platforms
    load('/data/vgg19_cache.mat');
    net = net_vgg19;
else
    net = vgg19();
end

fprintf('\n\n======================== PARSING VGG19 =======================\n');
nnvNet = matlab2nnv(net);

%% evaluate robustness

VT_star = zeros(P, M); % verification time of the approx-star method
VT_absdom = zeros(P, M); % verification time of the DeepPoly abstract domain method

r_star = zeros(P, M); % robustness percentage on an array of N tested input sets obtained by the approx-star method
r_absdom = zeros(P, M); % robustness percentage on an array of N tested input sets obtained by the DeepPoly abstract domain method


for i=1:P
    for j=1:M
        % Reachability analysis using approx-star
        tot_combo = 0;
        t = tic;
        for k=1:N
            res = nnvNet.verify_robustness(inputSetStar{i, j}(k), reachOpt_star, correct_labels{i, j}(k));
            tot_combo = tot_combo + res == 1;
        end
        r_star(i, j) = tot_combo;
        VT_star(i, j) = toc(t);

        % Reachability analysis using abs-dom
        tot_combo = 0;
        t = tic;
        for k=1:N
            res = nnvNet.verify_robustness(inputSetStar{i, j}(k), reachOpt_absdom, correct_labels{i, j}(k));
            tot_combo = tot_combo + res == 1;
        end
        r_absdom(i, j) = tot_combo;
        VT_absdom(i, j) = toc(t);
    end
end

% save VGG19_Results.mat r_star VT_star r_absdom VT_absdom;


%% Results - print to screen

fprintf('\n========================================================================================');
fprintf('\n          ROBUSTNESS VERIFICATION RESULTS (IN PERCENT) OF VGG19 UNDER DEEPFOOL ATTACK       ');
fprintf('\n========================================================================================\n\n');

for j=1:M
    fprintf("             delta = %.7f", delta(j));
end
fprintf("\n");
for j=1:M
    fprintf("           Polytope   ImageStar");
end
fprintf("\n");
for i=1:P
    fprintf("l = %.2f", l(i));
    for j=1:M
        fprintf("    %.2f        %.2f      ", 100*r_absdom(i, j), 100*r_star(i, j));
    end
    fprintf("\n");
end

fprintf('\n========================================================================================');
fprintf('\n                VERIFICATION TIMES (IN SECONDS) OF VGG19 UNDER DEEPFOOL ATTACK              ');
fprintf('\n========================================================================================\n\n');

for j=1:M
    fprintf("             delta = %.7f", delta(j));
end
fprintf("\n");
for j=1:M
    fprintf("           Polytope   ImageStar");
end
fprintf("\n");
for i=1:P
    fprintf("l = %.2f", l(i));
    for j=1:M
        fprintf("    %.2f       %.2f        ", VT_absdom(i, j), VT_star(i, j));
    end
    fprintf("\n");
end

