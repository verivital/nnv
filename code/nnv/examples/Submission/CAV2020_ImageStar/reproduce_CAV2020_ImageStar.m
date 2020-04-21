path_reproduce = pwd();

path_nnv_root = ['..', filesep, '..', filesep, '..', filesep]; % in cav2020 imagestar folder

cd(path_nnv_root);
% run installation if not on codeocean / not already set up
try
    is_codeocean();
catch
    install;
end

cd(path_reproduce);

path_out_mnist = [path_results(), filesep, 'MNIST', filesep];
path_out_vgg19 = [path_results(), filesep, 'vgg19', filesep];

mkdir();

cd MNIST_NETS/Small
plot_ranges
saveas(gcf, [path_out_mnist, 'figure8_mnist_small.png']);

% table 1
compare_star_absdom_short % ~ 5 min
%compare_star_absdom % full version: 10:41 total

% table 2
% next together: > 1.5 hours for short version
cd ../Medium
compare_star_absdom_short
%compare_star_absdom  % full version: 10:41 total

% table 3
cd ../Large
compare_star_absdom_short
%compare_star_absdom  % full version: 10:41 total


cd /code/nnv/examples/Submission/CAV2020_ImageStar/

% table 4, vgg16 part
cd VGG16/Compare_Polytope_ImageStar
verify_VGG16 % takes ~1:38 hours:min

% table 5, vgg16 part
cd ../Compare_Exact_vs_Approx
verify_robustness_delta_2e_07 % ~15 min

verify_robustness_delta_e_07 % ~15 min

cd /code/nnv/examples/Submission/CAV2020_ImageStar/

% table 4, vgg19 part
cd VGG19/Compare_Polytope_ImageStar
verify_VGG19 % 1:10

% table 5, vgg19 part
cd ../Compare_Exact_vs_Approx
verify_robustness_delta_2e_07 % ~34 min

verify_robustness_delta_e_07 % ~17 min

cd ../Plot_Figures

% 1:05 total for figs 9-12
% additionally, had an OOM, so checking
% figure 9
plot_vgg19_exact_range % ~3min
saveas(gcf, [path_out_vgg19, 'figure9_vgg19.png'])

% figure 10
%plot_vgg19_counter_example % ~3min, possible OOM
%saveas(gcf, [path_out_vgg19, 'figure10_vgg19.png'])

% figure 11
plot_vgg19_reachTime % ~4min
saveas(gcf, [path_out_vgg19, 'figure11_vgg19.png'])

% figure 12
plot_vgg19_inputSize_effect % ~55min
saveas(gcf, [path_out_vgg19, 'figure12_vgg19.png'])
