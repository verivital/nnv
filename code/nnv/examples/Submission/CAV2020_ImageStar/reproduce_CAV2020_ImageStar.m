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

% define output paths
path_out_mnist = [path_results(), filesep, 'MNIST', filesep];
path_out_vgg16 = [path_results(), filesep, 'vgg16', filesep];
path_out_vgg19 = [path_results(), filesep, 'vgg19', filesep];

mkdir(path_out_mnist);
mkdir(path_out_vgg16);
mkdir(path_out_vgg19);

% reproduce all computationally generated figures and tables
% each command indicates the figure/table number, with the details in the corresponding script that is called, e.g., plot_ranges for figure 8
% the diary commands save the matlab console interaction to a text file, which is easier to view than the full output log generated and displayed on screen

diary([path_out_mnist, 'figure8_mnist_small_log.txt'])
cd MNIST_NETS/Small
plot_ranges
saveas(gcf, [path_out_mnist, 'figure8_mnist_small.png']);
diary off;

% table 1
% for tables 1, 2, and 3, we reproduce a subset of the full tables due to runtime restrictions, which take ~1.5 hours for the short versions, versus ~11 hours for the full versions (run earlier in codeocean, see run #6830033)
diary([path_out_mnist, 'table1_mnist_small_log.txt'])
compare_star_absdom_short % ~ 5 min
compare_star_absdom % full version: ~10 hours, 41min total for small, medium, and large
diary off;

% table 2
% next together: > 1.5 hours for short version
diary([path_out_mnist, 'table2_mnist_medium_log.txt'])
cd ../Medium
compare_star_absdom_short
compare_star_absdom  % full version
diary off;

% table 3
diary([path_out_mnist, 'table3_mnist_large_log.txt'])
cd ../Large
compare_star_absdom_short
compare_star_absdom  % full version
diary off;

% figure 13 (appendix)
diary([path_out_mnist, 'figure13_mnist_architectures_log.txt'])
cd ../Architecture
plot_network_architectures;
saveas(gcf, [path_out_mnist, 'figure13_mnist_architectures.png']);
diary off;

cd(path_reproduce);

% table 4, vgg16 part
diary([path_out_vgg16, 'table4_vgg16_log.txt'])
cd VGG16/Compare_Polytope_ImageStar
verify_VGG16 % takes ~1:38 hours:min
diary off;

% table 5, vgg16 part
diary([path_out_vgg16, 'table5_vgg16_delta2e07_log.txt'])
cd ../Compare_Exact_vs_Approx
verify_robustness_delta_2e_07 % ~15 min
diary off;

diary([path_out_vgg16, 'table5_vgg16_delta1e07_log.txt'])
verify_robustness_delta_e_07 % ~15 min
diary off;

cd(path_reproduce);

% table 4, vgg19 part
diary([path_out_vgg19, 'table4_vgg19_log.txt']);
cd VGG19/Compare_Polytope_ImageStar
verify_VGG19 % 1:10
diary off;

% table 5, vgg19 part
diary([path_out_vgg19, 'table5_vgg19_delta2e07_log.txt'])
cd ../Compare_Exact_vs_Approx
verify_robustness_delta_2e_07 % ~34 min
diary off;
diary([path_out_vgg19, 'table5_vgg19_delta1e07_log.txt'])
verify_robustness_delta_e_07 % ~17 min
diary off;

cd ../Plot_Figures

% 1:05 total for figs 9-12
diary([path_out_vgg19, 'figure9_vgg19_log.txt'])
plot_vgg19_exact_range % ~3min
saveas(gcf, [path_out_vgg19, 'figure9_vgg19.png'])
diary off;

% figure 11
diary([path_out_vgg19, 'figure11_vgg19_log.txt'])
plot_vgg19_reachTime % ~4min
saveas(gcf, [path_out_vgg19, 'figure11_vgg19.png'])
diary off;

% figure 12
diary([path_out_vgg19, 'figure12_vgg19_log.txt'])
plot_vgg19_inputSize_effect % ~55min
saveas(gcf, [path_out_vgg19, 'figure12_vgg19.png'])
diary off;

% figure 10
diary([path_out_vgg19, 'figure10_vgg19_log.txt'])
plot_vgg19_counter_example % ~3min, possible OOM for some reason, doesn't seem to occur on other platforms, so maybe a platform-specific/codeocean problem, see workaround in the script file to prevent this, it appears to be a codeocean or matlab bug with respect to the job manager/scheduler
saveas(gcf, [path_out_vgg19, 'figure10_vgg19.png'])
diary off;