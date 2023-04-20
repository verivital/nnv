% check if on codeocean, error out if not
try
    if ~is_codeocean() % test codeocean detection
        'ERROR: run_codeocean.m should only be executed on the codeocean platform'
        return;
    end
catch
    'ERROR: path likely not set, run install.m'
end

matlabshared.supportpkg.setSupportPackageRoot('/usr/local/MATLAB/R2022b');
addpath(genpath('/usr/local/MATLAB'))

% default output path for path_results, ensure logs subdirectory there
mkdir('/results/logs/')

cd /code/nnv/examples/NNV2.0/Submission/CAV2023
%RE_cav23_short
run_cav23

% ARCH-COMP 2022
%cd /code/nnv/examples/Submission/ARCH-COMP2022/benchmarks/
%run_all

return 
cd /code/nnv/examples/Submission/FORMATS2022/

run_subset % run a subset of all FORMATS results, much faster
% run_all_nnv % run the complete results of the FORMATS paper

return
run_all_nnv

return

% individual examples
%return
cd ACC
verify_acc_orig % error with time, apparently is a bug in cora fixed in our fork
%verify_acc_plant
%verify_acc_tanhplant
return

% cav 2021 segmentation

cd /code/nnv/examples/Submission/CAV2021/
reproduce_CAV2021

% arch-comp 2021 next
return
cd /code/nnv/examples/Submission/ARCH-COMP2021
run_all

return
cd /code/nnv/examples/Submission/ARCH_COMP2020/benchmarks
run_all


return

cd /code/nnv/examples/Submission/NeurIPS2020/
produce_NeurIPS2020


return

% CAV 2020 tool paper reproducibility
cd /code/nnv/examples/Submission/CAV2020/

% examples to manually run nnv on acas-xu, these (and many more) are 
% run from run_scripts automatically
%verify_P0_N00_abs(1,1)
%verify_P0_N00_star_appr(1,1)
%verify_P0_N00_star(1,1)
%verify_P0_N00_zono(1,1)

cd ACC;
pwd

% run all closed-loop CPS examples
reproduce % will take ~32.5 minutes (see run 140515 or 152224)

% CAV 2020 ImageStar paper reproducibility
cd /code/nnv/examples/Submission/CAV2020_ImageStar/

% see details in this script for each figure/table recreation
reproduce_CAV2020_ImageStar;

% short test (uncomment and comment the above reproduce_CAV2020_ImageStar)
%cd MNIST_NETS/Small/
%plot_ranges
%saveas(gcf, '/results/figure8_mnist_small.png');


%return % end of CAV 2020 ImageStar reproducibility

% if interested, one can also look at these others, as well as how to manually run new examples in CodeOcean below


return; % stop here, comment/remove this to run all tests, other examples, etc.
% if you run all these, it will add another ~15 minutes of runtime

cd /code/nnv/tests/set/zono

% run a few tests before starting batch run (as examples to illustrate how to run individiual examples / tests)
pwd
test_zono_convexHull
pwd
saveas(gcf, '/results/results_fig1.png')

't1'
% star tests, all failed (yalmip)
test_star_plotBoxes_2D
pwd
saveas(gcf, '/results/results_fig2.png')

't2'
test_star_plotBoxes_3D
saveas(gcf, '/results/results_fig3.png')

't3'
pwd

dir_results = '/results/';

% global variables for batch test runs
% some tests need to be disabled for codeocean due to differences on running without a UI, etc.
i_d = 1; disabledTests = {'test_CNN_parse'} ; global i_d disabledTests;

cd '/code/nnv/tests';
pwd

% run all tests, recursively runs all tests in the code/nnv/tests directory
run_all_tests(dir_results, disabledTests, i_d) % add directory for results input and use dir_root/results

'after all batch tests run'

% Next, start some of the examples, these are not run inside of the 
% batch test execution, and can be extended/modified to run other 
% examples.
%
% We do not batch all of these as the runtime for everything would be on 
% the order of several days at least.
%'start nncs tests'
%cd /code/nnv/examples/NNCS/InvertedPendulum;
%reach_inverted_pendulum_control_sys;
%saveas(gcf, '/results/results_nncs_invpend_fig1.png')

%cd '/code/nnv/examples/NNCS/ACC/Verification/Scenarios 1'
%verify_controller_3_20;
%saveas(gcf, '/results/results_nncs_acc_s1_fig1.png')

%cd '/code/nnv/examples/NNCS/ACC/Verification/Scenarios 2'
%reach_nncACCsystem;
%saveas(gcf, '/results/results_nncs_acc_s2_fig1.png')

%'after nncs tests'