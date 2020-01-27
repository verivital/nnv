%install;

cd /code/nnv/examples/Submission/CAV2020/
is_codeocean

%verify_P0_N00_abs(1,1)
%verify_P0_N00_star_appr(1,1)
%verify_P0_N00_star(1,1)
%verify_P0_N00_zono(1,1)

cd ACC;
pwd

reproduce

return

% problem in plotting for this one on codeocean
%cd '../UUV Safety Mornitoring'
%pwd
%reach_run
%saveas(gcf, '/results/results_hscc2020re_uuv.png') % figure generation has a problem on codeocean for this one due to figure restrictions in codeoceans matlab integration

% return; % stop here, comment/remove this to run all tests, other examples, etc.

cd /code/nnv
cd tests;
pwd
cd set;
cd zono;

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
i_d = 1; disabledTests = {} ; global i_d disabledTests;

cd '/code/nnv/tests';
pwd

return

% run all tests, recursively runs all tests in the code/nnv/tests directory
run_all_tests(dir_results, disabledTests, i_d) % add directory for results input and use dir_root/results

'after all tests run'

% Next, start some of the examples, these are not run inside of the 
% batch test execution, and can be extended/modified to run other 
% examples.
%
% We do not batch all of these as the runtime for everything would be on 
% the order of several days at least.
'start nncs tests'
cd /code/nnv/examples/NNCS/InvertedPendulum;
reach_inverted_pendulum_control_sys;
saveas(gcf, '/results/results_nncs_ip_fig1.png')

return; % stop here

% some bugs in next, not sure if it's due to cora submodule and we need to upload (as codeocean does not support git submodules currently), or what
cd '/code/nnv/examples/NNCS/ACC/Verification/Scenarios 1'
verify_controller_3_20;
saveas(gcf, '/results/results_nncs_acc_s1_fig1.png')

cd '/code/nnv/examples/NNCS/ACC/Verification/Scenarios 2'
reach_nncACCsystem;
saveas(gcf, '/results/results_nncs_acc_s2_fig1.png')

'after nncs tests'