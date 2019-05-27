ls

cd nnv;

cd '/code/nnv/examples/NNCS/ACC/Verification/Scenarios 2'

pwd

cd /code/nnv;

pwd

install;

%dir_root = pwd; % note that install does a clear, must be defined here

'test after install'

cd /code/nnv

cd tests;

pwd

cd set;

%cd star;

cd zono;

% run a few tests
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

i_d = 1; disabledTests = {} ; global i_d disabledTests;

cd '/code/nnv/tests';

pwd

% try all test run
run_all_tests(dir_results, disabledTests, i_d) % add directory for results input and use dir_root/results

'after all tests run'

'start nncs tests'
cd /code/nnv/examples/NNCS/InvertedPendulum;
reach_inverted_pendulum_control_sys;
saveas(gcf, '/results/results_nncs_ip_fig1.png')

return;

% some bugs in next, not sure if it's due to cora submodule and we need to upload, or what
cd '/code/nnv/examples/NNCS/ACC/Verification/Scenarios 1'
verify_controller_3_20;
saveas(gcf, '/results/results_nncs_acc_s1_fig1.png')

cd '/code/nnv/examples/NNCS/ACC/Verification/Scenarios 2'
reach_nncACCsystem;
saveas(gcf, '/results/results_nncs_acc_s2_fig1.png')

'after nncs tests'