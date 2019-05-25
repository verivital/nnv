ls

cd nnv;

install;

%dir_root = pwd; % note that install does a clear, must be defined here

'test after install'

cd examples;

cd NN;


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

% run all tests doesn't work, I'm guessing what is happening is because some tests throw exceptions, even though they are handled, codeocean must be interpreting this as a failure and stopping the run, so we need to run only those that do not throw exceptions, or fix the ones with exceptions

'after all tests run'