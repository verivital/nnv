ls

cd nnv;

install;

'test after install'

cd examples;

cd NN;


cd /code/nnv

cd tests;

pwd

% run_all_tests % seems not to work, probably plotting or exception handling

cd set;

%cd star;

cd zono;

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
%test_star_affineMap
