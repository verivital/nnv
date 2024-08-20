%% Shape only data (3d)

disp("... adrenal ...");
verify_adrenal;

disp("... vessel ...");
verify_vessel;


%% Volume data (general 3D)


disp("... fracture ...");
verify_fracture;


disp("... nodule ...")
verify_nodule;


disp("... organ ...")
verify_organ;


disp("... synapse ...")
verify_synapse;


disp("... Creating plots...");
visualize_results_gen;
visualize_results_shape;
