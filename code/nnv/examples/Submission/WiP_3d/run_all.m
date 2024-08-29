% %% Shape only data (3d)
% 
% disp("... adrenal ...");
% verify_adrenal;
% 
% disp("... vessel ...");
% verify_vessel;
% 
% 
% %% Volume data (general 3D)
% 
% 
% disp("... fracture ...");
% verify_fracture;
% 
% 
% disp("... nodule ...")
% verify_nodule;
% 
% 
% disp("... organ ...")
% verify_organ;
% 
% 
% disp("... synapse ...")
% verify_synapse;


% disp("... Creating plots...");
% visualize_results_gen;
% visualize_results_shape;


%% GPU

disp("... adrenal ...");
verify_adrenal_gpu;

disp("... vessel ...");
verify_vessel_gpu;

disp("... fracture ...");
verify_fracture_gpu;


disp("... nodule ...")
verify_nodule_gpu;


disp("... organ ...")
verify_organ_gpu;


disp("... synapse ...")
verify_synapse_gpu;