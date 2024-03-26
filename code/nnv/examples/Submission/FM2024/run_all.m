%% Shape only data

% poolobj = gcp('nocreate');
% delete(poolobj);

% verify_adrenal;

% poolobj = gcp('nocreate');
% delete(poolobj);

% verify_vessel;

%% Volume data (general 3D)

% poolobj = gcp('nocreate');
% delete(poolobj);
% 
% disp("... fracture ...");
% try
%     verify_fracture;
% catch ME
%     warning(ME.message);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% poolobj = gcp('nocreate');
% delete(poolobj);
% 
% disp("... nodule ...")
% try
%     verify_nodule;
% catch ME
%     warning(ME.message);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% poolobj = gcp('nocreate');
% delete(poolobj);

% disp("... organ ...")
% try
%     verify_organ;
% catch ME
%     warning(ME.message);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% poolobj = gcp('nocreate');
% delete(poolobj);

disp("... synapse ...")
try
    verify_synapse;
catch ME
    warning(ME.message)
end