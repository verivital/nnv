function [status, tTime] = run_vnncomp_instance(category, onnx, vnnlib, outputfile)

% single script to run all instances (supported) from the vnncomp2025

status = 2; % unknown (to start with)

disp("We are running...")

tTime = 0;

fid = fopen(outputfile, 'w');
fprintf(fid, 'unknown \n');
fclose(fid);

end


