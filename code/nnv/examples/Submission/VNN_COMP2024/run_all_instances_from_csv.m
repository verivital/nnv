%% Run as many benchmarks as possible from 2023

vnncomp_path = "C:\Users\diego\Documents\Research\vnncomp2023_benchmarks\benchmarks\";

benchmarks = dir(vnncomp_path);

notSupported = {'test'; 'traffic_signs_recognition'; 'cctsdb_yolo'}; % skip for now

% for i=4:length(benchmarks)
for i=4 % only do acasxu

    if ~contains(notSupported, benchmarks(i).name) % skip evaluation of benchmarks not supported

        benchpath = vnncomp_path + benchmarks(i).name;
    
        instances = readmatrix(benchpath + "/instances.csv", "OutputType","string", "Delimiter",",");
    
        results_dir = "results_approx_" + benchmarks(i).name;
        mkdir(results_dir);
    
        for k=1:size(instances, 1)
            
            onnx = benchpath + filesep + instances(k,1);
            vnnlib = benchpath + filesep + instances(k,2);
    
            try
                run_vnncomp2024_instance(benchmarks(i).name,onnx,vnnlib, results_dir + "/instance_" + string(k)+".txt");
            catch ME
                fid = fopen(results_dir + "/instance_" + string(k)+".txt", 'w');
                fprintf(fid, [ME.message '\n']);
                fclose(fid);
            end
    
        end

    end

end