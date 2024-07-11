%% Run as many benchmarks as possible from 2024

vnncomp_path = "C:\Users\diego\Documents\Research\vnncomp2024_benchmarks\benchmarks\";
% vnncomp_path = "/home/manzand/Documents/MATLAB/vnncomp2024_benchmarks/benchmarks/";
% vnncomp_path = "/home/dieman95/Documents/MATLAB/vnncomp2024_benchmarks/benchmarks/";

benchmarks = dir(vnncomp_path);

notSupported = {'test'; 'traffic_signs_recognition'; 'cctsdb_yolo'; 'linearizenn'}; % skip for now, not even for falsification

regularTrack = {'acasxu'; 'nn4sys'; 'cora'; 'linearizenn'; 'safenlp'; 'dist_shift'; 'cifar100';...
    'tinyimagenet'; 'cgan'; 'metaroom'; 'tllverifybench'; 'collins_rul';
    };

extendedTrack = {'ml4acopf'; 'lsnc'; 'yolo'; 'cctsdb_yolo'; 'collins_aerospace';...
    'traffic_signs_recognition'; 'vggnet'; 'vit'; 
    }; % we don't really care much about this track, focus on the other one

% for i=3:length(benchmarks)
for i=13 % only do acasxu

    name_noyear = split(benchmarks(i).name, "_");
    if length(name_noyear) > 1
        name_noyear = strjoin(name_noyear(1:end-1), '_');
    else
        name_noyear = name_noyear{1};
    end

    if contains(name_noyear, regularTrack) % evaluate only regular track

    % if contains(extendedTrack, benchmarks(i).name) % evaluate only extended track

        if ~contains(name_noyear, notSupported) % skip evaluation of benchmarks not supported
    
            benchpath = vnncomp_path + benchmarks(i).name;
        
            instances = readmatrix(benchpath + "/instances.csv", "OutputType","string", "Delimiter",",");
        
            results_dir = "results_approx_" + benchmarks(i).name;
            mkdir(results_dir);
        
            for k=2:size(instances, 1)
                
                onnx = benchpath + filesep + instances(k,1);
                vnnlib = benchpath + filesep + instances(k,2);
        
                try
                    run_vnncomp2024_instance(name_noyear, onnx, vnnlib, results_dir + "/instance_" + string(k)+".txt");
                catch ME
                    fid = fopen(results_dir + "/instance_" + string(k)+".txt", 'w');
                    fprintf(fid, [ME.message '\n']);
                    fclose(fid);
                end
        
            end
    
        end

    end

end