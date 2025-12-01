
if exist('do_not_clear', 'var') && do_not_clear
    disp('Maintaining existing workspace...')
else
    clear
    
    % expt_path = 'results/MNIST_MLP';
    expt_path = 'results/test_file';
    
    n_layers_to_run_for_from_yaml_file = nan;
end

delete(gcp('nocreate'));
dbstop if error

expt = EXPT(expt_path);

file = '/var/conv_expt_any_layer.log';
log_file_size_lim_B = 100*2^20;
file_size_B = dir(file).bytes;
if file_size_B > log_file_size_lim_B
    error("Log file size " + file_size_B/2^20 + " MB larger than the limit " + log_file_size_lim_B/2^20 + " MB.")
end
diary(file)

% Load network
nn_name = char(expt.data.network);
switch nn_name
    case {"MNIST_MLP"}
        if strcmp(nn_name, "MNIST_MLP")
            folder = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, 'NN', filesep, 'MNIST', filesep, 'weightPerturb', filesep];
            mnist_model = load([folder 'mnist_model_fc.mat']);
        end
        matlabnet = mnist_model.net;
        classes = string(1:10);
    otherwise
        error(expt.data.network + " not supported yet.")
end

% Create NNV model
net = matlab2nnv(matlabnet);
net.Name = nn_name;
net.OutputSize = length(classes);
netbs = getByteStreamFromArray(net);    % bytestream to copy this network
    % for weight perturbations in parfor loop iterations
        
% Define the reachability options
reachOptions = struct; % initialize
reachOptions.reachMethod = 'approx-star'; % using approximate method
reachOptions.numCores = 1;
% reachOptions.device = 'gpu';
reachOptions.device = 'cpu';

% Load certain number of images from each class and set some options
switch expt.data.network
    
    case {"Very_Small_MNIST", "Small_MNIST", "Medium_MNIST", "Large_MNIST", "MNIST_MLP"}
        database = "mnist";
        
    otherwise
        error("Not supported yet.")
end
[images, labels] = load_images(database = database, ...
    n = expt.data.n_images, ...
    matlabnet = matlabnet);


fprintf("\n\n++++++++++++++++++++++++++++++++ Expt: %s ++++++++++++++++++++++++++++++++\n", expt.file);

if isnan(n_layers_to_run_for_from_yaml_file)
    n_layers_to_run_for_from_yaml_file = length(expt.data.layers);
end

for pert_layer_no = 1:n_layers_to_run_for_from_yaml_file
    fprintf("\n\n------------------------------ Layer: %s ------------------------------\n", expt.data.layers{pert_layer_no}.name)
    
    l = net.name2indx(expt.data.layers{pert_layer_no}.name);
        
    for frac_no = 1:length(expt.data.layers{pert_layer_no}.fracs)
        f = expt.data.layers{pert_layer_no}.fracs{frac_no}.frac;
        fprintf("\n\n############################# Fraction: %g #############################\n", f)
        
        percentages = cell2mat(expt.data.layers{pert_layer_no}.fracs{frac_no}.percents);
        times_y = cell2mat(expt.data.layers{pert_layer_no}.fracs{frac_no}.times);
        
        if all(~isnan([percentages times_y]), 'all')
            continue;   % if all results in this case are computed, skip.
        end
        
        net = getArrayFromByteStream(netbs);
        
        % perturb all weights in layer l's weight matrix by a fraction of the
        % range of weights in layer l
        p = f*WPutils.get_weights_range(net, l);
        net.Layers{l}.perturb_whole_layer(-p, p);
        
        % run robustness verification for chosen images in parallel
        no_of_images_verified = 0;
        no_of_images_verified_towards = 0;
        no_of_images_verified_formalizing = 0;
        
        t_modelstar = 0;
        t_towards = 0;
        t_formalizing = 0;
        
        t_total = tic;
        fprintf("\n\n")
        for image_no = 1:length(images)
            
            if isfield(reachOptions, 'dis_opt') && strcmp(reachOptions.dis_opt, 'display')
                fprintf("\n\n")
                disp("%%%%%%%%%%%%% Processing image no. " + image_no + " %%%%%%%%%%%%%");
                disp(' ')
            end
            
            img = images{image_no};
            target = labels{image_no};
            correct_output = net.evaluate(img);
            
            % ModelStar
            tic;
            result = WPutils.verify_robustness_for_3dim_img(net, reachOptions, input = img); % target not provided because
                % the images to be verified are selected for verification
                % in load_images only if the NN classifies them correctly
                % in the absence of perturbations.
            if result == 1
                no_of_images_verified = no_of_images_verified + 1;
                fprintf("+")
            else
                fprintf("-")
            end
            time_this_iter = toc;
            t_modelstar = t_modelstar + time_this_iter;
            
            if time_this_iter > 60
                fprintf("Image %d took %d s, finishing at %s.\n", image_no, round(time_this_iter), datetime('now'))
            end
            
            only_fc_relu_layers_after_l = 1;
            try
                WPutils.verify_has_only_fc_relu_placeholder_layers(net, l, []);
            catch ME
                if strcmp(ME.identifier, 'NNV:non_supported_layers_in_fc_relu_only_method')
                    only_fc_relu_layers_after_l = 0;
                else
                    rethrow(ME);
                end
            end
            
            if isa(net.Layers{l}, "FullyConnectedLayer") && only_fc_relu_layers_after_l && ~strcmp(expt.data.network, "Large_MNIST")
                
                % Towards Certificated ...
                tic;
                [towards_lb, towards_ub] = WPutils.compute_bounds_weng_2020_fc_layers_only(net, ...
                    first_fc_layer = l, ...
                    pert_mag = p);
                towards_ub_without_target_ub = towards_ub;
                towards_ub_without_target_ub(target) = NaN;
                if max(towards_ub_without_target_ub) < towards_lb(target)
                    no_of_images_verified_towards = no_of_images_verified_towards + 1;
                end
                t_towards = t_towards + toc;
                
                % Formalizing Robustness ...
                tic;
                image_verified = WPutils.formal_robust_verify_fc_layers_only(net, ...
                    first_fc_layer = l, ...
                    pert_mag = p, ...
                    correct_output = correct_output);
                no_of_images_verified_formalizing = no_of_images_verified_formalizing + image_verified;
                t_formalizing = t_formalizing + toc;
            end
        end
        t_total = toc(t_total);
        
        t_measured_in_parfor_total = t_modelstar + t_towards + t_formalizing;
        
        t_modelstar = t_modelstar/t_measured_in_parfor_total*t_total;
        t_towards = t_towards/t_measured_in_parfor_total*t_total;
        t_formalizing = t_formalizing/t_measured_in_parfor_total*t_total;
        
        fprintf('\n\n----------------------- Layer fraction results -------------------------------\n')
        disp('       With perturbation fraction (of weight range) f = ' + string(f) + ', ')
        
        percent_verif_modelstar = no_of_images_verified/length(images)*100;
        disp('       Percentage of images verified by Modelstar: ' + string(percent_verif_modelstar) + ' % in ' + string(t_modelstar) + ' s.')
        
        percent_verif_towards = no_of_images_verified_towards/length(images)*100;
        disp('       Percentage of images verified by "Towards": ' + string(percent_verif_towards) + ' % in ' + string(t_towards) + ' s.')
        
        percent_verif_formalizing = no_of_images_verified_formalizing/length(images)*100;
        % disp(['       Percentage of images verified by "Formalizing": ' + string(percent_verif_formalizing) + ' % in ' + string(t_formalizing) + ' s.'])
        
        expt = EXPT(expt_path);
        expt.data.layers{pert_layer_no}.fracs{frac_no}.percents = num2cell([percent_verif_modelstar, percent_verif_towards, percent_verif_formalizing]);
        expt.data.layers{pert_layer_no}.fracs{frac_no}.times = num2cell([t_modelstar, t_towards, t_formalizing]);
        expt.save;
    end
    
    expt.plot_results(layers_to_plot_for_from_yaml_file=1:n_layers_to_run_for_from_yaml_file);
    
end

diary off