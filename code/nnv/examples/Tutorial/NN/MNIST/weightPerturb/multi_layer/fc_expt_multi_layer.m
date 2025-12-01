%% Robustness verification of a NN (L infinity adversarial attack)
%  if f(x) = y, then forall x' in X s.t. ||x - x'||_{\infty} <= eps,
%  then f(x') = y = f(x)

clear

NPerClass = 1;   % test this many correctly classified images per class
pert_layer =   [7  7  9];     % layers to perturb one at a time
pert_layer_2 = [9 11 11];
% pert_layer = pert_layer(2:end);     % layers to perturb one at a time

figure(5);

til = tiledlayout(2, length(pert_layer), 'TileSpacing','compact');
percentage_axes = nan(1, length(pert_layer));
time_axes = percentage_axes;

for pert_layer_no = 1:length(pert_layer)
    % Load network 
    folder = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, 'NN', filesep, 'MNIST', filesep, 'weightPerturb', filesep];
    mnist_model = load([folder 'mnist_model_fc.mat']);
    
    times_modelstar = [];
    times_towards = [];
    percentages_verif_modelstar = [];
    percentages_verif_towards = [];
    
    
    % pert_fracs = round(logspace(-3, -1, 7), 1, 'significant');
    
    % pert_fracs = round(logspace(-3, -1, 10), 1, 'significant');
    
    % pert_fracs = [0.0001 linspace(0.002, 0.02, 10) 0.05];
    
    % pert_fracs = round(logspace(-3.7, -2.4, 6), 1, 'significant');
    
    % pert_fracs = round(logspace(-3.7, -1.7, 10), 1, 'significant');
    
    l = pert_layer(pert_layer_no);
    no_of_perts = 4;
    pert_fracs = round(logspace(-4, -3, no_of_perts), 1, 'significant');
    
    if l == 3
        pert_fracs = linspace(0.0002, 0.0005, no_of_perts);
    elseif l == 5
        pert_fracs = linspace(0.0004, 0.001, no_of_perts);
    elseif l == 7
        pert_fracs = [0.00001, 0.0001, 0.001, 0.002];
        if pert_layer_2(pert_layer_no) == 9
            pert_fracs = [0.00001, 0.0001, 0.0002, 0.0003];
        end
    elseif l == 9 || l == 11
        pert_fracs = linspace(0.0002, 0.0008, no_of_perts);
    elseif l == 13
        pert_fracs = linspace(0.005, 0.02, no_of_perts);
    end
    
    for f = pert_fracs
        % Create NNV model
        net = matlab2nnv(mnist_model.net);
        netbs = getByteStreamFromArray(net);    % bytestream to copy this network
            % for weight perturbations in parfor loop iterations
        
        % Load data (no download necessary)
        digitDatasetPath = fullfile(matlabroot,'toolbox','nnet','nndemos', ...
            'nndatasets','DigitDataset');
        % Images
        imds = imageDatastore(digitDatasetPath, ...
            'IncludeSubfolders',true,'LabelSource','foldernames');
        
        
        % Load certain number of images from each class
        % Evaluate robustness over them
        
        l1 = pert_layer(pert_layer_no);
        l12 = pert_layer_2(pert_layer_no);
        l2 = 13;
        
        p = f*net.get_weights_range(l1);
        
        fprintf("\n\n############--- Layers %s and %s frac %g pert mag. %g ---##############\n", net.Layers{l1}.Name,  net.Layers{l12}.Name, f, p)
        
        % subject entire weights matrices of layers l1 and l12 to
        % perturbation p
        net.Layers{l1}.perturb_whole_layer(-p, p);
        net.Layers{l12}.perturb_whole_layer(-p, p);
        
        % Define the reachability options
        reachOptions = struct; % initialize
        reachOptions.reachMethod = 'approx-star'; % using approximate method
        reachOptions.numCores = 32;
        reachOptions.device = 'cpu';
        % reachOptions.delete_old_sets = 1;
        reachOptions.free_mem_frac_for_LPs = 0.1;
        % reachOptions.dis_opt = 'display';
        % reachOptions.layer_specific_numCores = dictionary( ...
        %     ["start", "ReluLayer", "MaxPooling2DLayer", "end"], ...
        %     [  64   ,     -1     ,         -1         ,  -1  ]);
        reachOptions.layer_specific_numCores = dictionary( ...
            ["start", "ReluLayer", "MaxPooling2DLayer", "end"], ...
            [  64   ,     64     ,         64         ,  64  ]);
        % parpool(parcluster('local').NumWorkers);
        
        
        verify_by_sampling = 0;
        plot_range = 0;
        
        
        % choose the images to test
        numClasses = net.OutputSize;
        
        image_nos = NaN(1, NPerClass*numClasses);   % indices of images to use for
            % robustness verification
        image_nos_index = 1;
        img_no = 1;
        no_of_images_chosen_from_class = 0;
        class_no = 1;
        while class_no <= numClasses
            % Load one image in dataset
            [img, fileInfo] = readimage(imds, img_no);
            img = single(img); % change precision
            target = single(fileInfo.Label);
        
            Y_outputs = net.evaluate(img); 
            [~, yPred] = max(Y_outputs); % (expected: yPred = target)
        
            if target == yPred
                image_nos(image_nos_index) = img_no;
                image_nos_index = image_nos_index + 1;
                no_of_images_chosen_from_class = no_of_images_chosen_from_class + 1;
                if no_of_images_chosen_from_class == NPerClass
                    class_no = class_no + 1;
                    no_of_images_chosen_from_class = 0;
                    img_no = (class_no - 1)*1000;
                end
            end
            img_no = img_no + 1;
        end
        
        
        disturbance = zeros(size(img), 'single');
        
        % run robustness verification for chosen images in parallel
        no_of_images_verified = 0;
        no_of_images_verified_towards = 0;
        
        t_modelstar = 0;
        t_towards = 0;
        
        % t_total = tic;
        
        % parfor image_index = 1:length(image_nos)
        for image_index = 1:length(image_nos)
            % Load one image in dataset
            [img, fileInfo] = readimage(imds, image_nos(image_index));
            img = single(img); % change precision
            target = single(fileInfo.Label);
            
            tic;
            
            result = net.verify_robustness_for_3dim_img(reachOptions, input = img); % target not provided because
                % the images to be verified are selected for verification
                % in load_images only if the NN classifies them correctly
                % in the absence of perturbations.
            
            if result == 1
                no_of_images_verified = no_of_images_verified + 1;
                fprintf("+")
            else
                fprintf("-")
            end
            
            t_modelstar = t_modelstar + toc;
            
            
            % (existing approach) Towards Certificated ...
            
            tic;
            
            Y_outputs = net.evaluate(img); 
            [~, yPred] = max(Y_outputs); % (expected: yPred = target)
            
            xN = net.input_sets{l1}.V(:,:,:,1);
            xN = xN(:);
            
            pert_vec = p*ones(1, length(xN));
            
            % get bounds of the layer before l12 because l12 will be perturbed as well
            [pre_activation_lb_before_l12, pre_activation_ub_before_l12] = compute_bounds_weng_2020(net, pert_vec, inf, l1, l12 - 2);
            lb_before_l12 = max(0, pre_activation_lb_before_l12);
            ub_before_l12 = max(0, pre_activation_ub_before_l12);
            
            [towards_lb, towards_ub] = propagate_bounds_weng_2020(net, p, l12, l2, lb_before_l12, ub_before_l12);
            
            towards_ub_without_yPred = towards_ub;
            towards_ub_without_yPred(yPred) = NaN;
            
            if max(towards_ub_without_yPred) < towards_lb(yPred)
                no_of_images_verified_towards = no_of_images_verified_towards + 1;
            end
            
            t_towards = t_towards + toc;
        end
        
        disp(' ')
        disp(['With perturbation fraction of weight range f = ' + string(f) + ', '])
        
        percent_verif_modelstar = no_of_images_verified/length(image_nos)*100;
        disp(['Percentage of images verified by Modelstar: ' + string(percent_verif_modelstar) + ' % in ' + string(t_modelstar) + ' s.'])
        
        percent_verif_towards = no_of_images_verified_towards/length(image_nos)*100;
        disp(['Percentage of images verified by "Towards": ' + string(percent_verif_towards) + ' % in ' + string(t_towards) + ' s.'])
    
        times_modelstar = [times_modelstar t_modelstar];
        times_towards = [times_towards t_towards];
        percentages_verif_modelstar = [percentages_verif_modelstar percent_verif_modelstar];
        percentages_verif_towards = [percentages_verif_towards percent_verif_towards];
    end
    
    %% Plot results
    
    layers = [];
    for k = 1:length(net.Layers)
        if isa(net.Layers{k}, 'FullyConnectedLayer')
            layers = [layers k];
        end
    end
    l = find(layers == l1);
    
    percentages = [percentages_verif_modelstar' percentages_verif_towards'];
    times_y = [times_modelstar' times_towards'];
    
    percentage_axes(pert_layer_no) = nexttile(pert_layer_no);
    % scatter(fs, y)
    % plot(pert_fracs, percentages_verif_modelstar, '-o')
    % set(gca, 'XScale', 'log')
    bar(string(pert_fracs), percentages)
    set(gca, 'Xticklabel', [])
    if pert_layer_no == 1
        ylabel('%age of Images Verified')
    % else
    %     set(gca, 'YTick', [])
    end
    title(['FC Layers ' num2str(find(layers == pert_layer(pert_layer_no))) ' & ' num2str(find(layers == pert_layer_2(pert_layer_no)))])
    
    time_axes(pert_layer_no) = nexttile(length(pert_layer) + pert_layer_no);
    n_images = NPerClass*numClasses;
    bar(string(pert_fracs*100), times_y/n_images)
    % set(gca, 'YScale', 'log')
    if pert_layer_no == 1
        ylabel({"Average Execution", "Time per Image (s)"})
    % else
    %     set(gca, 'YTick', [])
    end
    
end

save("results/fc_expt_multi_layer_results.mat")

%% Plot results

% linkaxes(percentage_axes, 'xy')
% linkaxes(time_axes, 'xy')
linkaxes(percentage_axes, 'y')
% linkaxes(time_axes, 'y')
whole_or_half = 'Entire';
% title(til, [whole_or_half ' Weights Matrix of Single Layer Perturbed'])
% t = xlabel(til, "Fraction of Weights' Range used as Perturbation");
% t.VerticalAlignment = 'cap';
xlabel(time_axes(2), "Perturbation Magnitude as Percentage of First Perturbed Layer's Weights Range")
% lg = legend({'ModelStar', 'Certificated-Robust', 'Formal-Robust'}, 'Orientation', 'horizontal');
lg = legend({'ModelStar', 'Certificated-Robust'}, 'Orientation', 'horizontal');
lg.Layout.Tile = 'North';
lg.Box = 'off';
