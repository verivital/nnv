classdef EXPT < handle
% stores info of experiment configuration as
% well as the experimental result and
% manages saving and loading of it to and from
% yaml files
    
    properties
        data = [];
        file = "";
        loaded_from_file = 0;
    end
    
    
    methods
        
        % load data from yaml file
        function obj = EXPT(file)
            
            if nargin < 1
                obj = EXPT.create_test_expt();
                return
            end
            
            if strcmp(file, "empty_expt")
                return
            end
            
            obj.file = file;
            obj.adjust_extension;
            obj.data = yaml.loadFile(obj.file);
            obj.loaded_from_file = 1;
        end
        
        
        function save(obj)
            file_exists = isfile(obj.file);
            if file_exists && ~obj.loaded_from_file && ~strcmp(obj.file, "results/test_file.yaml")
                error("File " + obj.file + " already exists; not overwriting.")
            end
            if file_exists
                movefile(obj.file, obj.file + ".bak");
            elseif obj.loaded_from_file
                error("Current configuration has been loaded from file " + file + " but the file does not exist.")
            end
            yaml.dumpFile(obj.file, obj.data);
            if file_exists
                delete(obj.file + ".bak");
            end
        end
        
        
        function adjust_extension(obj)
            
            if isa(obj.file, "string")
                obj.file = char(obj.file);
            end
            
            if length(obj.file) < 5 || ~strcmp(obj.file(end - 4:end), ".yaml")
                obj.file = [obj.file '.yaml'];
            end
            obj.file = string(obj.file);
            
        end
        
        
        function plot_results(obj, args)
            
            arguments
                obj
                args.figure_no = 5
                args.reverse_order = 1
                args.disp_legend = 1
                args.x_axis_factor = 0.01
                args.per_image_execution_time = 1
                args.layers_to_plot_for_from_yaml_file = nan;
            end
            
            pert_layer = obj.data.layers;
            if ~isnan(args.layers_to_plot_for_from_yaml_file)
                pert_layer = pert_layer(args.layers_to_plot_for_from_yaml_file);
            end
            
            try
                clf(args.figure_no);
            catch ME
                if ~strcmp(ME.identifier, 'MATLAB:clf:InvalidFigureHandle')
                    rethrow(ME)
                end
            end
            
            figure(args.figure_no);
            
            til = tiledlayout(2, length(pert_layer), 'TileSpacing','compact');
            percentage_axes = nan(1, length(pert_layer));
            time_axes = percentage_axes;
            
            if args.reverse_order
                new_pert_layer = cell(size(pert_layer));
                for k = 1:length(pert_layer)
                    new_pert_layer{k} = pert_layer{end - k + 1};
                end
                pert_layer = new_pert_layer;
            end
            
            mid_layer_no = length(pert_layer)/2;
            if mod(length(pert_layer), 2) == 0
                xlabel_offset = 4;
            else
                mid_layer_no = ceil(mid_layer_no);
                xlabel_offset = 0;
            end
            
            for pert_layer_no = 1:length(pert_layer)    
                percentages = [];
                times_y = [];
                pert_fracs = [];
                for frac_no = 1:length(pert_layer{pert_layer_no}.fracs)
                    percentages = [percentages; cell2mat(pert_layer{pert_layer_no}.fracs{frac_no}.percents)];
                    times_y = [times_y; cell2mat(pert_layer{pert_layer_no}.fracs{frac_no}.times)];
                    pert_fracs = [pert_fracs pert_layer{pert_layer_no}.fracs{frac_no}.frac];
                end
                if any(isnan([percentages times_y]), 'all')
                    continue;
                end
                
                if args.per_image_execution_time
                    times_y = times_y/obj.data.n_images;
                end
                
                if sum(percentages(:, 2:3) + times_y(:, 2:3), 'all') == 0
                    percentages = percentages(:, 1);
                    times_y = times_y(:, 1);
                end
                
                percentage_axes(pert_layer_no) = nexttile(pert_layer_no);
                % set(gca, 'XScale', 'log')
                
                % bar_plot_xticks = string(pert_fracs);
                pert_fracs_for_xticks = pert_fracs;
                if args.x_axis_factor == 0.01
                    pert_fracs_for_xticks = pert_fracs_for_xticks*100;
                end
                bar_plot_xticks = string(compose('%.2g', pert_fracs_for_xticks));
                
                bar(bar_plot_xticks, percentages)
                set(gca, 'Xticklabel', [])
                if pert_layer_no == 1
                    ylabel('%age of Images Verified')
                % else
                %     set(gca, 'YTick', [])
                end
                
                n_neurons_color = 'black';
                n_weights_color = '0,0,0';
                if isfield(pert_layer{pert_layer_no}, 'split')
                    if pert_layer{pert_layer_no}.split && pert_layer{pert_layer_no}.n_neurons < pert_layer{pert_layer_no}.n_weights
                        n_neurons_color = 'red';
                    elseif ~pert_layer{pert_layer_no}.split
                        n_weights_color = '0,0.6,0';
                    end
                    
                end
                
                title({[char(strrep(pert_layer{pert_layer_no}.name, '_', '\_'))], ...
                    ['\color[rgb]{' n_weights_color '}' num2str(pert_layer{pert_layer_no}.n_weights) '\color{black}'], ...
                    ['\color{' n_neurons_color '}' num2str(pert_layer{pert_layer_no}.n_neurons) '\color{black}']})
                
                time_axes(pert_layer_no) = nexttile(length(pert_layer) + pert_layer_no);
                bar(bar_plot_xticks, times_y);
                
                if pert_layer_no == mid_layer_no
                    fraction_or_percentage = 'Fraction';
                    if args.x_axis_factor == 0.01
                        fraction_or_percentage = 'Percentage';
                    end
                    t = xlabel("Perturbation Magnitude as " + string(fraction_or_percentage) + " of Weights' Range", 'FontSize', 13);
                    t.VerticalAlignment = 'top';
                    
                    set(t, 'Position', get(t, 'Position') + [xlabel_offset 0 0]);
                    
                    % position = get(t, 'Position');
                    % position(2) = position(2) - 0.1;
                    % set(t, 'Position', position);
                end
                
                
                % % Displaying something next to the x-axis
                % 
                % % Get the current axes
                % ax = gca;
                % 
                % % Get the position of the axes in the figure
                % ax_position = get(ax, 'Position');
                % 
                % % Get the x-axis and y-axis tick values
                % x_ticks = double(get(ax, 'XTick'));
                % y_ticks = double(get(ax, 'YTick'));
                % 
                % % Convert axis limits to double
                % x_lim = double(ax.XLim);
                % y_lim = double(ax.YLim);
                % 
                % % Calculate the normalized position of the last x-tick
                % x_tick_pos = ax_position(1) + (x_ticks(end) - x_lim(1)) / diff(x_lim) * ax_position(3);
                % y_tick_pos = ax_position(2) + (y_ticks(1) - y_lim(1)) / diff(y_lim) * ax_position(4);
                % 
                % % Display the position
                % disp(['X-tick position: ', num2str(x_tick_pos)]);
                % disp(['Y-tick position: ', num2str(y_tick_pos)]);
                % 
                % % Add text next to the last x-tick
                % text(x_ticks(end) + 0.8 * diff(x_lim), y_lim(1) - 0.1 * diff(y_lim), 'Last X-Tick', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
                
                
                % set(gca, 'YScale', 'log')
                if pert_layer_no == 1
                    ylabel({'Average Execution', 'Time per image (s)'});
                    % ylabel('Execution time per image (s)')
                % else
                %     set(gca, 'YTick', [])
                end
                
            end
            
            % linkaxes(percentage_axes, 'xy')
            % linkaxes(time_axes, 'xy')
            linkaxes(percentage_axes, 'y')
            % linkaxes(time_axes, 'y')
            if ~isa(obj.data.title, 'yaml.Null')
                title(til, obj.data.title)
            end
            
            % % x-axis label for the whole figure
            % fraction_or_percentage = 'Fraction';
            % if args.x_axis_factor == 0.01
            %     fraction_or_percentage = 'Percentage';
            % end
            % t = xlabel(til, string(fraction_or_percentage) + " of Weights' Range used as Perturbation");
            % t.VerticalAlignment = 'top';
            % 
            % % position = get(t, 'Position');
            % % position(2) = position(2) - 0.1;
            % % set(t, 'Position', position);
            
            if args.disp_legend
                lg = legend({'ModelStar', 'Certificated-Robust', 'Formal-Robust'}, 'Orientation', 'horizontal');
                lg.Layout.Tile = 'North';
                lg.Box = 'off';
            end
            
            if isfield(obj.data, "position")
                fig = gcf;
                fig.Position = cell2mat(obj.data.position);
            end
            
            
        end
        
        
        function plot_large_mnist_results(obj)
            if strcmp(obj.file, 'results/Large_MNIST_with_splitting.yaml')
                obj.plot_results(disp_legend=0, reverse_order=1)
            else
                error("Wrong network.")
            end
        end
        
    end
    
    
    methods(Static)
        
        function e = create_and_save_expt_config_using_params(args)
            
            arguments
                args.plot_title
                args.matlabnet
                args.net_name
                args.n_images
                args.filename = []
                args.layer_names = []
                args.layer_indices = []     % not used unless args.layer_names is empty
                args.general_fracs = [0.0001 0.001 0.01]
                args.fracs_individual_layers = {}
            end
            
            file = args.filename;
            title = args.plot_title;
            net = matlab2nnv(args.matlabnet);
            net.Name = args.net_name;
            n_images = args.n_images;
            layer_names = args.layer_names;
            layer_indices = args.layer_indices;
            general_fracs = args.general_fracs;
            fracs_individual_layers = args.fracs_individual_layers;
            nn_info = analyzeNetwork(args.matlabnet);
            
            if isempty(file)
                file = "results/" + nn_name;
            end
            
            if isempty(layer_names)
                if isempty(layer_indices)
                    layer_indices = net.get_layer_indices(["Conv2DLayer", "FullyConnectedLayer"]);
                    layer_indices = layer_indices(end:-1:1);
                end
                for l = layer_indices
                    layer_names = [layer_names, string(net.Layers{l}.Name)];
                end
            else
                layer_indices = [];
                for l = 1:length(layer_names)
                    layer_indices = [layer_indices, net.name2indx(layer_names(l))];
                end
            end
            
            e = EXPT("empty_expt");
            
            e.file = file;
            e.adjust_extension;
            
            % prepare data object
            e.data.title = title;
            e.data.network = net.Name;
            e.data.n_images = n_images;
            e.data.layers = {};
            
            for k = 1:length(layer_indices)
                
                l = layer_indices(k);
                e.data.layers{k}.name = net.Layers{l}.Name;
                e.data.layers{k}.n_weights = numel(net.Layers{l}.Weights);
                tmp = nn_info.LayerInfo(l, 3);
                tmp = tmp{1, 1};
                e.data.layers{k}.n_neurons = prod(cell2mat(tmp(1)));
                
                if length(fracs_individual_layers) >= k && ~isempty(fracs_individual_layers{k})
                    fracs = fracs_individual_layers{k};
                else
                    fracs = general_fracs;
                end
                e.data.layers{k}.fracs = {};
                for m = 1:length(fracs)
                    e.data.layers{k}.fracs{m}.frac = fracs(m);
                    e.data.layers{k}.fracs{m}.percents = nan(3,1);
                    e.data.layers{k}.fracs{m}.times = nan(3,1);
                end
            end
            
            % save the object on disk
            e.save;
        end
        
        function expt = create_test_expt()
            % small NN on MNIST
            nn_name = 'Very_Small_MNIST';
            folder = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, 'NN', filesep, 'MNIST', filesep];
            mnist_model = load([folder 'mnist_model.mat']);
            matlabnet = mnist_model.net;
            
            % create data object
            expt = EXPT.create_and_save_expt_config_using_params( ...
                filename = "results/test_file", ...
                plot_title = "Entire Weights Tensor of Single Layer Perturbed", ...
                matlabnet = matlabnet, ...
                net_name = nn_name, ...
                n_images = 10);
        end
        
        function expt = create_vgg19_expt()
            % experiment for VGG19
            nn_name = "vgg19";
            matlabnet = imagePretrainedNetwork(nn_name);
            
            % create data object
            expt = EXPT.create_and_save_expt_config_using_params( ...
                filename = "results/vgg19", ...
                plot_title = "Entire Weights Matrix of Single Layer Perturbed", ...
                matlabnet = matlabnet, ...
                net_name = nn_name, ...
                n_images = 10);
        end
        
        function expt = create_squeezenet_expt()
            % experiment for squeezenet
            nn_name = "squeezenet";
            matlabnet = imagePretrainedNetwork(nn_name);
            
            % create data object
            expt = EXPT.create_and_save_expt_config_using_params( ...
                filename = "results/" + nn_name, ...
                plot_title = "Entire Weights Matrix of Single Layer Perturbed", ...
                matlabnet = matlabnet, ...
                net_name = nn_name, ...
                n_images = 10);
        end
        
        function expt = create_resnet18_expt()
            % experiment for ResNet-18
            nn_name = "resnet18";
            matlabnet = imagePretrainedNetwork(nn_name);
            
            % create data object
            expt = EXPT.create_and_save_expt_config_using_params( ...
                filename = "results/" + nn_name, ...
                plot_title = "Entire Weights Matrix of Single Layer Perturbed", ...
                matlabnet = matlabnet, ...
                net_name = nn_name, ...
                n_images = 10);
        end
        
        function expt = create_MNIST_expt(nn_size)
            % Load an MNIST network from the ImageStar paper
            folder = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Submission', filesep, 'CAV2020_ImageStar', filesep, 'MNIST_NETS', filesep, 'Architecture', filesep];
            nn_size = lower(char(nn_size));
            nn_size(1) = upper(nn_size(1));
            nn_name = [nn_size '_MNIST'];
            model_mat_file = load([folder nn_size '_ConvNet.mat']);
            matlabnet = model_mat_file.net;
            
            % create data object
            expt = EXPT.create_and_save_expt_config_using_params( ...
                filename = "results/" + nn_name, ...
                plot_title = "Entire Weights Matrix of Single Layer Perturbed", ...
                matlabnet = matlabnet, ...
                net_name = nn_name, ...
                n_images = 50);
        end
        
        function expt = create_MNIST_MLP_expt()
            % /data/usama/star_stuff/nnv_stuff/nnv_weight_perturb/code/nnv/examples/Tutorial/NN/my_MNIST_v12_comparison_with_others_bigger_architecture
            folder = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, 'NN', filesep, 'my_MNIST_v12_comparison_with_others_bigger_architecture', filesep];
            model_mat_file = load([folder 'mnist_model_fc.mat']);
            matlabnet = model_mat_file.net;
            nn_name = "MNIST_MLP";
            
            % create data object
            expt = EXPT.create_and_save_expt_config_using_params( ...
                filename = "results/" + nn_name, ...
                plot_title = "Entire Weights Matrix of Single Layer Perturbed", ...
                matlabnet = matlabnet, ...
                net_name = nn_name, ...
                n_images = 200);
            
        end
        
        
    end
    
    
end