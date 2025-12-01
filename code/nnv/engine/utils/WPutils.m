classdef WPutils
    % Utility functions for weight perturbation analysis
    
    % Author: Muhammad Usama Zubair
    
    methods (Static)
        
        function print_layers_info(net)
            nn_info = analyzeNetwork(net.matlabnet);
            for l_no = 1:length(net.Layers)
                l = net.Layers{l_no};
                fprintf("Layer %2d: %-30s %-25s ", l_no, class(l), l.Name)
                try
                    fprintf("num_Weights: %-10d ", numel(l.Weights));
                    tmp = nn_info.LayerInfo(l_no, 3);
                    tmp = tmp{1, 1};
                    n_neurons = prod(cell2mat(tmp(1)));
                    fprintf("num_Neurons: %-10d", n_neurons);
                catch ME
                    if ~strcmp(ME.identifier, 'MATLAB:noSuchMethodOrField')
                        rethrow(ME)
                    end
                end
                disp(' ')
            end
        end
        
        % for comparison with another approach in ModelStar experiments. Can be removed.
        function image_verified = formal_robust_verify_fc_layers_only(net, args)
            
            arguments
                net
                args.first_fc_layer
                args.pert_mag
                args.correct_output
                args.last_fc_layer = []
            end
            
            p = args.pert_mag;
            
            Y_outputs = args.correct_output;
            [~, yPred] = max(Y_outputs);
            
            N = args.first_fc_layer;
            K = args.last_fc_layer;
            fc_layers = WPutils.verify_has_only_fc_relu_placeholder_layers(net, N, K);
            K = fc_layers(end);
            
            xN = net.input_sets{N}.V(:,:,:,1);
            xN = xN(:);     % input to perturbed layer
            
            W = net.Layers{K}.Weights;      % last layer's weight matrix
            
            image_verified = 1;
            for j = 1:length(Y_outputs)
                if j ~= yPred
                    t1 = norm(W(yPred,:) - W(j,:), 1);
                    if N == K
                        t1 = 2;
                    end
                    
                    t2 = norm(xN, 1);                   % input to the only perturbed layer
                    
                    t3 = 1;
                    for l = fc_layers(2:end - 1)   % all fully-connected layers between the perturbed layer and the last layer
                        t3 = t3*max(norms(net.Layers{l}.Weights, 1));
                    end
                    
                    ub_on_decrease_in_margin = p*t1*t2*t3;
                    original_margin = Y_outputs(yPred) - Y_outputs(j);
                    if ub_on_decrease_in_margin >= original_margin
                        image_verified = 0;
                        break
                    end
                end
            end
        end
        
        % for comparison with another approach in ModelStar experiments. Can be removed.
        function [lb, ub] = compute_bounds_weng_2020_fc_layers_only(net, args)
            % call this function after performing reachability analysis
            % using star sets.
            
            arguments
                net
                args.first_fc_layer
                args.last_fc_layer = []
                args.pert_vec = []
                args.pert_mag = []
                args.pert_norm = Inf
            end
            
            p = args.pert_norm;
            N = args.first_fc_layer;
            K = args.last_fc_layer;
            pert_vec = args.pert_vec;
            % N is the perturbed layer.
            % pert_vec is the perturbation applied to every row of
            % layer N's weight matrix.
            
            fc_layers = WPutils.verify_has_only_fc_relu_placeholder_layers(net, N, K);
            
            W = @(k) net.Layers{fc_layers(k)}.Weights;
            b = @(k) net.Layers{fc_layers(k)}.Bias;
            
            l = cell(length(fc_layers), 1);
            u = l;
            alphaL = l;
            betaL = l;
            alphaU = l;
            betaU = l;
            delta = l;
            theta = l;
            
            xN = net.input_sets{N}.V(:,:,:,1);
            xN = xN(:);
            if isempty(pert_vec)
                pert_vec = args.pert_mag*ones(1, length(xN));
            end
            yN = W(1)*xN + b(1);
            pert = pert_vec*abs(xN);
            l{1} = yN - pert;
            u{1} = yN + pert;
            
            e = norm(pert_vec, p);
            q = 1/(1 - 1/p);
            
            % compute bounds layer-wise i.e., 2nd layer, then 3rd layer, and so on
            for K = 2:length(fc_layers)
                numOutputs = size(W(K), 1);
                
                A = cell(length(fc_layers), 1);
                O = A;
                lambda = A;
                omega = A;
                
                m = K - 1;
                
                alphaU{m} = u{m}./(u{m} - l{m});
                alphaL{m} = alphaU{m};
                 betaU{m} = -l{m};
                
                alphaU{m}(u{m} <= 0) = 0;
                alphaL{m}(u{m} <= 0) = 0;
                 betaU{m}(u{m} <= 0) = 0;
                
                alphaU{m}(l{m} >= 0) = 1;
                alphaL{m}(l{m} >= 0) = 1;
                 betaU{m}(l{m} >= 0) = 0;
                
                u{K} = nan(numOutputs, 1);
                l{K} = u{K};
                
                for j = 1:numOutputs
                    ej = zeros(numOutputs, 1);
                    ej(j) = 1;
                    
                    A{K}(j,:) = ej';
                    O{K}(j,:) = ej';
                    
                    for k = K-1 : -1 : 1
                        A{k} = nan(numOutputs, size(W(k+1), 2));
                        O{k} = A{k};
                        lambda{k} = A{k};
                        omega{k} = A{k};
                        delta{k} = A{k}';
                        theta{k} = A{k}';
                        
                        Wkp1 = W(k+1);
                        for i = 1:size(Wkp1, 2)
                            if A{k+1}(j,:) * Wkp1(:,i) >= 0
                                lambda{k}(j,i) = alphaU{k}(i);
                                 delta{k}(i,j) =  betaU{k}(i);
                            else
                                lambda{k}(j,i) = alphaL{k}(i);
                                 delta{k}(i,j) = 0;     % betaL is always 0
                            end
                            if O{k+1}(j,:) * Wkp1(:,i) >= 0
                                 omega{k}(j,i) = alphaL{k}(i);
                                 theta{k}(i,j) = 0;     % betaL is always 0
                            else
                                 omega{k}(j,i) = alphaU{k}(i);
                                 theta{k}(i,j) =  betaU{k}(i);
                            end
                        end
                        A{k}(j,:) = (A{k+1}(j,:) * W(k+1)).*lambda{k}(j,:);
                        O{k}(j,:) = (O{k+1}(j,:) * W(k+1)).* omega{k}(j,:);
                    end
                    
                    t4 = A{K}(j,:) * b(K);    % delta{K} is always zero
                    t5 = O{K}(j,:) * b(K);    % theta{K} is always zero
                    for k = 1:K-1
                        t4 = t4 + A{k}(j,:) * (b(k) + delta{k}(:,j));
                        t5 = t5 + O{k}(j,:) * (b(k) + theta{k}(:,j));
                    end
                    
                    u{K}(j) =  e * norm(A{1}(j,:), 1) * norm(xN, q) + A{1}(j,:) * W(1)*xN + t4;
                    l{K}(j) = -e * norm(O{1}(j,:), 1) * norm(xN, q) + O{1}(j,:) * W(1)*xN + t5;
                end
            end
            
            lb = l{end};
            ub = u{end};
        end
        
        function fc_layers = verify_has_only_fc_relu_placeholder_layers(net, l_start, l_end)
            
            arguments
                net
                l_start 
                l_end = []
            end
            
            if isempty(l_end)
                fc_layers = net.get_layer_indices("FullyConnectedLayer");
                l_end = fc_layers(end);
            end
            
            fc_layers = [];
            relu_encountered = 0;
            for l = l_start:l_end
                if isa(net.Layers{l}, "FullyConnectedLayer")
                    fc_layers = [fc_layers l];
                    relu_encountered = 0;
                elseif isa(net.Layers{l}, "ReluLayer")
                    if relu_encountered
                        error("Two ReLU layers without FullyConnectedLayer between them not supported.")
                    end
                    relu_encountered = 1;
                elseif isa(net.Layers{l}, "PlaceholderLayer")
                else
                    error_struct.message = ['Unsupported layer ' class(net.Layers{l}) '.'];
                    error_struct.identifier = "NNV:non_supported_layers_in_fc_relu_only_method";
                    error(error_struct);
                end
            end
            
            if ~isa(net.Layers{l_start}, "FullyConnectedLayer")
                error(['First layer should be a fully connected layer.'])
            end
            if ~isa(net.Layers{l_end}, "FullyConnectedLayer")
                error(['Last layer should be a fully connected layer.'])
            end
            
        end
        
        % Verify robustness of a neural network given a single input
        function [robustness_result_by_intersection, lb, ub, compare_bounds, mismatching_bounds] = verify_robustness_for_3dim_img(net, reachOptions, args)
            
            arguments
                net
                reachOptions
                args.input
                args.target = []
                args.matlabnet = []
                args.correct_output = []
                args.plot = false
                args.compare_with_sampling = false
                args.tolerance_for_comparison_with_sampling_bounds = 1e-4
                args.samples_per_pert = 10
                args.netbs = []
                % args.numWorkers_for_par_sampl = net.sampling_workers
            end
            
            if isempty(args.matlabnet)
                args.matlabnet = net.matlabnet;
            end
            target = args.target;
            [~, pred] = max(predict(args.matlabnet, args.input));
            if ~isempty(target) && pred ~= target
                error('Predicted class different from supplied target class!')
            end
            target = pred;
            
            V = args.input;
            V(:,:,:,2) = zeros(size(args.input));
            
            C = 0;
            d = 0;
            pred_lb = -1;
            pred_ub = 1;
            IS = ImageStar(V, C, d, pred_lb, pred_ub);
            
            robustness_result_by_intersection = net.verify_robustness(IS, reachOptions, target);
            
            if nargout > 1 % && args.compare_with_sampling
                [lb, ub, robustness_result_by_overlap] = WPutils.get_ranges_and_plot(net, ...
                    correct_output = args.correct_output, ...
                    plot = args.plot, ...
                    target = target);
                
                if isfield(reachOptions, 'disp_overlap_result')
                    if robustness_result_by_overlap == 1
                        disp('Neural network is ROBUST by overlap.')
                    else
                        disp('Neural network robustness UNKNOWN by overlap.')
                    end
                end
                
                if args.compare_with_sampling
                    [compare_bounds, mismatching_bounds] = WPutils.sample_weight_perturbed_nns(net, ...
                                                    input = args.input, ...
                                                    samples_per_pert = args.samples_per_pert, ...
                                                    lb_out = lb, ...
                                                    ub_out = ub, ...
                                                    robustness_result = robustness_result_by_intersection, ...
                                                    matlabnet = args.matlabnet, ...
                                                    netbs = args.netbs, ...
                                                    tolerance_for_comparison_of_bounds = args.tolerance_for_comparison_with_sampling_bounds);
                                                    % , ...
                                                    % numWorkers_for_par_sampl = args.numWorkers_for_par_sampl);
                end
            end
        end
        
        % get output ranges of the reachable set(s) of the neural network
        function [lb, ub, robustness_result_by_overlap] = get_ranges_and_plot(net, args)
            
            arguments
                net
                args.correct_output = []
                args.plot = false
                args.target = []
                % args.x_labels = []
                args.lp_solver = net.lp_solver
            end
            
            robustness_result_by_overlap = nan;

            % Get output reachable set
            R = net.reachSet{end};
            for k = 1:length(R)
                R(k).changeDevice('cpu');
            end
            
            [lb, ub] = R(1).getRanges(args.lp_solver, 'parallel');
            lb = squeeze(lb);
            ub = squeeze(ub);
            if ~isempty(args.target)
                robustness_result_by_overlap = net.is_robust_by_overlap(lb, ub, args.target);
            end
            
            % obj.start_pool;
            progress = zeros(length(R), 1);
            progress(1) = 1;
            parfor k = 2:length(R)
                % Get (overapproximate) ranges for each output index
                [lb_out_k, ub_out_k] = R(k).getRanges(args.lp_solver);
                lb_out_k = squeeze(lb_out_k);
                ub_out_k = squeeze(ub_out_k);
                if ~isempty(args.target)
                    robustness_result_by_overlap = robustness_result_by_overlap & net.is_robust_by_overlap(lb_out_k, ub_out_k, args.target);
                end
                lb = min(lb, lb_out_k);
                ub = max(ub, ub_out_k);
                progress(k) = 1;
            end
            
            if args.plot
                % Get middle point for each output and range sizes
                mid_range = (lb + ub)/2;
                range_size = ub - mid_range;
                
                % Label for x-axis
                x = [1:length(lb)]';
                
                % Visualize set ranges and evaluation points
                hold off;
                errorbar(x, mid_range, range_size, '.');
                xlim([(x(1) - 0.5) (x(end) + 0.5)]);
                if ~isempty(args.correct_output)
                    hold on;
                    scatter(x, args.correct_output, 'x', 'MarkerEdgeColor', 'r');
                    hold off;
                end
            end
            
        end
        
        % sample a neural network's perturbation space by plugging the
        % perturbed weights into the network and testing the accuracy over
        % the supplied input. It is assumed that the supplied input is
        % classified correctly by the unperturbed neural network, and
        % predictions different from that of the unperturbed neural network
        % are considered misclassifications.
        function [compare_bounds, mismatching_bounds] = sample_weight_perturbed_nns(net, args)
            
            arguments
                net                             % the NN object
                args.input                      % the input to the NN
                args.lb_out                     % lower bounds on the
                    % output of the neural network obtained using NNV
                args.ub_out                     % upper bounds on the
                    % output of the neural network obtained using NNV
                args.robustness_result          % robustness result
                    % obtained using NNV
                args.samples_per_pert = 10      % samples per perturbation
                args.display = 1                % whether to display results
                args.netbs = []                 % bytestream based copy of 
                    % the NNV NN object
                % args.numWorkers_for_par_sampl = 5
                args.tolerance_for_comparison_of_bounds = 1e-4
            end
            
            samples_per_pert = args.samples_per_pert;    % no. of samples per weight
            
            % gather the weight perturbations from the various layers of
            % the neural network
            weightPerturbSpecs = [];
            for l = 1:length(net.Layers)
                if isa(net.Layers{l}, 'FullyConnectedLayer') || isa(net.Layers{l}, 'Conv2DLayer')
                    weightPerturbSpecsThisLayer = net.Layers{l}.weightPerturb;
                    sz = size(weightPerturbSpecsThisLayer);
                    if sz(2) > 0
                        weightPerturbSpecs = [weightPerturbSpecs; l*ones(sz(1),1) weightPerturbSpecsThisLayer];
                    end
                elseif isprop(net.Layers{l}, 'weightPerturb')
                    error("Layer " + net.Layers{l}.Name + " not supported for sampling yet.")
                end
            end
            
            numPerts = size(weightPerturbSpecs, 1);
            weightUB = weightPerturbSpecs(:, end);
            weightLB = weightPerturbSpecs(:, end - 1);
            deltasWeightPerturb = (weightUB - weightLB)/(samples_per_pert - 1);
                % deltas of weight perturbations
            
            powersNSamplesPerPert = samples_per_pert.^(numPerts - 1 :-1: 0).';   % powers of n
            
            correct_output = net.evaluate(args.input);
            correct_output = correct_output(:);
            [~, target] = max(correct_output);
            
            lb_netOut = ones(size(correct_output))*Inf;
            ub_netOut = -lb_netOut;
            
            incorrectly_classified_samples = 0;
            
            total_combs = samples_per_pert^numPerts;
            % old_numCores = obj.numCores;
            % obj.numCores = args.numWorkers_for_par_sampl;
            % obj.start_pool;
            parfor m = 1:total_combs      % perturbation combination number
            % for m = 1:total_combs      % perturbation combination number
                multDeltasWeightPerturb = mod(ceil(m./powersNSamplesPerPert) - 1, samples_per_pert);  % multiples of deltas of weight perturbations
                weightPerturb = multDeltasWeightPerturb.*deltasWeightPerturb + weightLB;   % weight perturbations
                
                tmp_net = getArrayFromByteStream(args.netbs);
                
                for k = 1 : numPerts
                    l = weightPerturbSpecs(k, 1);
                    ind = weightPerturbSpecs(k, 2);
                    if isa(net.Layers{l}, 'FullyConnectedLayer') || isa(net.Layers{l}, 'Conv2DLayer')
                        n = numel(net.Layers{l}.Weights);
                        if ind <= n     % perturb weights matrix
                            tmp_net.Layers{l}.Weights(ind) = tmp_net.Layers{l}.Weights(ind) + weightPerturb(k);
                        else            % perturb bias matrix
                            ind = ind - n;
                            tmp_net.Layers{l}.Bias(ind) = tmp_net.Layers{l}.Bias(ind) + weightPerturb(k);
                        end
                    else
                        error(['Layer ' class(net.Layers{l}) ' not supported in sampling!'])
                    end
                end
                
                netOut = tmp_net.evaluate(args.input);
                netOut = netOut(:);
                [~, predThisComb] = max(netOut);
                if predThisComb ~= target
                    incorrectly_classified_samples = incorrectly_classified_samples + 1;
                end
                
                lb_netOut = min(lb_netOut, netOut);
                ub_netOut = max(ub_netOut, netOut);
            end
            
            % debug_disp
            lb_out = args.lb_out;
            ub_out = args.ub_out;
            compare_bounds = [lb_netOut, ub_netOut, lb_netOut - lb_out, ub_out - ub_netOut, lb_out, ub_out];
                % bounds from sampling, overapproximation of bounds
                % from star sets over bounds from sampling, bounds from
                % star sets
            
            percentage_of_incorrect_classification = incorrectly_classified_samples/total_combs*100;
            fprintf('\n\n------------------------- Comparison with Sampling -------------------------\n\n')
            disp(['Samples incorrectly classified: ' num2str(percentage_of_incorrect_classification) ' %'])
            
            if any(lb_netOut < 0, 'all')
                shift = abs(min(lb_netOut, [], 'all'));
            else
                shift = 0;
            end
            shifted_lb_netOut = lb_netOut + shift;
            shifted_ub_netOut = ub_netOut + shift;
            
            zero_thresh = args.tolerance_for_comparison_of_bounds*shifted_lb_netOut(target);
            
            matching_lbs = abs(compare_bounds(:,3)) < zero_thresh;
            matching_ubs = abs(compare_bounds(:,4)) < zero_thresh;
            matching_bounds = matching_lbs & matching_ubs;
            fprintf('Tolerance: %g \t ', args.tolerance_for_comparison_of_bounds)
            mismatching_bounds = [];
            if all(matching_bounds)
                fprintf('Bounds match.\n')
            else
                fprintf('Bounds MISMATCH! Mismatching bounds:\n\n')
                mismatching_bounds = [find(~matching_bounds), compare_bounds(~matching_bounds, :)];
                disp(mismatching_bounds)
            end
            
            if (incorrectly_classified_samples > 0) && (args.robustness_result == 1)
                % one or more samples resulted in incorrect classification, even though
                % neural network was determined to be robust by nnv
                error('Robustness mismatch!')
            end
            
            % net.numCores = old_numCores;
        end
        
        function the_range = get_weights_range(net, layer_no)
            the_range = range(net.Layers{layer_no}.Weights, 'all');
        end
        
        % specify a weight perturbation to be added to a layer's
        % weightPerturb specification matrix
        function add_pert_call_this_function_from_layer(layer, pert_indices, lb, ub)
            % @layer: layer to specify perturbation for
            % @pert_indices: indices of the weight to be perturbed
            % @lb: lower bound of the perturbation
            % @ub: upper bound of the perturbation
            
            sz = size(layer.Weights);
            if length(pert_indices) ~= length(sz)
                error(['Must specify location of perturbation in layer ' layer.Name ' as indices of its weights matrix which has size ' sz])
            end
            
            pert_indices = num2cell(pert_indices);
            linear_ind = sub2ind(sz, pert_indices{:});
            
            layer.weightPerturb = [layer.weightPerturb; linear_ind lb ub];
        end
        
        % specify a perturbation to be applied to the whole layer
        function perturb_whole_layer_call_this_function_from_layer(layer, lb, ub)
            % @lb: lower bound of the perturbation
            % @ub: upper bound of the perturbation
            
            n = numel(layer.Weights);
            
            temp_pert = nan(n, 3);
            temp_pert(:, 1) = [1:n]';
            temp_pert(:, 2) = lb;
            temp_pert(:, 3) = ub;
            
            layer.weightPerturb = [layer.weightPerturb; temp_pert];
        end
        
        % specify a symmetric perturbation to be applied to the whole layer
        % as a fraction of the range of the weights in the layer
        function pert_whole_layer_given_fraction_of_weights_range_call_fromLayer(layer, frac)
            % @frac: fraction of the weights' range to be applied as
            % perturbation to the whole layer
            
            p = frac*range(layer.Weights, 'all');
            NN.perturb_whole_layer_call_this_function_from_layer(layer, -p, p);
        end
        
        function mem = get_free_mem_B()
            [~, mem_str] = system("vmstat -s -S K | grep free\ memory | awk '{print $1}'");
            mem = str2double(mem_str)*2^10;
        end
        
        function idle_cpu = get_idle_cpu()
            [~, idle_cpu] = system("mpstat 1 1 | sed -n '4p' | awk '{print $NF}'");
            idle_cpu = str2double(idle_cpu)/100;
        end
        
        function mem = get_total_mem_B()
            [~, mem_str] = system("vmstat -s -S K | grep total\ memory | awk '{print $1}'");
            mem = str2double(mem_str)*2^10;
        end
        
        function frac = get_free_mem_frac()
            frac = NN.get_free_mem_B/NN.get_total_mem_B;
        end
        
        function mem = get_free_swap_B()
            [~, mem_str] = system("vmstat -s -S K | grep free\ swap | awk '{print $1}'");
            mem = str2double(mem_str)*2^10;
        end
        
        function mem_B = get_required_mem_B(sz, dataType)
            % sz is the size vector of the matrix
            
            n = prod(sz);
            
            switch dataType
                case 'double'
                    bytesPerElement = 8;
                case 'single'
                    bytesPerElement = 4;
                case 'int32'
                    bytesPerElement = 4;
                case 'int16'
                    bytesPerElement = 2;
                case 'int8'
                    bytesPerElement = 1;
                otherwise
                    error('Unsupported data type');
            end
            
            % Calculate the total memory required in bytes
            mem_B = n*bytesPerElement;
        end
        
        function [frac_used, free_mem_B] = get_nvidia_gpu_used_memory_frac()
            [~, gpu_str] = system("nvidia-smi | sed -n '10p'");
            gpu_str = strsplit(gpu_str);
            used_mem_MB = gpu_str{9};
            MiB_loc = strfind(used_mem_MB, "MiB");
            if isempty(MiB_loc)
                error("Check processing of output of nvidia-smi")
            end
            used_mem_B = str2double(used_mem_MB(1:MiB_loc - 1))*2^20;
            total_mem_MB = gpu_str{11};
            MiB_loc = strfind(total_mem_MB, "MiB");
            if isempty(MiB_loc)
                error("Check processing of output of nvidia-smi")
            end
            total_mem_B = str2double(total_mem_MB(1:MiB_loc - 1))*2^20;
            free_mem_B = total_mem_B - used_mem_B;
            frac_used = used_mem_B/total_mem_B;
        end
        
        % for ModelStar experiments. Can be removed.
        function net = acas_combine_affine_layer_with_fc(net, args)
            
            arguments
                net
                args.delete_affine_layers = 0
            end
            
            wrong_network = 0;
            % wrong_network = wrong_network || ~(isempty(obj.reachSet) && isempty(obj.input_sets));
            
            wrong_network = wrong_network || ...
                            ~isa(net.Layers{1},  'ImageInputLayer') || ~strcmp(net.Layers{1}.Normalization, 'none') || ...
                            ~isa(net.Layers{2}, 'PlaceholderLayer') || ...
                                 net.Layers{3}.DoOffset || ...
                            ~isa(net.Layers{4},     'FlattenLayer');
            
            acas_fc_layers = 5:3:23;
            wrong_network = wrong_network || norm(find(arrayfun(@(l) isa(l{1}, 'FullyConnectedLayer'), net.Layers)).' - acas_fc_layers, 1) ~= 0;
            
            acas_elementwise_affine_layers = 3:3:24;
            wrong_network = wrong_network || norm(find(arrayfun(@(l) isa(l{1}, 'ElementwiseAffineLayer'), net.Layers)).' - acas_elementwise_affine_layers, 1) ~= 0;
            wrong_network = wrong_network || any(arrayfun(@(l) isa(l{1}, 'ElementwiseAffineLayer') && l{1}.DoScale, net.Layers));
            
            acas_relu_layers = 7:3:22;
            wrong_network = wrong_network || norm(find(arrayfun(@(l) isa(l{1}, 'ReluLayer'), net.Layers)).' - acas_relu_layers, 1) ~= 0;
            
            wrong_network = wrong_network || ~isa(net.Layers{25}, 'PlaceholderLayer') || length(net.Layers) > 25;
            
            if wrong_network
                error("Network did not satisfy the expected ACAS network structure.")
            end
            
            for fc_layer = acas_fc_layers
                if (fc_layer < acas_fc_layers(end) && any(size(net.Layers{fc_layer + 1}.Offset) ~= [50 1])) || (fc_layer == acas_fc_layers(end) && any(size(net.Layers{fc_layer + 1}.Offset) ~= [5 1]))
                    error("Wrong size of offset.")
                end
                net.Layers{fc_layer}.Bias = net.Layers{fc_layer}.Bias + net.Layers{fc_layer + 1}.Offset;
                net.Layers{fc_layer + 1}.Offset = zeros(size(net.Layers{fc_layer + 1}.Offset));
                net.Layers{fc_layer + 1}.DoOffset = 0;
            end
            
            if any(arrayfun(@(l) isa(l{1}, 'ElementwiseAffineLayer') && (l{1}.DoScale || l{1}.DoOffset), net.Layers))
                error("There is some ElementwiseAffineLayer doing some scaling or offsetting.")
            end
            
            if args.delete_affine_layers
                if ~isempty(net.reachSet) || ~isempty(net.input_sets)
                    if any(arrayfun(@(l) net.input_sets{l} ~= net.reachSet{l}, acas_elementwise_affine_layers))
                        error("ElementwiseAffineLayer's input and output should have been identical i.e., it should have been combined with the preceding fully connected layer before reachability analysis.")
                    end
                end
                net.Layers(acas_elementwise_affine_layers) = [];
                if ~isempty(net.reachSet) || ~isempty(net.input_sets)
                    net.reachSet(acas_elementwise_affine_layers) = [];
                    net.input_sets(acas_elementwise_affine_layers) = [];
                end
            end
            
            net = net;
            
        end
        
    end
end