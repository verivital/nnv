classdef FullyConnectedLayer < handle
    % The FullyConnectedLayer layer class in CNN
    %   Contain constructor and reachability analysis methods
    % Main references:
    % 1) An intuitive explanation of convolutional neural networks: 
    %    https://ujjwalkarn.me/2016/08/11/intuitive-explanation-convnets/
    % 2) More detail about mathematical background of CNN
    %    http://cs231n.github.io/convolutional-networks/
    %    http://cs231n.github.io/convolutional-networks/#pool
    % 3) Matlab implementation of Convolution2DLayer and MaxPooling (for training and evaluating purpose)
    %    https://www.mathworks.com/help/deeplearning/ug/layers-of-a-convolutional-neural-network.html
    %    https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.fullyconnectedlayer.html
    
    %   Dung Tran: 6/26/2019
    
    properties
        Name = 'fully_connected_layer';
        % Hyperparameters
        InputSize = 0;  % number of input
        OutputSize = 0; % number of output
        Weights = []; % weight matrix
        Bias  = []; % bias vector     
        weightPerturb = [];   % perturb specified weights of this fully
            % connected layer in the specified range. The first element of
            % each row is the index of the weight to be perturbed in linear
            % format, so use sub2ind() on multi-dimensional indices before
            % specifying it here. The next two elements are the lower and
            % upper bounds of the perturbation.
            % If the index exceeds the number of elements in the weights
            % matrix, the perturbation is considered to be in the bias,
            % and the index of the perturbation in the bias is obtained by
            % substracting the number of elements in the weights matrix
            % from the supplied index.
    end
    
    
    methods % main methods
        
        % constructor of the class
        function obj = FullyConnectedLayer(varargin)
            % author: Dung Tran
            % date: 6/26/2019    
            % update: 
            
            switch nargin
                
                case 3
                    name = varargin{1};
                    W = varargin{2};
                    b = varargin{3};
                    if ~ischar(name)
                        error('Name is not char');
                    else
                        obj.Name = name;
                    end
                    if size(W, 1) ~= size(b, 1)
                        error('Inconsistent dimension between the weight matrix and bias vector');
                    end
                    if size(b,2) ~= 1
                        error('Bias vector should have one column');
                    end
                    
                    obj.InputSize = size(W,2);
                    obj.OutputSize = size(W,1);
                    obj.Weights = W;
                    obj.Bias = b;
                    
                case 2
                    
                    W = varargin{1};
                    b = varargin{2};
                    obj.Name = 'fully_connected_layer';
                    if size(W, 1) ~= size(b, 1)
                        error('Inconsistent dimension between the weight matrix and bias vector');
                    end
                    if size(b,2) ~= 1
                        error('Bias vector should have one column');
                    end
                    obj.InputSize = size(W,2);
                    obj.OutputSize = size(W,1);
                    obj.Weights = W;
                    obj.Bias = b;
                           
                case 0
                    
                    obj.Name = 'fully_connected_layer';
                    % Hyperparameters
                    obj.InputSize = 0;
                    obj.OutputSize = 0; % step size for traversing input
                    obj.Weights = [];
                    obj.Bias = [];
                    
                otherwise
                    error('Invalid number of inputs (should be 0, 2 or 3)');
            end 
             
        end
        
        % evaluation method
        function y = evaluate(obj, x)
            % @input: input
            % @y: output
            
            % author: Dung Tran
            % date: 6/26/2019
            % 
            % update: add support for sequence evaluation (neuralode, RNN)
            %   -date: 03/17/2023 (Diego Manzanas)
            
            n = size(x);
            if length(n) == 2
                x = reshape(x, [n(1)*n(2) 1]);
            elseif length(n) == 3
                x = reshape(x, [n(1)*n(2)*n(3) 1]);
            elseif length(n) == 4
                x = reshape(x, [n(1)*n(2)*n(3)*n(4) 1]);
            else
                error('Invalid input');
            end

            if size(x, 1) ~= size(obj.Weights, 2)
                error('Inconsistent input vector')
            end

            y1 = obj.Weights * x;
            if size(x, 2) ~= 1
                n = size(x, 2);
                for i=1:n
                    y1(:,i) = y1(:,i) + obj.Bias;
                end
                y = y1;
            else
                y = y1 + obj.Bias;
            end   
             
        end 
        
        function y = evaluateSequence(obj, input)
            % @input: input
            % @y: output
            % author: Neelanjana Pal
            % date: 1/6/2023

            y = obj.Weights*input + obj.Bias;
        end

        % main reachability analysis function
        function IS = reach(varargin)
            % @in_image: an input imagestar
            % @image: output set
            % @option: = 'single' or 'parallel' 
            
            % author: Dung Tran
            % date: 6/26/2019
            % update: 1/6/2020  update reason: add zonotope method
             
            switch nargin
                
                 case 7
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    % relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                    % lp_solver = varargin{7}; do not use
                
                case 6
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                
                case 4
                    obj = varargin{1};
                    in_images = varargin{2}; 
                    method = varargin{3};
                    option = varargin{4}; % computation option

                case 3
                    obj = varargin{1};
                    in_images = varargin{2}; % don't care the rest inputs
                    method = varargin{3};
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 2, 3, 4, 5, or 6)');
            end
            
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || contains(method, "relax-star")
                IS = obj.reach_star_multipleInputs(in_images, option);
            elseif strcmp(method, 'approx-zono')
                IS = obj.reach_zono_multipleInputs(in_images, option);
            else
                error('Unknown reachability method');
            end
            
        end
        
        function IS = reachSequence(varargin)
            % @in_image: an input imagestar
            % @image: output set
            % @option: = 'single' or 'parallel' 
             
            switch nargin
                
                 case 7
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                
                case 6
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                
                case 4
                    obj = varargin{1};
                    in_images = varargin{2}; 
                    method = varargin{3};
                    option = varargin{4}; % computation option

                case 3
                    obj = varargin{1};
                    in_images = varargin{2}; % don't care the rest inputs
                    method = varargin{3};
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 2, 3, 4, 5, or 6)');
            end
            
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || contains(method, "relax-star")
                IS = obj.reachSequence_star(in_images, option);
            elseif strcmp(method, 'approx-zono')
                IS = obj.reachSequence_zono(in_images, option);
            else
                error('Unknown reachability method');
            end
        end

        % change params to gpuArrays
        function obj = toGPU(obj)
            obj.Weights = gpuArray(obj.Weights);
            obj.Bias = gpuArray(obj.Bias);
        end

        % change params to cpuArrays
        function obj = toCPU(obj)
            obj.Weights = gather(obj.Weights);
            obj.Bias = gather(obj.Bias);
        end

        % Change params precision
        function obj = changeParamsPrecision(obj, precision)
            if strcmp(precision, "double")
                obj.Weights = double(obj.Weights);
                obj.Bias = double(obj.Bias);
            elseif strcmp(precision, "single")
                obj.Weights = single(obj.Weights);
                obj.Bias = single(obj.Bias);
            else
                error("Parameter numerical precision must be 'single' or 'double'");
            end
        end
    
    end   
     
    methods % reachability methods
        
        % reachability analysis using imagestar
        function image = reach_star_single_input(obj, in_image)
            % @in_image: input imagestar
            % @image: output set
            
            if ~isa(in_image, 'ImageStar') && ~isa(in_image, 'Star')
                error('Input set is not an ImageStar or Star');
            end
            
            if isa(in_image, 'ImageStar')
                % reach using ImageStar
                N = in_image.height*in_image.width*in_image.numChannel;
                if N~= obj.InputSize
                    error('Inconsistency between the size of the input image and the InputSize of the network');
                end
                           
                n = in_image.numPred;
                V(1, 1, :, :) = obj.Weights*reshape(in_image.V, N, n + 1);
                V(1, 1, :, 1) = reshape(V(1, 1, :, 1), obj.OutputSize, 1) + obj.Bias;
                C = in_image.C;
                d = in_image.d;
                pred_lb = in_image.pred_lb;
                pred_ub = in_image.pred_ub;
                
                if ~isempty(obj.weightPerturb)

                    % make sure that the weight perturbation matrix adheres
                    % to the 3 column format described above
                    sz = size(obj.weightPerturb);
                    if length(sz) > 2 || sz(2) > 3
                        error("Weight perturbation specification matrix must have only 3 elements in each row: [linear index of perturbed weight, lower bound of perturbation, upper bound of perturbation]")
                    end
                    
                    % split perturbation ranges that have positive as
                    % well as negative parts into only positive or only
                    % negative parts
                    
                    % if signs of upper and lower bounds are different
                    pert_with_pos_and_neg_range = (obj.weightPerturb(:, end - 1) < 0) & (obj.weightPerturb(:, end) > 0);

                    % copy these perturbations' entries, then set ub to 0
                    % in original entries and lb to 0 in copies
                    % (may use this non-positiveness or non-negativeness
                    % of a given predicate in determining constraints on
                    % products of predicates)
                    pert_copies = obj.weightPerturb(pert_with_pos_and_neg_range, :);
                    obj.weightPerturb(pert_with_pos_and_neg_range, end) = 0;
                    pert_copies(:, end - 1) = 0;
                    obj.weightPerturb = [obj.weightPerturb; pert_copies];
                    
                    % Combine weight perturbations that result in 
                    % perturbation in same neuron into a single
                    % perturbation.

                    % make a framework to hold these.
                    dim_wise_pred_lb_ub = zeros(size(obj.Weights, 1), 2);
                        % if both bounds are zero, no predicate along that
                        % dimension.

                    % break only those existing generators whose 
                    % predicate variables have no constraints on them,
                    % and accumulate these in framework.
                    preds_to_remove = [];
                    for pred_no = 1:size(V,4) - 1
                        if C(:, pred_no) == 0
                            for dim = 1:size(V,3)
                                if V(1, 1, dim, pred_no + 1) > 0
                                    dim_wise_pred_lb_ub(dim, :) = dim_wise_pred_lb_ub(dim, :) + V(1, 1, dim, pred_no + 1)*[pred_lb(pred_no) pred_ub(pred_no)];
                                elseif V(1, 1, dim, pred_no + 1) < 0
                                    dim_wise_pred_lb_ub(dim, :) = dim_wise_pred_lb_ub(dim, :) + V(1, 1, dim, pred_no + 1)*[pred_ub(pred_no) pred_lb(pred_no)];
                                end
                            end
                            preds_to_remove = [preds_to_remove, pred_no];
                        else
                            % error('Constraints encountered! not an error')
                        end
                    end

                    % remove the vs that have been added to framework,
                    % from V, thus completing transfer from ImageStar to
                    % framework
                    V(:, :, :, preds_to_remove + 1) = [];
                    C(:, preds_to_remove) = [];
                    pred_lb(preds_to_remove) = [];
                    pred_ub(preds_to_remove) = [];
                    
                    new_C = C;
                    new_d = d;
    
                    % accumulate the new vs in framework
                    in_vector = reshape(in_image.V(:, :, :, 1), N, 1);
                    for k = 1:size(obj.weightPerturb, 1)
                        weights_sz = size(obj.Weights);
                        ind = obj.weightPerturb(k, 1);
                        [row, col] = ind2sub(weights_sz, ind);
                        lb  = obj.weightPerturb(k, 2);
                        ub  = obj.weightPerturb(k, 3);
    
                        if col <= weights_sz(2)
                            
                            % for first term, without predicate products
                            mul_factor = in_vector(col);
                                % value of element of input vector being
                                % multiplied by the weight
                            
                            % multiplication with existing predicates
                            for pred_no = 1:size(in_image.V, 4) - 1
                                v = reshape(in_image.V(:, :, :, pred_no + 1), N, 1);
                                
                                if v(col) ~= 0
        
                                    % get limits of product of predicates
                                    lb1 = lb;
                                    ub1 = ub;
                                    lb2 = in_image.pred_lb(pred_no);
                                    ub2 = in_image.pred_ub(pred_no);
                                    
                                    b_prods = [lb1*lb2 lb1*ub2 ub1*lb2 ub1*ub2];
                                    pred_prod_lb = min(b_prods);
                                    pred_prod_ub = max(b_prods);
                                    
                                    % ignore the new
                                    % constraints on the products; just
                                    % accumulate the products of
                                    % predicates in the "framework".
    
                                    % integrate v's element into
                                    % limits and accumulate into 
                                    % framework.
                                    if v(col) > 0
                                        dim_wise_pred_lb_ub(row, :) = dim_wise_pred_lb_ub(row, :) + v(col)*[pred_prod_lb pred_prod_ub];
                                    else  % v(col) <= 0
                                        dim_wise_pred_lb_ub(row, :) = dim_wise_pred_lb_ub(row, :) + v(col)*[pred_prod_ub pred_prod_lb];
                                    end
                                end
                            end
                        else
                            % perturbation is in bias vector, which is
                            % simple addition, so the multiplicative
                            % factor is 1
                            mul_factor = 1;
                        end
    
                        % for first term, without predicate products
                        if mul_factor > 0
                            dim_wise_pred_lb_ub(row, :) = dim_wise_pred_lb_ub(row, :) + mul_factor*[lb ub];
                        elseif mul_factor < 0
                            dim_wise_pred_lb_ub(row, :) = dim_wise_pred_lb_ub(row, :) + mul_factor*[ub lb];
                        end
                    end

                    C = new_C;
                    d = new_d;

                    % transfer from framework to imagestar
                    V_cols = eye(size(dim_wise_pred_lb_ub, 1));
                    dims_to_skip = find((dim_wise_pred_lb_ub(:, 1) == 0) & (dim_wise_pred_lb_ub(:, 2) == 0));
                    V_cols(:, dims_to_skip) = [];
                    V(1, 1, :, size(V, 4) + (1:size(V_cols, 2))) = V_cols;

                    dim_wise_pred_lb_ub(dims_to_skip, :) = [];
                    C = [C zeros(size(C, 1), size(dim_wise_pred_lb_ub, 1))];
                    if isempty(pred_lb) && isempty(pred_ub)
                        dim_wise_pred_lb_ub = cast(dim_wise_pred_lb_ub, 'like', pred_lb);
                        pred_lb = [];
                        pred_ub = [];
                    end
                    pred_lb = [pred_lb; dim_wise_pred_lb_ub(:, 1)];
                    pred_ub = [pred_ub; dim_wise_pred_lb_ub(:, 2)];

                end
                
                if isempty(C)
                    V(1, 1, :, size(V, 4) + 1) = zeros(size(V, 3), 1);
                    C = 0;
                    pred_lb = 0;
                    pred_ub = 0;
                end

                % output set
                image = ImageStar(V, C, d, pred_lb, pred_ub);
            else % reach Star set
                image = in_image.affineMap(obj.Weights, obj.Bias);
            end
            
        end
        
        % handle multiple inputs
        function S = reach_star_multipleInputs(obj, inputs, option)
            % @inputs: an array of ImageStars
            % @option: = 'parallel' or 'single'
            % @S: output ImageStar
            
            n = length(inputs);
            if isa(inputs, "ImageStar")
                S(n) = ImageStar;
            elseif isa(inputs, "Star")
                S(n) = Star;
            else
                error("Input must be ImageStar or Star");
            end
            
            if strcmp(option, 'parallel')
                parfor i=1:n
                    S(i) = obj.reach_star_single_input(inputs(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    S(i) = obj.reach_star_single_input(inputs(i));
                end
            else
                error('Unknown computation option, should be parallel or single');
            end
            
        end
        
        %(reachability analysis using imagezono)
        function image = reach_zono(obj, in_image)
            % @in_image: input imagezono
            % @image: output set

            if ~isa(in_image, 'ImageZono') && ~isa(in_image, 'Zono')
                error('Input set is not an ImageZono or Zono');
            end

            if isa(in_image, 'ImageZono')
            
                N = in_image.height*in_image.width*in_image.numChannels;
                if N~= obj.InputSize
                    error('Inconsistency between the size of the input image and the InputSize of the network');
                end
                           
                n = in_image.numPreds;
                V(1, 1, :, in_image.numPreds + 1) = zeros(obj.OutputSize, 1);        
                for i=1:n+1
                    I = in_image.V(:,:,:,i);
                    I = reshape(I,N,1);
                    if i==1
                        V(1, 1,:,i) = double(obj.Weights)*I + double(obj.Bias);
                    else
                        V(1, 1,:,i) = double(obj.Weights)*I;
                    end
                end
                
                image = ImageZono(V);
            else
                image = in_image.affineMap(obj.Weights, obj.Bias);
            end
            
        end
        
        % handle mulitples inputs
        function S = reach_zono_multipleInputs(obj, inputs, option)
            % @inputs: an array of ImageZonos
            % @option: = 'parallel' or 'single'
            % @S: output ImageZono
            
            n = length(inputs);
            if isa(inputs, 'ImageZono')
                S(n) = ImageZono;
            elseif isa(inputs, 'Zono')
                S(n) = Zono;
            else
                error('Wrong input set. It must be ImageZono or Zono')
            end

            if strcmp(option, 'parallel')
                parfor i=1:n
                    S(i) = obj.reach_zono(inputs(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    S(i) = obj.reach_zono(inputs(i));
                end
            else
                error('Unknown computation option, should be parallel or single');
            end
            
        end
        
        function image = reachSequence_star(obj, in_image, ~)
            % @in_image: input imagestar
            % @image: output set            
            
            if ~isa(in_image, 'ImageStar')
                error('Input set is not an ImageStar');
            end
                       
            n = in_image.numPred;
            V(:, :, 1, in_image.numPred + 1) = zeros(obj.OutputSize, in_image.width);
            for i=1:n+1
                I = in_image.V(:,:,:,i);
                if i==1
                    V(:, :,1,i) = double(obj.Weights)*I + double(obj.Bias);
                else
                    V(:, :,1,i) = double(obj.Weights)*I;
                end
            end
            
            image = ImageStar(V, in_image.C, in_image.d, in_image.pred_lb, in_image.pred_ub);
            
        end
        
        % specify a weight perturbation to be added to the 
        % obj.weightPerturb specification matrix
        function add_pert(obj, pert_indices, lb, ub)
            % @pert_indices: indices of the weight to be perturbed
            % @lb: lower bound of the perturbation
            % @ub: upper bound of the perturbation
            
            NN.add_pert_call_this_function_from_layer(obj, pert_indices, lb, ub);
        end
        
        % specify a perturbation to be applied to the whole layer
        function perturb_whole_layer(obj, lb, ub)
            % @lb: lower bound of the perturbation
            % @ub: upper bound of the perturbation
            
            NN.perturb_whole_layer_call_this_function_from_layer(obj, lb, ub);
        end
        
        % specify a symmetric perturbation to be applied to the whole layer
        % as a fraction of the range of the weights in the layer
        function perturb_whole_layer_given_fraction_of_weights_range(obj, frac)
            % @frac: fraction of the weights' range to be applied as
            % perturbation to the whole layer
            
            NN.pert_whole_layer_given_fraction_of_weights_range_call_fromLayer(obj, frac);
        end
    
    end
    
    
    methods(Static)
        
        % parse a trained FullyConnectedLayer from matlab
        function L = parse(fully_connected_layer)
            % @fully_connecteted_Layer: a fully connected layer from matlab deep
            % neural network tool box
            % @L : a FullyConnectedLayer obj for reachability analysis purpose
            
            if ~isa(fully_connected_layer, 'nnet.cnn.layer.FullyConnectedLayer')
                error('Input is not a Matlab nnet.cnn.layer.FullyConnectedLayer class');
            end
            L = FullyConnectedLayer(fully_connected_layer.Name, fully_connected_layer.Weights, fully_connected_layer.Bias);         
            
        end
        
    end
    
end

