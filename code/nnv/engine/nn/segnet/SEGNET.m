classdef SEGNET < handle
    % Segmentation neural network class 
    % Author: Dung Tran
    % Date: 4/14/2020
    
    properties
        
        Name = 'segnet'; % name of the network
        Layers = {}; % An array of Layers, eg, Layers = [L1 L2 ...Ln]
        Connections = []; % table of connection
        
        InputNames = [];
        OutputNames = [];
        
        numLayers = 0; % number of Layers
        InputSize = 0; % number of Inputs
        OutputSize = 0; % number of Outputs
        
        % properties for reachability analysis        
        reachMethod = 'approx-star';    % reachable set computation scheme, default - 'star'
        reachOption = []; % parallel option, default - non-parallel computing
        numCores = 0; % number of cores (workers) using in computation
        reachSet = [];  % reachable set for each layers
        reachTime = []; % computation time for each layers
        totalReachTime = 0; % total computation time
        
        verifiedOutputSet = []; % verified output reach set
        groundTruthSegIms = []; % ground truth segmentation images
        % used for plotting verified output set
        numClasses = 0;
        
        % used for plot robustness information
        rb_val = []; % robustness value in percentage
        numMisPixels = []; % number of missclassified pixels
        numRbPixels = [];  % number of correctly classified pixels
        numPixels = []; % total number of pixels
    end
    
    
    methods % constructor, evaluation, sampling, print methods
        
        % constructor
        function obj = SEGNET(varargin)
            
            switch nargin
                case 5
                    name = varargin{1};
                    layers = varargin{2};
                    connections = varargin{3};
                    inputNames = varargin{4};
                    outputNames = varargin{5};
                                        
                    % check input arguments
                    if ~ischar(name)
                        error('Invalid type of network name, should be char');
                    end
                    
                    if ~ismatrix(layers)
                        error('Layers should be an array of layers');
                    end
                    
                    if ~isempty(connections) && ~istable(connections)
                        error('Connection should be a table');
                    end
                    
                    if ~iscell(inputNames)
                        error('InputNames should be a cell');
                    end
                    
                    if ~iscell(outputNames)
                        error('OutputNames should be a cell');
                    end
                    
                    
                    nL = length(layers); % number of Layers                                       
                    obj.Name = name;
                    obj.Layers = layers;
                    obj.numLayers = nL;    % number of layers
                    obj.Connections = connections;
                    obj.InputNames = inputNames;
                    obj.OutputNames = outputNames;
                    if isprop(layers{1}, 'InputSize')
                        obj.InputSize = layers{1}.InputSize; % input size
                    end
                    if isprop(layers{nL}, 'OutputSize')
                        obj.OutputSize = layers{nL}.OutputSize; % output size
                    end
                
                case 3
                    name = varargin{1};
                    layers = varargin{2};
                    connections = varargin{3};
                    % check input arguments
                    if ~ischar(name)
                        error('Invalid type of network name, should be char');
                    end
                    
                    if ~ismatrix(layers)
                        error('Layers should be an array of layers');
                    end
                    
                    if ~istable(connections)
                        error('Connection should be a table');
                    end
                                       
                    
                    nL = length(layers); % number of Layers                                       
                    obj.Name = name;
                    obj.Layers = layers;
                    obj.numLayers = nL;    % number of layers
                    obj.Connections = connections;
                    if isprop(layers(1), 'InputSize')
                        obj.InputSize = layers(1).InputSize; % input size
                    end
                    if isprop(layers(nL), 'OutputSize')
                        obj.OutputSize = layers(nL).OutputSize; % output size
                    end
                
                case 2
                    
                    layers = varargin{1};
                    connections = varargin{2};
                    % check input arguments
                    if ~ismatrix(layers)
                        error('Layers should be an array of layers');
                    end
                    
                    if ~istable(connections)
                        error('Connection should be a table');
                    end
                                       
                    
                    nL = length(layers); % number of Layers                                       
                    obj.Layers = layers;
                    obj.numLayers = nL;    % number of layers
                    obj.Connections = connections;
                    if isprop(layers(1), 'InputSize')
                        obj.InputSize = layers(1).InputSize; % input size
                    end
                    if isprop(layers(nL), 'OutputSize')
                        obj.OutputSize = layers(nL).OutputSize; % output size
                    end
                
                    
                case 0
                    
                    obj.Name = 'emptySEGNET';
                    obj.Layers = [];
                    obj.numLayers = 0;    % number of layers
                    obj.Connections = [];
                    obj.InputNames = {};
                    obj.OutputNames = {};
                    
                otherwise
                    error('Invalid number of inputs, should be 0, 2, 3 or 5');
            end
            
            obj.numClasses = length(obj.Layers{obj.numLayers}.Classes); % we introduces 2 more classes for analysis, unknown and misclass
            
                      
        end
                
        
        % Evaluation of a CNN
        function y = evaluate(varargin)
            % Evaluation of this FFNN
            % @x: input vector x
            % @y: output vector y
            % @layer_id: layer index
            
            switch nargin
                case 3
                    obj = varargin{1};
                    x = varargin{2};
                    layer_id = varargin{3};
                case 2
                    obj = varargin{1};
                    x = varargin{2};
                    layer_id = obj.numLayers;
                otherwise
                    error("Invalid number of input arguments, should be 1 or 2");
            end
            
            if layer_id < 1
                error('Invalid layer index');
            end
            
            y = x;
            for i=1:layer_id
                
                if ~isa(obj.Layers{i}, 'MaxUnpooling2DLayer')
                    if isa(obj.Layers{i}, 'MaxPooling2DLayer')
                        
                        if isempty(obj.Connections)
                            y = obj.Layers{i}.evaluate(y);
                        else
                            y = obj.Layers{i}.evaluate(y, 'segnet');
                        end
                        
                    else
                        y = obj.Layers{i}.evaluate(y);
                    end
                else
                    [maxIndx, outputSize] = obj.getMaxIndxAndInputSize(obj.Layers{i}.Name);
                    y = obj.Layers{i}.evaluate(y, maxIndx, outputSize); %
                end
            end
                
        end
        
        % evaluate parallel
        function seg_ims = evaluate_parallel(obj, images, nCores)
            % @images: a set of images
            % @nCores: number of cores used for computation
            % @seg_ims: segmentation images
            
            % author: Dung Tran
            % date:4/22/2020
            
            n = length(images);
            seg_ims = cell(1, n);
            
            if nCores < 1
                error("Invalid number of Cores");
            elseif nCores == 1
                for i=1:n
                    seg_ims{i} = obj.evaluate(images{i});
                end
            else
                obj.numCores = nCores;
                obj.start_pool();
                parfor i=1:n
                    seg_ims{i} = obj.evaluate(images{i});
                end
            end

        end
        
        
        % get the max indices and inputsize for max unpooling layer
        function [MaxIndx, InputSize] = getMaxIndxAndInputSize(obj, unpooling_layer_name)
            % @dest_name: destination name
            % @MaxIndx: max point index
            % @InputSize: the outputSize for max unpooling
            
            
            % author: Dung Tran
            % date: 4/20/2020
            
            
            if isempty(obj.Connections)
                error('No connection table');
            end
            
            if ~ischar(unpooling_layer_name)
                error('Invalid unpooling_layer_name');
            else
                dest_name = sprintf("%s/indices", unpooling_layer_name);
            end
            
            n = size(obj.Connections, 1);
            
            source_name = [];
            for i=1:n                
                if strcmp(obj.Connections.Destination(i), dest_name)
                    source_name = obj.Connections.Source(i);
                    break;
                end
            end
            
            if isempty(source_name)
                error('Unknown destination name');
            end
            
            maxpooling_layer_name = erase(source_name{1}, "/indices");
            MaxIndx = [];
            InputSize = [];
            for i=1:obj.numLayers
                if strcmp(obj.Layers{i}.Name, maxpooling_layer_name)
                    MaxIndx = obj.Layers{i}.MaxIndx;
                    InputSize = obj.Layers{i}.InputSize;
                    break;
                end
            end
            
        end
        
        
        % start parallel pool for computing
        function start_pool(obj)
            
            if obj.numCores > 1
                poolobj = gcp('nocreate'); % If no pool, do not create new one.
                if isempty(poolobj)
                    parpool('local', obj.numCores); 
                else
                    if poolobj.NumWorkers ~= obj.numCores
                        delete(poolobj); % delete the old poolobj
                        parpool('local', obj.numCores); % start the new one with new number of cores
                    end                    
                end
            end   
            
        end
        
        
         
    end
    
    
    methods % reachability analysis method
        
        function [IS, reachTime] = reach(varargin)            
            % @I: input set, a star set          
            % @method: = 'exact-star' or 'approx-star' -> compute reach set using stars
            %            'abs-dom' -> compute reach set using abstract
            %            domain (support in the future)
            %            'face-latice' -> compute reach set using
            %            face-latice (may support in the future)
            
            % @IS: output set is an ImageStar
            % @reachTime : reachable set computation time
            
            % author: Dung Tran
            % date:4/15/2020
            
            
            switch nargin 
                
                case 2
                    
                    obj = varargin{1};
                    inputSet = varargin{2};
                    obj.reachMethod = 'approx-star';
                    obj.numCores = 1; 
                    
                case 3 
                    
                    obj = varargin{1};
                    inputSet = varargin{2};
                    obj.reachMethod = varargin{3};
                    obj.numCores = 1; 
                    
                case 4
                    
                    obj = varargin{1};
                    inputSet = varargin{2};
                    obj.reachMethod = varargin{3};
                    obj.numCores = varargin{4}; 
                 
                otherwise 
                    
                    error('Invalid number of input arguments, the number should be 1, 2 or 3');
                
            end       
            
            if  obj.numCores > 1
                obj.start_pool;
                obj.reachOption = 'parallel';
            else
                obj.reachOption = [];
            end
            
            if ~isempty(obj.Connections)
                error("NNV have not yet support reachability for DAG networks");
            end
            
            obj.reachSet = cell(1, obj.numLayers+1);
            obj.reachTime = zeros(1, obj.numLayers);
            fprintf('\nPerform reachability analysis for the network %s...', obj.Name);
            obj.reachSet{1} = inputSet;
            for i=2:obj.numLayers+1
                fprintf('\nPerforming analysis for Layer %d (%s)...', i-1, obj.Layers{i-1}.Name);
                start_time = tic;
                if ~isa(obj.Layers{i-1}, 'PixelClassificationLayer')
                    obj.reachSet{i} = obj.Layers{i-1}.reach(obj.reachSet{i-1}, obj.reachMethod, obj.reachOption);
                else
                    [obj.reachSet{i}, ~, ~] = obj.Layers{i-1}.reach(obj.reachSet{i-1}, obj.reachMethod, obj.reachOption);
                end
                
                obj.reachTime(i-1) = toc(start_time);
                fprintf('\nReachability analysis for Layer %d (%s) is done in %.5f seconds', i-1, obj.Layers{i-1}.Name, obj.reachTime(i-1));
                fprintf('\nThe number of reachable sets at Layer %d (%s) is: %d', i-1, obj.Layers{i-1}.Name, length(obj.reachSet{i}));
            end
            fprintf('\nReachability analysis for the network %s is done in %.5f seconds', obj.Name, sum(obj.reachTime));
            fprintf('\nThe number ImageStar in the output sets is: %d', length(obj.reachSet{obj.numLayers+1}));
            obj.totalReachTime = sum(obj.reachTime);
            IS = obj.reachSet{obj.numLayers+1};
            reachTime = obj.totalReachTime;
        end
        
        
    end
    
    methods % verifier
        
        function [rb, n_mis, n_rb, ver_rs] = verify(varargin)
            % @in_images: an array of input set
            % @ground_truths: an array of ground truth images (without attack)
            % @method: reachability method
            % @nCores: number of cores used for computation
            % @rb: percentage of correctly classified pixels
            % @n_mis: number of misclassified pixels
            % @n_rb: number of robust pixels
            % @ver_rs: verified output reachable set, used for plot
                        
            % author: Dung Tran
            % date: 4/22/2020
            
            switch nargin
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    ground_truths = varargin{3};
                    method = varargin{4};
                    nCores = varargin{5};
                case 4
                    obj = varargin{1};
                    in_images = varargin{2};
                    ground_truths = varargin{3};
                    method = varargin{4};
                    nCores = 1;
                case 3
                    obj = varargin{1};
                    in_images = varargin{2};
                    ground_truths = varargin{3};
                    method = 'approx-star';
                    nCores = 1;
                otherwise
                    error("Invalid number of input arguments");
            end
            
            n1 = length(ground_truths);
            n2 = length(in_images);
            if n1 ~= n2
                error("Inconsistent number of ground truth images and input sets");
            end
            
            rb = zeros(1, n1);
            n_rb = zeros(1, n1);
            n_mis = zeros(1, n1);
            ver_rs = cell(1, n1);
                       
            % compute reachable set
            [seg_im_ids, ~] = obj.reach(in_images, method, nCores);
            
            % compute ground truth output segmentation image
            gr_seg_ims = obj.evaluate_parallel(ground_truths, nCores);
            
            % compute number of correctly classified pixels
            n_pixels = obj.OutputSize(1) * obj.OutputSize(2);
            
            for i=1:n1
                gr_im = gr_seg_ims{i};
                seg_im = seg_im_ids{i};
                c = (gr_im == seg_im);
                n_rb(i) = sum(c, 'all');
                n_mis(i) = n_pixels - n_rb(i);
                rb(i) = n_rb(i)/(n_pixels);
            end
            
            % obtain verify rearch set
            % 'unknown': the label of the pixel is unknown, may be correct may be not
            % 'misclass': the label of the pixel is misclassified 
            
            % we introduce 2 more classes for the reachable set
            % 1) unknown class
            % 2) misclassification class
                      
            for i=1:n1
                seg_im = seg_im_ids{i};
                gr_seg_im = gr_seg_ims{i};
                n = size(seg_im);
                ver_im = seg_im;
                for j=1:n(1)
                    for k=1:n(2)
                        % unk_id = obj.mis_id - 1
                        if (seg_im(j, k) ~= gr_seg_im(j,k)) && (seg_im(j,k) ~= obj.numClasses - 1)
                            ver_im(j,k) = obj.numClasses;
                        end
                    end
                end                
                ver_rs{i} = ver_im;
            end
            
            obj.verifiedOutputSet = ver_rs;
            obj.groundTruthSegIms = gr_seg_ims;
            obj.rb_val = rb;
            obj.numMisPixels = n_mis;
            obj.numPixels = n_pixels;
            obj.numRbPixels = n_rb;
            
        end
        
    end
    
    
    methods %plot
        
        % get color map corresponding to the RGB image output of label2rgb()
        % used to plot colorbar of a segmentation image
        % For example, RGB = label2rgb(C)
        % We want to find the color map of the RGB corresponding to the
        % idexes of C
        %       
        % We do: 
        %   1)  n = length(unique(C));
        %   2)  [IND, in_map] = rgb2ind(C,n)
        %   3) color_map = obj.getColorMap(C, IND, in_map);
        % the return color_map is corresponding to the label index in C
        
        function map = getColorMap(~, C, IND, in_map)
            % @C: is the label index matrix
            % @IND:
            % @in_map:
            % @map: return color_map
            
            % author: Dung Tran
            % date: 4/23/2020
            
            C_unique = unique(C);
            n = length(C_unique);
            [nC, mC] = size(C);
            map = zeros(n, 3);
            for i=1:n
                cat = C_unique(i);
                flag = 0;
                for j=1:nC
                    for k=1:mC
                        if C(j,k) == cat
                            id = IND(j, k) + 1;
                            map(i, :) = in_map(id, :);
                            flag = 1;
                            break;
                        end
                    end
                    if flag == 1
                        break;
                    end
                end

            end
 
        end
        
        % get classes corresponding to an class index array
        function classes = getClasses(obj, idxs)
            % @idxs: an array of class index
            % @classes: a string array of class name
            
            % author: Dung Tran
            % date: 4/23/2020
            
            classes  = obj.Layers{obj.numLayers}.getClasses(idxs);
        end
        
        
        % plot a segmentation image from an input image
        function plotSegmentationImage(varargin)
            % @im: in image
            % @figSize: figSize
            
            % author: Dung Tran
            % date: 4/23/2020
            
            
            switch nargin
                case 4
                    obj = varargin{1};
                    im = varargin{2};
                    figSize = varargin{3};
                    save_name = varargin{4};
                case 3
                    obj = varargin{1};
                    im = varargin{2};
                    figSize = varargin{3};
                    save_name = 'SegmentationImage.pdf';
                case 2
                    obj = varargin{1};
                    im = varargin{2};
                    figSize = [400 400];
                    save_name = 'SegmentationImage.pdf';
                
                otherwise
                    error('Invalid number of input arguments, should be 1, 2 or 3');
            end
            
            
            seg_im = obj.evaluate(im);
            RGB = label2rgb(seg_im); % get RGB image from label image
            seg_im_unique = unique(seg_im);
            m = length(seg_im_unique);
            [IND,in_map] = rgb2ind(RGB, m);
            map = obj.getColorMap(seg_im, IND, in_map);
            classes = obj.getClasses(seg_im_unique);
                       
            
            figure;
            subplot(1,2,1);
            imshow(im);
            
            ax = subplot(1,2,2);
            imshow(RGB);
            colormap(ax,map);
            cbh = colorbar(ax);
            xtick = 1/(2*m):1/m:1;
            cbh.Ticks = xtick;               
            cbh.TickLabels = classes;
            title("Segmentation image");
            truesize(figSize);
            saveas(gcf,save_name);
            
        end
        
        function plotPixelClassificationReachSet(varargin)
            % @ind: index of the reachable set (= index of the input set)
            
            % author: Dung Tran
            % date: 4/22/2020
            
            
            switch nargin
                case 4
                    obj = varargin{1};
                    im = varargin{2};
                    figSize = varargin{3};
                    save_name = varargin{4};
                case 3
                    obj = varargin{1};
                    ind = varargin{2};
                    figSize = varargin{3};
                    save_name = 'PixelClassificationReachSet.pdf';
                case 2
                    obj = varargin{1};
                    ind = varargin{2};
                    figSize = [400 400];
                    save_name = 'PixelClassificationReachSet.pdf';
                case 1
                    obj = varargin{1};
                    ind = 1; % plot the first reachable set
                    figSize = [400 400];
                    save_name = 'PixelClassificationReachSet.pdf';
                otherwise
                    error('Invalid number of input arguments, should be 0, 1, 2 or 3');
            end
            
            if isempty(obj.reachSet)
                error("Reachable set is empty, please perform reachability first");
            end
            
            if ~isa(obj.Layers{obj.numLayers}, 'PixelClassificationLayer')
                error("The last layer in the network is not a pixel classification layer, i.e., the network is not a segmentation network");
            end
            
            RS = obj.reachSet{obj.numLayers + 1};
                        
            if ind > length(RS) || ind < 1
                error("Invalid index");
            end

            seg_im = RS{ind};
            RGB = label2rgb(seg_im); % get RGB image from label image
            seg_im_unique = unique(seg_im);
            m = length(seg_im_unique);
            [IND,in_map] = rgb2ind(RGB, m);
            map = obj.getColorMap(seg_im, IND, in_map);
            classes = obj.getClasses(seg_im_unique);                       
            
            figure;
            imshow(RGB);
            colormap(gca,map);
            cbh = colorbar(gca);
            m = size(map, 1);
            xtick = 1/(2*m):1/m:1;
            cbh.Ticks = xtick;               
            cbh.TickLabels = classes;
            str = sprintf("%d^{th} segmentation reachable set", ind);
            title(str);
            truesize(figSize);
            saveas(gcf,save_name);

        end
        
        
        % plot verified output set
        function plotVerifiedOutputSet(varargin)
            % @ind: index of the reachable set (= index of the input set)
            
            % author: Dung Tran
            % date: 4/22/2020
            
            
            switch nargin
                
                case 4
                    obj = varargin{1};
                    ind = varargin{2};
                    figSize = varargin{3};
                    save_name = varargin{4};
                
                case 3
                    obj = varargin{1};
                    ind = varargin{2};
                    figSize = varargin{3};
                    save_name = 'VerifiedOutputSet.pdf';
                case 2
                    obj = varargin{1};
                    ind = varargin{2};
                    figSize = [400 400];
                    save_name = 'VerifiedOutputSet.pdf';
                case 1
                    obj = varargin{1};
                    ind = 1; % plot the first reachable set
                    figSize = [400 400];
                    save_name = 'VerifiedOutputSet.pdf';
                otherwise
                    error('Invalid number of input arguments, should be 0, 1, or 2');
            end
            
            if isempty(obj.verifiedOutputSet)
                error("Verified Output Reachable set is empty, please perform verify method first");
            end
            
            rs = obj.verifiedOutputSet;
            gr = obj.groundTruthSegIms;
            
            if ind > length(rs) || ind < 1
                error("Invalid index");
            end

            pl_rs = rs{ind};
            pl_gr = gr{ind};         
            gr_RGB = label2rgb(pl_gr);
            rs_RGB = label2rgb(pl_rs);
            
            gr_unique = unique(pl_gr);
            m1 = length(gr_unique);
            [IND1,in_map1] = rgb2ind(gr_RGB, m1);
            map1 = obj.getColorMap(pl_gr, IND1, in_map1);
            classes1 = obj.getClasses(gr_unique); 
            
            rs_unique = unique(pl_rs);
            m2 = length(rs_unique);
            [IND2,in_map2] = rgb2ind(rs_RGB, m2);
            map2 = obj.getColorMap(pl_rs, IND2, in_map2);
            classes2 = obj.getClasses(rs_unique); 
            
            figure;
            ax1 = subplot(1,2,1);
            imshow(gr_RGB);
            colormap(ax1,map1);
            cbh1 = colorbar(ax1);
            xtick = 1/(2*m1):1/m1:1;
            cbh1.Ticks = xtick;               
            cbh1.TickLabels = classes1;
            str = sprintf("%d^{th} segmentation image without attack", ind);
            title(str);

            ax2 = subplot(1,2,2);
            imshow(rs_RGB);
            colormap(ax2,map2);
            cbh2 = colorbar(ax2);
            xtick = 1/(2*m2):1/m2:1;
            cbh2.Ticks = xtick;               
            cbh2.TickLabels = classes2;
            str = sprintf("%d^{th} verified output segmentation image after attack", ind);
            title(str);

            truesize(figSize);
            saveas(gcf,save_name);
            
        end
        
        % plot robustness statistics
        function plotRobustnessStatistics(varargin)
            % @option: 'individual' or ''
            % @x: value of x axes
            
            % author: Dung Tran
            % date: 4/23/2020
            
            switch nargin
                case 1
                    obj = varargin{1};
                    x = length(obj.rb_val);
                    x_label = "Index";
                case 3
                    obj = varargin{1};
                    x = varargin{2};
                    x_label = varargin{3};
                otherwise
                    error("Invalid number of input arguments, should be 1 or 3");
            end
            
            if isempty(obj.rb_val)
                error("Robustness Statistics is empty, please perform verify method first");
            end
            
            n = length(obj.rb_val);
            
            if length(x) ~= n
                error('Inconsistency between the length of x and y');
            end
            
            figure;
            b = bar(x, obj.rb_val, 'y');
            xtips1 = b(1).XEndPoints;
            ytips1 = b(1).YEndPoints;
            labels1 = string(b(1).YData);
            text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom')
            xlabel(x_label);
            ylim([0 1.1]);
            title("Robustness Value");
            saveas(gcf, 'RobustValue.pdf');
            
            
            figure;
            b = bar(x, obj.numMisPixels, 'r');
            xtips1 = b(1).XEndPoints;
            ytips1 = b(1).YEndPoints;
            labels1 = string(b(1).YData);
            text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom');
            xlabel(x_label);
            ylim([0 max(obj.numMisPixels)+20]);
            title("Number of Misclassified Pixels");
            saveas(gcf, 'numMisPixels.pdf');
            
            figure;
            b = bar(x, obj.numRbPixels, 'b');
            xtips1 = b(1).XEndPoints;
            ytips1 = b(1).YEndPoints;
            labels1 = string(b(1).YData);
            text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom');
            xlabel(x_label);
            ylim([0 max(obj.numRbPixels)+500]);
            title("Number of Robust Pixels (Correctly Classified)");
            saveas(gcf, 'numRbPixels.pdf');
        end
        
    end
    

    methods(Static)
       
        % parse a dag neural network from matlab for reachability analysis
        function segnet = parse(varargin)
            % @net: input network, should be a dagnetwork or lgraph object
            % or a SeriesNetwork
                       
            % the constructed DAGNN for reachability analysis get rid of the
            % these following layers:
            % 2) Dropout Layer (is not used for prediction phase)
            
            % author: Dung Tran
            % date: 4/14/2020
            % update: 4/10/2020, 4/20/2020
            
            
            
            switch nargin
                case 1
                    net = varargin{1};
                    name = 'parsed_net';
                case 2
                    net = varargin{1};
                    name = varargin{2};
                otherwise
                    error('Invalid number of input arguments, should be 1 or 2');
            end
            
            
            
            if ~isa(net, 'DAGNetwork') && ~isa(net, 'LayerGraph') && ~isa(net, 'SeriesNetwork')
                error('the parsed object is not a Matlab DAGNetwork object, a LayerGraph object or a SeriesNetwork object');
            end
            
            n = length(net.Layers); % number of layers
                        
            Ls = {};
            j = 0;
            for i=1:n
                L = net.Layers(i);
                fprintf('\nParsing Layer %d...', i);
                
                if isa(L, 'nnet.cnn.layer.DropoutLayer') || isa(L, 'nnet.cnn.layer.ClassificationOutputLayer')                   
                    fprintf('\nLayer %d is a %s class which is neglected in the analysis phase', i, class(L));
                else
                    
                    if isa(L, 'nnet.cnn.layer.ImageInputLayer')
                        Li = ImageInputLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.Convolution2DLayer') 
                        Li = Conv2DLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.ReLULayer')
                        Li = ReluLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.BatchNormalizationLayer')
                        Li = BatchNormalizationLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.MaxPooling2DLayer')
                        Li = MaxPooling2DLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.AveragePooling2DLayer')
                        Li = AveragePooling2DLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.FullyConnectedLayer')
                        Li = FullyConnectedLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.MaxUnpooling2DLayer')
                        Li = MaxUnpooling2DLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.PixelClassificationLayer')
                        Li = PixelClassificationLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.SoftmaxLayer')
                        Li = SoftmaxLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.TransposedConvolution2DLayer')
                        Li = TransposedConv2DLayer.parse(L);
                    else                     
                        fprintf('\nLayer %d is a %s which have not supported yet in nnv, please consider removing this layer for the analysis', i, class(L));
                        error('\nUnsupported Class of Layer');                     
                    end
                    j = j + 1;
                    Ls{j} = Li;       
                end
                            
            end
            
            % update connection table
            if isprop(net, 'Connections')
                connects = net.Connections;
            else
                connects = [];
            end
            segnet = SEGNET(name, Ls, connects, net.InputNames, net.OutputNames);
            fprintf('\nParsing network is done successfully and %d Layers are neglected in the analysis phase', n - j);
            
            if ~isa(segnet.Layers{segnet.numLayers}, 'PixelClassificationLayer')
                error("The last layer in the network is not a pixel classification layer, i.e., the network is not a segmentation network");
            end
            
        end
        
             
        %
        
        
    end
    
    
    
end

