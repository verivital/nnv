classdef SequenceInputLayer < handle
    % The Sequence input layer class in CNN
    %   Contain constructor and reachability analysis methods   

    %   Author: Neelanjana Pal
    %   Date: 10/25/2021     
    
    properties
        Name = 'SequenceInputLayer';
        InputSize = [];
        Normalization = 'none'; %default
        MinLength = 1; %default
        Mean = []; % in 
        StandardDeviation = [];
        %Min = [];
        %Max = [];
        %NormalizationDimension = 'auto';

        % layer properties
        NumInputs = 0;
        InputNames = {};
        NumOutputs = 1;
        OutputNames = {'out'};
                
    end
    
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = SequenceInputLayer(varargin)           
            % author: Neelanjana Pal
            % date: 10/25/2021
            
            
            if mod(nargin, 2) ~= 0
                error('Invalid number of input arguments');
            end
            
            for i=1:nargin-1
                
                if mod(i, 2) ~= 0
                    
                    if strcmp(varargin{i}, 'Name')
                        obj.Name = varargin{i+1};
                    elseif strcmp(varargin{i}, 'InputSize')
                        obj.InputSize = varargin{i+1};
                    elseif strcmp(varargin{i}, 'Normalization')
                        obj.Normalization = varargin{i+1};
                    elseif strcmp(varargin{i}, 'MinLength')
                        obj.MinLength = varargin{i+1};
                    elseif strcmp(varargin{i}, 'Mean')
                        obj.Mean = double(varargin{i+1});
                    elseif strcmp(varargin{i}, 'StandardDeviation')
                        obj.StandardDeviation = double(varargin{i+1});
%                     elseif strcmp(varargin{i}, 'Min')
%                         obj.Mean = double(varargin{i+1});
%                     elseif strcmp(varargin{i}, 'Max')
%                         obj.Mean = double(varargin{i+1});
%                     elseif strcmp(varargin{i}, 'NormalizationDimension')
%                         obj.Mean = double(varargin{i+1});
                    elseif strcmp(varargin{i}, 'NumInputs')
                        obj.NumInputs = double(varargin{i+1});
                    elseif strcmp(varargin{i}, 'InputNames')
                        obj.InputNames = varargin{i+1};
                    elseif strcmp(varargin{i}, 'NumOutputs')
                        obj.NumOutputs = double(varargin{i+1});
                    elseif strcmp(varargin{i}, 'OutputNames')
                        obj.OutputNames = varargin{i+1};
                    end
                    
                end
                
            end
             
        end
        
        
    end
        
    % evaluation method
    methods
        
        function y = evaluateSequence(obj, input)
            % @input: input image
            % @y: output image with normalization
                             
            if isempty(obj.Mean) || isequal(obj.Normalization, 'none')
                y = double(input);
            elseif isequal(obj.Normalization, 'zerocenter')
                y = double(input) - double(obj.Mean); % zerocenter nomalization
            elseif isequal(obj.Normalization, 'zscore')
                y = (double(input) - double(obj.Mean))./double(obj.StandardDeviation); % zscore nomalization
            end
                               
        end

        function y = evaluate(obj, input)
            y = obj.evaluateSequence(input);
        end
         
    end
        
    
    methods % reachability methods
        
        function image = reach_star_single_input(obj, in_image)
            % @in_image: an input ImageStar
            % @image: an output ImageStar
            
            if ~isa(in_image, 'ImageStar')
                error('Input is not an ImageStar');
            end
            if isempty(obj.Mean) || isequal(obj.Normalization, 'none')
                image = in_image;
            elseif isequal(obj.Normalization, 'zerocenter')
                image = in_image.affineMap([], -obj.Mean); % zerocenter nomalization
            elseif isequal(obj.Normalization, 'zscore')
                image = in_image.affineMap(1./obj.StandardDeviation, -obj.Mean./obj.StandardDeviation); % zscore nomalization
            end
        end
        
        % handling multiple inputs
        function images = reach_star_multipleInputs(obj, in_images, option)
            % @in_images: an array of ImageStars
            % @method: = 'exact-star' or 'approx-star' or 'abs-dom'
            % @option: = 'parallel' or 'single' or empty
            % @images: an array of ImageStar (if we use 'exact-star' method)
            %         or a single ImageStar set
            
            n = length(in_images);
            images(n) = ImageStar;            
            if strcmp(option, 'parallel')
                parfor i=1:n
                    images(i) = obj.reach_star_single_input(in_images(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    images(i) = obj.reach_star_single_input(in_images(i));
                end
            else
                error('Unknown computation option');

            end
        end
        
        function image = reach_zono(obj, in_image)
            % @in_image: an input ImageZono
            % @image: an output ImageZono
            
            if ~isa(in_image, 'ImageZono')
                error('Input is not an ImageZono');
            end
            if isempty(obj.Mean) || isequal(obj.Normalization, 'none')
                image = in_image;
            elseif isequal(obj.Normalization, 'zerocenter')
                image = in_image.affineMap([], -obj.Mean); % zerocenter nomalization
            elseif isequal(obj.Normalization, 'zscore')
                image = in_image.affineMap(1./obj.StandardDeviation, -obj.Mean./obj.StandardDeviation); % zscore nomalization
            end
        end
        
        % handling multiple inputs
%         function images = reach_zono_multipleInputs(obj, in_images, option)
%             % @in_images: an array of ImageZonos
%             % @option: = 'parallel' or 'single' or empty
%             % @images: an array of ImageZono 
%             
%             % author: Dung Tran
%             % date: 1/7/2020
%             
%             n = length(in_images);
%             images(n) = ImageZono;            
%             if strcmp(option, 'parallel')
%                 parfor i=1:n
%                     images(i) = obj.reach_zono(in_images(i));
%                 end
%             elseif strcmp(option, 'single') || isempty(option)
%                 for i=1:n
%                     images(i) = obj.reach_zono(in_images(i));
%                 end
%             else
%                 error('Unknown computation option');
% 
%             end
%         end
        
        function signal = reach_star_Signal_single_input(obj, in_signal)
            % @in_Signal: an input SignalStar
            % @Signal: an output SignalStar
            
            if ~isa(in_signal, 'SignalStar')
                error('Input is not an SignalStar');
            end
            if isempty(obj.Mean) || isequal(obj.Normalization, 'none')
                signal = in_signal;
            elseif is(obj.Normalization == 'zerocenter')
                signal = in_signal.affineMap([], -obj.Mean); % zerocenter nomalization
%             elseif is(obj.Normalization == 'zscore')
%                 image = in_image.affineMap([], -obj.Mean);
%                 y = (double(input) - double(obj.Mean))./double(obj.StandardDeviation); % zscore nomalization
            end
        end
        
        % handling multiple inputs
%         function images = reach_star_Signal_multipleInputs(obj, in_signals, option)
%             % @in_signals: an array of ImageStars
%             % @method: = 'exact-star' or 'approx-star' or 'abs-dom'
%             % @option: = 'parallel' or 'single' or empty
%             % @signals: an array of SignalStar (if we use 'exact-star' method)
%             %         or a single SignalStar set
%             
%             % author: Dung Tran
%             % date: 1/7/2020
%             
%             n = length(in_signals);
%             images(n) = SignalStar;            
%             if strcmp(option, 'parallel')
%                 parfor i=1:n
%                     images(i) = obj.reach_star_Signal_single_input(in_signals(i));
%                 end
%             elseif strcmp(option, 'single') || isempty(option)
%                 for i=1:n
%                     images(i) = obj.reach_star_Signal_single_input(in_signals(i));
%                 end
%             else
%                 error('Unknown computation option');
% 
%             end
%         end

function sequences = reachSequence(varargin)
            % @in_image: an input imagestar
            % @image: output set
            % @option: = 'single' or 'parallel' 
             
            switch nargin
                
                 case 7
                    obj = varargin{1};
                    in_seqs = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    % relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                    % lp_solver = varargin{7}; do not use
                
                case 6
                    obj = varargin{1};
                    in_seqs = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                
                case 5
                    obj = varargin{1};
                    in_seqs = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                
                case 4
                    obj = varargin{1};
                    in_seqs = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                
                case 3
                    obj = varargin{1};
                    in_seqs = varargin{2};
                    method = varargin{3};
                    option = 'single';
                
                otherwise
                    error('Invalid number of input arguments (should be 2, 3, 4, 5, or 6)');
            end      
      
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || contains(method, "relax-star")
                sequences = obj.reach_star_multipleInputs(in_seqs, option);
            elseif strcmp(method, 'approx-zono')
                sequences = obj.reach_zono_multipleInputs(in_seqs, option);
            else
                error('Unknown reachability method');
            end   
            
                      
        end
        
        
    end
    
    
    methods(Static)
         % parse a trained input image layer from matlab
         function L = parse(input_sequence_layer)
            % @input_image_layer: input layer
            % @L: constructed layer
            
            
            if ~isa(input_sequence_layer, 'nnet.cnn.layer.SequenceInputLayer')
                error('Input is not a Matlab nnet.cnn.layer.SequenceInputLayer class');
            end
            
            if isprop(input_sequence_layer, 'Mean')
                L = SequenceInputLayer('Name', input_sequence_layer.Name, 'InputSize', input_sequence_layer.InputSize, 'Mean', input_sequence_layer.Mean,'StandardDeviation',input_sequence_layer.StandardDeviation,'Normalization', input_sequence_layer.Normalization,'MinLength',input_sequence_layer.MinLength);
            elseif isprop(input_image_layer, 'AverageImage')
                L = SequenceInputLayer('Name', input_sequence_layer.Name, 'InputSize', input_sequence_layer.InputSize, 'Mean', input_sequence_layer.AverageImage);
            else
                error('Mean or AverageImage property does not exist in the Input Image Layer');
            end
            
            fprintf('\nParsing a Matlab sequence input layer is done successfully');
            
        end
        
    end
    
    
    
end

