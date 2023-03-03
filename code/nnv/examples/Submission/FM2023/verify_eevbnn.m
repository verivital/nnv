% X_mnist = load("/home/crescendo/Workspace/CAV_2022/eevbnn_comparison/test_data/X_mnist").data;
% y_mnist = load("/home/crescendo/Workspace/CAV_2022/eevbnn_comparison/test_data/y_mnist").data;
% 
% X_cifar = load("/home/crescendo/Workspace/CAV_2022/eevbnn_comparison/test_data/X_cifar").data;
% y_cifar = load("/home/crescendo/Workspace/CAV_2022/eevbnn_comparison/test_data/y_cifar").data;
% 
% a1 = reshape(X_mnist(1,:,:), [28,28]);
% a2 = reshape(X_cifar(1,:,:), [32,32,3]);
% 
% imwrite(uint8(a1 * 255), '/home/crescendo/Workspace/CAV_2022/eevbnn_comparison/test_data/test_img.jpeg');
% imwrite(uint8(a2 * 255), '/home/crescendo/Workspace/CAV_2022/eevbnn_comparison/test_data/test_img1.jpeg');
% 
% 
% k = 0;

format long

root_directory = "/home/crescendo/Workspace/CAV_2022";

models_folder = "eevbnn_comparison/saved_models";

names = [
       "mnist-s-adv0.1", ...
       "mnist-s-adv0.3", ...
       "mnist-l-adv0.1-cbd3-vrf", ...
       "mnist-l-adv0.3-cbd3-vrf", ...
%        "cifar10-s-adv2", ...
%        "cifar10-s-adv8", ...
%        "cifar10-l-adv2-cbd3-vrf", ...
%        "cifar10-l-adv8-cbd3-vrf"
    ];

quantizing_steps = [
       0.61, ...
       0.61, ...
       0.61, ...
       0.61, ...
%        0.064, ...
%        0.064, ...
%        0.064, ...
%        0.064
    ];

deltas = [
        0.3, ...
        0.1, ...
        0.3, ...
        0.1, ...
%         8/255, ...
%         2/255, ...
%         8/255, ...
%         2/255
    ];

models_layers = {};

layers = {};

X_mnist = load("/home/crescendo/Workspace/CAV_2022/eevbnn_comparison/test_data/X_mnist").data;
y_mnist = load("/home/crescendo/Workspace/CAV_2022/eevbnn_comparison/test_data/y_mnist").data;

X_cifar = load("/home/crescendo/Workspace/CAV_2022/eevbnn_comparison/test_data/X_cifar").data;
y_cifar = load("/home/crescendo/Workspace/CAV_2022/eevbnn_comparison/test_data/y_cifar").data;

X = { 
    X_mnist, ...
      X_mnist, ...
      X_mnist, ...
      X_mnist, ...
%       X_cifar, ...
%       X_cifar, ...
%       X_cifar, ...
%       X_cifar 
      };
  
y = { 
    y_mnist, ...
      y_mnist, ...
      y_mnist, ...
      y_mnist, ...
%       y_cifar, ...
%       y_cifar, ...
%       y_cifar, ...
%       y_cifar 
      };
  
for i=1:length(names)
    current_dir = strcat(root_directory, '/', models_folder, '/', names(i));
    files = dir(current_dir);
    
    files = files(3:length(files));
    
    dir_names = {};
    
    for j=1:length(files)
        dir_names{j} = files(j).name;
    end
     
    [~, sorting_order] = natsortfiles(dir_names);
    sorted_directories = dir_names(sorting_order);
    
    for j=1:length(sorted_directories)
        disp(sorted_directories{j});
        layers{j} = load_and_parse_layer(convertCharsToStrings(sorted_directories{j}), current_dir);
    end
    
    models_layers{i} = layers;
end

verification_results = {};
counterex_ids = {};
counterexs = {};
imgs = {};

for i=1:length(models_layers)
	[verification_results{i}, counterex_ids{i}, counterexs{i}, imgs{i}] = verify_model(BNN(models_layers{i}), X{i}, y{i}, deltas(i), quantizing_steps(i));
    
    current_result = verification_results{i};
    current_counterex_ids = counterex_ids{i};
    current_counterexs = counterexs{i};
    current_imgs = imgs{i};
    
    current_name = strcat(names(i), "__vrf_result.mat");
    current_cexids_name = strcat(names(i), "__cexids_result.mat");
    current_cexs_name = strcat(names(i), "__cexs_result.mat");
    current_imgs_name = strcat(names(i), "__imgs_result.mat");
    
    save(current_name, 'current_result');
    save(current_cexids_name, 'current_counterex_ids');
    save(current_cexs_name, 'current_counterexs');
    save(current_imgs_name, 'current_imgs');
end

function [result, counterex_ids, counterexs, imgs] = verify_model(model, examples, labels, delta, quantizing_step)

    ex_num = size(examples, 1);

    result = zeros(ex_num, 6); % img_id, delta, UNSAT(1)/SAT(0), verif_time, valid_label, verified_label

    counterex_ids = [];
    counterexs = {};
    imgs = {};
    
    for i=1:ex_num
        if i == 8
            k = 0;
        end
        if length(size(examples)) == 3
            current_image = reshape(examples(i,:,:), [28,28]);
        else
            current_image = reshape(examples(i,:,:,:), [32,32,3]);
        end
                
        ev_res = model.evaluate(quantize_input(current_image, quantizing_step));
        [~, ev_class] = max(ev_res);

        valid_class = labels(i) + 1;
        
        lb = current_image - delta;
        ub = current_image + delta;

        lb(lb < 0) = 0;
        lb(lb > 255) = 255;
        ub(ub < 0) = 0;
        ub(ub > 255) = 255;

        lb_quant = quantize_input(lb, quantizing_step);
        ub_quant = quantize_input(ub, quantizing_step);

        S = ImageStar(lb_quant, ub_quant);

        %% Verify
        disp(i);
        [~, reachT{i}, signReachT{i}] = model.reach(S, 'approx-star', 1);

        [val, ver_class] = max(model.outputSet.V(:, 1));
        
        result(i, 1) = i;
        result(i, 2) = delta;
        
        if valid_class == ver_class
            result(i, 3) = 1;
        else
            [verdict, counterEx] = findCounterEx(S, model, 10, valid_class);
            if verdict == 0
                counterex_ids = [counterex_ids i];
                counterexs{length(counterexs) + 1} = counterEx;
                imgs{length(imgs) + 1} = current_image;
            end
            
            result(i, 3) = verdict;
        end
        
        result(i, 4) = reachT{i};
        result(i, 5) = valid_class;
        result(i, 6) = ver_class;
    end

end

function layer = load_and_parse_layer(unparsed_layer, current_directory)
    split_name = split(unparsed_layer, "_");
    
    if strcmp(split_name(2), 'BinConv2d')
        weight = double(load(strcat(current_directory, '/', unparsed_layer, '/weight.mat')).layer);
        stride = double(load(strcat(current_directory, '/', unparsed_layer, '/stride.mat')).layer);
        padding = double(load(strcat(current_directory, '/', unparsed_layer, '/padding.mat')).layer);
        dilation = double(load(strcat(current_directory, '/', unparsed_layer, '/dilation.mat')).layer);
        
        weight = fix_conv_kernel(weight);
        
        padding = fix_conv_padding(padding);
        
        layer = Conv2DLayer(weight, zeros(1, 1, size(weight, 4)), padding, stride, dilation);
    elseif strcmp(split_name(2), 'BatchNormStatsCallbak')
        num_channels = double(load(strcat(current_directory, '/', unparsed_layer, '/num_channels.mat')).layer);
%         var = double(load(strcat(current_directory, '/', unparsed_layer, '/var.mat')).layer);
        scale = double(load(strcat(current_directory, '/', unparsed_layer, '/scale.mat')).layer);
%         mean = double(load(strcat(current_directory, '/', unparsed_layer, '/mean.mat')).layer);
        bias = double(load(strcat(current_directory, '/', unparsed_layer, '/bias.mat')).layer);
        
%         var = reshape(var, [1 size(var)]);
        if length(size(scale)) < 3
            scale = reshape(scale, [1 size(scale)]);
        end
%         mean = reshape(mean, [1 size(mean)]);
        
        if length(size(bias)) < 3
            bias = reshape(bias, [1 size(bias)]);
        end
            
        layer = BatchNormalizationLayer('NumChannels', num_channels, 'Offset', bias, 'Scale', scale);

    elseif strcmp(split_name(2), 'Binarize01Act')
        layer = SignLayer(0, 'nonnegative_zero_to_pos_one');
    elseif strcmp(split_name(2), 'BinConv2dPos')
        weight = double(load(strcat(current_directory, '/', unparsed_layer, '/weight.mat')).layer);
        bias = double(load(strcat(current_directory, '/', unparsed_layer, '/bias.mat')).layer);
        stride = double(load(strcat(current_directory, '/', unparsed_layer, '/stride.mat')).layer);
        padding = double(load(strcat(current_directory, '/', unparsed_layer, '/padding.mat')).layer);
        dilation = double(load(strcat(current_directory, '/', unparsed_layer, '/dilation.mat')).layer);
        
        weight = fix_conv_kernel(weight);
        
        padding = fix_conv_padding(padding);
        
        layer = Conv2DLayer(weight, reshape(bias, [1,1,size(bias,2)]), padding, stride, dilation);
    elseif strcmp(split_name(2), 'BinLinearPos')
        weight = double(load(strcat(current_directory, '/', unparsed_layer, '/weight.mat')).layer);
        % TODO: might need to transpose the weights instead
        bias = double(load(strcat(current_directory, '/', unparsed_layer, '/bias.mat')).layer)';
        
        layer = FullyConnectedLayer(weight, bias); %, 'purelin'
    elseif strcmp(split_name(2), 'Flatten')
        layer = FlattenLayer('nnet.keras.layer.FlattenCStyleLayer', 'row');
    end
end

function kernel = fix_conv_kernel(weight)
    kernel = zeros(flip(size(weight)));
    
    weight_size = size(weight);
    
    buffer = zeros(weight_size(3), weight_size(4));
    
    for i = 1:weight_size(1)
        for j = 1:weight_size(2)
            buffer(:,:) = weight(i,j,:,:);
            
            kernel(:,:,j,i) = buffer;
        end
    end
end

function fixed_padding = fix_conv_padding(padding)
    fixed_padding = [padding ones(1, 4 - size(padding,2))];
end

function quantized_input = quantize_input(input, step)
    quantized_input = round(input / step);
    
    quantized_input = quantized_input * step;
end

function [flag, counterEx] = findCounterEx(In, bnn, samplesNum, groundLabel)
    flag = -1;
    counterEx = [];

    Out = bnn.outputSet;
    v_size = size(In.V);
    
    cL = 1;
    
    for i=1:length(v_size)-1
        cL = cL * v_size(i);
    end
    
    l = v_size(length(v_size) - 1);
    
    inputs = zeros(cL, samplesNum);
    
    for i = 1:samplesNum
        predicateVecs = zeros(cL, l);
        inputs(:, i) = reshape(In.V(:, :, :, 1), [cL, 1]);
        
        for j = 2:l
            lb = Out.pred_lb(j);
            ub = Out.pred_ub(j);
            
            predicateVecs(:, j-1) = lb + (ub - lb) .* rand(cL,1);
            
            current_pred_vector = reshape(In.V(:, :, :, j + 1), [cL, 1]);
            
            inputs(:, i) = inputs(:, i) + current_pred_vector .* predicateVecs(:, j-1);
        end
         
        current_input = reshape(inputs(:, i), v_size(1:length(v_size) - 1));
        
        [res, pos] = max(bnn.evaluate(current_input));
        
        if pos ~= groundLabel
            flag = 0;
            counterEx = current_input;
            return;
        end
    end
end

function S = constructInput(data, disturbances)
    S = [];

    for i=1:size(data, 2)
        lb = data(:, i) - disturbances(i, 1);
        ub = data(:, i) + disturbances(i, 1);

        lb(lb < 0) = 0;
        lb(lb > 255) = 255;
        ub(ub < 0) = 0;
        ub (ub > 255) = 255;

        S = [S Star(lb, ub)];
    end
end

function v = constructVec(data)
    v = zeros(size(data,3),1);
    
    for i = 1:length(data)
        v(i) = data(1,1,i);
    end
end