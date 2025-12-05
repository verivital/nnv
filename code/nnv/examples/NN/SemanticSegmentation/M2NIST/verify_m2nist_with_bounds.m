function verify_m2nist_with_bounds()
    % Run M2NIST segmentation verification with pixel-wise bounds visualization
    % Similar to Weizmann Horse example but for M2NIST digit segmentation

    fprintf('=== M2NIST Segmentation Verification with Bounds Visualization ===\n\n');

    % Load network
    fprintf('Loading M2NIST segmentation network...\n');
    net = load("models/m2nist_dilated_72iou_24layer.mat");
    net = matlab2nnv(net.net);

    % Load images
    fprintf('Loading test images...\n');
    images = load('m2nist_6484_test_images.mat');
    im_data = images.im_data;

    % Create example input set
    Nmax = 50;   % maximum allowable number of attacked pixels
    de = 0.0001; % disturbance
    Nt = 150;    % threshold value

    % Randomly select 1 image to verify
    rng(0);
    img_idx = randperm(1000,1);

    fprintf('Selected image index: %d\n', img_idx);
    fprintf('Image size: 64x84\n\n');

    % Create input set from adversarial perturbation
    fprintf('Creating adversarial input set...\n');
    fprintf('  Max attacked pixels: %d\n', Nmax);
    fprintf('  Disturbance: %.4f\n', de);
    fprintf('  Brightness threshold: %d\n\n', Nt);

    % Initialize vars
    ct = 0;
    flag = 0;
    im = im_data(:,:,img_idx);
    at_im = im;
    for i=1:64
        for j=1:84
            if im(i,j) > Nt
                at_im(i,j) = 0;
                ct = ct + 1;
                if ct == Nmax
                    flag = 1;
                    break;
                end
            end
        end
        if flag == 1
            break;
        end
    end

    fprintf('Attacked %d pixels\n', ct);

    % Define input set as ImageStar
    dif_im = im - at_im;
    noise = -dif_im;
    V(:,:,:,1) = double(im);
    V(:,:,:,2) = double(noise);
    C = [1; -1];
    d = [1; de-1];
    IS = ImageStar(V, C, d, 1-de, 1);
    GrTruth = {im};

    %% Evaluate network on original and adversarial images

    fprintf('\nEvaluating segmentation...\n');
    seg_original = net.evaluate(double(im));
    seg_adversarial = net.evaluate(double(at_im));

    fprintf('  Original output size: [%s]\n', num2str(size(seg_original)));
    fprintf('  Adversarial output size: [%s]\n', num2str(size(seg_adversarial)));

    %% Verify network

    fprintf('\n--- Reachability-based Segmentation Verification ---\n');
    reachOptions.reachMethod = 'approx-star';
    t = tic;
    [riou, rv, rs, n_rb, n_mis, n_unk, n_att, ver_rs, eval_seg_ims] = ...
        net.verify_segmentation(IS, GrTruth, reachOptions);
    verification_time = toc(t);

    fprintf('Verification completed in %.2f seconds\n\n', verification_time);

    fprintf('Results:\n');
    fprintf('  Robust IoU: %.2f%%\n', riou * 100);
    fprintf('  Robust verification: %d pixels\n', rv);
    fprintf('  Robust segmentation: %d pixels\n', rs);
    fprintf('  Robust pixels (both): %d\n', n_rb);
    fprintf('  Misclassified pixels: %d\n', n_mis);
    fprintf('  Unknown pixels: %d\n', n_unk);
    fprintf('  Attacked pixels: %d\n', n_att);

    %% Visualize traditional segmentation output

    fprintf('\nVisualizing segmentation results...\n');
    figure('Name', 'M2NIST Segmentation Verification');
    net.plot_segmentation_output_set(ver_rs{1}, eval_seg_ims{1});

    %% Visualize Pixel-wise Bounds (MNIST-style interval plot)

    fprintf('\n--- Visualizing Pixel-wise Reachability Bounds ---\n');

    % Get actual network output
    % For M2NIST, we need the softmax probabilities
    if length(size(seg_original)) == 3
        actual_output = seg_original;
        fprintf('  Using actual softmax output\n');
    else
        fprintf('  Using midrange of bounds as proxy\n');
        actual_output = [];
    end

    % Get the reachability sets from the network
    fprintf('  Network has %d layers with reach sets\n', length(net.reachSet));

    % Use softmax layer output (one before pixel classification)
    if length(net.reachSet) >= 2
        final_reach_sets = net.reachSet{end-1}; % Softmax layer
    else
        final_reach_sets = net.reachSet{end};
    end

    fprintf('  Using reach set from layer %d\n', length(net.reachSet) - 1);

    % Determine number of classes from output
    if ~isempty(actual_output) && length(size(actual_output)) == 3
        num_classes = size(actual_output, 3);
    else
        % M2NIST typically has 11 classes (background + digits 0-9)
        num_classes = 11;
    end

    class_names = cell(1, num_classes);
    class_names{1} = 'background';
    for i = 2:num_classes
        class_names{i} = sprintf('digit_%d', i-2);
    end

    fprintf('  Classes: %d\n', num_classes);

    % Image size
    img_size = [64, 84];

    % Call visualization function (from Weizmann Horse example)
    % Add path to function if needed
    addpath('../../../Tutorial/NN/WeizmannHorse');
    plot_segmentation_bounds(final_reach_sets, actual_output, class_names, img_size);

    %% Summary

    fprintf('\n=== Verification Summary ===\n');
    fprintf('Dataset: M2NIST (digit segmentation)\n');
    fprintf('Image size: 64x84\n');
    fprintf('Classes: %d (%s)\n', num_classes, strjoin(class_names, ', '));
    fprintf('Attack: Darkening %d pixels\n', ct);
    fprintf('Perturbation: %.4f\n', de);
    fprintf('Robust IoU: %.2f%%\n', riou * 100);
    fprintf('Robust pixels: %d/%d (%.1f%%)\n', n_rb, numel(im), n_rb/numel(im)*100);
    fprintf('Verification time: %.2f seconds\n', verification_time);

end
