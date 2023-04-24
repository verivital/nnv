function example_transposed()
    % Run a segmentation example using a segnet with transposed convolution
    
    % Load network
    net = load("models/m2nist_75iou_transposedcnn_avgpool.mat");
    net = matlab2nnv(net.net);
    
    % Load images
    images = load('data/M2NIST/m2nist_6484_test_images.mat');
    im_data = images.im_data;
    
    % Create example input set
    Nmax = 50; % maximum allowable number of attacked pixels
    de = 0.0001; % disturbance
    Nt = 150; % threshold value
    % Randomly select 1 images to verify
    rng(0);
    img_idx = randperm(1000,1);
    
    % Create input set from adversarial perturbation

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
    
    % Define input set as ImageStar
    dif_im = im - at_im;
    noise = -dif_im;
    V(:,:,:,1) = double(im);
    V(:,:,:,2) = double(noise);
    C = [1; -1];
    d = [1; de-1];
    IS = ImageStar(V, C, d, 1-de, 1);
    GrTruth = {im};
    
    % Verify network
    reachOptions.reachMethod = 'approx-star';
    t = tic;
    [riou, rv, rs, ~, ~, ~, ~, ver_rs, eval_seg_ims] = net.verify_segmentation(IS, GrTruth, reachOptions);
    t = toc(t);

    % Visualize results
    [f1, f3] = net.plot_segmentation_output_set(ver_rs{1}, eval_seg_ims{1});

    % Save results
    save("transposed_results_"+string(de)+".mat", 't', 'riou', 'rs', 'rv')
    if is_codeocean
        exportgraphics(f1,'/results/logs/transposed_results_0.0001_1.pdf', 'ContentType', 'vector');
        exportgraphics(f3,'/results/logs/transposed_results_0.0001_3.pdf', 'ContentType', 'vector');
        saveas(f1,'/results/logs/transposed_results_0.0001_1.png');
        saveas(f3,'/results/logs/transposed_results_0.0001_3.png');
    else
        exportgraphics(f1,'../transposed_results_0.0001_1.pdf','ContentType', 'vector');
        exportgraphics(f3,'../transposed_results_0.0001_3.pdf','ContentType', 'vector');
    end

end