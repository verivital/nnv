function example_transposed()
    % Run a segmentation example using a segnet with transposed convolution
    
    % Load network
    net = load("m2nist_75iou_transposedcnn_avgpool.mat");
    net = SEGNET.parse(net.net, 'm2nist_75iou_transposedcnn_avgpool');
    
    % Load images
    images = load('m2nist_6484_test_images.mat');
    im_data = images.im_data;
    
    % Create example input set
    Nmax = 50; % maximum allowable number of attacked pixels
    de = 0.0001;
    Nt = 150;
    % Randomly select 1 images to verify
    rng(0);
    img_idxs = randperm(1000,1);
    
    %% create input set
    N1 = length(de);  
    
    IS(N1) = ImageStar;
    GrTruth = cell(1,N1);
    for l=1:N1
        ct = 0;
        flag = 0;
        im = im_data(:,:,img_idxs(l));
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
    
        dif_im = im - at_im;
        noise = -dif_im;
        V(:,:,:,1) = double(im);
        V(:,:,:,2) = double(noise);
        C = [1; -1];
        d = [1; de(l)-1];
        S = ImageStar(V, C, d, 1-de(l), 1);
        IS(l) = S; 
        GrTruth{l} = {im};
    end
    
    for i=1:N1
        % Verify network
        t = tic;
        [riou, rv, rs, ~, ~, ~, ~, ~] = net.verify(IS(i), GrTruth{1,i}, "approx-star");
        t = toc(t);
        % Visualize results
        plot_segmentation_output_set(net, 1, "transposed_results_"+string(de(i)))
        % Save results
        save("transposed_results_"+string(de(i))+".mat", 't', 'riou', 'rs', 'rv')
    end

end