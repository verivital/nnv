function run_mnist_inf()

    % Parameters
    pix = 784; % pixels per image to attack
    numT = 5;
    noise = 0.5/255; % (L_inf norm)
    rng(2022); % Set random seed

    % Load all test images
    Xall = processMNISTimages('t10k-images.idx3-ubyte');
    Yall = processMNISTlabels('t10k-labels.idx1-ubyte');
    Xall = single(extractdata(Xall));
    Yall = single(Yall);
    
    % Add images to run
    XTest = zeros(28,28,1,numT, 'single');
    YTest = zeros(numT,1, 'single');
    for ck=1:numT
        YTest(ck) = Yall(ck);
        XTest(:,:,:,ck) = Xall(:,:,:,ck);
    end

    % Reach computation
    for noiseT = noise
        % Run medium network
        reach_cnn_medium(pix,numT,noiseT,XTest,YTest,'inf')
        % Run tiny network
        reach_cnn_tiny(pix,numT,noiseT,XTest,YTest,'inf');
    end

end
