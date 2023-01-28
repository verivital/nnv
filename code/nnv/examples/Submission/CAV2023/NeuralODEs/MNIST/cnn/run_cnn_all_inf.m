function run_cnn_all_inf()

    % Run all 
    pix = 784; % pixels per image to attack
    numT = 50;
    noise = [0.5/255,1/255,2/255]; % (L_inf norm)
    rng(2022); % Set random seed
    % Load all test images
    Xall = processMNISTimages('t10k-images.idx3-ubyte');
    Yall = processMNISTlabels('t10k-labels.idx1-ubyte');
    Xall = extractdata(Xall);
    Yall = double(Yall);
    
    % For a fair comparison, let's evaluate an equal number of image categories
    cat_max = numT/10; % Max images per category
    cnt_cat = zeros(10,1); % Keep track of images added
    XTest = zeros(28,28,1,numT);
    YTest = zeros(numT,1);
    ck = 1; % Image check
    adImg = 1; % Number of images added

    % Add images to run
    while sum(cnt_cat) < numT
        idx = Yall(ck);
        if cnt_cat(idx) < cat_max
            YTest(adImg) = Yall(ck);
            XTest(:,:,:,adImg) = double(Xall(:,:,:,ck));
            cnt_cat(idx) = cnt_cat(idx)+1;
            adImg = adImg + 1;
        end
        ck = ck+1;
    end
    
    cora = false; % No reach with cora, only nnv

    % Reach computation
    for noiseT = noise
        % Run smaller network
        reach_cnn_small(pix,numT,noiseT,XTest,YTest,cora,'inf');
        % Run medium network
        reach_cnn_medium(pix,numT,noiseT,XTest,YTest,cora,'inf')
        % Run tiny network
        reach_cnn_tiny(pix,numT,noiseT,XTest,YTest,cora,'inf');
    end
end
