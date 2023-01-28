function eval_cnn_all()

    % Test all CNNs
    % Load all test images
    Xall = processMNISTimages('t10k-images.idx3-ubyte');
    Yall = processMNISTlabels('t10k-labels.idx1-ubyte');
    Xall = extractdata(Xall);
    Yall = double(Yall);
    
    % Run smaller network
    acc_small = eval_cnn_small(Xall,Yall);
    % Run medium network
    acc_medium = eval_cnn_medium(Xall,Yall);
    % Run larger network
    acc_tiny = eval_cnn_tiny(Xall,Yall);
    
    % Save results
    save('eval.mat','acc_small','acc_medium','acc_tiny');

end
