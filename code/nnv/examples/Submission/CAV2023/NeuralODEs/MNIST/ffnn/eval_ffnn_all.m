function eval_ffnn_all()

    % Test all ffnns neural ODEs 

    % Load all test images
    Xall = processMNISTimages('t10k-images.idx3-ubyte');
    Yall = processMNISTlabels('t10k-labels.idx1-ubyte');
    Xall = extractdata(Xall);
    Yall = double(Yall);

    % Run smaller network
    acc_small = eval_ffnn_small(Xall,Yall);
    % Run medium network
    acc_medium = eval_ffnn_mid(Xall,Yall);
    % Run larger network
    acc_large = eval_ffnn_large(Xall,Yall);

    % Save results
    save('eval.mat','acc_small','acc_medium','acc_large');

end
