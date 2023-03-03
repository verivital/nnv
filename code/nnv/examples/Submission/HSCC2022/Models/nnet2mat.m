function nnet2mat(filefun,nntype)

% Converts NN in .nnet format to .mat



NN = GetNN(filefun);
numofNN = length(NN);

for i = 1:numofNN
    
    act_fcns = 'linear';
    diff =  length(act_fcns)- length(nntype);
    
    for k = 1:diff
        nntype = [nntype ' '];
    end
    
    W = NN{i}.W;
    b = NN{i}.B;
    numoflay = NN{i}.nlayers;
    
    for j = 1:numoflay-1
        
    act_fcns = [nntype ; act_fcns];
    
    end
    
    savedir = 'VCAS/nnv_format/';
    save(fullfile(savedir,sprintf('nn_vcas_ReLu_%d',i)),'W','b','act_fcns');
    
end

end