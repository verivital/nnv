function [res, time] = verify_tllverify_matlab(onnxF, vnnlibF)
    % Verify vnnlib property (vnnlibF) of onnx network (onnxF) using MATLAB
    % results = verify_tllvverify_matlab(onnxF, vnnlibF);

    % Load network
    net = importONNXNetwork(onnxF, InputDataFormats="BC");
    % Remove ElementWise and output layers
    Layers = net.Layers;
    n = length(Layers);
    good_idxs = [];
    for i=1:(n-1)
        if isa(Layers(i), 'nnet.onnx.layer.ElementwiseAffineLayer')
            Layers(i-1).Bias = Layers(i).Offset;
        else
            good_idxs = [good_idxs i];
        end
    end
    Layers = Layers(good_idxs);
    net = dlnetwork(Layers); % ready for verification
    
    % Load vnnlib property
    [XLower, XUpper, output] = load_vnnlib_matlab(vnnlibF);
    XLower = dlarray(XLower, "CB");
    XUpper = dlarray(XUpper, "CB");

    % Reachability Computation
    t = tic;
    [lb, ub] = estimateNetworkOutputBounds(net, XLower, XUpper);
    res = verifyMAT(lb, ub, output);
    time = toc(t);
end

function res = verifyMAT(lb, ub, output)
    result = ones(length(output),1);
    for i = 1:length(output)
        result(i) = eval(output{i}{1}); % Verification result (1: sat, 0: unsat or unknown)
        if ~result(i) % check it property is sat
            result(i) = eval(output{i}{2}); % if result = 1, this means we prove unsat
            if ~ result(i)
                result(i) = 2; % unknown
            else
                result(i) = 0; % unsat
                break;
            end
        end
    end % end while loop
    % check the global result of the specifications
    if all(result == 1)
        res = 1;
    elseif any(result== 0)
        res = 0;
    else
        res = 2;
    end
end

