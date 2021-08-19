function preProcessing(varargin)
    switch nargin
        case 3
          onnxfile =  varargin{1};
          vnnlibfile = varargin{2};
          ip_shape =  varargin{3};
        case 2
          onnxfile =  varargin{1};
          vnnlibfile = varargin{2};
    end

    [~,vnnfile,~] = fileparts(vnnlibfile);
    vnnfile = "../nnv/code/nnv/examples/Submission/VNN_COMP2021/intermediateFiles/"+vnnfile + ".mat";
    load(vnnfile);

    if nargin == 3
        net = parse_onnx_to_nnvnet(onnxfile,ip_shape);
    else
        net = parse_onnx_to_nnvnet(onnxfile);
    end

    imagestar_set = prepare_input(ip_bounds, net);

    [~,netfilename,~] = fileparts(onnxfile);
    netfilename ="../nnv/code/nnv/examples/Submission/VNN_COMP2021/intermediateFiles/"+ netfilename + ".mat";

    save(netfilename,'net');
    save(vnnfile,'imagestar_set' ,'-append');
end

%% creates star set from ip bounds for reachability analysis
function imagestar_set = prepare_input(ip_bounds,net)
    lb_i = ip_bounds(:,1);
    ub_i = ip_bounds(:,2);
    inputSize = net.Layers{1,1}.InputSize;
    
     % need to change for cifar10
    if length(inputSize)==3 && inputSize(3)==1
        lb = reshape(lb_i,inputSize)';
        ub = reshape(ub_i,inputSize)'; 
    elseif length(inputSize)==3 && inputSize(3)==3
        lb_i = reshape(lb_i,[inputSize(1)*inputSize(2),inputSize(3)]);
        ub_i = reshape(ub_i,[inputSize(1)*inputSize(2),inputSize(3)]);
        for i = 1:inputSize(3)
            lb(:,:,i) = reshape(lb_i(:,i),inputSize(1),inputSize(2))';
            ub(:,:,i) = reshape(ub_i(:,i),inputSize(1),inputSize(2))';
        end
    end
    
    imagezono_set = ImageZono(lb,ub);
    imagestar_set = imagezono_set.toImageStar();
end