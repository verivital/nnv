function preProcessing(onnxfile, vnnlibfile)
[~,vnnfile,~] = fileparts(vnnlibfile);
vnnfile = vnnfile + ".mat";
load(vnnfile);

net = parse_onnx_to_nnvnet(onnxfile);

imagestar_set = prepare_input(ip_bounds, net);

[~,netfilename,~] = fileparts(onnxfile);
netfilename = netfilename + ".mat";

save(netfilename,'net');
save(vnnfile,'imagestar_set' ,'-append');
end

%% creates star set from ip bounds for reachability analysis
function imagestar_set = prepare_input(ip_bounds, nnv_net)
    lb_i = ip_bounds(:,1);
    ub_i = ip_bounds(:,2);
    inputSize = nnv_net.Layers{1,1}.InputSize;
    
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