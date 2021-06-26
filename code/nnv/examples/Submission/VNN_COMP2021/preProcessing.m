function preProcessing(onnxfile, vnnlibfile)
[~,vnnfile,~] = fileparts(vnnlibfile);
vnnfile = vnnfile + ".mat";
load(vnnfile);

net = parse_onnx_to_nnvnet(onnxfile);

imagestar_set = prepare_input(ip_bounds, net);

[~,netfilename,~] = fileparts(onnxfile);
netfilename = netfilename + ".mat";

save(netfilename,net,imagestar_set,ip_bounds,op;
save(vnnfile,'imagestar_set' ,'-append');
end

%% creates star set from ip bounds for reachability analysis
function imagestar_set = prepare_input(ip_bounds, nnv_net)
    lb = ip_bounds(:,1);
    ub = ip_bounds(:,2);
    inputSize = nnv_net.Layers{1,1}.InputSize;
    
     % need to change for cifar10
    if length(inputSize)==3 && inputSize(3)==1
        lb = reshape(lb,inputSize)';
        ub = reshape(ub,inputSize)'; 
    end
    
    imagezono_set = ImageZono(lb,ub);
    imagestar_set = imagezono_set.toImageStar();
end