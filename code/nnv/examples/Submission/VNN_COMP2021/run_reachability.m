function [status, total_time] = run_reachability(onnxfile,vnnlibfile,category)%)%ip_shape(net,imagestar_set,op_specs_mat, op_specs_vec
    
    status = 2 ; %intially Unknown

    c = parcluster('local'); % build the 'local' cluster object
    numCores = c.NumWorkers;
    
    %variables: imagestar_set,op_specs_mat, op_specs_vec,ip_bounds 
    [~,vnnfile,~] = fileparts(vnnlibfile);
    vnnfile = "../nnv/code/nnv/examples/Submission/VNN_COMP2021/intermediateFiles/"+ vnnfile + ".mat";
    load(vnnfile);
    
    %variables: net for nnv format 
    [~,netfilename,~] = fileparts(onnxfile);
    netfilename = "../nnv/code/nnv/examples/Submission/VNN_COMP2021/intermediateFiles/"+netfilename + ".mat";
    load(netfilename);
    
%     if category == "acasxu" || category == "test"
%         method = 'exact-star';
%     else
%         method = 'approx-star';
%     end
    % Op Imagrstar
   [op_bounds, total_time] = net.reach(imagestar_set,'approx-star');

   for i = 1: size(op_bounds,2)
       if isempty(op_bounds(1,i).toStar.intersectHalfSpace(op_specs_mat, op_specs_vec(:,1)))%(:,j)(i,j) for exact star
            status = 1; % Safe if no intersection with the op_spec/unsafe region--- violates the counter-ex spec
       else
            status = 0; %Unknown if intersection with the op_spec/unsafe region; for approx-star-- holds/unknown for the counter-ex spec
            break
       end
   end
end