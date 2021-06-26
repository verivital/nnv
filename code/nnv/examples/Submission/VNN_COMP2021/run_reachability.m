function [status, total_time] = run_reachability(onnxfile,vnnlibfile)%ip_shape
    status = 0 ; %intially Unknown

    c = parcluster('local'); % build the 'local' cluster object
    numCores = c.NumWorkers;
    
    %variables: imagestar_set,op_specs_mat, op_specs_vec,ip_bounds 
    [~,vnnfile,~] = fileparts(vnnlibfile);
    vnnfile = vnnfile + ".mat";
    load(vnnfile); 
    
    %variables: net for nnv format 
    [~,netfilename,~] = fileparts(onnxfile);
    netfilename = netfilename + ".mat";
    load(netfilename); 
    
    % Op Imagrstar
    [op_bounds, total_time] = net.reach(imagestar_set);

    if isempty(op_bounds.toStar.intersectHalfSpace(op_specs_mat, op_specs_vec))%(:,j)(i,j) for exact star
        status = 1; % Safe if no intersection with the op_spec/unsafe region
    else
        status = 0; %Unknown if intersection with the op_spec/unsafe region; for approx-star)
        %break
    end
endfunction [status, total_time] = run_reachability(onnxfile,vnnlibfile)%ip_shape
    status = 0 ; %intially Unknown

    c = parcluster('local'); % build the 'local' cluster object
    numCores = c.NumWorkers;
    
    %variables: imagestar_set,op_specs_mat, op_specs_vec,ip_bounds 
    [~,vnnfile,~] = fileparts(vnnlibfile);
    vnnfile = vnnfile + ".mat";
    load(vnnfile); 
    
    %variables: net for nnv format 
    [~,netfilename,~] = fileparts(onnxfile);
    netfilename = netfilename + ".mat";
    load(netfilename); 
    
    % Op Imagrstar
    [op_bounds, total_time] = net.reach(imagestar_set);

    if isempty(op_bounds.toStar.intersectHalfSpace(op_specs_mat, op_specs_vec))%(:,j)(i,j) for exact star
        status = 1; % Safe if no intersection with the op_spec/unsafe region
    else
        status = 0; %Unknown if intersection with the op_spec/unsafe region; for approx-star)
        %break
    end
end