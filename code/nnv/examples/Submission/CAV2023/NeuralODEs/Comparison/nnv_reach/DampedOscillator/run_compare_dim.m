function run_compare_dim()

    % Set parameters
    unc = 0.01;
    nnpath = "../../benchmark_dynamics/DampedOsc/ilnode(";
    models = {@node0,@node1};
    for dim=0:1
        path_node = nnpath + string(dim) + ")/5/model.mat";
        reach_ilnode(path_node,models{dim+1},dim, unc, true);
    end

end