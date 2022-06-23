% Set parameters
unc = 0.01;
nnpath = "../../benchmark_dynamics/DampedOsc/ilnode(";
models = {@node0,@node1,@node2,@node3,@node4,@node5};
for dim=2:2
    path_node = nnpath + string(dim) + ")/5/model.mat";
    reach_ilnode(path_node,models{dim+1},dim, unc, true);
    disp("Model with augmented dimensions "+string(dim)+" has finished");
end