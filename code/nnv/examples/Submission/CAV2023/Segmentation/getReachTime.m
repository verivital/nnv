function [relu_reachTime, others_reachTime, total_reachTime] = getReachTime(net)
    total_reachTime = sum(net.reachTime);
    n = length(net.Layers);
    relu_reachTime = 0;
    others_reachTime = 0; 
    for i=1:n
        fprintf("\nHi %d", i);
        if isa(net.Layers{i}, 'ReluLayer')
            fprintf("\nHi there");
            relu_reachTime = relu_reachTime + net.reachTime(i);
        else
            others_reachTime = others_reachTime + net.reachTime(i);
        end
    end
end

