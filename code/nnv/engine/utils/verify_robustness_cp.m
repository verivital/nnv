function result = verify_robustness_cp(net,inputSet, reachOptions, target, numClasses)
    % net: dlnetwork
    % Compute reachable set
    R = Prob_reach(net, inputSet, reachOptions);
    % Create spec
    spec = create_robustness_spec(target, numClasses);
    % Check robustness
    result = verify_specification(R,spec);
end

function spec = create_robustness_spec(target, outputSize)
    % Define HalfSpace Matrix and vector
    G = ones(outputSize,1);
    G = diag(G);
    G(target, :) = [];
    G = -G;
    G(:, target) = 1;


    % Create HalfSapce to define robustness specification
    Hs = [];
    for i=1:height(G)
        Hs = [Hs; HalfSpace(G(i,:), 0)];
    end

    spec = Hs;

end

