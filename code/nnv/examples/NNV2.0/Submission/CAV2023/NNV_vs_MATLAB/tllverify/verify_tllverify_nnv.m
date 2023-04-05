function [res, time] = verify_tllverify_nnv(onnxF, vnnlibF, reachOpt)
    % Verification of tllverify using NNV
    
    % load network
    loadOpt.InputDataFormat = "BC";
    nn = onnx2nnv(onnxF, loadOpt);
    
    % load property
    property = load_vnnlib(vnnlibF);
    IS = ImageStar(property.lb, property.ub);

    % Reach
    t = tic;
    R = nn.reach(IS, reachOpt);

    % Verify
    res = verifyNNV(R, property.prop);
    time = toc(t);
end % Main function


%% Helper functions
function res = verifyNNV(R, prop)
%     nc = length(prop); % should be 1
    nr = length(R);    % number of output sets (for approx should be 1)
    prop = prop{1};
    np = length(prop);
    if np == 1 % only one halfspace
        for k = 1:nr
            Set = R(k).toStar;
            S = Set.intersectHalfSpace(prop.Hg.G, prop.Hg.g); % 
            if isempty(S)
%                 res = categorical("verified"); % unsat
                res = 0;
            elseif isempty(Set.intersectHalfSpace(-prop.Hg.G, -prop.Hg.g))
%                 res = categorical("violated"); % sat
                res = 1;
                break;
            else
%                 res = categorical("nope");    % unknown if approx, unsat if exact
                res = 2;
                break;
            end
        end
    else
        cp = 1; % current halfspace we are looking at
        res = 0; % start assuming property is unsat (no intersection)
        while cp < np % multiple halfspaces, which means OR assertion
            for k = 1:nr % check every reach set vs OR property
                Set = R(k).toStar;
                S = Set.intersectHalfSpace(prop.Hg.G, prop.Hg.g); % 
                if isempty(S)
%                     res = categorical("verified"); % unsat
                    continue; % does nothing, just need an statement, wanted to make this clear
                elseif isempty(Set.intersectHalfSpace(-prop.Hg.G, -prop.Hg.g))
%                     res = categorical("violated"); % sat
                    res = 1; % sat, it violates the property
                    cp = np;
                else
%                     res = categorical("nope");    % unknown if approx, unsat if exact
                    res = 2; %  unknown if approx, sat if exact
                    break;
                end
            end
            cp = cp+1;
        end
    end

end