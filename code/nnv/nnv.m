% main entry point for nnv when used in 
% command-line like fashion
%
% alternative execution can be built with custom 
% matlab scripts, refer to the examples folder
% for illustration
%
% general call:
% nnv(path_nn, method)
%
% example call (using default optinos):
% nnv('examples\Manual\Engine_Toy_Tansig_net.mat')
%
% example call (specifying method):
% nnv('examples\Manual\Engine_Toy_Tansig_net.mat', 'exact-star')
%
function [R,t] = nnv(varargin)
    if length(varargin) < 1
        'error: need input arg'
        return;
    end

    % base call: no method specified, no cores, use defaults
    if length(varargin) >= 1
        fp = varargin{1};
        if length(varargin) == 1
            method = 'approx-star';
            numCores = 1;
        end
    end

    % file and method    
    if length(varargin) >= 2
        method = varargin{2};
        if length(varargin) == 2
            numCores = 1;
        end
    end
    
    if length(varargin) >= 3
        if length(varargin) == 3
            numCores = varargin{3};
        end
    end

    if length(varargin) >= 4
        'error: unsupported number of inputs'
        return;
    end

    if isfile(fp)
        [filepath,filename,ext] = fileparts(fp);

        if ext == '.mat'
            load(fp);
                
            % NOTE / TODO: we assume the network is in the variable named net
            % it would be good to generalize if possible
            nnvNet = FFNNS.parse(net);
        else
            'error: unsupported network type (need to call nnmt automatically)'
        end
    else
        'error: invalid file or could not locate'
        return;
    end
    
    % TODO: add parsing input set
    lb = [-100; -200]; % lower bound vector
    ub = [100; 200]; % upper bound vector
    I = Star(lb, ub); % star input set

    [R,t] = nnvNet.reach(I, method, numCores);

end
