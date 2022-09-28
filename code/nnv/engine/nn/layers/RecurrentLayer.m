classdef RecurrentLayer
    % RecurrentLayer is a vanilla recurrent layer class that contains reachability analysis method using
    % stars
    % author: Dung Tran
    % date: 5/27/2021
    
    properties
        % hidden states
        Wh; % weights_mat for hidden states
        bh; % bias vector for hidden states
        fh; % activation function for hidden nodes;
        gh = 0; % used only for LeakyReLU
        
        % output states
        Wo; % weights_mat for output states
        bo; % bias vector for output states
        fo; % activation function for output nodes
        go = 0; % used only for LeakyReLU
        
        % inputs
        Wi; % weights_mat for input states
        
        % h[t] = fh(Wh * h[t-1] + Wi * x[t] + bh)
        % y[t] = fo(Wo * h[t] + bo)
        
        nH; % number of hidden nodes
        nI; % number of inputs
        nO; % number of outputs
        
        option = []; % parallel option, 'parallel' or []
        dis_opt = []; % display option, 'display' or []
        lp_solver = 'glpk'; % lp solver option, 'linprog' or 'glpk'
        relaxFactor = 0; % use only for approx-star method
    end
    
    methods % constructor - evaluation - sampling
        % Constructor
        function obj = RecurrentLayer(varargin)
            
            if isstruct(varargin{1})
                % hidden states
                if isfield(varargin{1}, 'Wh')
                    obj.Wh = varargin{1}.Wh;
                else
                    error('Input should have Wh field');
                end
                if isfield(varargin{1}, 'bh')
                    obj.bh = varargin{1}.bh;
                else
                    error('Input should have bh field');
                end
                if isfield(varargin{1}, 'fh')
                    obj.fh = varargin{1}.fh;
                else
                    error('Input should have fh field');
                end
                if isfield(varargin{1}, 'gh')
                    obj.gh = varargin{1}.gh;
                end
                
                % output states
                if isfield(varargin{1}, 'Wo')
                    obj.Wo = varargin{1}.Wo;
                else
                    error('Input should have Wo field');
                end
                if isfield(varargin{1}, 'bo')
                    obj.bo = varargin{1}.bo;
                else
                    error('Input should have bo field');
                end
                if isfield(varargin{1}, 'fo')
                    obj.fo = varargin{1}.fo;
                else
                    error('Input should have fo field');
                end
                if isfield(varargin{1}, 'go')
                    obj.go = varargin{1}.go;
                end
                
                % input
                if isfield(varargin{1}, 'Wi')
                    obj.Wi = varargin{1}.Wi;
                else
                    error('Input should have Wi field');
                end
                
            else
                error('Input should be a struct array');
            end
            
            if size(obj.Wh, 1) ~= size(obj.bh, 1) || size(obj.Wh, 1) ~= size(obj.Wi, 1)
                error('Inconsistent dimensions between weights matrix of hidden states and its bias vector or the weight matrix of the input')
            else
                obj.nH = size(obj.Wh, 1);
            end
            
            if size(obj.Wo, 1) ~= size(obj.bo, 1)
                error('Inconsistent dimensions between weights matrix and bias vector for output states');
            else
                obj.nO = size(obj.Wo, 1);
            end
            
            if size(obj.Wo, 2) ~= size(obj.Wh, 1)
                error('Inconsistent dimension between output weight matrix and number of hidden states');
            end
            
            if (obj.gh >= 1) || (obj.go >= 1)
                error('Invalid parameter for leakyReLU, gamma should be <= 1');
            end
            
            obj.nI = size(obj.Wi, 2);
                           
        end
        
        
        
        % Evaluation method
        function y = evaluate(obj, x)  % evaluation of this layer with a specific vector
            % @x: an input sequence of length n, x[i] is the input at time
            % step i
            % @y: output sequence of length n
            % @h: hidden states sequence of length n
            
            % author: Dung Tran
            % date: 5/27/2021
            
            [m, n] = size(x); 
            if m~= obj.nI
                error('Inconsistent dimension of the input vector and the network input matrix Wi');
            end
            
            if n <= 0
                error('Invalid input sequence');
            end
            
            % hidden state
            h1 = zeros(obj.nH, n);            
            Wx = obj.Wi*x;
            for i=1:n
                if i==1
                    h1(:,i) = Wx(:,i) + obj.bh;
                    if strcmp(obj.fh, 'poslin')
                        h1(:,i) = poslin(h1(:,i));
                    elseif strcmp(obj.fh, 'purelin')
                    elseif strcmp(obj.fh, 'satlin')
                        h1(:,i) = satlin(h1(:,i));
                    elseif strcmp(obj.fh, 'satlins')
                        h1(:,i) = satlins(h1(:,i));
                    elseif strcmp(obj.fh, 'tansig')
                        h1(:,i) = tansig(h1(:,i));
                    elseif strcmp(obj.fh, 'logsig')
                        h1(:,i) = logsig(h1(:,i));
                    elseif strcmp(obj.fh, 'softmax')
                        h1(:,i) = softmax(h1(:,i));
                    elseif strcmp(obj.fh, 'leakyrelu')
                        h2 = h1(:,i);
                        h2(find(h2 < 0)) = obj.gh*h2(find(h2<0));
                        h1(:,i) = h2;
                    else
                        error('Unknown or unsupported activation function');
                    end
                else
                    h1(:,i) = obj.Wh*h1(:,i-1) + Wx(:,i) + obj.bh;
                    if strcmp(obj.fh, 'poslin')
                        h1(:,i) = poslin(h1(:,i));
                    elseif strcmp(obj.fh, 'purelin')
                    elseif strcmp(obj.fh, 'satlin')
                        h1(:,i) = satlin(h1(:,i));
                    elseif strcmp(obj.fh, 'satlins')
                        h1(:,i) = satlins(h1(:,i));
                    elseif strcmp(obj.fh, 'tansig')
                        h1(:,i) = tansig(h1(:,i));
                    elseif strcmp(obj.fh, 'logsig')
                        h1(:,i) = logsig(h1(:,i));
                    elseif strcmp(obj.fh, 'softmax')
                        h1(:,i) = softmax(h1(:,i));
                    elseif strcmp(obj.fh, 'leakyrelu')
                        h2 = h1(:,i);
                        h2(find(h2 < 0)) = obj.gh*h2(find(h2<0));
                        h1(:,i) = h2;
                    else
                        error('Unknown or unsupported activation function');
                    end
                end              
                
            end
                      
            
            y1 = obj.Wo*h1 + obj.bo;
            
            if strcmp(obj.fo, 'poslin')
                y = poslin(y1);
            elseif strcmp(obj.fo, 'purelin')
                y = y1;
            elseif strcmp(obj.fo, 'satlin')
                y = satlin(y1);
            elseif strcmp(obj.fo, 'satlins')
                y = satlins(y1);
            elseif strcmp(obj.fo, 'tansig')
                y = tansig(y1);
            elseif strcmp(obj.fo, 'logsig')
                y = logsig(y1);
            elseif strcmp(obj.fo, 'softmax')
                y = softmax(y1);
            elseif strcmp(obj.fo, 'leakyrelu')
                y = y1;
                y(find(y < 0)) = obj.go*y(find(y<0));
            else
                error('Unknown or unsupported activation function');
            end
        end
        
            
    end
    
    
    methods % reachability analysis method
        
        function O = reach(varargin)
            % @I: an array of inputs set sequence
            % @method: 'exact-star' or 'approx-star'
            % @option:  'parallel' use parallel computing
            %           '[]' or not declared -> don't use parallel
            %           computing
            % @O: a cell of output sets sequence, length(O) = length(I)
            % @H: a cell of hidden sets sequence, length(H) = length(I)
            
            % author: Dung Tran
            % date: 5/27/2021
            
            % parse inputs 
            switch nargin
                
                
                case 7
                    obj = varargin{1};
                    I = varargin{2};
                    method = varargin{3};
                    obj.option = varargin{4};
                    obj.relaxFactor = varargin{5}; % only use for approx-star method
                    obj.dis_opt = varargin{6};
                    obj.lp_solver = varargin{7};
                case 6
                    obj = varargin{1};
                    I = varargin{2};
                    method = varargin{3};
                    obj.option = varargin{4};
                    obj.relaxFactor = varargin{5}; % only use for approx-star method
                    obj.dis_opt = varargin{6};
                case 5
                    obj = varargin{1};
                    I = varargin{2};
                    method = varargin{3};
                    obj.option = varargin{4};
                    obj.relaxFactor = varargin{5}; % only use for approx-star method
                case 4
                    obj = varargin{1};
                    I = varargin{2};
                    method = varargin{3};
                    obj.option = varargin{4};
                case 3
                    obj = varargin{1};
                    I = varargin{2};
                    method = varargin{3};
                case 2 
                    obj = varargin{1};
                    I = varargin{2};
                    method = 'approx-star';
                otherwise
                    error('Invalid number of input arguments (should be 1, 2, 3, 4, 5, or 6)');
            end
            
            if ~strcmp(method, 'exact-star') && ~strcmp(method, 'approx-star') && ~contains(method, 'relax-star') && ~strcmp(method, 'abs-dom')
                error('Unknown reachability analysis method');
            end
            
            if strcmp(method, 'exact-star') && (~strcmp(obj.fh, 'purelin') && ~strcmp(obj.fh, 'leakyrelu') && ~strcmp(obj.fh, 'poslin') && ~strcmp(obj.fh, 'satlin') && ~strcmp(obj.fh, 'satlins'))
                method = 'approx-star';
                fprintf('\nThe current layer has %s activation function for hidden nodes -> cannot compute exact reachable set for the current layer, we use approx-star method instead', obj.fh);
            end
                                  
            n = length(I);
            O = cell(1, n); % output reachable set sequence
            Wh1 = obj.Wh;
            bh1 = obj.bh;
            fh1 = obj.fh;
            gh1 = obj.gh;
            Wo1 = obj.Wo;
            bo1 = obj.bo;
            fo1 = obj.fo;
            go1 = obj.go;
            Wi1 = obj.Wi;
            
            rF = obj.relaxFactor;
            dis = obj.dis_opt;
            lps = obj.lp_solver;
            WI = []; % mapped input set: WI = Wi*I + bh
            
            if strcmp(obj.option, 'parallel') % reachability analysis using star set
                
                parfor i=1:n               
                    % affine mapping y = Wx + b;
                    if isa(I(i), 'Polyhedron')
                        error('Do not accept polyhedron input set, please convert to star set');
                    else
                        WI = [WI I(i).affineMap(Wi1, bh1)];
                    end
                end
                
                H = cell(1, n);
                for i=1:n
                    if i==1
                        S1 = [];
                        if strcmp(fh1, 'purelin')
                            S1 = [S1 WI(1)];
                        elseif strcmp(fh1, 'poslin')
                            S1 = [S1 PosLin.reach(WI(1), method, [], rF, dis, lps)];
                        elseif strcmp(fh1, 'satlin')
                            S1 = [S1 SatLin.reach(WI(1), method, [], dis, lps)];
                        elseif strcmp(fh1, 'satlins')
                            S1 = [S1 SatLins.reach(WI(1), method)];
                        elseif strcmp(fh1, 'leakyrelu')
                            S1 = [S1 LeakyReLU.reach(WI(1), gh1, method, [], rF, dis, lps)];
                        elseif strcmp(fh1, 'logsig')
                            S1 = [S1 LogSig.reach(WI(1), method,[], rF, dis, lps)];
                        elseif strcmp(fh1, 'tansig')
                            S1 = [S1 TanSig.reach(WI(1), method, [], rF, dis, lps)];
                        elseif strcmp(fh1, 'softmax')
                            fprintf("\nSoftmax reachability is neglected in verification");
                            S1 = [S1 WI(1)];
                        else
                            error('Unsupported activation function, currently support purelin, poslin(ReLU), satlin, logsig, tansig');
                        end
                        H{i} = S1; % hidden state reach set
                        
                        m = length(S1);
                        M2 = [];
                        parfor j=1:m
                            M1 = S1(j).affineMap(Wo1, bo1);
                            if strcmp(fo1, 'purelin')
                                M2 = [M2 M1];
                            elseif strcmp(fo1, 'poslin')
                                M2 = [M2 PosLin.reach(M1, method, [], rF, dis, lps)];
                            elseif strcmp(fo1, 'satlin')
                                M2 = [M2 SatLin.reach(M1, method, [], dis, lps)];
                            elseif strcmp(fo1, 'satlins')
                                M2 = [M2 SatLins.reach(M1, method)];
                            elseif strcmp(fo1, 'leakyrelu')
                                M2 = [M2 LeakyReLU.reach(M1, go1, method, [], rF, dis, lps)];
                            elseif strcmp(fo1, 'logsig')
                                M2 = [M2 LogSig.reach(M1, method,[], rF, dis, lps)];
                            elseif strcmp(fo1, 'tansig')
                                M2 = [M2 TanSig.reach(M1, method, [], rF, dis, lps)];
                            elseif strcmp(fo1, 'softmax')
                                fprintf("\nSoftmax reachability is neglected in verification");
                                M2 = [M2 M1];
                            else
                                error('Unsupported activation function, currently support purelin, poslin(ReLU), satlin, logsig, tansig');
                            end
                        end
                        O{i} = M2;

                    else
                        m = length(H{i-1});
                        S2 = H{i-1};
                        S3 = [];
                        WIi = WI(i);
                        parfor j=1:m
                            S4 = S2(j).affineMap(Wh1, []);
                            S5 = S4.Sum(WIi);
                            if strcmp(fh1, 'purelin')
                                S3 = [S3 S5];
                            elseif strcmp(fh1, 'poslin')
                                S3 = [S3 PosLin.reach(S5, method, [], rF, dis, lps)];
                            elseif strcmp(fh1, 'satlin')
                                S3 = [S3 SatLin.reach(S5, method, [], dis, lps)];
                            elseif strcmp(fh1, 'satlins')
                                S3 = [S3 SatLins.reach(S5, method)];
                            elseif strcmp(fh1, 'leakyrelu')
                                S3 = [S3 LeakyReLU.reach(S5, gh1, method, [], rF, dis, lps)];
                            elseif strcmp(fh1, 'logsig')
                                S3 = [S3 LogSig.reach(S5, method,[], rF, dis, lps)];
                            elseif strcmp(fh1, 'tansig')
                                S3 = [S3 TanSig.reach(S5, method, [], rF, dis, lps)];
                            elseif strcmp(fh1, 'softmax')
                                fprintf("\nSoftmax reachability is neglected in verification");
                                S3 = [S3 S5];
                            else
                                error('Unsupported activation function, currently support purelin, poslin(ReLU), satlin, logsig, tansig');
                            end
                        end
                        H{i} = S3;
                        
                        m = length(S3);
                        M2 = [];
                        parfor j=1:m
                            M1 = S3(j).affineMap(Wo1, bo1);
                            if strcmp(fo1, 'purelin')
                                M2 = [M2 M1];
                            elseif strcmp(fo1, 'poslin')
                                M2 = [M2 PosLin.reach(M1, method, [], rF, dis, lps)];
                            elseif strcmp(fo1, 'satlin')
                                M2 = [M2 SatLin.reach(M1, method, [], dis, lps)];
                            elseif strcmp(fo1, 'satlins')
                                M2 = [M2 SatLins.reach(M1, method)];
                            elseif strcmp(fo1, 'leakyrelu')
                                M2 = [M2 LeakyReLU.reach(M1, go1, method, [], rF, dis, lps)];
                            elseif strcmp(fo1, 'logsig')
                                M2 = [M2 LogSig.reach(M1, method,[], rF, dis, lps)];
                            elseif strcmp(fo1, 'tansig')
                                M2 = [M2 TanSig.reach(M1, method, [], rF, dis, lps)];
                            elseif strcmp(fo1, 'softmax')
                                fprintf("\nSoftmax reachability is neglected in verification");
                                M2 = [M2 M1];
                            else
                                error('Unsupported activation function, currently support purelin, poslin(ReLU), satlin, logsig, tansig');
                            end
                        end
                        O{i} = M2;
                        
                        
                    end
                end
                
 
                
            else
                
                for i=1:n
                    
                    % affine mapping y = Wx + b;
                    if isa(I(i), 'Polyhedron')
                        error('Do not accept polyhedron input set, please convert to star set');
                    else
                        WI = [WI I(i).affineMap(Wi1, bh1)];
                    end
                end
                
                H = cell(1, n);
                for i=1:n
                    if i==1
                        S1 = [];
                        if strcmp(fh1, 'purelin')
                            S1 = [S1 WI(1)];
                        elseif strcmp(fh1, 'poslin')
                            S1 = [S1 PosLin.reach(WI(1), method, [], rF, dis, lps)];
                        elseif strcmp(fh1, 'satlin')
                            S1 = [S1 SatLin.reach(WI(1), method, [], dis, lps)];
                        elseif strcmp(fh1, 'satlins')
                            S1 = [S1 SatLins.reach(WI(1), method)];
                        elseif strcmp(fh1, 'leakyrelu')
                            S1 = [S1 LeakyReLU.reach(WI(1), gh1, method, [], rF, dis, lps)];
                        elseif strcmp(fh1, 'logsig')
                            S1 = [S1 LogSig.reach(WI(1), method,[], rF, dis, lps)];
                        elseif strcmp(fh1, 'tansig')
                            S1 = [S1 TanSig.reach(WI(1), method, [], rF, dis, lps)];
                        elseif strcmp(fh1, 'softmax')
                            fprintf("\nSoftmax reachability is neglected in verification");
                            S1 = [S1 WI(1)];
                        else
                            error('Unsupported activation function, currently support purelin, poslin(ReLU), satlin, logsig, tansig');
                        end
                        H{i} = S1;                       
                        m = length(S1);
                        M2 = [];
                        for j=1:m
                            M1 = S1(j).affineMap(Wo1, bo1);
                            if strcmp(fo1, 'purelin')
                                M2 = [M2 M1];
                            elseif strcmp(fo1, 'poslin')
                                M2 = [M2 PosLin.reach(M1, method, [], rF, dis, lps)];
                            elseif strcmp(fo1, 'satlin')
                                M2 = [M2 SatLin.reach(M1, method, [], dis, lps)];
                            elseif strcmp(fo1, 'satlins')
                                M2 = [M2 SatLins.reach(M1, method)];
                            elseif strcmp(fo1, 'leakyrelu')
                                M2 = [M2 LeakyReLU.reach(M1, go1, method, [], rF, dis, lps)];
                            elseif strcmp(fo1, 'logsig')
                                M2 = [M2 LogSig.reach(M1, method,[], rF, dis, lps)];
                            elseif strcmp(fo1, 'tansig')
                                M2 = [M2 TanSig.reach(M1, method, [], rF, dis, lps)];
                            elseif strcmp(fo1, 'softmax')
                                fprintf("\nSoftmax reachability is neglected in verification");
                                M2 = [M2 M1];
                            else
                                error('Unsupported activation function, currently support purelin, poslin(ReLU), satlin, logsig, tansig');
                            end
                        end
                        O{i} = M2;

                    else
                        m = length(H{i-1});
                        S2 = H{i-1};
                        S3 = [];
                        WIi = WI(i);
                        for j=1:m
                            S4 = S2(j).affineMap(Wh1, []);
                            S5 = S4.Sum(WIi);
                            if strcmp(fh1, 'purelin')
                                S3 = [S3 S5];
                            elseif strcmp(fh1, 'poslin')
                                S3 = [S3 PosLin.reach(S5, method, [], rF, dis, lps)];
                            elseif strcmp(fh1, 'satlin')
                                S3 = [S3 SatLin.reach(S5, method, [], dis, lps)];
                            elseif strcmp(fh1, 'satlins')
                                S3 = [S3 SatLins.reach(S5, method)];
                            elseif strcmp(fh1, 'leakyrelu')
                                S3 = [S3 LeakyReLU.reach(S5, gh1, method, [], rF, dis, lps)];
                            elseif strcmp(fh1, 'logsig')
                                S3 = [S3 LogSig.reach(S5, method,[], rF, dis, lps)];
                            elseif strcmp(fh1, 'tansig')
                                S3 = [S3 TanSig.reach(S5, method, [], rF, dis, lps)];
                            elseif strcmp(fh1, 'softmax')
                                fprintf("\nSoftmax reachability is neglected in verification");
                                S3 = [S3 S5];
                            else
                                error('Unsupported activation function, currently support purelin, poslin(ReLU), satlin, logsig, tansig');
                            end
                        end
                        H{i} = S3;
                        
                        m = length(S3);
                        M2 = [];
                        for j=1:m
                            M1 = S3(j).affineMap(Wo1, bo1);
                            if strcmp(fo1, 'purelin')
                                M2 = [M2 M1];
                            elseif strcmp(fo1, 'poslin')
                                M2 = [M2 PosLin.reach(M1, method, [], rF, dis, lps)];
                            elseif strcmp(fo1, 'satlin')
                                M2 = [M2 SatLin.reach(M1, method, [], dis, lps)];
                            elseif strcmp(fo1, 'satlins')
                                M2 = [M2 SatLins.reach(M1, method)];
                            elseif strcmp(fo1, 'leakyrelu')
                                M2 = [M2 LeakyReLU.reach(M1, go1, method, [], rF, dis, lps)];
                            elseif strcmp(fo1, 'logsig')
                                M2 = [M2 LogSig.reach(M1, method,[], rF, dis, lps)];
                            elseif strcmp(fo1, 'tansig')
                                M2 = [M2 TanSig.reach(M1, method, [], rF, dis, lps)];
                            elseif strcmp(fo1, 'softmax')
                                fprintf("\nSoftmax reachability is neglected in verification");
                                M2 = [M2 M1];
                            else
                                error('Unsupported activation function, currently support purelin, poslin(ReLU), satlin, logsig, tansig');
                            end
                        end
                        O{i} = M2;                        
                        
                    end
                end

            end
            
            if strcmp(method, 'approx-star') || contains(method, 'relax-star')
                O1 = [];
                for i=1:n
                    O1 = [O1 O{i}];
                end
                O = O1;
            end
            
       
        end
        
        
        
    end
    
    
    methods % flattening a layer into a sequence of operation
        
       
    end
    
    methods(Static) 
        
        function L = rand(nHiddens, nOut, nIn, actfunc)
            % @nHiddens: number of hidden nodes
            % @nOut: number of outputs
            % @nIn: number of inputs
            % @nMem: number of memory units
            % @actfunc: activation function for hidden states and output
            % states
            
            % author: Dung Tran
            % date: 5/27/2021
            
            if nHiddens <= 0 || nOut <= 0 || nIn <= 0 
                error('Invalid number of hidden nodes or outputs or inputs');
            end
            
            if ~iscell(actfunc) || length(actfunc) ~= 2
                error('Invalid activation function array');
            end
            
            S.Wh = rand(nHiddens, nHiddens);
            S.bh = -rand(nHiddens, 1);
            S.fh = actfunc{1};
            
            S.Wo = rand(nOut, nHiddens);
            S.bo = -rand(nOut, 1);
            S.fo = actfunc{2};
            
            S.Wi = rand(nHiddens, nIn);
            
            L = RecurrentLayer(S);
        end
        
        
    end
    
    
    
    
end

