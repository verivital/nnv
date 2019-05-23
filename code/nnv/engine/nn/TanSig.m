classdef TanSig
    % TanSig class contains method for reachability analysis for Layer with
    % Tanh activation function (Matlab called TanSig)
    % author: Dung Tran
    % date: 1/3/2019
    
    properties
    end
    
    methods(Static) % evaluate method and over-approximate reachability analysis with stars
    
        % evaluation
        function y = evaluate(x)
            y = tansig(x);
        end
        
                
        % reachability analysis with star
        function S = reach_star_approx(I)
            % @I: input star
            % @S: output star
            
            % author: Dung Tran
            % date: 1/3/2019
            
            % method: approximate sigmoid function by a zonotope
            % reference: Fast and Effective Robustness Certification,
            % Gagandeep Singh, NIPS, 2018
            
            if ~isa(I, 'Star')
                error('Input set is not a star');
            end
            
            B = I.getBox;
            if isempty(B)
                S = [];
            else
                lb = B.lb;
                ub = B.ub;
                Z = [tansig('dn', lb) tansig('dn', ub)];
                gamma_opt = min(Z, [], 2);
                gamma_mat = diag(gamma_opt);
                mu1 = 0.5 * (tansig(ub) + tansig(lb) - gamma_mat * (ub + lb));
                mu2 = 0.5 * (tansig(ub) - tansig(lb) - gamma_mat * (ub - lb));
                S1 = I.affineMap(gamma_mat, mu1);
                n = I.dim;
                new_V = diag(mu2);
                new_C = [eye(n); -eye(n)];
                new_d = ones(2*n, 1);
                
                V = [S1.V new_V];
                C = blkdiag(S1.C, new_C);
                d = [S1.d; new_d];
                
                S = Star(V, C, d);
                            
            end
                  
        end
        
    end
    
    
    methods(Static) % over-approximate reachability analysis with zonotope
        
        
        % reachability analysis with zonotope
        function Z = reach_zono_approx(I)
            % @I: input star
            % @Z: output star
            
            % author: Dung Tran
            % date: 5/3/2019
            
            % method: approximate sigmoid function by a zonotope
            % reference: Fast and Effective Robustness Certification,
            % Gagandeep Singh, NIPS, 2018
            
            if ~isa(I, 'Zono')
                error('Input set is not a Zonotope');
            end
            
            B = I.getBox;
            
            lb = B.lb;
            ub = B.ub;
            G = [tansig('dn', lb) tansig('dn', ub)];
            gamma_opt = min(G, [], 2);
            gamma_mat = diag(gamma_opt);
            mu1 = 0.5 * (tansig(ub) + tansig(lb) - gamma_mat * (ub + lb));
            mu2 = 0.5 * (tansig(ub) - tansig(lb) - gamma_mat * (ub - lb));
            Z1 = I.affineMap(gamma_mat, mu1);
            new_V = diag(mu2);
            
            V = [Z1.V new_V];
            Z = Zono(Z1.c, V);
                  
        end
        
    end
    
    
end

