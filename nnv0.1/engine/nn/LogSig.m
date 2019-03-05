classdef LogSig
    % LOGSIG Class contains methods for reachability analysis of layer with
    % Sigmoid activation function.
    % Reference: https://www.mathworks.com/help/deeplearning/ref/logsig.html
    % Author: Dung Tran
    % Date: 28/2/2019
    
    properties
        
    end
    
    methods(Static)  % evaluate method and over-approximate reachability analysis with stars
        
        % evaluation
        function y = evaluate(x)
            y = logsig(x);
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
                Z = [logsig('dn', lb) logsig('dn', ub)];
                gamma_opt = min(Z, [], 2);
                gamma_mat = diag(gamma_opt);
                mu1 = 0.5 * (logsig(ub) + logsig(lb) - gamma_mat * (ub + lb));
                mu2 = 0.5 * (logsig(ub) - logsig(lb) - gamma_mat * (ub - lb));
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
    
    
    methods(Static) % over-approximate reachability analysis using Zonotope
        
        function Z = reach_zono_approx(I)
            % @I: zonotope input set
            % @Z: zonotope output set
            
            % author: Dung Tran
            % date: 5/3/2019
            
            % reference: Fast and Effective Robustness Certification,
            % Gagandeep Singh, NIPS, 2018
            
            if ~isa(I, 'Zono')
                error('Input set is not a Zonotope');
            end
            
            B = I.getBox;
            
            lb = B.lb;
            ub = B.ub;
            G = [logsig('dn', lb) logsig('dn', ub)];
            gamma_opt = min(G, [], 2);
            gamma_mat = diag(gamma_opt);
            mu1 = 0.5 * (logsig(ub) + logsig(lb) - gamma_mat * (ub + lb));
            mu2 = 0.5 * (logsig(ub) - logsig(lb) - gamma_mat * (ub - lb));
            Z1 = I.affineMap(gamma_mat, mu1);
            new_V = diag(mu2);
            V = [Z1.V new_V];
            Z = Zono(Z1.c, V);
            
        end
        
    end
    
    
end

