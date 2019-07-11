function [A,D1,D2] = mpt_scale_matrix(A)
% scales matrix A by finding D1 and D2 in An = D1*A*D2  such that infinity
% norm of each row and column belongs approaches 1
%
% find D1, D2  (in the sense minimize ||D1*A*D2|| )
% 
%  s.t. 
%     An = D1*A*D2 
%
% Scaling matrix is used in solving A*x=b for badly scaled matrix A as
% follows:
%
%                  A*x = b      
%               D1*A*x = D1*b    / multiply from left by D1
%    D1*A*D2*inv(D2)*x = D1*b    / insert D2*inv(D2) 
%    D1*A*D2*y = D1*b            / substitue y = inv(D2)*x
%         An*y = bn              / substitue An = D1*A*D2, bn = D1*b
%
% First solve An*y = bn, then obtain x = D2*y
%
% Details of the method are in file drRAL2001034.ps.gz at
% http://www.numerical.rl.ac.uk/reports/reports.html.

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% check if A has any row or column equal zero
if ~all(any(A,1)) || ~all(any(A,2))
   error('Given matrix contains a row or column with all elements equal zero.');
end

validate_realmatrix(A);

% get dimensions
[m,n] = size(A);

% initialize 
D1 = eye(m);
D2 = eye(n);
row_norm = 1;
column_norm = 1;
k=0;

% tolerance
tol = MPTOPTIONS.rel_tol;

while (row_norm>tol) && (column_norm>tol)
    nA = max(abs(A), [], 2);
    r = sqrt(nA);
    rn = 1 - nA;
    DR = diag(1./r);

    nA = max(abs(A), [], 1)';
    c = sqrt(nA);
    cn = 1-nA;
    DC = diag(1./c);
    
    % matrix updates
    A = DR * A * DC;
    D1 = D1*DR;
    D2 = D2*DC;
    
    row_norm = norm(rn, Inf);
    column_norm = norm(cn, Inf);
    
    k = k+1;
end
