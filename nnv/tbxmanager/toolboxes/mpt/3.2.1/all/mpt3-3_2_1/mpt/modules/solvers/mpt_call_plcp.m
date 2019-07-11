function ret = mpt_call_plcp(opt)
%
% a gateway routine to parametric LCP solver
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

if strcmpi(opt.problem_type,'LCP')
% LCP 
    if any(eig(opt.M)<-MPTOPTIONS.abs_tol);
        error('mpt_call_plcp: PLCP solver does not solve indefinite LCP problems. Use ENUMPLCP solver instead.');    
    end

    % obtain solution
    ret = mpt_plcp(opt);

%     % exitflag
%     ret.exitflag = ret.lcpSol.exitflag;
%     % regions
%     ret.regions = ret.lcpSol.regions;
% 
%     % separate z (primal), w (dual) variables
%     ret.primal  = PolySet;
%     ret.dual    = PolySet;
% 
%     % set primal, dual variables
%     for i=1:ret.regions.numP
%         %I = ret.lcpSol.x.P(i).Data.I;
%         H = ret.lcpSol.x.P(i).user.H;
% 
%         % primal variables n+1:2*n
%         ret.primal.add(  PolyAffineFunc(ret.regions.P(i),'funcDat', H(opt.n+1:2*opt.n,:)) );
%         % dual variables 1:n
%         ret.dual.add( PolyAffineFunc(ret.regions.P(i), 'funcDat', H(1:opt.n,:)) );
%     end

    
elseif any(strcmpi(opt.problem_type,{'LP','QP'}))
% LP or QP
    % 1. Convert to LCP format
    lc = opt.copy;
    lc.qp2lcp;
    
    % 2. Solve the LCP
    ret = mpt_plcp(lc);
    
%     % 3. Convert back to original format
%     % The parameters did not change, so all the CRs are fine
%     % Just need the primal / dual values from z and w
%     regions = lcpSol.regions;
%     primal  = PolySet;
%     dual    = PolySet;
%     
%     if norm(opt.pF) < MPTOPTIONS.abs_tol
%         obj = PolySet('isConvex', true);
%     else
%         obj = PolySet;
%     end
%     
%     % Re-order variables if this came from YALMIP
%     P = speye(opt.n);
%     if ~isempty(opt.varOrder)
%         P = P(opt.varOrder.requested_variables,:);
%     end
%     for i=1:regions.numP
%         % Compute affine mapping from parameter to primal
%         x = lcpSol.x.P(i).user.H;
%         
%         T = lc.recover.uX*x + lc.recover.uTh;
%         
%         % Compute cost [x;1]'*Q*[x;1]
%         if ~isempty(opt.H)
%             Q = 0.5*T'*opt.H*T + [opt.pF opt.f]'*T;
%         else
%             Q = [opt.pF opt.f]'*T;
%         end
%         obj.add(   PolyFunc(      regions.P(i), 'func', @(x) [x(:)', 1]*Q*[x(:); 1], 'user', struct('Q', Q)));
%         
%         % Compute primal/dual variables requested by user
%         primal.add(PolyAffineFunc(regions.P(i), 'funcDat', P*T));
%         dual.add(  PolyAffineFunc(regions.P(i), 'funcDat', lc.recover.lambdaX*x + lc.recover.lambdaTh));
%         
% %         % compute primal variables
% %         Lprimal = P*T;
% %         regions.P(i).addFunction(AffFunction(Lprimal(:,1:end-1),Lprimal(:,end)),'primal');
% %         
% %         % compute dual variables
% %         Ldual = lc.recover.lambdaX*x + lc.recover.lambdaTh;
% %         regions.P(i).addFunction(AffFunction(Ldual(:,1:end-1),Ldual(:,end)),'dual');
% %         
% %         % compute the objective value
% %         Y = T(:,1:end-1);
% %         R = T(:,end);
% % 
% %         if ~isempty(opt.H)
% %             qt = 0.5*Y'*opt.H*Y + opt.pF'*Y;
% %             lt = R'*opt.H*Y + R'*opt.pF + opt.f'*Y;
% %             at = 0.5*R'*opt.H*R + opt.f'*R;
% %             regions.P(i).addFunction(QuadFunction(qt,lt,at),'obj');
% %         else
% %             lt = R'*opt.pF + opt.f'*Y;
% %             at = opt.f'*R;
% %             regions.P(i).addFunction(AffFunction(lt,at),'obj');
% %         end
%         
%     end
%     
%     ret.regions = regions;
%     ret.primal  = primal;
%     ret.dual    = dual;
%     ret.obj     = obj;
%     
%     ret.lcpSol = lcpSol;
%     ret.exitflag = lcpSol.exitflag;
    
else
    error('mpt_call_plcp: PLCP solver does not solve %s problems.',opt.problem_type);    
end