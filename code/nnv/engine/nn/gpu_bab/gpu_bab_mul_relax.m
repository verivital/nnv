function [aL, bL, cL, aU, bU, cU] = gpu_bab_mul_relax(xl, xu, yl, yu, alphaL, alphaU, precision)
% GPU_BAB_MUL_RELAX  McCormick linear relaxation of an elementwise product z = x.*y over the
%   box [xl,xu] x [yl,yu]. Returns per-element plane coefficients such that, for every x,y in
%   the box,
%       aL.*x + bL.*y + cL  <=  x.*y  <=  aU.*x + bU.*y + cU.
%   (GenBaB / Shi et al. 2019; see research/IMPL_ALGORITHMS_2026-06-17.md s.GenBaB.) The two
%   LOWER planes pass through corners (xl,yl) and (xu,yu); the two UPPER planes through (xl,yu)
%   and (xu,yl). CROWN must pick ONE plane per side; GenBaB interpolates the two with an
%   optimizable alpha in [0,1] (alpha=1 -> corner-solution-1, alpha=0 -> corner-solution-2).
%
%   alphaL, alphaU : optional per-element interpolation weights in [0,1] (clipped here for
%                    soundness). If EMPTY, use the Shi-2019 fixed rule: pick the tighter corner
%                    per element by comparing the two planes at the box CENTER (the standard
%                    McCormick choice) -- a sound, division-free default.
%
%   Fully ELEMENTWISE: xl/xu/yl/yu may be n-by-1 (single box) or n-by-B (batched subdomains);
%   the outputs match their shape. No division (unlike ReLU's u/(u-l)) -> no eps guard / NaN
%   risk; degenerate widths (xl==xu) are valid (the two corner planes coincide).
%
%   SOUNDNESS: each corner plane is a supporting under/over-estimator of x*y (Shi-2019 Eq.9/10);
%   a convex combination (alpha in [0,1]) of two valid under-estimators is a valid under-estimator
%   (the envelope's feasible set in (a,b,c) is convex), and likewise for over. The CALLER must
%   apply the plane sign-awarely (positive output-coeff -> LOWER plane for a spec lower bound;
%   negative -> UPPER plane), exactly like the ReLU al/au rule.

    if nargin < 5, alphaL = []; end
    if nargin < 6, alphaU = []; end
    if nargin < 7 || isempty(precision), precision = 'single'; end
    xl = cast(xl, precision); xu = cast(xu, precision);
    yl = cast(yl, precision); yu = cast(yu, precision);

    if isempty(alphaL)
        % Shi-2019 fixed lower corner: compare L1 (corner xl,yl) vs L2 (corner xu,yu) at the
        % box center; use whichever is the LARGER (tighter) under-estimator. alphaL=1 -> L1.
        xc = (xl + xu) / 2; yc = (yl + yu) / 2;
        L1 = yl .* xc + xl .* yc - xl .* yl;
        L2 = yu .* xc + xu .* yc - xu .* yu;
        alphaL = cast(L1 >= L2, precision);
    else
        alphaL = min(max(cast(alphaL, precision), 0), 1);   % CLIP to [0,1] (soundness)
    end
    if isempty(alphaU)
        xc = (xl + xu) / 2; yc = (yl + yu) / 2;
        U1 = yu .* xc + xl .* yc - xl .* yu;
        U2 = yl .* xc + xu .* yc - xu .* yl;
        alphaU = cast(U1 <= U2, precision);                 % pick the SMALLER (tighter) over-estimator
    else
        alphaU = min(max(cast(alphaU, precision), 0), 1);
    end

    % LOWER plane: convex combo of corner (xl,yl) [alpha=1] and (xu,yu) [alpha=0]
    aL = alphaL .* yl + (1 - alphaL) .* yu;
    bL = alphaL .* xl + (1 - alphaL) .* xu;
    cL = -alphaL .* (xl .* yl) - (1 - alphaL) .* (xu .* yu);
    % UPPER plane: convex combo of corner (xl,yu) [alpha=1] and (xu,yl) [alpha=0]
    aU = alphaU .* yu + (1 - alphaU) .* yl;
    bU = alphaU .* xl + (1 - alphaU) .* xu;
    cU = -alphaU .* (xl .* yu) - (1 - alphaU) .* (xu .* yl);
end
