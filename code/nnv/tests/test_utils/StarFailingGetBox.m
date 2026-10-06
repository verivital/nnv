classdef StarFailingGetBox < Star
    % StarFailingGetBox: a Star whose getBox reports an LP failure.
    %
    % Star.getBox returns [] when lpsolver ends with an exit flag other
    % than "l1" or "g5" (for example glpk "g2", a feasible but not proven
    % optimal point, on a non-empty set). This subclass reproduces that
    % outcome deterministically, so callers' handling of it can be tested
    % without depending on a particular solver failing.
    %
    % Usage: S = StarFailingGetBox(V, C, d, pred_lb, pred_ub)

    methods

        function obj = StarFailingGetBox(varargin)
            obj@Star(varargin{:});
        end

        function B = getBox(~)
            B = [];
        end

    end

end
