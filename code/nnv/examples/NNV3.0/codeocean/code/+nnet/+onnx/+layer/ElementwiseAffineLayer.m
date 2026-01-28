classdef ElementwiseAffineLayer
    % Stub class to satisfy class references in NNV when ONNX support is not installed
    % This is never actually used - it just prevents class resolution errors

    properties
        Name
        Scale
        Offset
        DoScale
        DoOffset
    end

    methods
        function obj = ElementwiseAffineLayer()
            obj.Name = '';
            obj.Scale = 1;
            obj.Offset = 0;
            obj.DoScale = false;
            obj.DoOffset = false;
        end
    end
end
