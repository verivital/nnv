classdef FeatureInputLayer
    % Stub class to satisfy class references in NNV when ONNX support is not installed

    properties
        Name
        InputSize
    end

    methods
        function obj = FeatureInputLayer()
            obj.Name = '';
            obj.InputSize = [];
        end
    end
end
