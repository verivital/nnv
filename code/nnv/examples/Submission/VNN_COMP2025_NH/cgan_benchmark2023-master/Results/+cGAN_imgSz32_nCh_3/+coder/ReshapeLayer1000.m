classdef ReshapeLayer1000 < nnet.layer.Layer & nnet.layer.Formattable
    % A custom layer auto-generated while importing an ONNX network.
    %#codegen

    %#ok<*PROPLC>
    %#ok<*NBRAK>
    %#ok<*INUSL>
    %#ok<*VARARG>
    properties (Learnable)
    end

    properties (State)
    end

    properties
        Vars
        NumDims
    end

    methods(Static, Hidden)
        % Specify the properties of the class that will not be modified
        % after the first assignment.
        function p = matlabCodegenNontunableProperties(~)
            p = {
                % Constants, i.e., Vars, NumDims and all learnables and states
                'Vars'
                'NumDims'
                };
        end
    end


    methods(Static, Hidden)
        % Instantiate a codegenable layer instance from a MATLAB layer instance
        function this_cg = matlabCodegenToRedirected(mlInstance)
            this_cg = cGAN_imgSz32_nCh_3.coder.ReshapeLayer1000(mlInstance);
        end
        function this_ml = matlabCodegenFromRedirected(cgInstance)
            this_ml = cGAN_imgSz32_nCh_3.ReshapeLayer1000(cgInstance.Name);
            if isstruct(cgInstance.Vars)
                names = fieldnames(cgInstance.Vars);
                for i=1:numel(names)
                    fieldname = names{i};
                    this_ml.Vars.(fieldname) = dlarray(cgInstance.Vars.(fieldname));
                end
            else
                this_ml.Vars = [];
            end
            this_ml.NumDims = cgInstance.NumDims;
        end
    end

    methods
        function this = ReshapeLayer1000(mlInstance)
            this.Name = mlInstance.Name;
            this.OutputNames = {'input'};
            if isstruct(mlInstance.Vars)
                names = fieldnames(mlInstance.Vars);
                for i=1:numel(names)
                    fieldname = names{i};
                    this.Vars.(fieldname) = cGAN_imgSz32_nCh_3.coder.ops.extractIfDlarray(mlInstance.Vars.(fieldname));
                end
            else
                this.Vars = [];
            end

            this.NumDims = mlInstance.NumDims;
        end

        function [input] = predict(this, onnx__Reshape_58__)
            if isdlarray(onnx__Reshape_58__)
                onnx__Reshape_58_ = stripdims(onnx__Reshape_58__);
            else
                onnx__Reshape_58_ = onnx__Reshape_58__;
            end
            onnx__Reshape_58NumDims = 2;
            onnx__Reshape_58 = cGAN_imgSz32_nCh_3.coder.ops.permuteInputVar(onnx__Reshape_58_, [2 1], 2);

            [input__, inputNumDims__] = ReshapeGraph1000(this, onnx__Reshape_58, onnx__Reshape_58NumDims, false);
            input_ = cGAN_imgSz32_nCh_3.coder.ops.permuteOutputVar(input__, [3 4 2 1], 4);

            input = dlarray(single(input_), 'SSCB');
        end

        function [input, inputNumDims1001] = ReshapeGraph1000(this, onnx__Reshape_58, onnx__Reshape_58NumDims, Training)

            % Execute the operators:
            % Reshape:
            [shape1000, inputNumDims] = cGAN_imgSz32_nCh_3.coder.ops.prepareReshapeArgs(onnx__Reshape_58, this.Vars.onnx__Reshape_59, onnx__Reshape_58NumDims, 0);
            input = reshape(onnx__Reshape_58, shape1000{:});

            % Set graph output arguments
            inputNumDims1001 = inputNumDims;

        end

    end

end