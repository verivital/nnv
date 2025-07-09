classdef ReshapeLayer1000 < nnet.layer.Layer & nnet.layer.Formattable
    % A custom layer auto-generated while importing an ONNX network.

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
        % Specify the path to the class that will be used for codegen
        function name = matlabCodegenRedirect(~)
            name = 'cGAN_imgSz32_nCh_3.coder.ReshapeLayer1000';
        end
    end


    methods
        function this = ReshapeLayer1000(name)
            this.Name = name;
            this.OutputNames = {'input'};
        end

        function [input] = predict(this, onnx__Reshape_58)
            if isdlarray(onnx__Reshape_58)
                onnx__Reshape_58 = stripdims(onnx__Reshape_58);
            end
            onnx__Reshape_58NumDims = 2;
            onnx__Reshape_58 = cGAN_imgSz32_nCh_3.ops.permuteInputVar(onnx__Reshape_58, [2 1], 2);

            [input, inputNumDims] = ReshapeGraph1000(this, onnx__Reshape_58, onnx__Reshape_58NumDims, false);
            input = cGAN_imgSz32_nCh_3.ops.permuteOutputVar(input, [3 4 2 1], 4);

            input = dlarray(single(input), 'SSCB');
        end

        function [input] = forward(this, onnx__Reshape_58)
            if isdlarray(onnx__Reshape_58)
                onnx__Reshape_58 = stripdims(onnx__Reshape_58);
            end
            onnx__Reshape_58NumDims = 2;
            onnx__Reshape_58 = cGAN_imgSz32_nCh_3.ops.permuteInputVar(onnx__Reshape_58, [2 1], 2);

            [input, inputNumDims] = ReshapeGraph1000(this, onnx__Reshape_58, onnx__Reshape_58NumDims, true);
            input = cGAN_imgSz32_nCh_3.ops.permuteOutputVar(input, [3 4 2 1], 4);

            input = dlarray(single(input), 'SSCB');
        end

        function [input, inputNumDims1001] = ReshapeGraph1000(this, onnx__Reshape_58, onnx__Reshape_58NumDims, Training)

            % Execute the operators:
            % Reshape:
            [shape, inputNumDims] = cGAN_imgSz32_nCh_3.ops.prepareReshapeArgs(onnx__Reshape_58, this.Vars.onnx__Reshape_59, onnx__Reshape_58NumDims, 0);
            input = reshape(onnx__Reshape_58, shape{:});

            % Set graph output arguments
            inputNumDims1001 = inputNumDims;

        end

    end

end