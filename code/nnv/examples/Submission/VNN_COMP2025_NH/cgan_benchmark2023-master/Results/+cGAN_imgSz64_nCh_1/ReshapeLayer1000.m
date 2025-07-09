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
            name = 'cGAN_imgSz64_nCh_1.coder.ReshapeLayer1000';
        end
    end


    methods
        function this = ReshapeLayer1000(name)
            this.Name = name;
            this.OutputNames = {'input'};
        end

        function [input] = predict(this, onnx__Reshape_72)
            if isdlarray(onnx__Reshape_72)
                onnx__Reshape_72 = stripdims(onnx__Reshape_72);
            end
            onnx__Reshape_72NumDims = 2;
            onnx__Reshape_72 = cGAN_imgSz64_nCh_1.ops.permuteInputVar(onnx__Reshape_72, [2 1], 2);

            [input, inputNumDims] = ReshapeGraph1000(this, onnx__Reshape_72, onnx__Reshape_72NumDims, false);
            input = cGAN_imgSz64_nCh_1.ops.permuteOutputVar(input, [3 4 2 1], 4);

            input = dlarray(single(input), 'SSCB');
        end

        function [input] = forward(this, onnx__Reshape_72)
            if isdlarray(onnx__Reshape_72)
                onnx__Reshape_72 = stripdims(onnx__Reshape_72);
            end
            onnx__Reshape_72NumDims = 2;
            onnx__Reshape_72 = cGAN_imgSz64_nCh_1.ops.permuteInputVar(onnx__Reshape_72, [2 1], 2);

            [input, inputNumDims] = ReshapeGraph1000(this, onnx__Reshape_72, onnx__Reshape_72NumDims, true);
            input = cGAN_imgSz64_nCh_1.ops.permuteOutputVar(input, [3 4 2 1], 4);

            input = dlarray(single(input), 'SSCB');
        end

        function [input, inputNumDims1001] = ReshapeGraph1000(this, onnx__Reshape_72, onnx__Reshape_72NumDims, Training)

            % Execute the operators:
            % Reshape:
            [shape, inputNumDims] = cGAN_imgSz64_nCh_1.ops.prepareReshapeArgs(onnx__Reshape_72, this.Vars.onnx__Reshape_122, onnx__Reshape_72NumDims, 0);
            input = reshape(onnx__Reshape_72, shape{:});

            % Set graph output arguments
            inputNumDims1001 = inputNumDims;

        end

    end

end