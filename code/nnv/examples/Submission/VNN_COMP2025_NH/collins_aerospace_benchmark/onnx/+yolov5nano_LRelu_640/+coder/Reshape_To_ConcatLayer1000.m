classdef Reshape_To_ConcatLayer1000 < nnet.layer.Layer & nnet.layer.Formattable
    % A custom layer auto-generated while importing an ONNX network.
    %#codegen

    %#ok<*PROPLC>
    %#ok<*NBRAK>
    %#ok<*INUSL>
    %#ok<*VARARG>
    properties (Learnable)
        x_model_24_Consta_18
        x_model_24_Consta_3
        x_model_24_Consta_12
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
                'x_model_24_Consta_18'
                'x_model_24_Consta_3'
                'x_model_24_Consta_12'
                };
        end
    end


    methods(Static, Hidden)
        % Instantiate a codegenable layer instance from a MATLAB layer instance
        function this_cg = matlabCodegenToRedirected(mlInstance)
            this_cg = yolov5nano_LRelu_640.coder.Reshape_To_ConcatLayer1000(mlInstance);
        end
        function this_ml = matlabCodegenFromRedirected(cgInstance)
            this_ml = yolov5nano_LRelu_640.Reshape_To_ConcatLayer1000(cgInstance.Name);
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
            this_ml.x_model_24_Consta_18 = cgInstance.x_model_24_Consta_18;
            this_ml.x_model_24_Consta_3 = cgInstance.x_model_24_Consta_3;
            this_ml.x_model_24_Consta_12 = cgInstance.x_model_24_Consta_12;
        end
    end

    methods
        function this = Reshape_To_ConcatLayer1000(mlInstance)
            this.Name = mlInstance.Name;
            this.NumInputs = 3;
            this.OutputNames = {'output0'};
            if isstruct(mlInstance.Vars)
                names = fieldnames(mlInstance.Vars);
                for i=1:numel(names)
                    fieldname = names{i};
                    this.Vars.(fieldname) = yolov5nano_LRelu_640.coder.ops.extractIfDlarray(mlInstance.Vars.(fieldname));
                end
            else
                this.Vars = [];
            end

            this.NumDims = mlInstance.NumDims;
            this.x_model_24_Consta_18 = mlInstance.x_model_24_Consta_18;
            this.x_model_24_Consta_3 = mlInstance.x_model_24_Consta_3;
            this.x_model_24_Consta_12 = mlInstance.x_model_24_Consta_12;
        end

        function [output0] = predict(this, x_model_24_m_0_Conv___, x_model_24_m_1_Conv___, x_model_24_m_2_Conv___)
            if isdlarray(x_model_24_m_0_Conv___)
                x_model_24_m_0_Conv__ = stripdims(x_model_24_m_0_Conv___);
            else
                x_model_24_m_0_Conv__ = x_model_24_m_0_Conv___;
            end
            if isdlarray(x_model_24_m_1_Conv___)
                x_model_24_m_1_Conv__ = stripdims(x_model_24_m_1_Conv___);
            else
                x_model_24_m_1_Conv__ = x_model_24_m_1_Conv___;
            end
            if isdlarray(x_model_24_m_2_Conv___)
                x_model_24_m_2_Conv__ = stripdims(x_model_24_m_2_Conv___);
            else
                x_model_24_m_2_Conv__ = x_model_24_m_2_Conv___;
            end
            x_model_24_m_0_Conv_NumDims = 4;
            x_model_24_m_1_Conv_NumDims = 4;
            x_model_24_m_2_Conv_NumDims = 4;
            x_model_24_m_0_Conv_ = yolov5nano_LRelu_640.coder.ops.permuteInputVar(x_model_24_m_0_Conv__, [4 3 1 2], 4);
            x_model_24_m_1_Conv_ = yolov5nano_LRelu_640.coder.ops.permuteInputVar(x_model_24_m_1_Conv__, [4 3 1 2], 4);
            x_model_24_m_2_Conv_ = yolov5nano_LRelu_640.coder.ops.permuteInputVar(x_model_24_m_2_Conv__, [4 3 1 2], 4);

            [output0__, output0NumDims__] = Reshape_To_ConcatGraph1000(this, x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_, x_model_24_m_0_Conv_NumDims, x_model_24_m_1_Conv_NumDims, x_model_24_m_2_Conv_NumDims, false);
            output0_ = yolov5nano_LRelu_640.coder.ops.permuteOutputVar(output0__, ['as-is'], 3);

            output0 = dlarray(single(output0_), repmat('U', 1, max(2, output0NumDims__)));
        end

        function [output0, output0NumDims1007] = Reshape_To_ConcatGraph1000(this, x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_, x_model_24_m_0_Conv_NumDims, x_model_24_m_1_Conv_NumDims, x_model_24_m_2_Conv_NumDims, Training)

            % Execute the operators:
            % Reshape:
            [shape1000, x_model_24_Reshape_oNumDims] = yolov5nano_LRelu_640.coder.ops.prepareReshapeArgs(x_model_24_m_0_Conv_, this.Vars.x_model_24_Consta_23, x_model_24_m_0_Conv_NumDims, 0);
            x_model_24_Reshape_o = reshape(x_model_24_m_0_Conv_, shape1000{:});

            % Transpose:
            [perm1001, x_model_24_Transpo_2NumDims] = yolov5nano_LRelu_640.coder.ops.prepareTransposeArgs(this.Vars.TransposePerm1001, x_model_24_Reshape_oNumDims);
            if isempty(perm1001)
                x_model_24_Transpo_2 = x_model_24_Reshape_o;
            else
                x_model_24_Transpo_2 = permute(yolov5nano_LRelu_640.coder.ops.extractIfDlarray(x_model_24_Reshape_o), perm1001);
            end

            % Sigmoid:
            X1002 = dlarray(yolov5nano_LRelu_640.coder.ops.extractIfDlarray(x_model_24_Transpo_2));
            Y1003 = sigmoid(X1002);
            x_model_24_Sigmoid_o = yolov5nano_LRelu_640.coder.ops.extractIfDlarray(Y1003);
            x_model_24_Sigmoid_oNumDims = x_model_24_Transpo_2NumDims;

            % Split:
            [x_model_24_Split_out, x_model_24_Split_o_1, x_model_24_Split_o_2, x_model_24_Split_outNumDims, x_model_24_Split_o_1NumDims, x_model_24_Split_o_2NumDims] = yolov5nano_LRelu_640.coder.ops.onnxSplit(x_model_24_Sigmoid_o, 4, this.Vars.SplitSplit1002, 0, x_model_24_Sigmoid_oNumDims);

            % Mul:
            x_model_24_Mul_outpu = x_model_24_Split_out .* this.Vars.x_model_24_Consta_10;
            x_model_24_Mul_outpuNumDims = max(x_model_24_Split_outNumDims, this.NumDims.x_model_24_Consta_10);

            % Add:
            x_model_24_Add_outpu = x_model_24_Mul_outpu + this.Vars.x_model_24_Consta_15;
            x_model_24_Add_outpuNumDims = max(x_model_24_Mul_outpuNumDims, this.NumDims.x_model_24_Consta_15);

            % Mul:
            x_model_24_Mul_1_out = x_model_24_Add_outpu .* this.Vars.x_model_24_Consta_16;
            x_model_24_Mul_1_outNumDims = max(x_model_24_Add_outpuNumDims, this.NumDims.x_model_24_Consta_16);

            % Mul:
            x_model_24_Mul_2_out = x_model_24_Split_o_1 .* this.Vars.x_model_24_Consta_17;
            x_model_24_Mul_2_outNumDims = max(x_model_24_Split_o_1NumDims, this.NumDims.x_model_24_Consta_17);

            % Pow:
            x_model_24_Pow_outpu = power(x_model_24_Mul_2_out, this.x_model_24_Consta_18);
            x_model_24_Pow_outpuNumDims = max(x_model_24_Mul_2_outNumDims, this.NumDims.x_model_24_Consta_18);

            % Mul:
            x_model_24_Mul_3_out = x_model_24_Pow_outpu .* this.Vars.x_model_24_Consta_19;
            x_model_24_Mul_3_outNumDims = max(x_model_24_Pow_outpuNumDims, this.NumDims.x_model_24_Consta_19);

            % Concat:
            [x_model_24_Concat_ou, x_model_24_Concat_ouNumDims] = yolov5nano_LRelu_640.coder.ops.onnxConcat(4, {x_model_24_Mul_1_out, x_model_24_Mul_3_out, x_model_24_Split_o_2}, [x_model_24_Mul_1_outNumDims, x_model_24_Mul_3_outNumDims, x_model_24_Split_o_2NumDims]);

            % Reshape:
            [shape1004, x_model_24_Reshape_1NumDims] = yolov5nano_LRelu_640.coder.ops.prepareReshapeArgs(x_model_24_Concat_ou, this.Vars.x_model_24_Consta_20, x_model_24_Concat_ouNumDims, 0);
            x_model_24_Reshape_1 = reshape(x_model_24_Concat_ou, shape1004{:});

            % Reshape:
            [shape1005, x_model_24_Reshape_2NumDims] = yolov5nano_LRelu_640.coder.ops.prepareReshapeArgs(x_model_24_m_1_Conv_, this.Vars.x_model_24_Consta_21, x_model_24_m_1_Conv_NumDims, 0);
            x_model_24_Reshape_2 = reshape(x_model_24_m_1_Conv_, shape1005{:});

            % Transpose:
            [perm1006, x_model_24_TransposeNumDims] = yolov5nano_LRelu_640.coder.ops.prepareTransposeArgs(this.Vars.TransposePerm1003, x_model_24_Reshape_2NumDims);
            if isempty(perm1006)
                x_model_24_Transpose = x_model_24_Reshape_2;
            else
                x_model_24_Transpose = permute(yolov5nano_LRelu_640.coder.ops.extractIfDlarray(x_model_24_Reshape_2), perm1006);
            end

            % Sigmoid:
            X1007 = dlarray(yolov5nano_LRelu_640.coder.ops.extractIfDlarray(x_model_24_Transpose));
            Y1008 = sigmoid(X1007);
            x_model_24_Sigmoid_1 = yolov5nano_LRelu_640.coder.ops.extractIfDlarray(Y1008);
            x_model_24_Sigmoid_1NumDims = x_model_24_TransposeNumDims;

            % Split:
            [x_model_24_Split_1_o, x_model_24_Split_1_1, x_model_24_Split_1_2, x_model_24_Split_1_oNumDims, x_model_24_Split_1_1NumDims, x_model_24_Split_1_2NumDims] = yolov5nano_LRelu_640.coder.ops.onnxSplit(x_model_24_Sigmoid_1, 4, this.Vars.SplitSplit1004, 0, x_model_24_Sigmoid_1NumDims);

            % Mul:
            x_model_24_Mul_4_out = x_model_24_Split_1_o .* this.Vars.x_model_24_Consta_22;
            x_model_24_Mul_4_outNumDims = max(x_model_24_Split_1_oNumDims, this.NumDims.x_model_24_Consta_22);

            % Add:
            x_model_24_Add_1_out = x_model_24_Mul_4_out + this.Vars.x_model_24_Constant_;
            x_model_24_Add_1_outNumDims = max(x_model_24_Mul_4_outNumDims, this.NumDims.x_model_24_Constant_);

            % Mul:
            x_model_24_Mul_5_out = x_model_24_Add_1_out .* this.Vars.x_model_24_Consta_1;
            x_model_24_Mul_5_outNumDims = max(x_model_24_Add_1_outNumDims, this.NumDims.x_model_24_Consta_1);

            % Mul:
            x_model_24_Mul_6_out = x_model_24_Split_1_1 .* this.Vars.x_model_24_Consta_2;
            x_model_24_Mul_6_outNumDims = max(x_model_24_Split_1_1NumDims, this.NumDims.x_model_24_Consta_2);

            % Pow:
            x_model_24_Pow_1_out = power(x_model_24_Mul_6_out, this.x_model_24_Consta_3);
            x_model_24_Pow_1_outNumDims = max(x_model_24_Mul_6_outNumDims, this.NumDims.x_model_24_Consta_3);

            % Mul:
            x_model_24_Mul_7_out = x_model_24_Pow_1_out .* this.Vars.x_model_24_Consta_4;
            x_model_24_Mul_7_outNumDims = max(x_model_24_Pow_1_outNumDims, this.NumDims.x_model_24_Consta_4);

            % Concat:
            [x_model_24_Concat_1_, x_model_24_Concat_1_NumDims] = yolov5nano_LRelu_640.coder.ops.onnxConcat(4, {x_model_24_Mul_5_out, x_model_24_Mul_7_out, x_model_24_Split_1_2}, [x_model_24_Mul_5_outNumDims, x_model_24_Mul_7_outNumDims, x_model_24_Split_1_2NumDims]);

            % Reshape:
            [shape1009, x_model_24_Reshape_3NumDims] = yolov5nano_LRelu_640.coder.ops.prepareReshapeArgs(x_model_24_Concat_1_, this.Vars.x_model_24_Consta_5, x_model_24_Concat_1_NumDims, 0);
            x_model_24_Reshape_3 = reshape(x_model_24_Concat_1_, shape1009{:});

            % Reshape:
            [shape1010, x_model_24_Reshape_4NumDims] = yolov5nano_LRelu_640.coder.ops.prepareReshapeArgs(x_model_24_m_2_Conv_, this.Vars.x_model_24_Consta_6, x_model_24_m_2_Conv_NumDims, 0);
            x_model_24_Reshape_4 = reshape(x_model_24_m_2_Conv_, shape1010{:});

            % Transpose:
            [perm1011, x_model_24_Transpo_1NumDims] = yolov5nano_LRelu_640.coder.ops.prepareTransposeArgs(this.Vars.TransposePerm1005, x_model_24_Reshape_4NumDims);
            if isempty(perm1011)
                x_model_24_Transpo_1 = x_model_24_Reshape_4;
            else
                x_model_24_Transpo_1 = permute(yolov5nano_LRelu_640.coder.ops.extractIfDlarray(x_model_24_Reshape_4), perm1011);
            end

            % Sigmoid:
            X1012 = dlarray(yolov5nano_LRelu_640.coder.ops.extractIfDlarray(x_model_24_Transpo_1));
            Y1013 = sigmoid(X1012);
            x_model_24_Sigmoid_2 = yolov5nano_LRelu_640.coder.ops.extractIfDlarray(Y1013);
            x_model_24_Sigmoid_2NumDims = x_model_24_Transpo_1NumDims;

            % Split:
            [x_model_24_Split_2_o, x_model_24_Split_2_1, x_model_24_Split_2_2, x_model_24_Split_2_oNumDims, x_model_24_Split_2_1NumDims, x_model_24_Split_2_2NumDims] = yolov5nano_LRelu_640.coder.ops.onnxSplit(x_model_24_Sigmoid_2, 4, this.Vars.SplitSplit1006, 0, x_model_24_Sigmoid_2NumDims);

            % Mul:
            x_model_24_Mul_8_out = x_model_24_Split_2_o .* this.Vars.x_model_24_Consta_7;
            x_model_24_Mul_8_outNumDims = max(x_model_24_Split_2_oNumDims, this.NumDims.x_model_24_Consta_7);

            % Add:
            x_model_24_Add_2_out = x_model_24_Mul_8_out + this.Vars.x_model_24_Consta_8;
            x_model_24_Add_2_outNumDims = max(x_model_24_Mul_8_outNumDims, this.NumDims.x_model_24_Consta_8);

            % Mul:
            x_model_24_Mul_9_out = x_model_24_Add_2_out .* this.Vars.x_model_24_Consta_9;
            x_model_24_Mul_9_outNumDims = max(x_model_24_Add_2_outNumDims, this.NumDims.x_model_24_Consta_9);

            % Mul:
            x_model_24_Mul_10_ou = x_model_24_Split_2_1 .* this.Vars.x_model_24_Consta_11;
            x_model_24_Mul_10_ouNumDims = max(x_model_24_Split_2_1NumDims, this.NumDims.x_model_24_Consta_11);

            % Pow:
            x_model_24_Pow_2_out = power(x_model_24_Mul_10_ou, this.x_model_24_Consta_12);
            x_model_24_Pow_2_outNumDims = max(x_model_24_Mul_10_ouNumDims, this.NumDims.x_model_24_Consta_12);

            % Mul:
            x_model_24_Mul_11_ou = x_model_24_Pow_2_out .* this.Vars.x_model_24_Consta_13;
            x_model_24_Mul_11_ouNumDims = max(x_model_24_Pow_2_outNumDims, this.NumDims.x_model_24_Consta_13);

            % Concat:
            [x_model_24_Concat_2_, x_model_24_Concat_2_NumDims] = yolov5nano_LRelu_640.coder.ops.onnxConcat(4, {x_model_24_Mul_9_out, x_model_24_Mul_11_ou, x_model_24_Split_2_2}, [x_model_24_Mul_9_outNumDims, x_model_24_Mul_11_ouNumDims, x_model_24_Split_2_2NumDims]);

            % Reshape:
            [shape1014, x_model_24_Reshape_5NumDims] = yolov5nano_LRelu_640.coder.ops.prepareReshapeArgs(x_model_24_Concat_2_, this.Vars.x_model_24_Consta_14, x_model_24_Concat_2_NumDims, 0);
            x_model_24_Reshape_5 = reshape(x_model_24_Concat_2_, shape1014{:});

            % Concat:
            [output0, output0NumDims] = yolov5nano_LRelu_640.coder.ops.onnxConcat(1, {x_model_24_Reshape_1, x_model_24_Reshape_3, x_model_24_Reshape_5}, [x_model_24_Reshape_1NumDims, x_model_24_Reshape_3NumDims, x_model_24_Reshape_5NumDims]);

            % Set graph output arguments
            output0NumDims1007 = output0NumDims;

        end

    end

end