classdef Reshape_To_ConcatLayer1159 < nnet.layer.Layer & nnet.layer.Formattable
    % A custom layer auto-generated while importing an ONNX network.
    
    %#codegen
    %#ok<*PROPLC>
    %#ok<*NBRAK>
    %#ok<*INUSL>
    %#ok<*VARARG>
    
    properties (Learnable)
        x_model_24_Consta_12
        x_model_24_Consta_18
        x_model_24_Consta_3
    end
    
    properties
        ONNXParams         % An ONNXParameters object containing parameters used by this layer.
    end
    
    methods
        function this = Reshape_To_ConcatLayer1159(name, onnxParams)
            this.Name = name;
            this.NumInputs = 3;
            this.OutputNames = {'output0'};
            this.ONNXParams = onnxParams;
            this.x_model_24_Consta_12 = onnxParams.Learnables.x_model_24_Consta_12;
            this.x_model_24_Consta_18 = onnxParams.Learnables.x_model_24_Consta_18;
            this.x_model_24_Consta_3 = onnxParams.Learnables.x_model_24_Consta_3;
        end
        
        function [output0] = predict(this, x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_)
            if isdlarray(x_model_24_m_0_Conv_)
                x_model_24_m_0_Conv_ = stripdims(x_model_24_m_0_Conv_);
            end
            if isdlarray(x_model_24_m_1_Conv_)
                x_model_24_m_1_Conv_ = stripdims(x_model_24_m_1_Conv_);
            end
            if isdlarray(x_model_24_m_2_Conv_)
                x_model_24_m_2_Conv_ = stripdims(x_model_24_m_2_Conv_);
            end
            x_model_24_m_0_Conv_NumDims = 4;
            x_model_24_m_1_Conv_NumDims = 4;
            x_model_24_m_2_Conv_NumDims = 4;
            onnxParams = this.ONNXParams;
            onnxParams.Learnables.x_model_24_Consta_12 = this.x_model_24_Consta_12;
            onnxParams.Learnables.x_model_24_Consta_18 = this.x_model_24_Consta_18;
            onnxParams.Learnables.x_model_24_Consta_3 = this.x_model_24_Consta_3;
            [output0, output0NumDims] = Reshape_To_ConcatFcn(x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_, x_model_24_m_0_Conv_NumDims, x_model_24_m_1_Conv_NumDims, x_model_24_m_2_Conv_NumDims, onnxParams, 'Training', false, ...
                'InputDataPermutation', {[4 3 1 2], [4 3 1 2], [4 3 1 2], ['as-is'], ['as-is'], ['as-is']}, ...
                'OutputDataPermutation', {['as-is'], ['as-is']});
            if any(cellfun(@(A)~isnumeric(A), {output0}))
                fprintf('Runtime error in network. The custom layer ''%s'' output a non-numeric value.\n', 'Reshape_To_ConcatLayer1159');
                error(message('nnet_cnn_onnx:onnx:BadCustomLayerRuntimeOutput', 'Reshape_To_ConcatLayer1159'));
            end
            output0 = dlarray(single(output0), repmat('U', 1, max(2, output0NumDims)));
            if ~coder.target('MATLAB')
                output0 = extractdata(output0);
            end
        end
        
        function [output0] = forward(this, x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_)
            if isdlarray(x_model_24_m_0_Conv_)
                x_model_24_m_0_Conv_ = stripdims(x_model_24_m_0_Conv_);
            end
            if isdlarray(x_model_24_m_1_Conv_)
                x_model_24_m_1_Conv_ = stripdims(x_model_24_m_1_Conv_);
            end
            if isdlarray(x_model_24_m_2_Conv_)
                x_model_24_m_2_Conv_ = stripdims(x_model_24_m_2_Conv_);
            end
            x_model_24_m_0_Conv_NumDims = 4;
            x_model_24_m_1_Conv_NumDims = 4;
            x_model_24_m_2_Conv_NumDims = 4;
            onnxParams = this.ONNXParams;
            onnxParams.Learnables.x_model_24_Consta_12 = this.x_model_24_Consta_12;
            onnxParams.Learnables.x_model_24_Consta_18 = this.x_model_24_Consta_18;
            onnxParams.Learnables.x_model_24_Consta_3 = this.x_model_24_Consta_3;
            [output0, output0NumDims] = Reshape_To_ConcatFcn(x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_, x_model_24_m_0_Conv_NumDims, x_model_24_m_1_Conv_NumDims, x_model_24_m_2_Conv_NumDims, onnxParams, 'Training', true, ...
                'InputDataPermutation', {[4 3 1 2], [4 3 1 2], [4 3 1 2], ['as-is'], ['as-is'], ['as-is']}, ...
                'OutputDataPermutation', {['as-is'], ['as-is']});
            if any(cellfun(@(A)~isnumeric(A), {output0}))
                fprintf('Runtime error in network. The custom layer ''%s'' output a non-numeric value.\n', 'Reshape_To_ConcatLayer1159');
                error(message('nnet_cnn_onnx:onnx:BadCustomLayerRuntimeOutput', 'Reshape_To_ConcatLayer1159'));
            end
            output0 = dlarray(single(output0), repmat('U', 1, max(2, output0NumDims)));
            if ~coder.target('MATLAB')
                output0 = extractdata(output0);
            end
        end
    end
end

function [output0, output0NumDims, state] = Reshape_To_ConcatFcn(x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_, x_model_24_m_0_Conv_NumDims, x_model_24_m_1_Conv_NumDims, x_model_24_m_2_Conv_NumDims, params, varargin)
%RESHAPE_TO_CONCATFCN Function implementing an imported ONNX network.
%
% THIS FILE WAS AUTO-GENERATED BY importONNXFunction.
% ONNX Operator Set Version: 12
%
% Variable names in this function are taken from the original ONNX file.
%
% [OUTPUT0] = Reshape_To_ConcatFcn(X_MODEL_24_M_0_CONV_, X_MODEL_24_M_1_CONV_, X_MODEL_24_M_2_CONV_, PARAMS)
%			- Evaluates the imported ONNX network RESHAPE_TO_CONCATFCN with input(s)
%			X_MODEL_24_M_0_CONV_, X_MODEL_24_M_1_CONV_, X_MODEL_24_M_2_CONV_ and the imported network parameters in PARAMS. Returns
%			network output(s) in OUTPUT0.
%
% [OUTPUT0, STATE] = Reshape_To_ConcatFcn(X_MODEL_24_M_0_CONV_, X_MODEL_24_M_1_CONV_, X_MODEL_24_M_2_CONV_, PARAMS)
%			- Additionally returns state variables in STATE. When training,
%			use this form and set TRAINING to true.
%
% [__] = Reshape_To_ConcatFcn(X_MODEL_24_M_0_CONV_, X_MODEL_24_M_1_CONV_, X_MODEL_24_M_2_CONV_, PARAMS, 'NAME1', VAL1, 'NAME2', VAL2, ...)
%			- Specifies additional name-value pairs described below:
%
% 'Training'
% 			Boolean indicating whether the network is being evaluated for
%			prediction or training. If TRAINING is true, state variables
%			will be updated.
%
% 'InputDataPermutation'
%			'auto' - Automatically attempt to determine the permutation
%			 between the dimensions of the input data and the dimensions of
%			the ONNX model input. For example, the permutation from HWCN
%			(MATLAB standard) to NCHW (ONNX standard) uses the vector
%			[4 3 1 2]. See the documentation for IMPORTONNXFUNCTION for
%			more information about automatic permutation.
%
%			'none' - Input(s) are passed in the ONNX model format. See 'Inputs'.
%
%			numeric vector - The permutation vector describing the
%			transformation between input data dimensions and the expected
%			ONNX input dimensions.%
%			cell array - If the network has multiple inputs, each cell
%			contains 'auto', 'none', or a numeric vector.
%
% 'OutputDataPermutation'
%			'auto' - Automatically attempt to determine the permutation
%			between the dimensions of the output and a conventional MATLAB
%			dimension ordering. For example, the permutation from NC (ONNX
%			standard) to CN (MATLAB standard) uses the vector [2 1]. See
%			the documentation for IMPORTONNXFUNCTION for more information
%			about automatic permutation.
%
%			'none' - Return output(s) as given by the ONNX model. See 'Outputs'.
%
%			numeric vector - The permutation vector describing the
%			transformation between the ONNX output dimensions and the
%			desired output dimensions.%
%			cell array - If the network has multiple outputs, each cell
%			contains 'auto', 'none' or a numeric vector.
%
% Inputs:
% -------
% X_MODEL_24_M_0_CONV_, X_MODEL_24_M_1_CONV_, X_MODEL_24_M_2_CONV_
%			- Input(s) to the ONNX network.
%			  The input size(s) expected by the ONNX file are:
%				  X_MODEL_24_M_0_CONV_:		[Unknown, Unknown, Unknown, Unknown]				Type: FLOAT
%				  X_MODEL_24_M_1_CONV_:		[Unknown, Unknown, Unknown, Unknown]				Type: FLOAT
%				  X_MODEL_24_M_2_CONV_:		[Unknown, Unknown, Unknown, Unknown]				Type: FLOAT
%			  By default, the function will try to permute the input(s)
%			  into this dimension ordering. If the default is incorrect,
%			  use the 'InputDataPermutation' argument to control the
%			  permutation.
%
%
% PARAMS	- Network parameters returned by 'importONNXFunction'.
%
%
% Outputs:
% --------
% OUTPUT0
%			- Output(s) of the ONNX network.
%			  Without permutation, the size(s) of the outputs are:
%				  OUTPUT0:		[1, 25200, 11]				Type: FLOAT
%			  By default, the function will try to permute the output(s)
%			  from this dimension ordering into a conventional MATLAB
%			  ordering. If the default is incorrect, use the
%			  'OutputDataPermutation' argument to control the permutation.
%
% STATE		- (Optional) State variables. When TRAINING is true, these will
% 			  have been updated from the original values in PARAMS.State.
%
%
%  See also importONNXFunction

% Preprocess the input data and arguments:
[x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_, Training, outputDataPerms, anyDlarrayInputs] = preprocessInput(x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_, params, varargin{:});
% Put all variables into a single struct to implement dynamic scoping:
[Vars, NumDims] = packageVariables(params, {'x_model_24_m_0_Conv_', 'x_model_24_m_1_Conv_', 'x_model_24_m_2_Conv_'}, {x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_}, [x_model_24_m_0_Conv_NumDims x_model_24_m_1_Conv_NumDims x_model_24_m_2_Conv_NumDims]);
% Call the top-level graph function:
[output0, output0NumDims, state] = Reshape_To_ConcatGraph1148(x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_, NumDims.x_model_24_m_0_Conv_, NumDims.x_model_24_m_1_Conv_, NumDims.x_model_24_m_2_Conv_, Vars, NumDims, Training, params.State);
% Postprocess the output data
[output0] = postprocessOutput(output0, outputDataPerms, anyDlarrayInputs, Training, varargin{:});
end

function [output0, output0NumDims1158, state] = Reshape_To_ConcatGraph1148(x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_, x_model_24_m_0_Conv_NumDims1155, x_model_24_m_1_Conv_NumDims1156, x_model_24_m_2_Conv_NumDims1157, Vars, NumDims, Training, state)
% Function implementing the graph 'Reshape_To_ConcatGraph1148'
% Update Vars and NumDims from the graph's formal input parameters. Note that state variables are already in Vars.
Vars.x_model_24_m_0_Conv_ = x_model_24_m_0_Conv_;
NumDims.x_model_24_m_0_Conv_ = x_model_24_m_0_Conv_NumDims1155;
Vars.x_model_24_m_1_Conv_ = x_model_24_m_1_Conv_;
NumDims.x_model_24_m_1_Conv_ = x_model_24_m_1_Conv_NumDims1156;
Vars.x_model_24_m_2_Conv_ = x_model_24_m_2_Conv_;
NumDims.x_model_24_m_2_Conv_ = x_model_24_m_2_Conv_NumDims1157;

% Execute the operators:
% Reshape:
[shape, NumDims.x_model_24_Reshape_o] = prepareReshapeArgs(Vars.x_model_24_m_0_Conv_, Vars.x_model_24_Consta_23, NumDims.x_model_24_m_0_Conv_, 0);
Vars.x_model_24_Reshape_o = reshape(Vars.x_model_24_m_0_Conv_, shape{:});

% Transpose:
[perm, NumDims.x_model_24_Transpo_2] = prepareTransposeArgs(Vars.TransposePerm1149, NumDims.x_model_24_Reshape_o);
if ~isempty(perm)
    Vars.x_model_24_Transpo_2 = permute(Vars.x_model_24_Reshape_o, perm);
end

% Sigmoid:
Vars.x_model_24_Sigmoid_o = sigmoid(Vars.x_model_24_Transpo_2);
NumDims.x_model_24_Sigmoid_o = NumDims.x_model_24_Transpo_2;

% Split:
[Vars.x_model_24_Split_out, Vars.x_model_24_Split_o_1, Vars.x_model_24_Split_o_2, NumDims.x_model_24_Split_out, NumDims.x_model_24_Split_o_1, NumDims.x_model_24_Split_o_2] = onnxSplit(Vars.x_model_24_Sigmoid_o, 4, Vars.SplitSplit1150, 0, NumDims.x_model_24_Sigmoid_o);

% Mul:
Vars.x_model_24_Mul_outpu = Vars.x_model_24_Split_out .* Vars.x_model_24_Consta_10;
NumDims.x_model_24_Mul_outpu = max(NumDims.x_model_24_Split_out, NumDims.x_model_24_Consta_10);

% Add:
Vars.x_model_24_Add_outpu = Vars.x_model_24_Mul_outpu + Vars.x_model_24_Consta_15;
NumDims.x_model_24_Add_outpu = max(NumDims.x_model_24_Mul_outpu, NumDims.x_model_24_Consta_15);

% Mul:
Vars.x_model_24_Mul_1_out = Vars.x_model_24_Add_outpu .* Vars.x_model_24_Consta_16;
NumDims.x_model_24_Mul_1_out = max(NumDims.x_model_24_Add_outpu, NumDims.x_model_24_Consta_16);

% Mul:
Vars.x_model_24_Mul_2_out = Vars.x_model_24_Split_o_1 .* Vars.x_model_24_Consta_17;
NumDims.x_model_24_Mul_2_out = max(NumDims.x_model_24_Split_o_1, NumDims.x_model_24_Consta_17);

% Pow:
Vars.x_model_24_Pow_outpu = power(Vars.x_model_24_Mul_2_out, Vars.x_model_24_Consta_18);
NumDims.x_model_24_Pow_outpu = max(NumDims.x_model_24_Mul_2_out, NumDims.x_model_24_Consta_18);

% Mul:
Vars.x_model_24_Mul_3_out = Vars.x_model_24_Pow_outpu .* Vars.x_model_24_Consta_19;
NumDims.x_model_24_Mul_3_out = max(NumDims.x_model_24_Pow_outpu, NumDims.x_model_24_Consta_19);

% Concat:
[Vars.x_model_24_Concat_ou, NumDims.x_model_24_Concat_ou] = onnxConcat(4, {Vars.x_model_24_Mul_1_out, Vars.x_model_24_Mul_3_out, Vars.x_model_24_Split_o_2}, [NumDims.x_model_24_Mul_1_out, NumDims.x_model_24_Mul_3_out, NumDims.x_model_24_Split_o_2]);

% Reshape:
[shape, NumDims.x_model_24_Reshape_1] = prepareReshapeArgs(Vars.x_model_24_Concat_ou, Vars.x_model_24_Consta_20, NumDims.x_model_24_Concat_ou, 0);
Vars.x_model_24_Reshape_1 = reshape(Vars.x_model_24_Concat_ou, shape{:});

% Reshape:
[shape, NumDims.x_model_24_Reshape_2] = prepareReshapeArgs(Vars.x_model_24_m_1_Conv_, Vars.x_model_24_Consta_21, NumDims.x_model_24_m_1_Conv_, 0);
Vars.x_model_24_Reshape_2 = reshape(Vars.x_model_24_m_1_Conv_, shape{:});

% Transpose:
[perm, NumDims.x_model_24_Transpose] = prepareTransposeArgs(Vars.TransposePerm1151, NumDims.x_model_24_Reshape_2);
if ~isempty(perm)
    Vars.x_model_24_Transpose = permute(Vars.x_model_24_Reshape_2, perm);
end

% Sigmoid:
Vars.x_model_24_Sigmoid_1 = sigmoid(Vars.x_model_24_Transpose);
NumDims.x_model_24_Sigmoid_1 = NumDims.x_model_24_Transpose;

% Split:
[Vars.x_model_24_Split_1_o, Vars.x_model_24_Split_1_1, Vars.x_model_24_Split_1_2, NumDims.x_model_24_Split_1_o, NumDims.x_model_24_Split_1_1, NumDims.x_model_24_Split_1_2] = onnxSplit(Vars.x_model_24_Sigmoid_1, 4, Vars.SplitSplit1152, 0, NumDims.x_model_24_Sigmoid_1);

% Mul:
Vars.x_model_24_Mul_4_out = Vars.x_model_24_Split_1_o .* Vars.x_model_24_Consta_22;
NumDims.x_model_24_Mul_4_out = max(NumDims.x_model_24_Split_1_o, NumDims.x_model_24_Consta_22);

% Add:
Vars.x_model_24_Add_1_out = Vars.x_model_24_Mul_4_out + Vars.x_model_24_Constant_;
NumDims.x_model_24_Add_1_out = max(NumDims.x_model_24_Mul_4_out, NumDims.x_model_24_Constant_);

% Mul:
Vars.x_model_24_Mul_5_out = Vars.x_model_24_Add_1_out .* Vars.x_model_24_Consta_1;
NumDims.x_model_24_Mul_5_out = max(NumDims.x_model_24_Add_1_out, NumDims.x_model_24_Consta_1);

% Mul:
Vars.x_model_24_Mul_6_out = Vars.x_model_24_Split_1_1 .* Vars.x_model_24_Consta_2;
NumDims.x_model_24_Mul_6_out = max(NumDims.x_model_24_Split_1_1, NumDims.x_model_24_Consta_2);

% Pow:
Vars.x_model_24_Pow_1_out = power(Vars.x_model_24_Mul_6_out, Vars.x_model_24_Consta_3);
NumDims.x_model_24_Pow_1_out = max(NumDims.x_model_24_Mul_6_out, NumDims.x_model_24_Consta_3);

% Mul:
Vars.x_model_24_Mul_7_out = Vars.x_model_24_Pow_1_out .* Vars.x_model_24_Consta_4;
NumDims.x_model_24_Mul_7_out = max(NumDims.x_model_24_Pow_1_out, NumDims.x_model_24_Consta_4);

% Concat:
[Vars.x_model_24_Concat_1_, NumDims.x_model_24_Concat_1_] = onnxConcat(4, {Vars.x_model_24_Mul_5_out, Vars.x_model_24_Mul_7_out, Vars.x_model_24_Split_1_2}, [NumDims.x_model_24_Mul_5_out, NumDims.x_model_24_Mul_7_out, NumDims.x_model_24_Split_1_2]);

% Reshape:
[shape, NumDims.x_model_24_Reshape_3] = prepareReshapeArgs(Vars.x_model_24_Concat_1_, Vars.x_model_24_Consta_5, NumDims.x_model_24_Concat_1_, 0);
Vars.x_model_24_Reshape_3 = reshape(Vars.x_model_24_Concat_1_, shape{:});

% Reshape:
[shape, NumDims.x_model_24_Reshape_4] = prepareReshapeArgs(Vars.x_model_24_m_2_Conv_, Vars.x_model_24_Consta_6, NumDims.x_model_24_m_2_Conv_, 0);
Vars.x_model_24_Reshape_4 = reshape(Vars.x_model_24_m_2_Conv_, shape{:});

% Transpose:
[perm, NumDims.x_model_24_Transpo_1] = prepareTransposeArgs(Vars.TransposePerm1153, NumDims.x_model_24_Reshape_4);
if ~isempty(perm)
    Vars.x_model_24_Transpo_1 = permute(Vars.x_model_24_Reshape_4, perm);
end

% Sigmoid:
Vars.x_model_24_Sigmoid_2 = sigmoid(Vars.x_model_24_Transpo_1);
NumDims.x_model_24_Sigmoid_2 = NumDims.x_model_24_Transpo_1;

% Split:
[Vars.x_model_24_Split_2_o, Vars.x_model_24_Split_2_1, Vars.x_model_24_Split_2_2, NumDims.x_model_24_Split_2_o, NumDims.x_model_24_Split_2_1, NumDims.x_model_24_Split_2_2] = onnxSplit(Vars.x_model_24_Sigmoid_2, 4, Vars.SplitSplit1154, 0, NumDims.x_model_24_Sigmoid_2);

% Mul:
Vars.x_model_24_Mul_8_out = Vars.x_model_24_Split_2_o .* Vars.x_model_24_Consta_7;
NumDims.x_model_24_Mul_8_out = max(NumDims.x_model_24_Split_2_o, NumDims.x_model_24_Consta_7);

% Add:
Vars.x_model_24_Add_2_out = Vars.x_model_24_Mul_8_out + Vars.x_model_24_Consta_8;
NumDims.x_model_24_Add_2_out = max(NumDims.x_model_24_Mul_8_out, NumDims.x_model_24_Consta_8);

% Mul:
Vars.x_model_24_Mul_9_out = Vars.x_model_24_Add_2_out .* Vars.x_model_24_Consta_9;
NumDims.x_model_24_Mul_9_out = max(NumDims.x_model_24_Add_2_out, NumDims.x_model_24_Consta_9);

% Mul:
Vars.x_model_24_Mul_10_ou = Vars.x_model_24_Split_2_1 .* Vars.x_model_24_Consta_11;
NumDims.x_model_24_Mul_10_ou = max(NumDims.x_model_24_Split_2_1, NumDims.x_model_24_Consta_11);

% Pow:
Vars.x_model_24_Pow_2_out = power(Vars.x_model_24_Mul_10_ou, Vars.x_model_24_Consta_12);
NumDims.x_model_24_Pow_2_out = max(NumDims.x_model_24_Mul_10_ou, NumDims.x_model_24_Consta_12);

% Mul:
Vars.x_model_24_Mul_11_ou = Vars.x_model_24_Pow_2_out .* Vars.x_model_24_Consta_13;
NumDims.x_model_24_Mul_11_ou = max(NumDims.x_model_24_Pow_2_out, NumDims.x_model_24_Consta_13);

% Concat:
[Vars.x_model_24_Concat_2_, NumDims.x_model_24_Concat_2_] = onnxConcat(4, {Vars.x_model_24_Mul_9_out, Vars.x_model_24_Mul_11_ou, Vars.x_model_24_Split_2_2}, [NumDims.x_model_24_Mul_9_out, NumDims.x_model_24_Mul_11_ou, NumDims.x_model_24_Split_2_2]);

% Reshape:
[shape, NumDims.x_model_24_Reshape_5] = prepareReshapeArgs(Vars.x_model_24_Concat_2_, Vars.x_model_24_Consta_14, NumDims.x_model_24_Concat_2_, 0);
Vars.x_model_24_Reshape_5 = reshape(Vars.x_model_24_Concat_2_, shape{:});

% Concat:
[Vars.output0, NumDims.output0] = onnxConcat(1, {Vars.x_model_24_Reshape_1, Vars.x_model_24_Reshape_3, Vars.x_model_24_Reshape_5}, [NumDims.x_model_24_Reshape_1, NumDims.x_model_24_Reshape_3, NumDims.x_model_24_Reshape_5]);

% Set graph output arguments from Vars and NumDims:
output0 = Vars.output0;
output0NumDims1158 = NumDims.output0;
% Set output state from Vars:
state = updateStruct(state, Vars);
end

function [inputDataPerms, outputDataPerms, Training] = parseInputs(x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_, numDataOutputs, params, varargin)
% Function to validate inputs to Reshape_To_ConcatFcn:
p = inputParser;
isValidArrayInput = @(x)isnumeric(x) || isstring(x);
isValidONNXParameters = @(x)isa(x, 'ONNXParameters');
addRequired(p, 'x_model_24_m_0_Conv_', isValidArrayInput);
addRequired(p, 'x_model_24_m_1_Conv_', isValidArrayInput);
addRequired(p, 'x_model_24_m_2_Conv_', isValidArrayInput);
addRequired(p, 'params', isValidONNXParameters);
addParameter(p, 'InputDataPermutation', 'auto');
addParameter(p, 'OutputDataPermutation', 'auto');
addParameter(p, 'Training', false);
parse(p, x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_, params, varargin{:});
inputDataPerms = p.Results.InputDataPermutation;
outputDataPerms = p.Results.OutputDataPermutation;
Training = p.Results.Training;
if isnumeric(inputDataPerms)
    inputDataPerms = {inputDataPerms};
end
if isstring(inputDataPerms) && isscalar(inputDataPerms) || ischar(inputDataPerms)
    inputDataPerms = repmat({inputDataPerms},1,3);
end
if isnumeric(outputDataPerms)
    outputDataPerms = {outputDataPerms};
end
if isstring(outputDataPerms) && isscalar(outputDataPerms) || ischar(outputDataPerms)
    outputDataPerms = repmat({outputDataPerms},1,numDataOutputs);
end
end

function [x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_, Training, outputDataPerms, anyDlarrayInputs] = preprocessInput(x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_, params, varargin)
% Parse input arguments
[inputDataPerms, outputDataPerms, Training] = parseInputs(x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_, 1, params, varargin{:});
anyDlarrayInputs = any(cellfun(@(x)isa(x, 'dlarray'), {x_model_24_m_0_Conv_, x_model_24_m_1_Conv_, x_model_24_m_2_Conv_}));
% Make the input variables into unlabelled dlarrays:
x_model_24_m_0_Conv_ = makeUnlabeledDlarray(x_model_24_m_0_Conv_);
x_model_24_m_1_Conv_ = makeUnlabeledDlarray(x_model_24_m_1_Conv_);
x_model_24_m_2_Conv_ = makeUnlabeledDlarray(x_model_24_m_2_Conv_);
% Permute inputs if requested:
x_model_24_m_0_Conv_ = permuteInputVar(x_model_24_m_0_Conv_, inputDataPerms{1}, 4);
x_model_24_m_1_Conv_ = permuteInputVar(x_model_24_m_1_Conv_, inputDataPerms{2}, 4);
x_model_24_m_2_Conv_ = permuteInputVar(x_model_24_m_2_Conv_, inputDataPerms{3}, 4);
end

function [output0] = postprocessOutput(output0, outputDataPerms, anyDlarrayInputs, Training, varargin)
% Set output type:
if ~anyDlarrayInputs && ~Training
    if isdlarray(output0)
        output0 = extractdata(output0);
    end
end
% Permute outputs if requested:
output0 = permuteOutputVar(output0, outputDataPerms{1}, 3);
end


%% dlarray functions implementing ONNX operators:

function [Y, numDimsY] = onnxConcat(ONNXAxis, XCell, numDimsXArray)
% Concatentation that treats all empties the same. Necessary because
% dlarray.cat does not allow, for example, cat(1, 1x1, 1x0) because the
% second dimension sizes do not match.
numDimsY = numDimsXArray(1);
XCell(cellfun(@isempty, XCell)) = [];
if isempty(XCell)
    Y = dlarray([]);
else
    if ONNXAxis<0
        ONNXAxis = ONNXAxis + numDimsY;
    end
    DLTAxis = numDimsY - ONNXAxis;
    Y = cat(DLTAxis, XCell{:});
end
end

function varargout = onnxSplit(X, ONNXaxis, splits, numSplits, numDimsX)
% Implements the ONNX Split operator

% ONNXaxis is origin 0. splits is a vector of the lengths of each segment.
% If numSplits is nonzero, instead split into segments of equal length.
if ONNXaxis<0
    ONNXaxis = ONNXaxis + numDimsX;
end
DLTAxis = numDimsX - ONNXaxis;
if numSplits > 0
    C       = size(X, DLTAxis);
    sz      = floor(C/numSplits);
    splits	= repmat(sz, 1, numSplits);
else
    splits = extractdata(splits);
end
S      = struct;
S.type = '()';
S.subs = repmat({':'}, 1, ndims(X));
splitIndices = [0 cumsum(splits(:)')];
numY = numel(splitIndices)-1;
for i = 1:numY
    from            = splitIndices(i) + 1;
    to              = splitIndices(i+1);
    S.subs{DLTAxis}	= from:to;
    % The first numY outputs are the Y's. The second numY outputs are their
    % numDims. We assume all the outputs of Split have the same numDims as
    % the input.
    varargout{i}        = subsref(X, S);
    varargout{i + numY} = numDimsX;
end
end

function [DLTShape, numDimsY] = prepareReshapeArgs(X, ONNXShape, numDimsX, allowzero)
% Prepares arguments for implementing the ONNX Reshape operator
ONNXShape = flip(extractdata(ONNXShape));            % First flip the shape to make it correspond to the dimensions of X.
% In ONNX, 0 means "unchanged" if allowzero is false, and -1 means "infer". In DLT, there is no
% "unchanged", and [] means "infer".
DLTShape = num2cell(ONNXShape);                      % Make a cell array so we can include [].
% Replace zeros with the actual size if allowzero is true
if any(ONNXShape==0) && allowzero==0
    i0 = find(ONNXShape==0);
    DLTShape(i0) = num2cell(size(X, numDimsX - numel(ONNXShape) + i0));  % right-align the shape vector and dims
end
if any(ONNXShape == -1)
    % Replace -1 with []
    i = ONNXShape == -1;
    DLTShape{i} = [];
end
if numel(DLTShape)==1
    DLTShape = [DLTShape 1];
end
numDimsY = numel(ONNXShape);
end

function [perm, numDimsA] = prepareTransposeArgs(ONNXPerm, numDimsA)
% Prepares arguments for implementing the ONNX Transpose operator
if numDimsA <= 1        % Tensors of numDims 0 or 1 are unchanged by ONNX Transpose.
    perm = [];
else
    if isempty(ONNXPerm)        % Empty ONNXPerm means reverse the dimensions.
        perm = numDimsA:-1:1;
    else
        perm = numDimsA-flip(ONNXPerm);
    end
end
end

%% Utility functions:

function s = appendStructs(varargin)
% s = appendStructs(s1, s2,...). Assign all fields in s1, s2,... into s.
if isempty(varargin)
    s = struct;
else
    s = varargin{1};
    for i = 2:numel(varargin)
        fromstr = varargin{i};
        fs = fieldnames(fromstr);
        for j = 1:numel(fs)
            s.(fs{j}) = fromstr.(fs{j});
        end
    end
end
end

function checkInputSize(inputShape, expectedShape, inputName)

if numel(expectedShape)==0
    % The input is a scalar
    if ~isequal(inputShape, [1 1])
        inputSizeStr = makeSizeString(inputShape);
        error(message('nnet_cnn_onnx:onnx:InputNeedsResize',inputName, "[1,1]", inputSizeStr));
    end
elseif numel(expectedShape)==1
    % The input is a vector
    if ~shapeIsColumnVector(inputShape) || ~iSizesMatch({inputShape(1)}, expectedShape)
        expectedShape{2} = 1;
        expectedSizeStr = makeSizeString(expectedShape);
        inputSizeStr = makeSizeString(inputShape);
        error(message('nnet_cnn_onnx:onnx:InputNeedsResize',inputName, expectedSizeStr, inputSizeStr));
    end
else
    % The input has 2 dimensions or more
    
    % The input dimensions have been reversed; flip them back to compare to the
    % expected ONNX shape.
    inputShape = fliplr(inputShape);
    
    % If the expected shape has fewer dims than the input shape, error.
    if numel(expectedShape) < numel(inputShape)
        expectedSizeStr = strjoin(["[", strjoin(string(expectedShape), ","), "]"], "");
        error(message('nnet_cnn_onnx:onnx:InputHasGreaterNDims', inputName, expectedSizeStr));
    end
    
    % Prepad the input shape with trailing ones up to the number of elements in
    % expectedShape
    inputShape = num2cell([ones(1, numel(expectedShape) - length(inputShape)) inputShape]);
    
    % Find the number of variable size dimensions in the expected shape
    numVariableInputs = sum(cellfun(@(x) isa(x, 'char') || isa(x, 'string'), expectedShape));
    
    % Find the number of input dimensions that are not in the expected shape
    % and cannot be represented by a variable dimension
    nonMatchingInputDims = setdiff(string(inputShape), string(expectedShape));
    numNonMatchingInputDims  = numel(nonMatchingInputDims) - numVariableInputs;
    
    expectedSizeStr = makeSizeString(expectedShape);
    inputSizeStr = makeSizeString(inputShape);
    if numNonMatchingInputDims == 0 && ~iSizesMatch(inputShape, expectedShape)
        % The actual and expected input dimensions match, but in
        % a different order. The input needs to be permuted.
        error(message('nnet_cnn_onnx:onnx:InputNeedsPermute',inputName, expectedSizeStr, inputSizeStr));
    elseif numNonMatchingInputDims > 0
        % The actual and expected input sizes do not match.
        error(message('nnet_cnn_onnx:onnx:InputNeedsResize',inputName, expectedSizeStr, inputSizeStr));
    end
end
end

function doesMatch = iSizesMatch(inputShape, expectedShape)
% Check whether the input and expected shapes match, in order.
% Size elements match if (1) the elements are equal, or (2) the expected
% size element is a variable (represented by a character vector or string)
doesMatch = true;
for i=1:numel(inputShape)
    if ~(isequal(inputShape{i},expectedShape{i}) || ischar(expectedShape{i}) || isstring(expectedShape{i}))
        doesMatch = false;
        return
    end
end
end

function sizeStr = makeSizeString(shape)
sizeStr = strjoin(["[", strjoin(string(shape), ","), "]"], "");
end

function isVec = shapeIsColumnVector(shape)
if numel(shape) == 2 && shape(2) == 1
    isVec = true;
else
    isVec = false;
end
end
function X = makeUnlabeledDlarray(X)
% Make numeric X into an unlabelled dlarray
if isa(X, 'dlarray')
    X = stripdims(X);
elseif isnumeric(X)
    if isinteger(X)
        % Make ints double so they can combine with anything without
        % reducing precision
        X = double(X);
    end
    X = dlarray(X);
end
end

function [Vars, NumDims] = packageVariables(params, inputNames, inputValues, inputNumDims)
% inputNames, inputValues are cell arrays. inputRanks is a numeric vector.
Vars = appendStructs(params.Learnables, params.Nonlearnables, params.State);
NumDims = params.NumDimensions;
% Add graph inputs
for i = 1:numel(inputNames)
    Vars.(inputNames{i}) = inputValues{i};
    NumDims.(inputNames{i}) = inputNumDims(i);
end
end

function X = permuteInputVar(X, userDataPerm, onnxNDims)
% Returns reverse-ONNX ordering
if onnxNDims == 0
    return;
elseif onnxNDims == 1 && isvector(X)
    X = X(:);
    return;
elseif isnumeric(userDataPerm)
    % Permute into reverse ONNX ordering
    if numel(userDataPerm) ~= onnxNDims
        error(message('nnet_cnn_onnx:onnx:InputPermutationSize', numel(userDataPerm), onnxNDims));
    end
    perm = fliplr(userDataPerm);
elseif isequal(userDataPerm, 'auto') && onnxNDims == 4
    % Permute MATLAB HWCN to reverse onnx (WHCN)
    perm = [2 1 3 4];
elseif isequal(userDataPerm, 'as-is')
    % Do not permute the input
    perm = 1:ndims(X);
else
    % userDataPerm is either 'none' or 'auto' with no default, which means
    % it's already in onnx ordering, so just make it reverse onnx
    perm = max(2,onnxNDims):-1:1;
end
X = permute(X, perm);
end

function Y = permuteOutputVar(Y, userDataPerm, onnxNDims)
switch onnxNDims
    case 0
        perm = [];
    case 1
        if isnumeric(userDataPerm)
            % Use the user's permutation because Y is a column vector which
            % already matches ONNX.
            perm = userDataPerm;
        elseif isequal(userDataPerm, 'auto')
            % Treat the 1D onnx vector as a 2D column and transpose it
            perm = [2 1];
        else
            % userDataPerm is 'none'. Leave Y alone because it already
            % matches onnx.
            perm = [];
        end
    otherwise
        % ndims >= 2
        if isnumeric(userDataPerm)
            % Use the inverse of the user's permutation. This is not just the
            % flip of the permutation vector.
            perm = onnxNDims + 1 - userDataPerm;
        elseif isequal(userDataPerm, 'auto')
            if onnxNDims == 2
                % Permute reverse ONNX CN to DLT CN (do nothing)
                perm = [];
            elseif onnxNDims == 4
                % Permute reverse onnx (WHCN) to MATLAB HWCN
                perm = [2 1 3 4];
            else
                % User wants the output in ONNX ordering, so just reverse it from
                % reverse onnx
                perm = onnxNDims:-1:1;
            end
        elseif isequal(userDataPerm, 'as-is')
            % Do not permute the input
            perm = 1:ndims(Y);
        else
            % userDataPerm is 'none', so just make it reverse onnx
            perm = onnxNDims:-1:1;
        end
end
if ~isempty(perm)
    Y = permute(Y, perm);
end
end

function s = updateStruct(s, t)
% Set all existing fields in s from fields in t, ignoring extra fields in t.
for name = transpose(fieldnames(s))
    s.(name{1}) = t.(name{1});
end
end
