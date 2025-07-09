classdef MatMul_To_ConcatLayer1000 < nnet.layer.Layer & nnet.layer.Formattable
    % A custom layer auto-generated while importing an ONNX network.
    %#codegen

    %#ok<*PROPLC>
    %#ok<*NBRAK>
    %#ok<*INUSL>
    %#ok<*VARARG>
    properties (Learnable)
        controller_net_net_1
        controller_net_net_0
        controller_net_net_3
        controller_net_net_2
        x_controller_net__2
        onnx__MatMul_124
        onnx__MatMul_125
        dynamics_dynamics__1
        dynamics_dynamics_re
        dynamics_dynamics__3
        dynamics_dynamics__2
        dynamics_dynamics__5
        dynamics_dynamics__4
        x_dynamics_relu_m_2
        onnx__MatMul_134
        onnx__MatMul_135
        onnx__MatMul_136
        onnx__MatMul_137
        lyapunov_R
        onnx__ReduceSum_106
        onnx__ReduceSum_112
        x_Constant_13_output
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
                'controller_net_net_1'
                'controller_net_net_0'
                'controller_net_net_3'
                'controller_net_net_2'
                'x_controller_net__2'
                'onnx__MatMul_124'
                'onnx__MatMul_125'
                'dynamics_dynamics__1'
                'dynamics_dynamics_re'
                'dynamics_dynamics__3'
                'dynamics_dynamics__2'
                'dynamics_dynamics__5'
                'dynamics_dynamics__4'
                'x_dynamics_relu_m_2'
                'onnx__MatMul_134'
                'onnx__MatMul_135'
                'onnx__MatMul_136'
                'onnx__MatMul_137'
                'lyapunov_R'
                'onnx__ReduceSum_106'
                'onnx__ReduceSum_112'
                'x_Constant_13_output'
                };
        end
    end


    methods(Static, Hidden)
        % Instantiate a codegenable layer instance from a MATLAB layer instance
        function this_cg = matlabCodegenToRedirected(mlInstance)
            this_cg = relu_quadrotor2d_state.coder.MatMul_To_ConcatLayer1000(mlInstance);
        end
        function this_ml = matlabCodegenFromRedirected(cgInstance)
            this_ml = relu_quadrotor2d_state.MatMul_To_ConcatLayer1000(cgInstance.Name);
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
            this_ml.controller_net_net_1 = cgInstance.controller_net_net_1;
            this_ml.controller_net_net_0 = cgInstance.controller_net_net_0;
            this_ml.controller_net_net_3 = cgInstance.controller_net_net_3;
            this_ml.controller_net_net_2 = cgInstance.controller_net_net_2;
            this_ml.x_controller_net__2 = cgInstance.x_controller_net__2;
            this_ml.onnx__MatMul_124 = cgInstance.onnx__MatMul_124;
            this_ml.onnx__MatMul_125 = cgInstance.onnx__MatMul_125;
            this_ml.dynamics_dynamics__1 = cgInstance.dynamics_dynamics__1;
            this_ml.dynamics_dynamics_re = cgInstance.dynamics_dynamics_re;
            this_ml.dynamics_dynamics__3 = cgInstance.dynamics_dynamics__3;
            this_ml.dynamics_dynamics__2 = cgInstance.dynamics_dynamics__2;
            this_ml.dynamics_dynamics__5 = cgInstance.dynamics_dynamics__5;
            this_ml.dynamics_dynamics__4 = cgInstance.dynamics_dynamics__4;
            this_ml.x_dynamics_relu_m_2 = cgInstance.x_dynamics_relu_m_2;
            this_ml.onnx__MatMul_134 = cgInstance.onnx__MatMul_134;
            this_ml.onnx__MatMul_135 = cgInstance.onnx__MatMul_135;
            this_ml.onnx__MatMul_136 = cgInstance.onnx__MatMul_136;
            this_ml.onnx__MatMul_137 = cgInstance.onnx__MatMul_137;
            this_ml.lyapunov_R = cgInstance.lyapunov_R;
            this_ml.onnx__ReduceSum_106 = cgInstance.onnx__ReduceSum_106;
            this_ml.onnx__ReduceSum_112 = cgInstance.onnx__ReduceSum_112;
            this_ml.x_Constant_13_output = cgInstance.x_Constant_13_output;
        end
    end

    methods
        function this = MatMul_To_ConcatLayer1000(mlInstance)
            this.Name = mlInstance.Name;
            this.NumInputs = 2;
            this.OutputNames = {'x123'};
            if isstruct(mlInstance.Vars)
                names = fieldnames(mlInstance.Vars);
                for i=1:numel(names)
                    fieldname = names{i};
                    this.Vars.(fieldname) = relu_quadrotor2d_state.coder.ops.extractIfDlarray(mlInstance.Vars.(fieldname));
                end
            else
                this.Vars = [];
            end

            this.NumDims = mlInstance.NumDims;
            this.controller_net_net_1 = mlInstance.controller_net_net_1;
            this.controller_net_net_0 = mlInstance.controller_net_net_0;
            this.controller_net_net_3 = mlInstance.controller_net_net_3;
            this.controller_net_net_2 = mlInstance.controller_net_net_2;
            this.x_controller_net__2 = mlInstance.x_controller_net__2;
            this.onnx__MatMul_124 = mlInstance.onnx__MatMul_124;
            this.onnx__MatMul_125 = mlInstance.onnx__MatMul_125;
            this.dynamics_dynamics__1 = mlInstance.dynamics_dynamics__1;
            this.dynamics_dynamics_re = mlInstance.dynamics_dynamics_re;
            this.dynamics_dynamics__3 = mlInstance.dynamics_dynamics__3;
            this.dynamics_dynamics__2 = mlInstance.dynamics_dynamics__2;
            this.dynamics_dynamics__5 = mlInstance.dynamics_dynamics__5;
            this.dynamics_dynamics__4 = mlInstance.dynamics_dynamics__4;
            this.x_dynamics_relu_m_2 = mlInstance.x_dynamics_relu_m_2;
            this.onnx__MatMul_134 = mlInstance.onnx__MatMul_134;
            this.onnx__MatMul_135 = mlInstance.onnx__MatMul_135;
            this.onnx__MatMul_136 = mlInstance.onnx__MatMul_136;
            this.onnx__MatMul_137 = mlInstance.onnx__MatMul_137;
            this.lyapunov_R = mlInstance.lyapunov_R;
            this.onnx__ReduceSum_106 = mlInstance.onnx__ReduceSum_106;
            this.onnx__ReduceSum_112 = mlInstance.onnx__ReduceSum_112;
            this.x_Constant_13_output = mlInstance.x_Constant_13_output;
        end

        function [x123] = predict(this, onnx__Mul_0__, onnx__Mul_0NumDims__)
            if isdlarray(onnx__Mul_0__)
                onnx__Mul_0_ = stripdims(onnx__Mul_0__);
            else
                onnx__Mul_0_ = onnx__Mul_0__;
            end
            onnx__Mul_0NumDims = numel(onnx__Mul_0NumDims__);
            onnx__Mul_0 = relu_quadrotor2d_state.coder.ops.permuteInputVar(onnx__Mul_0_, ['as-is'], 2);

            [x123__, x123NumDims__] = MatMul_To_ConcatGraph1000(this, onnx__Mul_0, onnx__Mul_0NumDims, false);
            x123_ = relu_quadrotor2d_state.coder.ops.permuteOutputVar(x123__, ['as-is'], 2);

            x123 = dlarray(single(x123_), repmat('U', 1, max(2, x123NumDims__)));
        end

        function [x123, x123NumDims1011] = MatMul_To_ConcatGraph1000(this, onnx__Mul_0, onnx__Mul_0NumDims, Training)

            % Execute the operators:
            % Mul:
            x_Mul_output_0 = onnx__Mul_0 .* this.Vars.x_Constant_output_0;
            x_Mul_output_0NumDims = max(onnx__Mul_0NumDims, this.NumDims.x_Constant_output_0);

            % Gemm:
            [A1000, B1001, C1002, alpha1003, beta1004, x_controller_net_netNumDims] = relu_quadrotor2d_state.coder.ops.prepareGemmArgs(x_Mul_output_0, this.controller_net_net_1, this.controller_net_net_0, this.Vars.Gemmalpha1001, this.Vars.Gemmbeta1002, 0, 1, this.NumDims.controller_net_net_0);
            x_controller_net_net = alpha1003*B1001*A1000 + beta1004*C1002;

            % Relu:
            X1005 = dlarray(relu_quadrotor2d_state.coder.ops.extractIfDlarray(x_controller_net_net));
            Y1006 = relu(X1005);
            x_controller_net__4 = relu_quadrotor2d_state.coder.ops.extractIfDlarray(Y1006);
            x_controller_net__4NumDims = x_controller_net_netNumDims;

            % Gemm:
            [A1007, B1008, C1009, alpha1010, beta1011, x_controller_net__6NumDims] = relu_quadrotor2d_state.coder.ops.prepareGemmArgs(x_controller_net__4, this.controller_net_net_3, this.controller_net_net_2, this.Vars.Gemmalpha1003, this.Vars.Gemmbeta1004, 0, 1, this.NumDims.controller_net_net_2);
            x_controller_net__6 = alpha1010*B1008*A1007 + beta1011*C1009;

            % MatMul:
            [x_controller_net__3, x_controller_net__3NumDims] = relu_quadrotor2d_state.coder.ops.onnxMatMul(this.x_controller_net__2, this.onnx__MatMul_124, this.NumDims.x_controller_net__2, this.NumDims.onnx__MatMul_124);

            % Add:
            x_controller_net__1 = this.controller_net_net_0 + x_controller_net__3;
            x_controller_net__1NumDims = max(this.NumDims.controller_net_net_0, x_controller_net__3NumDims);

            % Relu:
            X1012 = dlarray(relu_quadrotor2d_state.coder.ops.extractIfDlarray(x_controller_net__1));
            Y1013 = relu(X1012);
            x_controller_net__5 = relu_quadrotor2d_state.coder.ops.extractIfDlarray(Y1013);
            x_controller_net__5NumDims = x_controller_net__1NumDims;

            % MatMul:
            [x_controller_net__8, x_controller_net__8NumDims] = relu_quadrotor2d_state.coder.ops.onnxMatMul(x_controller_net__5, this.onnx__MatMul_125, x_controller_net__5NumDims, this.NumDims.onnx__MatMul_125);

            % Add:
            x_controller_net__7 = this.controller_net_net_2 + x_controller_net__8;
            x_controller_net__7NumDims = max(this.NumDims.controller_net_net_2, x_controller_net__8NumDims);

            % Sub:
            x_controller_Sub_out = x_controller_net__6 - x_controller_net__7;
            x_controller_Sub_outNumDims = max(x_controller_net__6NumDims, x_controller_net__7NumDims);

            % Add:
            x_controller_Add_out = x_controller_Sub_out + this.Vars.x_controller_Const_4;
            x_controller_Add_outNumDims = max(x_controller_Sub_outNumDims, this.NumDims.x_controller_Const_4);

            % Sub:
            x_controller_Sub_1_o = x_controller_Add_out - this.Vars.x_controller_Constan;
            x_controller_Sub_1_oNumDims = max(x_controller_Add_outNumDims, this.NumDims.x_controller_Constan);

            % Relu:
            X1014 = dlarray(relu_quadrotor2d_state.coder.ops.extractIfDlarray(x_controller_Sub_1_o));
            Y1015 = relu(X1014);
            x_controller_Relu_ou = relu_quadrotor2d_state.coder.ops.extractIfDlarray(Y1015);
            x_controller_Relu_ouNumDims = x_controller_Sub_1_oNumDims;

            % Add:
            x_controller_Add_1_o = x_controller_Relu_ou + this.Vars.x_controller_Const_1;
            x_controller_Add_1_oNumDims = max(x_controller_Relu_ouNumDims, this.NumDims.x_controller_Const_1);

            % Sub:
            x_controller_Sub_2_o = this.Vars.x_controller_Const_2 - x_controller_Add_1_o;
            x_controller_Sub_2_oNumDims = max(this.NumDims.x_controller_Const_2, x_controller_Add_1_oNumDims);

            % Relu:
            X1016 = dlarray(relu_quadrotor2d_state.coder.ops.extractIfDlarray(x_controller_Sub_2_o));
            Y1017 = relu(X1016);
            x_controller_Relu_1_ = relu_quadrotor2d_state.coder.ops.extractIfDlarray(Y1017);
            x_controller_Relu_1_NumDims = x_controller_Sub_2_oNumDims;

            % Sub:
            x_controller_Sub_3_o = x_controller_Relu_1_ - this.Vars.x_controller_Const_3;
            x_controller_Sub_3_oNumDims = max(x_controller_Relu_1_NumDims, this.NumDims.x_controller_Const_3);

            % Neg:
            x_controller_Neg_out = -(x_controller_Sub_3_o);
            x_controller_Neg_outNumDims = x_controller_Sub_3_oNumDims;

            % Div:
            x_Div_output_0 = x_Mul_output_0 ./ this.Vars.x_Constant_1_output_;
            x_Div_output_0NumDims = max(x_Mul_output_0NumDims, this.NumDims.x_Constant_1_output_);

            % Cast:
            x_Cast_output_0 = single(x_controller_Neg_out);
            x_Cast_output_0NumDims = x_controller_Neg_outNumDims;

            % Slice:
            [indices1018, x_Slice_output_0NumDims] = relu_quadrotor2d_state.coder.ops.prepareSliceArgs(x_Div_output_0, this.Vars.x_Constant_3_output_, this.Vars.x_Constant_4_output_, this.Vars.x_Constant_2_output_, this.Vars.x_Constant_5_output_, x_Div_output_0NumDims);
            x_Slice_output_0 = x_Div_output_0(indices1018{:});

            % Slice:
            [indices1019, x_Slice_1_output_0NumDims] = relu_quadrotor2d_state.coder.ops.prepareSliceArgs(x_Div_output_0, this.Vars.x_Constant_7_output_, this.Vars.x_Constant_8_output_, this.Vars.x_Constant_6_output_, this.Vars.x_Constant_9_output_, x_Div_output_0NumDims);
            x_Slice_1_output_0 = x_Div_output_0(indices1019{:});

            % Gather:
            [x_Gather_output_0, x_Gather_output_0NumDims] = relu_quadrotor2d_state.coder.ops.onnxGather(x_Div_output_0, this.Vars.onnx__Gather_65, 1, x_Div_output_0NumDims, this.NumDims.onnx__Gather_65);

            % Concat:
            [x_Concat_output_0, x_Concat_output_0NumDims] = relu_quadrotor2d_state.coder.ops.onnxConcat(-1, {x_Gather_output_0, x_Cast_output_0}, [x_Gather_output_0NumDims, x_Cast_output_0NumDims]);

            % Gemm:
            [A1020, B1021, C1022, alpha1023, beta1024, x_dynamics_relu_modeNumDims] = relu_quadrotor2d_state.coder.ops.prepareGemmArgs(x_Concat_output_0, this.dynamics_dynamics__1, this.dynamics_dynamics_re, this.Vars.Gemmalpha1005, this.Vars.Gemmbeta1006, 0, 1, this.NumDims.dynamics_dynamics_re);
            x_dynamics_relu_mode = alpha1023*B1021*A1020 + beta1024*C1022;

            % Relu:
            X1025 = dlarray(relu_quadrotor2d_state.coder.ops.extractIfDlarray(x_dynamics_relu_mode));
            Y1026 = relu(X1025);
            x_dynamics_relu_m_4 = relu_quadrotor2d_state.coder.ops.extractIfDlarray(Y1026);
            x_dynamics_relu_m_4NumDims = x_dynamics_relu_modeNumDims;

            % Gemm:
            [A1027, B1028, C1029, alpha1030, beta1031, x_dynamics_relu_m_6NumDims] = relu_quadrotor2d_state.coder.ops.prepareGemmArgs(x_dynamics_relu_m_4, this.dynamics_dynamics__3, this.dynamics_dynamics__2, this.Vars.Gemmalpha1007, this.Vars.Gemmbeta1008, 0, 1, this.NumDims.dynamics_dynamics__2);
            x_dynamics_relu_m_6 = alpha1030*B1028*A1027 + beta1031*C1029;

            % Relu:
            X1032 = dlarray(relu_quadrotor2d_state.coder.ops.extractIfDlarray(x_dynamics_relu_m_6));
            Y1033 = relu(X1032);
            x_dynamics_relu_m_9 = relu_quadrotor2d_state.coder.ops.extractIfDlarray(Y1033);
            x_dynamics_relu_m_9NumDims = x_dynamics_relu_m_6NumDims;

            % Gemm:
            [A1034, B1035, C1036, alpha1037, beta1038, x_dynamics_relu_m_11NumDims] = relu_quadrotor2d_state.coder.ops.prepareGemmArgs(x_dynamics_relu_m_9, this.dynamics_dynamics__5, this.dynamics_dynamics__4, this.Vars.Gemmalpha1009, this.Vars.Gemmbeta1010, 0, 1, this.NumDims.dynamics_dynamics__4);
            x_dynamics_relu_m_11 = alpha1037*B1035*A1034 + beta1038*C1036;

            % MatMul:
            [x_dynamics_relu_m_3, x_dynamics_relu_m_3NumDims] = relu_quadrotor2d_state.coder.ops.onnxMatMul(this.x_dynamics_relu_m_2, this.onnx__MatMul_134, this.NumDims.x_dynamics_relu_m_2, this.NumDims.onnx__MatMul_134);

            % Add:
            x_dynamics_relu_m_1 = this.dynamics_dynamics_re + x_dynamics_relu_m_3;
            x_dynamics_relu_m_1NumDims = max(this.NumDims.dynamics_dynamics_re, x_dynamics_relu_m_3NumDims);

            % Relu:
            X1039 = dlarray(relu_quadrotor2d_state.coder.ops.extractIfDlarray(x_dynamics_relu_m_1));
            Y1040 = relu(X1039);
            x_dynamics_relu_m_5 = relu_quadrotor2d_state.coder.ops.extractIfDlarray(Y1040);
            x_dynamics_relu_m_5NumDims = x_dynamics_relu_m_1NumDims;

            % MatMul:
            [x_dynamics_relu_m_8, x_dynamics_relu_m_8NumDims] = relu_quadrotor2d_state.coder.ops.onnxMatMul(x_dynamics_relu_m_5, this.onnx__MatMul_135, x_dynamics_relu_m_5NumDims, this.NumDims.onnx__MatMul_135);

            % Add:
            x_dynamics_relu_m_7 = this.dynamics_dynamics__2 + x_dynamics_relu_m_8;
            x_dynamics_relu_m_7NumDims = max(this.NumDims.dynamics_dynamics__2, x_dynamics_relu_m_8NumDims);

            % Relu:
            X1041 = dlarray(relu_quadrotor2d_state.coder.ops.extractIfDlarray(x_dynamics_relu_m_7));
            Y1042 = relu(X1041);
            x_dynamics_relu_m_10 = relu_quadrotor2d_state.coder.ops.extractIfDlarray(Y1042);
            x_dynamics_relu_m_10NumDims = x_dynamics_relu_m_7NumDims;

            % MatMul:
            [x_dynamics_relu_m_13, x_dynamics_relu_m_13NumDims] = relu_quadrotor2d_state.coder.ops.onnxMatMul(x_dynamics_relu_m_10, this.onnx__MatMul_136, x_dynamics_relu_m_10NumDims, this.NumDims.onnx__MatMul_136);

            % Add:
            x_dynamics_relu_m_12 = this.dynamics_dynamics__4 + x_dynamics_relu_m_13;
            x_dynamics_relu_m_12NumDims = max(this.NumDims.dynamics_dynamics__4, x_dynamics_relu_m_13NumDims);

            % Add:
            x_Add_output_0 = x_Slice_1_output_0 + x_dynamics_relu_m_11;
            x_Add_output_0NumDims = max(x_Slice_1_output_0NumDims, x_dynamics_relu_m_11NumDims);

            % Sub:
            x_Sub_output_0 = x_Add_output_0 - x_dynamics_relu_m_12;
            x_Sub_output_0NumDims = max(x_Add_output_0NumDims, x_dynamics_relu_m_12NumDims);

            % Add:
            x_Add_1_output_0 = x_Slice_1_output_0 + x_Sub_output_0;
            x_Add_1_output_0NumDims = max(x_Slice_1_output_0NumDims, x_Sub_output_0NumDims);

            % Mul:
            x_Mul_1_output_0 = x_Add_1_output_0 .* this.Vars.x_Constant_10_output;
            x_Mul_1_output_0NumDims = max(x_Add_1_output_0NumDims, this.NumDims.x_Constant_10_output);

            % Div:
            x_Div_1_output_0 = x_Mul_1_output_0 ./ this.Vars.x_Constant_11_output;
            x_Div_1_output_0NumDims = max(x_Mul_1_output_0NumDims, this.NumDims.x_Constant_11_output);

            % Add:
            x_Add_2_output_0 = x_Slice_output_0 + x_Div_1_output_0;
            x_Add_2_output_0NumDims = max(x_Slice_output_0NumDims, x_Div_1_output_0NumDims);

            % Concat:
            [x_Concat_1_output_0, x_Concat_1_output_0NumDims] = relu_quadrotor2d_state.coder.ops.onnxConcat(-1, {x_Add_2_output_0, x_Sub_output_0}, [x_Add_2_output_0NumDims, x_Sub_output_0NumDims]);

            % Sub:
            x_Sub_1_output_0 = x_Concat_1_output_0 - x_Div_output_0;
            x_Sub_1_output_0NumDims = max(x_Concat_1_output_0NumDims, x_Div_output_0NumDims);

            % Mul:
            x_Mul_2_output_0 = x_Sub_1_output_0 .* this.Vars.x_Constant_12_output;
            x_Mul_2_output_0NumDims = max(x_Sub_1_output_0NumDims, this.NumDims.x_Constant_12_output);

            % Add:
            x_Add_3_output_0 = x_Mul_output_0 + x_Mul_2_output_0;
            x_Add_3_output_0NumDims = max(x_Mul_output_0NumDims, x_Mul_2_output_0NumDims);

            % Sub:
            x_lyapunov_Sub_outpu = x_Mul_output_0 - this.Vars.x_lyapunov_Constan_1;
            x_lyapunov_Sub_outpuNumDims = max(x_Mul_output_0NumDims, this.NumDims.x_lyapunov_Constan_1);

            % MatMul:
            [x_lyapunov_MatMul_ou, x_lyapunov_MatMul_ouNumDims] = relu_quadrotor2d_state.coder.ops.onnxMatMul(this.onnx__MatMul_137, this.lyapunov_R, this.NumDims.onnx__MatMul_137, this.NumDims.lyapunov_R);

            % Add:
            x_lyapunov_Add_outpu = this.Vars.x_lyapunov_Constant_ + x_lyapunov_MatMul_ou;
            x_lyapunov_Add_outpuNumDims = max(this.NumDims.x_lyapunov_Constant_, x_lyapunov_MatMul_ouNumDims);

            % MatMul:
            [x_lyapunov_MatMul_1_, x_lyapunov_MatMul_1_NumDims] = relu_quadrotor2d_state.coder.ops.onnxMatMul(x_lyapunov_Sub_outpu, x_lyapunov_Add_outpu, x_lyapunov_Sub_outpuNumDims, x_lyapunov_Add_outpuNumDims);

            % Mul:
            x_lyapunov_Mul_outpu = x_lyapunov_Sub_outpu .* x_lyapunov_MatMul_1_;
            x_lyapunov_Mul_outpuNumDims = max(x_lyapunov_Sub_outpuNumDims, x_lyapunov_MatMul_1_NumDims);

            % ReduceSum:
            dims1043 = relu_quadrotor2d_state.coder.ops.prepareReduceArgs(this.onnx__ReduceSum_106, this.NumDims.onnx__ReduceSum_106);
            if isempty(dims1043)
                   dims1043 = 1; % or dims1043 = 'all'; depending on what the ONNX intended
            end
            xReduced1044 = sum(x_lyapunov_Mul_outpu, dims1043);
            x_lyapunov_ReduceSum = xReduced1044;
            x_lyapunov_ReduceSumNumDims = x_lyapunov_Mul_outpuNumDims;

            % Sub:
            x_Sub_2_output_0 = x_Add_3_output_0 - x_Mul_output_0;
            x_Sub_2_output_0NumDims = max(x_Add_3_output_0NumDims, x_Mul_output_0NumDims);

            % Add:
            x_Add_4_output_0 = x_Add_3_output_0 + x_Mul_output_0;
            x_Add_4_output_0NumDims = max(x_Add_3_output_0NumDims, x_Mul_output_0NumDims);

            % MatMul:
            [x_MatMul_output_0, x_MatMul_output_0NumDims] = relu_quadrotor2d_state.coder.ops.onnxMatMul(x_Add_4_output_0, x_lyapunov_Add_outpu, x_Add_4_output_0NumDims, x_lyapunov_Add_outpuNumDims);

            % Mul:
            x_Mul_3_output_0 = x_Sub_2_output_0 .* x_MatMul_output_0;
            x_Mul_3_output_0NumDims = max(x_Sub_2_output_0NumDims, x_MatMul_output_0NumDims);

            % ReduceSum:
            dims1045 = relu_quadrotor2d_state.coder.ops.prepareReduceArgs(this.onnx__ReduceSum_112, this.NumDims.onnx__ReduceSum_112);
            if isempty(dims1045)
                   dims1045 = 1; % or dims1043 = 'all'; depending on what the ONNX intended
            end
            xReduced1046 = sum(x_Mul_3_output_0, dims1045);
            x_ReduceSum_output_0 = xReduced1046;
            x_ReduceSum_output_0NumDims = x_Mul_3_output_0NumDims;

            % MatMul:
            [x_MatMul_1_output_0, x_MatMul_1_output_0NumDims] = relu_quadrotor2d_state.coder.ops.onnxMatMul(this.x_Constant_13_output, x_lyapunov_Add_outpu, this.NumDims.x_Constant_13_output, x_lyapunov_Add_outpuNumDims);

            % Mul:
            x_Mul_4_output_0 = x_Sub_2_output_0 .* x_MatMul_1_output_0;
            x_Mul_4_output_0NumDims = max(x_Sub_2_output_0NumDims, x_MatMul_1_output_0NumDims);

            % ReduceSum:
            dims1047 = relu_quadrotor2d_state.coder.ops.prepareReduceArgs(this.onnx__ReduceSum_112, this.NumDims.onnx__ReduceSum_112);
            if isempty(dims1047)
                   dims1047 = 1; % or dims1043 = 'all'; depending on what the ONNX intended
            end
            xReduced1048 = sum(x_Mul_4_output_0, dims1047);
            x_ReduceSum_1_output = xReduced1048;
            x_ReduceSum_1_outputNumDims = x_Mul_4_output_0NumDims;

            % Mul:
            x_Mul_5_output_0 = x_ReduceSum_1_output .* this.Vars.x_Constant_14_output;
            x_Mul_5_output_0NumDims = max(x_ReduceSum_1_outputNumDims, this.NumDims.x_Constant_14_output);

            % Sub:
            x_Sub_3_output_0 = x_ReduceSum_output_0 - x_Mul_5_output_0;
            x_Sub_3_output_0NumDims = max(x_ReduceSum_output_0NumDims, x_Mul_5_output_0NumDims);

            % Div:
            x_Div_2_output_0 = x_Add_3_output_0 ./ this.Vars.x_Constant_15_output;
            x_Div_2_output_0NumDims = max(x_Add_3_output_0NumDims, this.NumDims.x_Constant_15_output);

            % Concat:
            [x123, x123NumDims] = relu_quadrotor2d_state.coder.ops.onnxConcat(-1, {x_Sub_3_output_0, x_lyapunov_ReduceSum, x_Div_2_output_0}, [x_Sub_3_output_0NumDims, x_lyapunov_ReduceSumNumDims, x_Div_2_output_0NumDims]);

            % Set graph output arguments
            x123NumDims1011 = x123NumDims;

        end

    end

end