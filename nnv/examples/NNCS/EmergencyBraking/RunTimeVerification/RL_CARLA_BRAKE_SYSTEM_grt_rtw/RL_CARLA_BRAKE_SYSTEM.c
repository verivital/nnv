/*
 * RL_CARLA_BRAKE_SYSTEM.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "RL_CARLA_BRAKE_SYSTEM".
 *
 * Model version              : 1.41
 * Simulink Coder version : 9.0 (R2018b) 24-May-2018
 * C source code generated on : Tue Apr 23 14:37:57 2019
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "RL_CARLA_BRAKE_SYSTEM.h"
#include "RL_CARLA_BRAKE_SYSTEM_private.h"
#include "RL_CARLA_BRAKE_SYSTEM_dt.h"

/* Block signals (default storage) */
B_RL_CARLA_BRAKE_SYSTEM_T RL_CARLA_BRAKE_SYSTEM_B;

/* Block states (default storage) */
DW_RL_CARLA_BRAKE_SYSTEM_T RL_CARLA_BRAKE_SYSTEM_DW;

/* Real-time model */
RT_MODEL_RL_CARLA_BRAKE_SYSTE_T RL_CARLA_BRAKE_SYSTEM_M_;
RT_MODEL_RL_CARLA_BRAKE_SYSTE_T *const RL_CARLA_BRAKE_SYSTEM_M =
  &RL_CARLA_BRAKE_SYSTEM_M_;

/* Model step function */
void RL_CARLA_BRAKE_SYSTEM_step(void)
{
  /* local block i/o variables */
  real_T rtb_Gain2;
  real_T rtb_NormalizedSpeed;
  real_T tmp[50];
  real_T tmp_0[30];
  real_T tmp_1;
  real_T tmp_2;
  real_T tmp_3;
  real_T tmp_4;
  real_T tmp_5;
  real_T tmp_6;
  real_T tmp_7;
  real_T tmp_8;
  real_T tmp_9;
  real_T tmp_a;
  real_T tmp_b;
  real_T tmp_c;
  real_T tmp_d;
  real_T tmp_e;
  real_T tmp_f;
  real_T tmp_g;
  real_T tmp_h;
  real_T tmp_i;
  real_T tmp_j;
  real_T tmp_k;
  real_T tmp_l;
  real_T tmp_m;
  real_T tmp_n;
  real_T tmp_o;
  real_T tmp_p;
  real_T tmp_q;
  real_T tmp_r;
  real_T tmp_s;
  int32_T i;
  real_T rtb_TmpSignalConversionAtDotP_0;
  real_T rtb_TmpSignalConversionAtDotP_1;
  real_T u0;

  /* DiscreteStateSpace: '<Root>/Discrete State-Space' */
  {
    RL_CARLA_BRAKE_SYSTEM_B.DiscreteStateSpace[0] =
      (RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_C[0])*
      RL_CARLA_BRAKE_SYSTEM_DW.DiscreteStateSpace_DSTATE[0];
    RL_CARLA_BRAKE_SYSTEM_B.DiscreteStateSpace[1] =
      (RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_C[1])*
      RL_CARLA_BRAKE_SYSTEM_DW.DiscreteStateSpace_DSTATE[1];
    RL_CARLA_BRAKE_SYSTEM_B.DiscreteStateSpace[2] =
      (RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_C[2])*
      RL_CARLA_BRAKE_SYSTEM_DW.DiscreteStateSpace_DSTATE[2];
  }

  /* Gain: '<S2>/Gain' */
  rtb_NormalizedSpeed = RL_CARLA_BRAKE_SYSTEM_P.Gain_Gain *
    RL_CARLA_BRAKE_SYSTEM_B.DiscreteStateSpace[1];

  /* Stop: '<Root>/Stop Simulation' incorporates:
   *  Constant: '<S1>/Constant'
   *  RelationalOperator: '<S1>/Compare'
   */
  if (rtb_NormalizedSpeed <= RL_CARLA_BRAKE_SYSTEM_P.Constant_Value) {
    rtmSetStopRequested(RL_CARLA_BRAKE_SYSTEM_M, 1);
  }

  /* End of Stop: '<Root>/Stop Simulation' */

  /* SignalConversion: '<S13>/TmpSignal ConversionAtDot ProductInport2' incorporates:
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S2>/Gain2'
   */
  rtb_TmpSignalConversionAtDotP_0 = RL_CARLA_BRAKE_SYSTEM_P.Gain1_Gain *
    RL_CARLA_BRAKE_SYSTEM_B.DiscreteStateSpace[0];
  rtb_TmpSignalConversionAtDotP_1 = RL_CARLA_BRAKE_SYSTEM_P.Gain2_Gain *
    RL_CARLA_BRAKE_SYSTEM_B.DiscreteStateSpace[2];

  /* Sum: '<S5>/netsum' incorporates:
   *  Constant: '<S11>/IW{1,1}(1,:)''
   *  Constant: '<S11>/IW{1,1}(10,:)''
   *  Constant: '<S11>/IW{1,1}(11,:)''
   *  Constant: '<S11>/IW{1,1}(12,:)''
   *  Constant: '<S11>/IW{1,1}(13,:)''
   *  Constant: '<S11>/IW{1,1}(14,:)''
   *  Constant: '<S11>/IW{1,1}(15,:)''
   *  Constant: '<S11>/IW{1,1}(16,:)''
   *  Constant: '<S11>/IW{1,1}(17,:)''
   *  Constant: '<S11>/IW{1,1}(18,:)''
   *  Constant: '<S11>/IW{1,1}(19,:)''
   *  Constant: '<S11>/IW{1,1}(2,:)''
   *  Constant: '<S11>/IW{1,1}(20,:)''
   *  Constant: '<S11>/IW{1,1}(21,:)''
   *  Constant: '<S11>/IW{1,1}(22,:)''
   *  Constant: '<S11>/IW{1,1}(23,:)''
   *  Constant: '<S11>/IW{1,1}(24,:)''
   *  Constant: '<S11>/IW{1,1}(25,:)''
   *  Constant: '<S11>/IW{1,1}(26,:)''
   *  Constant: '<S11>/IW{1,1}(27,:)''
   *  Constant: '<S11>/IW{1,1}(28,:)''
   *  Constant: '<S11>/IW{1,1}(29,:)''
   *  Constant: '<S11>/IW{1,1}(3,:)''
   *  Constant: '<S11>/IW{1,1}(30,:)''
   *  Constant: '<S11>/IW{1,1}(31,:)''
   *  Constant: '<S11>/IW{1,1}(32,:)''
   *  Constant: '<S11>/IW{1,1}(33,:)''
   *  Constant: '<S11>/IW{1,1}(34,:)''
   *  Constant: '<S11>/IW{1,1}(35,:)''
   *  Constant: '<S11>/IW{1,1}(36,:)''
   *  Constant: '<S11>/IW{1,1}(37,:)''
   *  Constant: '<S11>/IW{1,1}(38,:)''
   *  Constant: '<S11>/IW{1,1}(39,:)''
   *  Constant: '<S11>/IW{1,1}(4,:)''
   *  Constant: '<S11>/IW{1,1}(40,:)''
   *  Constant: '<S11>/IW{1,1}(41,:)''
   *  Constant: '<S11>/IW{1,1}(42,:)''
   *  Constant: '<S11>/IW{1,1}(43,:)''
   *  Constant: '<S11>/IW{1,1}(44,:)''
   *  Constant: '<S11>/IW{1,1}(45,:)''
   *  Constant: '<S11>/IW{1,1}(46,:)''
   *  Constant: '<S11>/IW{1,1}(47,:)''
   *  Constant: '<S11>/IW{1,1}(48,:)''
   *  Constant: '<S11>/IW{1,1}(49,:)''
   *  Constant: '<S11>/IW{1,1}(5,:)''
   *  Constant: '<S11>/IW{1,1}(50,:)''
   *  Constant: '<S11>/IW{1,1}(6,:)''
   *  Constant: '<S11>/IW{1,1}(7,:)''
   *  Constant: '<S11>/IW{1,1}(8,:)''
   *  Constant: '<S11>/IW{1,1}(9,:)''
   *  DotProduct: '<S13>/Dot Product'
   *  DotProduct: '<S14>/Dot Product'
   *  DotProduct: '<S15>/Dot Product'
   *  DotProduct: '<S16>/Dot Product'
   *  DotProduct: '<S17>/Dot Product'
   *  DotProduct: '<S18>/Dot Product'
   *  DotProduct: '<S19>/Dot Product'
   *  DotProduct: '<S20>/Dot Product'
   *  DotProduct: '<S21>/Dot Product'
   *  DotProduct: '<S22>/Dot Product'
   *  DotProduct: '<S23>/Dot Product'
   *  DotProduct: '<S24>/Dot Product'
   *  DotProduct: '<S25>/Dot Product'
   *  DotProduct: '<S26>/Dot Product'
   *  DotProduct: '<S27>/Dot Product'
   *  DotProduct: '<S28>/Dot Product'
   *  DotProduct: '<S29>/Dot Product'
   *  DotProduct: '<S30>/Dot Product'
   *  DotProduct: '<S31>/Dot Product'
   *  DotProduct: '<S32>/Dot Product'
   *  DotProduct: '<S33>/Dot Product'
   *  DotProduct: '<S34>/Dot Product'
   *  DotProduct: '<S35>/Dot Product'
   *  DotProduct: '<S36>/Dot Product'
   *  DotProduct: '<S37>/Dot Product'
   *  DotProduct: '<S38>/Dot Product'
   *  DotProduct: '<S39>/Dot Product'
   *  DotProduct: '<S40>/Dot Product'
   *  DotProduct: '<S41>/Dot Product'
   *  DotProduct: '<S42>/Dot Product'
   *  DotProduct: '<S43>/Dot Product'
   *  DotProduct: '<S44>/Dot Product'
   *  DotProduct: '<S45>/Dot Product'
   *  DotProduct: '<S46>/Dot Product'
   *  DotProduct: '<S47>/Dot Product'
   *  DotProduct: '<S48>/Dot Product'
   *  DotProduct: '<S49>/Dot Product'
   *  DotProduct: '<S50>/Dot Product'
   *  DotProduct: '<S51>/Dot Product'
   *  DotProduct: '<S52>/Dot Product'
   *  DotProduct: '<S53>/Dot Product'
   *  DotProduct: '<S54>/Dot Product'
   *  DotProduct: '<S55>/Dot Product'
   *  DotProduct: '<S56>/Dot Product'
   *  DotProduct: '<S57>/Dot Product'
   *  DotProduct: '<S58>/Dot Product'
   *  DotProduct: '<S59>/Dot Product'
   *  DotProduct: '<S60>/Dot Product'
   *  DotProduct: '<S61>/Dot Product'
   *  DotProduct: '<S62>/Dot Product'
   *  SignalConversion: '<S13>/TmpSignal ConversionAtDot ProductInport2'
   */
  tmp[0] = (RL_CARLA_BRAKE_SYSTEM_P.IW111_Value[0] *
            rtb_TmpSignalConversionAtDotP_0 +
            RL_CARLA_BRAKE_SYSTEM_P.IW111_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW111_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[1] = (RL_CARLA_BRAKE_SYSTEM_P.IW112_Value[0] *
            rtb_TmpSignalConversionAtDotP_0 +
            RL_CARLA_BRAKE_SYSTEM_P.IW112_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW112_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[2] = (RL_CARLA_BRAKE_SYSTEM_P.IW113_Value[0] *
            rtb_TmpSignalConversionAtDotP_0 +
            RL_CARLA_BRAKE_SYSTEM_P.IW113_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW113_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[3] = (RL_CARLA_BRAKE_SYSTEM_P.IW114_Value[0] *
            rtb_TmpSignalConversionAtDotP_0 +
            RL_CARLA_BRAKE_SYSTEM_P.IW114_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW114_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[4] = (RL_CARLA_BRAKE_SYSTEM_P.IW115_Value[0] *
            rtb_TmpSignalConversionAtDotP_0 +
            RL_CARLA_BRAKE_SYSTEM_P.IW115_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW115_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[5] = (RL_CARLA_BRAKE_SYSTEM_P.IW116_Value[0] *
            rtb_TmpSignalConversionAtDotP_0 +
            RL_CARLA_BRAKE_SYSTEM_P.IW116_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW116_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[6] = (RL_CARLA_BRAKE_SYSTEM_P.IW117_Value[0] *
            rtb_TmpSignalConversionAtDotP_0 +
            RL_CARLA_BRAKE_SYSTEM_P.IW117_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW117_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[7] = (RL_CARLA_BRAKE_SYSTEM_P.IW118_Value[0] *
            rtb_TmpSignalConversionAtDotP_0 +
            RL_CARLA_BRAKE_SYSTEM_P.IW118_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW118_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[8] = (RL_CARLA_BRAKE_SYSTEM_P.IW119_Value[0] *
            rtb_TmpSignalConversionAtDotP_0 +
            RL_CARLA_BRAKE_SYSTEM_P.IW119_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW119_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[9] = (RL_CARLA_BRAKE_SYSTEM_P.IW1110_Value[0] *
            rtb_TmpSignalConversionAtDotP_0 +
            RL_CARLA_BRAKE_SYSTEM_P.IW1110_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1110_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[10] = (RL_CARLA_BRAKE_SYSTEM_P.IW1111_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1111_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1111_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[11] = (RL_CARLA_BRAKE_SYSTEM_P.IW1112_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1112_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1112_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[12] = (RL_CARLA_BRAKE_SYSTEM_P.IW1113_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1113_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1113_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[13] = (RL_CARLA_BRAKE_SYSTEM_P.IW1114_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1114_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1114_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[14] = (RL_CARLA_BRAKE_SYSTEM_P.IW1115_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1115_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1115_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[15] = (RL_CARLA_BRAKE_SYSTEM_P.IW1116_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1116_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1116_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[16] = (RL_CARLA_BRAKE_SYSTEM_P.IW1117_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1117_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1117_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[17] = (RL_CARLA_BRAKE_SYSTEM_P.IW1118_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1118_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1118_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[18] = (RL_CARLA_BRAKE_SYSTEM_P.IW1119_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1119_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1119_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[19] = (RL_CARLA_BRAKE_SYSTEM_P.IW1120_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1120_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1120_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[20] = (RL_CARLA_BRAKE_SYSTEM_P.IW1121_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1121_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1121_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[21] = (RL_CARLA_BRAKE_SYSTEM_P.IW1122_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1122_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1122_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[22] = (RL_CARLA_BRAKE_SYSTEM_P.IW1123_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1123_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1123_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[23] = (RL_CARLA_BRAKE_SYSTEM_P.IW1124_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1124_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1124_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[24] = (RL_CARLA_BRAKE_SYSTEM_P.IW1125_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1125_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1125_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[25] = (RL_CARLA_BRAKE_SYSTEM_P.IW1126_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1126_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1126_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[26] = (RL_CARLA_BRAKE_SYSTEM_P.IW1127_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1127_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1127_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[27] = (RL_CARLA_BRAKE_SYSTEM_P.IW1128_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1128_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1128_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[28] = (RL_CARLA_BRAKE_SYSTEM_P.IW1129_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1129_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1129_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[29] = (RL_CARLA_BRAKE_SYSTEM_P.IW1130_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1130_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1130_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[30] = (RL_CARLA_BRAKE_SYSTEM_P.IW1131_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1131_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1131_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[31] = (RL_CARLA_BRAKE_SYSTEM_P.IW1132_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1132_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1132_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[32] = (RL_CARLA_BRAKE_SYSTEM_P.IW1133_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1133_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1133_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[33] = (RL_CARLA_BRAKE_SYSTEM_P.IW1134_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1134_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1134_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[34] = (RL_CARLA_BRAKE_SYSTEM_P.IW1135_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1135_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1135_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[35] = (RL_CARLA_BRAKE_SYSTEM_P.IW1136_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1136_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1136_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[36] = (RL_CARLA_BRAKE_SYSTEM_P.IW1137_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1137_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1137_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[37] = (RL_CARLA_BRAKE_SYSTEM_P.IW1138_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1138_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1138_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[38] = (RL_CARLA_BRAKE_SYSTEM_P.IW1139_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1139_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1139_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[39] = (RL_CARLA_BRAKE_SYSTEM_P.IW1140_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1140_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1140_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[40] = (RL_CARLA_BRAKE_SYSTEM_P.IW1141_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1141_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1141_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[41] = (RL_CARLA_BRAKE_SYSTEM_P.IW1142_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1142_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1142_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[42] = (RL_CARLA_BRAKE_SYSTEM_P.IW1143_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1143_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1143_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[43] = (RL_CARLA_BRAKE_SYSTEM_P.IW1144_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1144_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1144_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[44] = (RL_CARLA_BRAKE_SYSTEM_P.IW1145_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1145_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1145_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[45] = (RL_CARLA_BRAKE_SYSTEM_P.IW1146_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1146_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1146_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[46] = (RL_CARLA_BRAKE_SYSTEM_P.IW1147_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1147_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1147_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[47] = (RL_CARLA_BRAKE_SYSTEM_P.IW1148_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1148_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1148_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[48] = (RL_CARLA_BRAKE_SYSTEM_P.IW1149_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1149_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1149_Value[2] * rtb_TmpSignalConversionAtDotP_1;
  tmp[49] = (RL_CARLA_BRAKE_SYSTEM_P.IW1150_Value[0] *
             rtb_TmpSignalConversionAtDotP_0 +
             RL_CARLA_BRAKE_SYSTEM_P.IW1150_Value[1] * rtb_NormalizedSpeed) +
    RL_CARLA_BRAKE_SYSTEM_P.IW1150_Value[2] * rtb_TmpSignalConversionAtDotP_1;

  /* DotProduct: '<S66>/Dot Product' */
  rtb_TmpSignalConversionAtDotP_0 = 0.0;

  /* DotProduct: '<S77>/Dot Product' */
  rtb_TmpSignalConversionAtDotP_1 = 0.0;

  /* DotProduct: '<S88>/Dot Product' */
  tmp_1 = 0.0;

  /* DotProduct: '<S90>/Dot Product' */
  tmp_2 = 0.0;

  /* DotProduct: '<S91>/Dot Product' */
  tmp_3 = 0.0;

  /* DotProduct: '<S92>/Dot Product' */
  tmp_4 = 0.0;

  /* DotProduct: '<S93>/Dot Product' */
  tmp_5 = 0.0;

  /* DotProduct: '<S94>/Dot Product' */
  tmp_6 = 0.0;

  /* DotProduct: '<S95>/Dot Product' */
  tmp_7 = 0.0;

  /* DotProduct: '<S67>/Dot Product' */
  tmp_8 = 0.0;

  /* DotProduct: '<S68>/Dot Product' */
  tmp_9 = 0.0;

  /* DotProduct: '<S69>/Dot Product' */
  tmp_a = 0.0;

  /* DotProduct: '<S70>/Dot Product' */
  tmp_b = 0.0;

  /* DotProduct: '<S71>/Dot Product' */
  tmp_c = 0.0;

  /* DotProduct: '<S72>/Dot Product' */
  tmp_d = 0.0;

  /* DotProduct: '<S73>/Dot Product' */
  tmp_e = 0.0;

  /* DotProduct: '<S74>/Dot Product' */
  tmp_f = 0.0;

  /* DotProduct: '<S75>/Dot Product' */
  tmp_g = 0.0;

  /* DotProduct: '<S76>/Dot Product' */
  tmp_h = 0.0;

  /* DotProduct: '<S78>/Dot Product' */
  tmp_i = 0.0;

  /* DotProduct: '<S79>/Dot Product' */
  tmp_j = 0.0;

  /* DotProduct: '<S80>/Dot Product' */
  tmp_k = 0.0;

  /* DotProduct: '<S81>/Dot Product' */
  tmp_l = 0.0;

  /* DotProduct: '<S82>/Dot Product' */
  tmp_m = 0.0;

  /* DotProduct: '<S83>/Dot Product' */
  tmp_n = 0.0;

  /* DotProduct: '<S84>/Dot Product' */
  tmp_o = 0.0;

  /* DotProduct: '<S85>/Dot Product' */
  tmp_p = 0.0;

  /* DotProduct: '<S86>/Dot Product' */
  tmp_q = 0.0;

  /* DotProduct: '<S87>/Dot Product' */
  tmp_r = 0.0;

  /* DotProduct: '<S89>/Dot Product' */
  tmp_s = 0.0;
  for (i = 0; i < 50; i++) {
    /* Saturate: '<S12>/Saturation' incorporates:
     *  Constant: '<S5>/b{1}'
     *  Sum: '<S5>/netsum'
     */
    u0 = tmp[i] + RL_CARLA_BRAKE_SYSTEM_P.b1_Value[i];
    if (u0 > RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat) {
      u0 = RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat;
    } else {
      if (u0 < RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat) {
        u0 = RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat;
      }
    }

    /* End of Saturate: '<S12>/Saturation' */

    /* DotProduct: '<S66>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(1,:)''
     */
    rtb_TmpSignalConversionAtDotP_0 += RL_CARLA_BRAKE_SYSTEM_P.IW211_Value[i] *
      u0;

    /* DotProduct: '<S77>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(2,:)''
     */
    rtb_TmpSignalConversionAtDotP_1 += RL_CARLA_BRAKE_SYSTEM_P.IW212_Value[i] *
      u0;

    /* DotProduct: '<S88>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(3,:)''
     */
    tmp_1 += RL_CARLA_BRAKE_SYSTEM_P.IW213_Value[i] * u0;

    /* DotProduct: '<S90>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(4,:)''
     */
    tmp_2 += RL_CARLA_BRAKE_SYSTEM_P.IW214_Value[i] * u0;

    /* DotProduct: '<S91>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(5,:)''
     */
    tmp_3 += RL_CARLA_BRAKE_SYSTEM_P.IW215_Value[i] * u0;

    /* DotProduct: '<S92>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(6,:)''
     */
    tmp_4 += RL_CARLA_BRAKE_SYSTEM_P.IW216_Value[i] * u0;

    /* DotProduct: '<S93>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(7,:)''
     */
    tmp_5 += RL_CARLA_BRAKE_SYSTEM_P.IW217_Value[i] * u0;

    /* DotProduct: '<S94>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(8,:)''
     */
    tmp_6 += RL_CARLA_BRAKE_SYSTEM_P.IW218_Value[i] * u0;

    /* DotProduct: '<S95>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(9,:)''
     */
    tmp_7 += RL_CARLA_BRAKE_SYSTEM_P.IW219_Value[i] * u0;

    /* DotProduct: '<S67>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(10,:)''
     */
    tmp_8 += RL_CARLA_BRAKE_SYSTEM_P.IW2110_Value[i] * u0;

    /* DotProduct: '<S68>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(11,:)''
     */
    tmp_9 += RL_CARLA_BRAKE_SYSTEM_P.IW2111_Value[i] * u0;

    /* DotProduct: '<S69>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(12,:)''
     */
    tmp_a += RL_CARLA_BRAKE_SYSTEM_P.IW2112_Value[i] * u0;

    /* DotProduct: '<S70>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(13,:)''
     */
    tmp_b += RL_CARLA_BRAKE_SYSTEM_P.IW2113_Value[i] * u0;

    /* DotProduct: '<S71>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(14,:)''
     */
    tmp_c += RL_CARLA_BRAKE_SYSTEM_P.IW2114_Value[i] * u0;

    /* DotProduct: '<S72>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(15,:)''
     */
    tmp_d += RL_CARLA_BRAKE_SYSTEM_P.IW2115_Value[i] * u0;

    /* DotProduct: '<S73>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(16,:)''
     */
    tmp_e += RL_CARLA_BRAKE_SYSTEM_P.IW2116_Value[i] * u0;

    /* DotProduct: '<S74>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(17,:)''
     */
    tmp_f += RL_CARLA_BRAKE_SYSTEM_P.IW2117_Value[i] * u0;

    /* DotProduct: '<S75>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(18,:)''
     */
    tmp_g += RL_CARLA_BRAKE_SYSTEM_P.IW2118_Value[i] * u0;

    /* DotProduct: '<S76>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(19,:)''
     */
    tmp_h += RL_CARLA_BRAKE_SYSTEM_P.IW2119_Value[i] * u0;

    /* DotProduct: '<S78>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(20,:)''
     */
    tmp_i += RL_CARLA_BRAKE_SYSTEM_P.IW2120_Value[i] * u0;

    /* DotProduct: '<S79>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(21,:)''
     */
    tmp_j += RL_CARLA_BRAKE_SYSTEM_P.IW2121_Value[i] * u0;

    /* DotProduct: '<S80>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(22,:)''
     */
    tmp_k += RL_CARLA_BRAKE_SYSTEM_P.IW2122_Value[i] * u0;

    /* DotProduct: '<S81>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(23,:)''
     */
    tmp_l += RL_CARLA_BRAKE_SYSTEM_P.IW2123_Value[i] * u0;

    /* DotProduct: '<S82>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(24,:)''
     */
    tmp_m += RL_CARLA_BRAKE_SYSTEM_P.IW2124_Value[i] * u0;

    /* DotProduct: '<S83>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(25,:)''
     */
    tmp_n += RL_CARLA_BRAKE_SYSTEM_P.IW2125_Value[i] * u0;

    /* DotProduct: '<S84>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(26,:)''
     */
    tmp_o += RL_CARLA_BRAKE_SYSTEM_P.IW2126_Value[i] * u0;

    /* DotProduct: '<S85>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(27,:)''
     */
    tmp_p += RL_CARLA_BRAKE_SYSTEM_P.IW2127_Value[i] * u0;

    /* DotProduct: '<S86>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(28,:)''
     */
    tmp_q += RL_CARLA_BRAKE_SYSTEM_P.IW2128_Value[i] * u0;

    /* DotProduct: '<S87>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(29,:)''
     */
    tmp_r += RL_CARLA_BRAKE_SYSTEM_P.IW2129_Value[i] * u0;

    /* DotProduct: '<S89>/Dot Product' incorporates:
     *  Constant: '<S64>/IW{2,1}(30,:)''
     */
    tmp_s += RL_CARLA_BRAKE_SYSTEM_P.IW2130_Value[i] * u0;
  }

  /* Sum: '<S6>/netsum' incorporates:
   *  DotProduct: '<S66>/Dot Product'
   *  DotProduct: '<S67>/Dot Product'
   *  DotProduct: '<S68>/Dot Product'
   *  DotProduct: '<S69>/Dot Product'
   *  DotProduct: '<S70>/Dot Product'
   *  DotProduct: '<S71>/Dot Product'
   *  DotProduct: '<S72>/Dot Product'
   *  DotProduct: '<S73>/Dot Product'
   *  DotProduct: '<S74>/Dot Product'
   *  DotProduct: '<S75>/Dot Product'
   *  DotProduct: '<S76>/Dot Product'
   *  DotProduct: '<S77>/Dot Product'
   *  DotProduct: '<S78>/Dot Product'
   *  DotProduct: '<S79>/Dot Product'
   *  DotProduct: '<S80>/Dot Product'
   *  DotProduct: '<S81>/Dot Product'
   *  DotProduct: '<S82>/Dot Product'
   *  DotProduct: '<S83>/Dot Product'
   *  DotProduct: '<S84>/Dot Product'
   *  DotProduct: '<S85>/Dot Product'
   *  DotProduct: '<S86>/Dot Product'
   *  DotProduct: '<S87>/Dot Product'
   *  DotProduct: '<S88>/Dot Product'
   *  DotProduct: '<S89>/Dot Product'
   *  DotProduct: '<S90>/Dot Product'
   *  DotProduct: '<S91>/Dot Product'
   *  DotProduct: '<S92>/Dot Product'
   *  DotProduct: '<S93>/Dot Product'
   *  DotProduct: '<S94>/Dot Product'
   *  DotProduct: '<S95>/Dot Product'
   */
  tmp_0[0] = rtb_TmpSignalConversionAtDotP_0;
  tmp_0[1] = rtb_TmpSignalConversionAtDotP_1;
  tmp_0[2] = tmp_1;
  tmp_0[3] = tmp_2;
  tmp_0[4] = tmp_3;
  tmp_0[5] = tmp_4;
  tmp_0[6] = tmp_5;
  tmp_0[7] = tmp_6;
  tmp_0[8] = tmp_7;
  tmp_0[9] = tmp_8;
  tmp_0[10] = tmp_9;
  tmp_0[11] = tmp_a;
  tmp_0[12] = tmp_b;
  tmp_0[13] = tmp_c;
  tmp_0[14] = tmp_d;
  tmp_0[15] = tmp_e;
  tmp_0[16] = tmp_f;
  tmp_0[17] = tmp_g;
  tmp_0[18] = tmp_h;
  tmp_0[19] = tmp_i;
  tmp_0[20] = tmp_j;
  tmp_0[21] = tmp_k;
  tmp_0[22] = tmp_l;
  tmp_0[23] = tmp_m;
  tmp_0[24] = tmp_n;
  tmp_0[25] = tmp_o;
  tmp_0[26] = tmp_p;
  tmp_0[27] = tmp_q;
  tmp_0[28] = tmp_r;
  tmp_0[29] = tmp_s;

  /* DotProduct: '<S99>/Dot Product' incorporates:
   *  Constant: '<S97>/IW{3,2}(1,:)''
   */
  rtb_TmpSignalConversionAtDotP_0 = 0.0;
  for (i = 0; i < 30; i++) {
    /* Saturate: '<S65>/Saturation' incorporates:
     *  Constant: '<S6>/b{2}'
     *  Constant: '<S97>/IW{3,2}(1,:)''
     *  Sum: '<S6>/netsum'
     */
    u0 = tmp_0[i] + RL_CARLA_BRAKE_SYSTEM_P.b2_Value[i];
    if (u0 > RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_e) {
      u0 = RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_e;
    } else {
      if (u0 < RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat_p) {
        u0 = RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat_p;
      }
    }

    /* End of Saturate: '<S65>/Saturation' */
    rtb_TmpSignalConversionAtDotP_0 += RL_CARLA_BRAKE_SYSTEM_P.IW321_Value[i] *
      u0;
  }

  /* Sum: '<S7>/netsum' incorporates:
   *  Constant: '<S7>/b{3}'
   *  DotProduct: '<S99>/Dot Product'
   */
  u0 = rtb_TmpSignalConversionAtDotP_0 + RL_CARLA_BRAKE_SYSTEM_P.b3_Value;

  /* Saturate: '<S98>/Saturation' */
  if (u0 > RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_j) {
    u0 = RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_j;
  } else {
    if (u0 < RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat_n) {
      u0 = RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat_n;
    }
  }

  /* Sum: '<S101>/netsum' incorporates:
   *  Constant: '<S107>/IW{1,1}(1,:)''
   *  Constant: '<S107>/IW{1,1}(10,:)''
   *  Constant: '<S107>/IW{1,1}(11,:)''
   *  Constant: '<S107>/IW{1,1}(12,:)''
   *  Constant: '<S107>/IW{1,1}(13,:)''
   *  Constant: '<S107>/IW{1,1}(14,:)''
   *  Constant: '<S107>/IW{1,1}(15,:)''
   *  Constant: '<S107>/IW{1,1}(16,:)''
   *  Constant: '<S107>/IW{1,1}(17,:)''
   *  Constant: '<S107>/IW{1,1}(18,:)''
   *  Constant: '<S107>/IW{1,1}(19,:)''
   *  Constant: '<S107>/IW{1,1}(2,:)''
   *  Constant: '<S107>/IW{1,1}(20,:)''
   *  Constant: '<S107>/IW{1,1}(21,:)''
   *  Constant: '<S107>/IW{1,1}(22,:)''
   *  Constant: '<S107>/IW{1,1}(23,:)''
   *  Constant: '<S107>/IW{1,1}(24,:)''
   *  Constant: '<S107>/IW{1,1}(25,:)''
   *  Constant: '<S107>/IW{1,1}(26,:)''
   *  Constant: '<S107>/IW{1,1}(27,:)''
   *  Constant: '<S107>/IW{1,1}(28,:)''
   *  Constant: '<S107>/IW{1,1}(29,:)''
   *  Constant: '<S107>/IW{1,1}(3,:)''
   *  Constant: '<S107>/IW{1,1}(30,:)''
   *  Constant: '<S107>/IW{1,1}(31,:)''
   *  Constant: '<S107>/IW{1,1}(32,:)''
   *  Constant: '<S107>/IW{1,1}(33,:)''
   *  Constant: '<S107>/IW{1,1}(34,:)''
   *  Constant: '<S107>/IW{1,1}(35,:)''
   *  Constant: '<S107>/IW{1,1}(36,:)''
   *  Constant: '<S107>/IW{1,1}(37,:)''
   *  Constant: '<S107>/IW{1,1}(38,:)''
   *  Constant: '<S107>/IW{1,1}(39,:)''
   *  Constant: '<S107>/IW{1,1}(4,:)''
   *  Constant: '<S107>/IW{1,1}(40,:)''
   *  Constant: '<S107>/IW{1,1}(41,:)''
   *  Constant: '<S107>/IW{1,1}(42,:)''
   *  Constant: '<S107>/IW{1,1}(43,:)''
   *  Constant: '<S107>/IW{1,1}(44,:)''
   *  Constant: '<S107>/IW{1,1}(45,:)''
   *  Constant: '<S107>/IW{1,1}(46,:)''
   *  Constant: '<S107>/IW{1,1}(47,:)''
   *  Constant: '<S107>/IW{1,1}(48,:)''
   *  Constant: '<S107>/IW{1,1}(49,:)''
   *  Constant: '<S107>/IW{1,1}(5,:)''
   *  Constant: '<S107>/IW{1,1}(50,:)''
   *  Constant: '<S107>/IW{1,1}(6,:)''
   *  Constant: '<S107>/IW{1,1}(7,:)''
   *  Constant: '<S107>/IW{1,1}(8,:)''
   *  Constant: '<S107>/IW{1,1}(9,:)''
   *  DotProduct: '<S109>/Dot Product'
   *  DotProduct: '<S110>/Dot Product'
   *  DotProduct: '<S111>/Dot Product'
   *  DotProduct: '<S112>/Dot Product'
   *  DotProduct: '<S113>/Dot Product'
   *  DotProduct: '<S114>/Dot Product'
   *  DotProduct: '<S115>/Dot Product'
   *  DotProduct: '<S116>/Dot Product'
   *  DotProduct: '<S117>/Dot Product'
   *  DotProduct: '<S118>/Dot Product'
   *  DotProduct: '<S119>/Dot Product'
   *  DotProduct: '<S120>/Dot Product'
   *  DotProduct: '<S121>/Dot Product'
   *  DotProduct: '<S122>/Dot Product'
   *  DotProduct: '<S123>/Dot Product'
   *  DotProduct: '<S124>/Dot Product'
   *  DotProduct: '<S125>/Dot Product'
   *  DotProduct: '<S126>/Dot Product'
   *  DotProduct: '<S127>/Dot Product'
   *  DotProduct: '<S128>/Dot Product'
   *  DotProduct: '<S129>/Dot Product'
   *  DotProduct: '<S130>/Dot Product'
   *  DotProduct: '<S131>/Dot Product'
   *  DotProduct: '<S132>/Dot Product'
   *  DotProduct: '<S133>/Dot Product'
   *  DotProduct: '<S134>/Dot Product'
   *  DotProduct: '<S135>/Dot Product'
   *  DotProduct: '<S136>/Dot Product'
   *  DotProduct: '<S137>/Dot Product'
   *  DotProduct: '<S138>/Dot Product'
   *  DotProduct: '<S139>/Dot Product'
   *  DotProduct: '<S140>/Dot Product'
   *  DotProduct: '<S141>/Dot Product'
   *  DotProduct: '<S142>/Dot Product'
   *  DotProduct: '<S143>/Dot Product'
   *  DotProduct: '<S144>/Dot Product'
   *  DotProduct: '<S145>/Dot Product'
   *  DotProduct: '<S146>/Dot Product'
   *  DotProduct: '<S147>/Dot Product'
   *  DotProduct: '<S148>/Dot Product'
   *  DotProduct: '<S149>/Dot Product'
   *  DotProduct: '<S150>/Dot Product'
   *  DotProduct: '<S151>/Dot Product'
   *  DotProduct: '<S152>/Dot Product'
   *  DotProduct: '<S153>/Dot Product'
   *  DotProduct: '<S154>/Dot Product'
   *  DotProduct: '<S155>/Dot Product'
   *  DotProduct: '<S156>/Dot Product'
   *  DotProduct: '<S157>/Dot Product'
   *  DotProduct: '<S158>/Dot Product'
   *  Saturate: '<S98>/Saturation'
   *  SignalConversion: '<S109>/TmpSignal ConversionAtDot ProductInport2'
   */
  tmp[0] = RL_CARLA_BRAKE_SYSTEM_P.IW111_Value_e[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW111_Value_e[1] * u0;
  tmp[1] = RL_CARLA_BRAKE_SYSTEM_P.IW112_Value_i[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW112_Value_i[1] * u0;
  tmp[2] = RL_CARLA_BRAKE_SYSTEM_P.IW113_Value_p[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW113_Value_p[1] * u0;
  tmp[3] = RL_CARLA_BRAKE_SYSTEM_P.IW114_Value_p[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW114_Value_p[1] * u0;
  tmp[4] = RL_CARLA_BRAKE_SYSTEM_P.IW115_Value_m[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW115_Value_m[1] * u0;
  tmp[5] = RL_CARLA_BRAKE_SYSTEM_P.IW116_Value_k[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW116_Value_k[1] * u0;
  tmp[6] = RL_CARLA_BRAKE_SYSTEM_P.IW117_Value_p[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW117_Value_p[1] * u0;
  tmp[7] = RL_CARLA_BRAKE_SYSTEM_P.IW118_Value_f[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW118_Value_f[1] * u0;
  tmp[8] = RL_CARLA_BRAKE_SYSTEM_P.IW119_Value_n[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW119_Value_n[1] * u0;
  tmp[9] = RL_CARLA_BRAKE_SYSTEM_P.IW1110_Value_d[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1110_Value_d[1] * u0;
  tmp[10] = RL_CARLA_BRAKE_SYSTEM_P.IW1111_Value_g[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1111_Value_g[1] * u0;
  tmp[11] = RL_CARLA_BRAKE_SYSTEM_P.IW1112_Value_b[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1112_Value_b[1] * u0;
  tmp[12] = RL_CARLA_BRAKE_SYSTEM_P.IW1113_Value_l[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1113_Value_l[1] * u0;
  tmp[13] = RL_CARLA_BRAKE_SYSTEM_P.IW1114_Value_b[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1114_Value_b[1] * u0;
  tmp[14] = RL_CARLA_BRAKE_SYSTEM_P.IW1115_Value_c[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1115_Value_c[1] * u0;
  tmp[15] = RL_CARLA_BRAKE_SYSTEM_P.IW1116_Value_p[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1116_Value_p[1] * u0;
  tmp[16] = RL_CARLA_BRAKE_SYSTEM_P.IW1117_Value_i[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1117_Value_i[1] * u0;
  tmp[17] = RL_CARLA_BRAKE_SYSTEM_P.IW1118_Value_c[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1118_Value_c[1] * u0;
  tmp[18] = RL_CARLA_BRAKE_SYSTEM_P.IW1119_Value_b[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1119_Value_b[1] * u0;
  tmp[19] = RL_CARLA_BRAKE_SYSTEM_P.IW1120_Value_m[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1120_Value_m[1] * u0;
  tmp[20] = RL_CARLA_BRAKE_SYSTEM_P.IW1121_Value_i[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1121_Value_i[1] * u0;
  tmp[21] = RL_CARLA_BRAKE_SYSTEM_P.IW1122_Value_k[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1122_Value_k[1] * u0;
  tmp[22] = RL_CARLA_BRAKE_SYSTEM_P.IW1123_Value_b[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1123_Value_b[1] * u0;
  tmp[23] = RL_CARLA_BRAKE_SYSTEM_P.IW1124_Value_i[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1124_Value_i[1] * u0;
  tmp[24] = RL_CARLA_BRAKE_SYSTEM_P.IW1125_Value_l[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1125_Value_l[1] * u0;
  tmp[25] = RL_CARLA_BRAKE_SYSTEM_P.IW1126_Value_h[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1126_Value_h[1] * u0;
  tmp[26] = RL_CARLA_BRAKE_SYSTEM_P.IW1127_Value_c[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1127_Value_c[1] * u0;
  tmp[27] = RL_CARLA_BRAKE_SYSTEM_P.IW1128_Value_a[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1128_Value_a[1] * u0;
  tmp[28] = RL_CARLA_BRAKE_SYSTEM_P.IW1129_Value_p[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1129_Value_p[1] * u0;
  tmp[29] = RL_CARLA_BRAKE_SYSTEM_P.IW1130_Value_h[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1130_Value_h[1] * u0;
  tmp[30] = RL_CARLA_BRAKE_SYSTEM_P.IW1131_Value_l[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1131_Value_l[1] * u0;
  tmp[31] = RL_CARLA_BRAKE_SYSTEM_P.IW1132_Value_i[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1132_Value_i[1] * u0;
  tmp[32] = RL_CARLA_BRAKE_SYSTEM_P.IW1133_Value_g[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1133_Value_g[1] * u0;
  tmp[33] = RL_CARLA_BRAKE_SYSTEM_P.IW1134_Value_d[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1134_Value_d[1] * u0;
  tmp[34] = RL_CARLA_BRAKE_SYSTEM_P.IW1135_Value_n[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1135_Value_n[1] * u0;
  tmp[35] = RL_CARLA_BRAKE_SYSTEM_P.IW1136_Value_k[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1136_Value_k[1] * u0;
  tmp[36] = RL_CARLA_BRAKE_SYSTEM_P.IW1137_Value_f[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1137_Value_f[1] * u0;
  tmp[37] = RL_CARLA_BRAKE_SYSTEM_P.IW1138_Value_o[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1138_Value_o[1] * u0;
  tmp[38] = RL_CARLA_BRAKE_SYSTEM_P.IW1139_Value_o[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1139_Value_o[1] * u0;
  tmp[39] = RL_CARLA_BRAKE_SYSTEM_P.IW1140_Value_l[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1140_Value_l[1] * u0;
  tmp[40] = RL_CARLA_BRAKE_SYSTEM_P.IW1141_Value_c[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1141_Value_c[1] * u0;
  tmp[41] = RL_CARLA_BRAKE_SYSTEM_P.IW1142_Value_l[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1142_Value_l[1] * u0;
  tmp[42] = RL_CARLA_BRAKE_SYSTEM_P.IW1143_Value_i[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1143_Value_i[1] * u0;
  tmp[43] = RL_CARLA_BRAKE_SYSTEM_P.IW1144_Value_d[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1144_Value_d[1] * u0;
  tmp[44] = RL_CARLA_BRAKE_SYSTEM_P.IW1145_Value_c[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1145_Value_c[1] * u0;
  tmp[45] = RL_CARLA_BRAKE_SYSTEM_P.IW1146_Value_p[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1146_Value_p[1] * u0;
  tmp[46] = RL_CARLA_BRAKE_SYSTEM_P.IW1147_Value_d[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1147_Value_d[1] * u0;
  tmp[47] = RL_CARLA_BRAKE_SYSTEM_P.IW1148_Value_l[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1148_Value_l[1] * u0;
  tmp[48] = RL_CARLA_BRAKE_SYSTEM_P.IW1149_Value_d[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1149_Value_d[1] * u0;
  tmp[49] = RL_CARLA_BRAKE_SYSTEM_P.IW1150_Value_f[0] * rtb_NormalizedSpeed +
    RL_CARLA_BRAKE_SYSTEM_P.IW1150_Value_f[1] * u0;

  /* DotProduct: '<S162>/Dot Product' */
  rtb_TmpSignalConversionAtDotP_0 = 0.0;

  /* DotProduct: '<S173>/Dot Product' */
  rtb_TmpSignalConversionAtDotP_1 = 0.0;

  /* DotProduct: '<S184>/Dot Product' */
  tmp_1 = 0.0;

  /* DotProduct: '<S186>/Dot Product' */
  tmp_2 = 0.0;

  /* DotProduct: '<S187>/Dot Product' */
  tmp_3 = 0.0;

  /* DotProduct: '<S188>/Dot Product' */
  tmp_4 = 0.0;

  /* DotProduct: '<S189>/Dot Product' */
  tmp_5 = 0.0;

  /* DotProduct: '<S190>/Dot Product' */
  tmp_6 = 0.0;

  /* DotProduct: '<S191>/Dot Product' */
  tmp_7 = 0.0;

  /* DotProduct: '<S163>/Dot Product' */
  tmp_8 = 0.0;

  /* DotProduct: '<S164>/Dot Product' */
  tmp_9 = 0.0;

  /* DotProduct: '<S165>/Dot Product' */
  tmp_a = 0.0;

  /* DotProduct: '<S166>/Dot Product' */
  tmp_b = 0.0;

  /* DotProduct: '<S167>/Dot Product' */
  tmp_c = 0.0;

  /* DotProduct: '<S168>/Dot Product' */
  tmp_d = 0.0;

  /* DotProduct: '<S169>/Dot Product' */
  tmp_e = 0.0;

  /* DotProduct: '<S170>/Dot Product' */
  tmp_f = 0.0;

  /* DotProduct: '<S171>/Dot Product' */
  tmp_g = 0.0;

  /* DotProduct: '<S172>/Dot Product' */
  tmp_h = 0.0;

  /* DotProduct: '<S174>/Dot Product' */
  tmp_i = 0.0;

  /* DotProduct: '<S175>/Dot Product' */
  tmp_j = 0.0;

  /* DotProduct: '<S176>/Dot Product' */
  tmp_k = 0.0;

  /* DotProduct: '<S177>/Dot Product' */
  tmp_l = 0.0;

  /* DotProduct: '<S178>/Dot Product' */
  tmp_m = 0.0;

  /* DotProduct: '<S179>/Dot Product' */
  tmp_n = 0.0;

  /* DotProduct: '<S180>/Dot Product' */
  tmp_o = 0.0;

  /* DotProduct: '<S181>/Dot Product' */
  tmp_p = 0.0;

  /* DotProduct: '<S182>/Dot Product' */
  tmp_q = 0.0;

  /* DotProduct: '<S183>/Dot Product' */
  tmp_r = 0.0;

  /* DotProduct: '<S185>/Dot Product' */
  tmp_s = 0.0;
  for (i = 0; i < 50; i++) {
    /* Saturate: '<S108>/Saturation' incorporates:
     *  Constant: '<S101>/b{1}'
     *  Sum: '<S101>/netsum'
     */
    u0 = tmp[i] + RL_CARLA_BRAKE_SYSTEM_P.b1_Value_i[i];
    if (u0 > RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_f) {
      u0 = RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_f;
    } else {
      if (u0 < RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat_b) {
        u0 = RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat_b;
      }
    }

    /* End of Saturate: '<S108>/Saturation' */

    /* DotProduct: '<S162>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(1,:)''
     */
    rtb_TmpSignalConversionAtDotP_0 += RL_CARLA_BRAKE_SYSTEM_P.IW211_Value_a[i] *
      u0;

    /* DotProduct: '<S173>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(2,:)''
     */
    rtb_TmpSignalConversionAtDotP_1 += RL_CARLA_BRAKE_SYSTEM_P.IW212_Value_j[i] *
      u0;

    /* DotProduct: '<S184>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(3,:)''
     */
    tmp_1 += RL_CARLA_BRAKE_SYSTEM_P.IW213_Value_b[i] * u0;

    /* DotProduct: '<S186>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(4,:)''
     */
    tmp_2 += RL_CARLA_BRAKE_SYSTEM_P.IW214_Value_a[i] * u0;

    /* DotProduct: '<S187>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(5,:)''
     */
    tmp_3 += RL_CARLA_BRAKE_SYSTEM_P.IW215_Value_d[i] * u0;

    /* DotProduct: '<S188>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(6,:)''
     */
    tmp_4 += RL_CARLA_BRAKE_SYSTEM_P.IW216_Value_p[i] * u0;

    /* DotProduct: '<S189>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(7,:)''
     */
    tmp_5 += RL_CARLA_BRAKE_SYSTEM_P.IW217_Value_d[i] * u0;

    /* DotProduct: '<S190>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(8,:)''
     */
    tmp_6 += RL_CARLA_BRAKE_SYSTEM_P.IW218_Value_j[i] * u0;

    /* DotProduct: '<S191>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(9,:)''
     */
    tmp_7 += RL_CARLA_BRAKE_SYSTEM_P.IW219_Value_f[i] * u0;

    /* DotProduct: '<S163>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(10,:)''
     */
    tmp_8 += RL_CARLA_BRAKE_SYSTEM_P.IW2110_Value_m[i] * u0;

    /* DotProduct: '<S164>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(11,:)''
     */
    tmp_9 += RL_CARLA_BRAKE_SYSTEM_P.IW2111_Value_a[i] * u0;

    /* DotProduct: '<S165>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(12,:)''
     */
    tmp_a += RL_CARLA_BRAKE_SYSTEM_P.IW2112_Value_a[i] * u0;

    /* DotProduct: '<S166>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(13,:)''
     */
    tmp_b += RL_CARLA_BRAKE_SYSTEM_P.IW2113_Value_c[i] * u0;

    /* DotProduct: '<S167>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(14,:)''
     */
    tmp_c += RL_CARLA_BRAKE_SYSTEM_P.IW2114_Value_e[i] * u0;

    /* DotProduct: '<S168>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(15,:)''
     */
    tmp_d += RL_CARLA_BRAKE_SYSTEM_P.IW2115_Value_j[i] * u0;

    /* DotProduct: '<S169>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(16,:)''
     */
    tmp_e += RL_CARLA_BRAKE_SYSTEM_P.IW2116_Value_j[i] * u0;

    /* DotProduct: '<S170>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(17,:)''
     */
    tmp_f += RL_CARLA_BRAKE_SYSTEM_P.IW2117_Value_e[i] * u0;

    /* DotProduct: '<S171>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(18,:)''
     */
    tmp_g += RL_CARLA_BRAKE_SYSTEM_P.IW2118_Value_n[i] * u0;

    /* DotProduct: '<S172>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(19,:)''
     */
    tmp_h += RL_CARLA_BRAKE_SYSTEM_P.IW2119_Value_h[i] * u0;

    /* DotProduct: '<S174>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(20,:)''
     */
    tmp_i += RL_CARLA_BRAKE_SYSTEM_P.IW2120_Value_j[i] * u0;

    /* DotProduct: '<S175>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(21,:)''
     */
    tmp_j += RL_CARLA_BRAKE_SYSTEM_P.IW2121_Value_f[i] * u0;

    /* DotProduct: '<S176>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(22,:)''
     */
    tmp_k += RL_CARLA_BRAKE_SYSTEM_P.IW2122_Value_k[i] * u0;

    /* DotProduct: '<S177>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(23,:)''
     */
    tmp_l += RL_CARLA_BRAKE_SYSTEM_P.IW2123_Value_j[i] * u0;

    /* DotProduct: '<S178>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(24,:)''
     */
    tmp_m += RL_CARLA_BRAKE_SYSTEM_P.IW2124_Value_a[i] * u0;

    /* DotProduct: '<S179>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(25,:)''
     */
    tmp_n += RL_CARLA_BRAKE_SYSTEM_P.IW2125_Value_i[i] * u0;

    /* DotProduct: '<S180>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(26,:)''
     */
    tmp_o += RL_CARLA_BRAKE_SYSTEM_P.IW2126_Value_g[i] * u0;

    /* DotProduct: '<S181>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(27,:)''
     */
    tmp_p += RL_CARLA_BRAKE_SYSTEM_P.IW2127_Value_p[i] * u0;

    /* DotProduct: '<S182>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(28,:)''
     */
    tmp_q += RL_CARLA_BRAKE_SYSTEM_P.IW2128_Value_j[i] * u0;

    /* DotProduct: '<S183>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(29,:)''
     */
    tmp_r += RL_CARLA_BRAKE_SYSTEM_P.IW2129_Value_h[i] * u0;

    /* DotProduct: '<S185>/Dot Product' incorporates:
     *  Constant: '<S160>/IW{2,1}(30,:)''
     */
    tmp_s += RL_CARLA_BRAKE_SYSTEM_P.IW2130_Value_f[i] * u0;
  }

  /* Sum: '<S102>/netsum' incorporates:
   *  DotProduct: '<S162>/Dot Product'
   *  DotProduct: '<S163>/Dot Product'
   *  DotProduct: '<S164>/Dot Product'
   *  DotProduct: '<S165>/Dot Product'
   *  DotProduct: '<S166>/Dot Product'
   *  DotProduct: '<S167>/Dot Product'
   *  DotProduct: '<S168>/Dot Product'
   *  DotProduct: '<S169>/Dot Product'
   *  DotProduct: '<S170>/Dot Product'
   *  DotProduct: '<S171>/Dot Product'
   *  DotProduct: '<S172>/Dot Product'
   *  DotProduct: '<S173>/Dot Product'
   *  DotProduct: '<S174>/Dot Product'
   *  DotProduct: '<S175>/Dot Product'
   *  DotProduct: '<S176>/Dot Product'
   *  DotProduct: '<S177>/Dot Product'
   *  DotProduct: '<S178>/Dot Product'
   *  DotProduct: '<S179>/Dot Product'
   *  DotProduct: '<S180>/Dot Product'
   *  DotProduct: '<S181>/Dot Product'
   *  DotProduct: '<S182>/Dot Product'
   *  DotProduct: '<S183>/Dot Product'
   *  DotProduct: '<S184>/Dot Product'
   *  DotProduct: '<S185>/Dot Product'
   *  DotProduct: '<S186>/Dot Product'
   *  DotProduct: '<S187>/Dot Product'
   *  DotProduct: '<S188>/Dot Product'
   *  DotProduct: '<S189>/Dot Product'
   *  DotProduct: '<S190>/Dot Product'
   *  DotProduct: '<S191>/Dot Product'
   */
  tmp_0[0] = rtb_TmpSignalConversionAtDotP_0;
  tmp_0[1] = rtb_TmpSignalConversionAtDotP_1;
  tmp_0[2] = tmp_1;
  tmp_0[3] = tmp_2;
  tmp_0[4] = tmp_3;
  tmp_0[5] = tmp_4;
  tmp_0[6] = tmp_5;
  tmp_0[7] = tmp_6;
  tmp_0[8] = tmp_7;
  tmp_0[9] = tmp_8;
  tmp_0[10] = tmp_9;
  tmp_0[11] = tmp_a;
  tmp_0[12] = tmp_b;
  tmp_0[13] = tmp_c;
  tmp_0[14] = tmp_d;
  tmp_0[15] = tmp_e;
  tmp_0[16] = tmp_f;
  tmp_0[17] = tmp_g;
  tmp_0[18] = tmp_h;
  tmp_0[19] = tmp_i;
  tmp_0[20] = tmp_j;
  tmp_0[21] = tmp_k;
  tmp_0[22] = tmp_l;
  tmp_0[23] = tmp_m;
  tmp_0[24] = tmp_n;
  tmp_0[25] = tmp_o;
  tmp_0[26] = tmp_p;
  tmp_0[27] = tmp_q;
  tmp_0[28] = tmp_r;
  tmp_0[29] = tmp_s;

  /* DotProduct: '<S195>/Dot Product' incorporates:
   *  Constant: '<S193>/IW{3,2}(1,:)''
   */
  rtb_TmpSignalConversionAtDotP_0 = 0.0;
  for (i = 0; i < 30; i++) {
    /* Saturate: '<S161>/Saturation' incorporates:
     *  Constant: '<S102>/b{2}'
     *  Constant: '<S193>/IW{3,2}(1,:)''
     *  Sum: '<S102>/netsum'
     */
    u0 = tmp_0[i] + RL_CARLA_BRAKE_SYSTEM_P.b2_Value_h[i];
    if (u0 > RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_c) {
      u0 = RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_c;
    } else {
      if (u0 < RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat_bl) {
        u0 = RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat_bl;
      }
    }

    /* End of Saturate: '<S161>/Saturation' */
    rtb_TmpSignalConversionAtDotP_0 += RL_CARLA_BRAKE_SYSTEM_P.IW321_Value_m[i] *
      u0;
  }

  /* Gain: '<S4>/Gain2' incorporates:
   *  Constant: '<S103>/b{3}'
   *  DotProduct: '<S195>/Dot Product'
   *  Gain: '<S4>/Gain'
   *  Gain: '<S4>/Gain1'
   *  Sum: '<S103>/netsum'
   *  Sum: '<S4>/Add'
   */
  rtb_Gain2 = ((rtb_TmpSignalConversionAtDotP_0 +
                RL_CARLA_BRAKE_SYSTEM_P.b3_Value_b) *
               RL_CARLA_BRAKE_SYSTEM_P.Gain_Gain_h -
               RL_CARLA_BRAKE_SYSTEM_P.Gain1_Gain_m * rtb_NormalizedSpeed) *
    RL_CARLA_BRAKE_SYSTEM_P.Gain2_Gain_o;

  /* Update for DiscreteStateSpace: '<Root>/Discrete State-Space' */
  {
    real_T xnew[3];
    xnew[0] = (RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_A[0])*
      RL_CARLA_BRAKE_SYSTEM_DW.DiscreteStateSpace_DSTATE[0]
      + (RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_A[1])*
      RL_CARLA_BRAKE_SYSTEM_DW.DiscreteStateSpace_DSTATE[1];
    xnew[1] = (RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_A[2])*
      RL_CARLA_BRAKE_SYSTEM_DW.DiscreteStateSpace_DSTATE[1];
    xnew[1] += (RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_B[0])*rtb_Gain2;
    xnew[2] = (RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_B[1])*rtb_Gain2;
    (void) memcpy(&RL_CARLA_BRAKE_SYSTEM_DW.DiscreteStateSpace_DSTATE[0], xnew,
                  sizeof(real_T)*3);
  }

  /* Matfile logging */
  rt_UpdateTXYLogVars(RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo,
                      (&RL_CARLA_BRAKE_SYSTEM_M->Timing.taskTime0));

  /* External mode */
  rtExtModeUploadCheckTrigger(1);

  {                                    /* Sample time: [0.066666666666666666s, 0.0s] */
    rtExtModeUpload(0, (real_T)RL_CARLA_BRAKE_SYSTEM_M->Timing.taskTime0);
  }

  /* signal main to stop simulation */
  {                                    /* Sample time: [0.066666666666666666s, 0.0s] */
    if ((rtmGetTFinal(RL_CARLA_BRAKE_SYSTEM_M)!=-1) &&
        !((rtmGetTFinal(RL_CARLA_BRAKE_SYSTEM_M)-
           RL_CARLA_BRAKE_SYSTEM_M->Timing.taskTime0) >
          RL_CARLA_BRAKE_SYSTEM_M->Timing.taskTime0 * (DBL_EPSILON))) {
      rtmSetErrorStatus(RL_CARLA_BRAKE_SYSTEM_M, "Simulation finished");
    }

    if (rtmGetStopRequested(RL_CARLA_BRAKE_SYSTEM_M)) {
      rtmSetErrorStatus(RL_CARLA_BRAKE_SYSTEM_M, "Simulation finished");
    }
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++RL_CARLA_BRAKE_SYSTEM_M->Timing.clockTick0)) {
    ++RL_CARLA_BRAKE_SYSTEM_M->Timing.clockTickH0;
  }

  RL_CARLA_BRAKE_SYSTEM_M->Timing.taskTime0 =
    RL_CARLA_BRAKE_SYSTEM_M->Timing.clockTick0 *
    RL_CARLA_BRAKE_SYSTEM_M->Timing.stepSize0 +
    RL_CARLA_BRAKE_SYSTEM_M->Timing.clockTickH0 *
    RL_CARLA_BRAKE_SYSTEM_M->Timing.stepSize0 * 4294967296.0;
}

/* Model initialize function */
void RL_CARLA_BRAKE_SYSTEM_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat = rtInf;
  RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_e = rtInf;
  RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_f = rtInf;
  RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_c = rtInf;

  /* initialize real-time model */
  (void) memset((void *)RL_CARLA_BRAKE_SYSTEM_M, 0,
                sizeof(RT_MODEL_RL_CARLA_BRAKE_SYSTE_T));
  rtmSetTFinal(RL_CARLA_BRAKE_SYSTEM_M, 8.0);
  RL_CARLA_BRAKE_SYSTEM_M->Timing.stepSize0 = 0.066666666666666666;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = NULL;
    RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo, (NULL));
    rtliSetLogT(RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo, "tout");
    rtliSetLogX(RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo, "");
    rtliSetLogXFinal(RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo, 0);
    rtliSetLogDecimation(RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo, 1);
    rtliSetLogY(RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo, (NULL));
  }

  /* External mode info */
  RL_CARLA_BRAKE_SYSTEM_M->Sizes.checksums[0] = (112044904U);
  RL_CARLA_BRAKE_SYSTEM_M->Sizes.checksums[1] = (1786873236U);
  RL_CARLA_BRAKE_SYSTEM_M->Sizes.checksums[2] = (1306172202U);
  RL_CARLA_BRAKE_SYSTEM_M->Sizes.checksums[3] = (2742050055U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    RL_CARLA_BRAKE_SYSTEM_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(RL_CARLA_BRAKE_SYSTEM_M->extModeInfo,
      &RL_CARLA_BRAKE_SYSTEM_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(RL_CARLA_BRAKE_SYSTEM_M->extModeInfo,
                        RL_CARLA_BRAKE_SYSTEM_M->Sizes.checksums);
    rteiSetTPtr(RL_CARLA_BRAKE_SYSTEM_M->extModeInfo, rtmGetTPtr
                (RL_CARLA_BRAKE_SYSTEM_M));
  }

  /* block I/O */
  (void) memset(((void *) &RL_CARLA_BRAKE_SYSTEM_B), 0,
                sizeof(B_RL_CARLA_BRAKE_SYSTEM_T));

  /* states (dwork) */
  (void) memset((void *)&RL_CARLA_BRAKE_SYSTEM_DW, 0,
                sizeof(DW_RL_CARLA_BRAKE_SYSTEM_T));

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    RL_CARLA_BRAKE_SYSTEM_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 14;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(RL_CARLA_BRAKE_SYSTEM_M->rtwLogInfo, 0.0,
    rtmGetTFinal(RL_CARLA_BRAKE_SYSTEM_M),
    RL_CARLA_BRAKE_SYSTEM_M->Timing.stepSize0, (&rtmGetErrorStatus
    (RL_CARLA_BRAKE_SYSTEM_M)));

  /* InitializeConditions for DiscreteStateSpace: '<Root>/Discrete State-Space' */
  RL_CARLA_BRAKE_SYSTEM_DW.DiscreteStateSpace_DSTATE[0] =
    (RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_InitialCondi[0]);
  RL_CARLA_BRAKE_SYSTEM_DW.DiscreteStateSpace_DSTATE[1] =
    (RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_InitialCondi[1]);
  RL_CARLA_BRAKE_SYSTEM_DW.DiscreteStateSpace_DSTATE[2] =
    (RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_InitialCondi[2]);
}

/* Model terminate function */
void RL_CARLA_BRAKE_SYSTEM_terminate(void)
{
  /* (no terminate code required) */
}
