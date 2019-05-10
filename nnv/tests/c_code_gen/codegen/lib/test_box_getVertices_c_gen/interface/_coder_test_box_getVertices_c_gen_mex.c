/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_test_box_getVertices_c_gen_mex.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 23-Apr-2019 15:56:18
 */

/* Include Files */
#include "_coder_test_box_getVertices_c_gen_api.h"
#include "_coder_test_box_getVertices_c_gen_mex.h"

/* Function Declarations */
static void c_test_box_getVertices_c_gen_me(int32_T nlhs, mxArray *plhs[1],
  int32_T nrhs, const mxArray *prhs[2]);

/* Function Definitions */

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[1]
 *                int32_T nrhs
 *                const mxArray *prhs[2]
 * Return Type  : void
 */
static void c_test_box_getVertices_c_gen_me(int32_T nlhs, mxArray *plhs[1],
  int32_T nrhs, const mxArray *prhs[2])
{
  const mxArray *outputs[1];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 2, 4,
                        26, "test_box_getVertices_c_gen");
  }

  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 26,
                        "test_box_getVertices_c_gen");
  }

  /* Call the function. */
  test_box_getVertices_c_gen_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  emlrtReturnArrays(1, plhs, outputs);
}

/*
 * Arguments    : int32_T nlhs
 *                mxArray * const plhs[]
 *                int32_T nrhs
 *                const mxArray * const prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(test_box_getVertices_c_gen_atexit);

  /* Module initialization. */
  test_box_getVertices_c_gen_initialize();

  /* Dispatch the entry-point. */
  c_test_box_getVertices_c_gen_me(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  test_box_getVertices_c_gen_terminate();
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/*
 * File trailer for _coder_test_box_getVertices_c_gen_mex.c
 *
 * [EOF]
 */
