#include "mex.h"
#include <math.h>

static float get_scalar_as_single(const mxArray *a) {
  if (!mxIsNumeric(a) || mxIsComplex(a) || mxIsEmpty(a)) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_mex:Type", "Scalar input must be real numeric.");
  }
  if (mxIsSingle(a)) {
    return mxGetSingles(a)[0];
  }
  return (float)mxGetScalar(a);
}

static void copy_vector_to_single(const mxArray *a, float *out, mwSize n) {
  mwSize m = mxGetM(a);
  mwSize k = mxGetN(a);
  mwSize len = m * k;
  mwSize i;
  if (!mxIsNumeric(a) || mxIsComplex(a)) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_mex:Type", "Vector input must be real numeric.");
  }
  if (len != n) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_mex:Size", "Vector length mismatch.");
  }
  if (mxIsSingle(a)) {
    const float *p = mxGetSingles(a);
    for (i = 0; i < n; ++i) out[i] = p[i];
  } else {
    const double *p = mxGetDoubles(a);
    for (i = 0; i < n; ++i) out[i] = (float)p[i];
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  const mxArray *invL_in, *oldN_in, *c2c1h_in, *c2c2h_in;
  const float *invL = NULL;
  mwSize ld, nmax, oldN, i, j, k;
  float *c, *l21, *inv_new_old;
  float c2c2h, s, l22, invl22;

  if (nrhs != 4) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_mex:Args", "Expected 4 inputs: invL_CCH, oldN, C2C1H, C2C2H.");
  }
  if (nlhs != 3) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_mex:Args", "Expected 3 outputs: L21_row, L22, invL_new_old.");
  }

  invL_in = prhs[0];
  oldN_in = prhs[1];
  c2c1h_in = prhs[2];
  c2c2h_in = prhs[3];

  if (!mxIsSingle(invL_in) || mxIsComplex(invL_in) || mxGetNumberOfDimensions(invL_in) != 2) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_mex:Type", "invL_CCH must be a real single 2D matrix.");
  }
  if (mxGetM(invL_in) != mxGetN(invL_in)) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_mex:Size", "invL_CCH must be square.");
  }

  oldN = (mwSize)mxGetScalar(oldN_in);
  ld = mxGetM(invL_in);
  nmax = mxGetN(invL_in);
  if (oldN > ld || oldN > nmax) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_mex:Size", "oldN exceeds matrix size.");
  }

  plhs[0] = mxCreateNumericMatrix(1, oldN, mxSINGLE_CLASS, mxREAL); /* L21_row */
  plhs[1] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);    /* L22 */
  plhs[2] = mxCreateNumericMatrix(1, oldN, mxSINGLE_CLASS, mxREAL); /* invL_new_old */
  l21 = mxGetSingles(plhs[0]);
  inv_new_old = mxGetSingles(plhs[2]);

  invL = mxGetSingles(invL_in);
  c = (float *)mxCalloc((oldN > 0 ? oldN : 1), sizeof(float));
  copy_vector_to_single(c2c1h_in, c, oldN);
  c2c2h = get_scalar_as_single(c2c2h_in);

  if (oldN == 0) {
    s = c2c2h;
    if (s <= 0.0f) {
      mxFree(c);
      mexErrMsgIdAndTxt("isdf_schur_rank1_mex:SPD", "Non-positive Schur complement (oldN=0).");
    }
    mxGetSingles(plhs[1])[0] = sqrtf(s);
    mxFree(c);
    return;
  }

  for (k = 0; k < oldN; ++k) {
    double acc = 0.0;
    for (j = 0; j < oldN; ++j) {
      /* invL11' access: invL(k, j) */
      acc += (double)c[j] * (double)invL[k + j * ld];
    }
    l21[k] = (float)acc;
  }

  s = c2c2h;
  for (k = 0; k < oldN; ++k) {
    s -= l21[k] * l21[k];
  }
  if (s <= 0.0f) {
    mxFree(c);
    mexErrMsgIdAndTxt("isdf_schur_rank1_mex:SPD", "Non-positive Schur complement in rank-1 update.");
  }

  l22 = sqrtf(s);
  mxGetSingles(plhs[1])[0] = l22;
  invl22 = 1.0f / l22;

  for (j = 0; j < oldN; ++j) {
    double acc = 0.0;
    for (i = 0; i < oldN; ++i) {
      acc += (double)(invl22 * l21[i]) * (double)invL[i + j * ld];
    }
    inv_new_old[j] = (float)(-acc);
  }

  mxFree(c);
}
