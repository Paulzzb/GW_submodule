#include "mex.h"
#include "blas.h"
#include <math.h>
#include <stddef.h>

static void require_single_complex_matrix(const mxArray *a, const char *name) {
  if (!mxIsSingle(a) || !mxIsComplex(a) || mxGetNumberOfDimensions(a) != 2) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:Type", "%s must be a single complex 2D matrix.", name);
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  const mxArray *invL_in, *oldN_in, *psi_new_in, *phi_new_in, *psi_all_in, *phi_all_in;
  const float *invL = NULL;
  const mxComplexSingle *psi_new = NULL, *phi_new = NULL, *psi_all = NULL, *phi_all = NULL;
  mwSize ld, nmax, oldN, nb1, nb2, i, j;
  ptrdiff_t m_, n_, k_, lda, ldb, ldc;
  char trans_n = 'N', trans_c = 'C';
  float alpha_c[2] = {1.0f, 0.0f};
  float beta_c[2] = {0.0f, 0.0f};
  mxComplexSingle *m1 = NULL, *m2 = NULL;
  float *l21, *inv_new_old;
  float c2c2h, s, l22, invl22;
  float *c = NULL;

  if (nrhs != 6) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:Args",
                      "Expected 6 inputs: invL_CCH, oldN, Psi_new, Phi_new, Psi_all, Phi_all.");
  }
  if (nlhs != 3) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:Args",
                      "Expected 3 outputs: L21_row, L22, invL_new_old.");
  }

  invL_in = prhs[0];
  oldN_in = prhs[1];
  psi_new_in = prhs[2];
  phi_new_in = prhs[3];
  psi_all_in = prhs[4];
  phi_all_in = prhs[5];

  if (!mxIsSingle(invL_in) || mxIsComplex(invL_in) || mxGetNumberOfDimensions(invL_in) != 2) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:Type", "invL_CCH must be a real single 2D matrix.");
  }
  if (mxGetM(invL_in) != mxGetN(invL_in)) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:Size", "invL_CCH must be square.");
  }
  require_single_complex_matrix(psi_new_in, "Psi_new");
  require_single_complex_matrix(phi_new_in, "Phi_new");
  require_single_complex_matrix(psi_all_in, "Psi_all");
  require_single_complex_matrix(phi_all_in, "Phi_all");

  if (mxGetM(psi_new_in) != 1 || mxGetM(phi_new_in) != 1) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:Size", "Psi_new and Phi_new must be 1-by-k row matrices (Nadd=1).");
  }

  oldN = (mwSize)mxGetScalar(oldN_in);
  ld = mxGetM(invL_in);
  nmax = mxGetN(invL_in);
  if (oldN > ld || oldN > nmax) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:Size", "oldN exceeds invL_CCH size.");
  }

  nb1 = mxGetN(psi_new_in);
  nb2 = mxGetN(phi_new_in);
  if (mxGetN(psi_all_in) != nb1 || mxGetN(phi_all_in) != nb2) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:Size", "Column counts of Psi/Phi new and all must match.");
  }
  if (mxGetM(psi_all_in) < oldN || mxGetM(phi_all_in) < oldN) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:Size", "Psi_all/Phi_all row count must be >= oldN.");
  }

  plhs[0] = mxCreateNumericMatrix(1, oldN, mxSINGLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
  plhs[2] = mxCreateNumericMatrix(1, oldN, mxSINGLE_CLASS, mxREAL);
  l21 = mxGetSingles(plhs[0]);
  inv_new_old = mxGetSingles(plhs[2]);

  invL = mxGetSingles(invL_in);
  psi_new = mxGetComplexSingles(psi_new_in);
  phi_new = mxGetComplexSingles(phi_new_in);
  psi_all = mxGetComplexSingles(psi_all_in);
  phi_all = mxGetComplexSingles(phi_all_in);

  {
    double npsi = 0.0;
    double nphi = 0.0;
    for (j = 0; j < nb1; ++j) {
      float xr = psi_new[j].real, xi = psi_new[j].imag;
      npsi += (double)xr * (double)xr + (double)xi * (double)xi;
    }
    for (j = 0; j < nb2; ++j) {
      float xr = phi_new[j].real, xi = phi_new[j].imag;
      nphi += (double)xr * (double)xr + (double)xi * (double)xi;
    }
    c2c2h = (float)(npsi * nphi);
  }

  if (oldN == 0) {
    s = c2c2h;
    if (s <= 0.0f) {
      mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:SPD", "Non-positive Schur complement (oldN=0).");
    }
    mxGetSingles(plhs[1])[0] = sqrtf(s);
    return;
  }

  m_ = 1;
  n_ = (ptrdiff_t)oldN;
  k_ = (ptrdiff_t)nb1;
  lda = 1;
  ldb = (ptrdiff_t)mxGetM(psi_all_in);
  ldc = 1;
  m1 = (mxComplexSingle *)mxCalloc((oldN > 0 ? oldN : 1), sizeof(mxComplexSingle));
  cgemm(&trans_n, &trans_c, &m_, &n_, &k_, alpha_c, (float *)psi_new, &lda, (float *)psi_all, &ldb, beta_c, (float *)m1, &ldc);

  k_ = (ptrdiff_t)nb2;
  ldb = (ptrdiff_t)mxGetM(phi_all_in);
  m2 = (mxComplexSingle *)mxCalloc((oldN > 0 ? oldN : 1), sizeof(mxComplexSingle));
  cgemm(&trans_n, &trans_c, &m_, &n_, &k_, alpha_c, (float *)phi_new, &lda, (float *)phi_all, &ldb, beta_c, (float *)m2, &ldc);

  c = (float *)mxCalloc((oldN > 0 ? oldN : 1), sizeof(float));
  for (i = 0; i < oldN; ++i) {
    float ar = m1[i].real;
    float ai = -m1[i].imag;
    float br = m2[i].real;
    float bi = m2[i].imag;
    {
      float cr = ar * br - ai * bi;
      float ci = ar * bi + ai * br;
      float tol = 1e-5f * (fabsf(cr) + 1.0f);
      if (fabsf(ci) > tol) {
        mxFree(c);
        mxFree(m1);
        mxFree(m2);
        mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:Complex",
                          "C2C1H has significant imaginary part; fused rank-1 kernel expects near-real values.");
      }
      c[i] = cr;
    }
  }

  for (i = 0; i < oldN; ++i) {
    double acc = 0.0;
    for (j = 0; j < oldN; ++j) {
      acc += (double)c[j] * (double)invL[i + j * ld];
    }
    l21[i] = (float)acc;
  }

  s = c2c2h;
  for (i = 0; i < oldN; ++i) {
    s -= l21[i] * l21[i];
  }
  if (s <= 0.0f) {
    mxFree(c);
    mxFree(m1);
    mxFree(m2);
    mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:SPD", "Non-positive Schur complement in rank-1 update.");
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
  mxFree(m1);
  mxFree(m2);
}
