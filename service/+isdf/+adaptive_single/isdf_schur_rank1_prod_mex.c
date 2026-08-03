// License-Identifier: BSD-3-Clause
//
// Copyright (C) 2026
//
// Authors (see AUTHORS file for details): ZZ
//
// Last modified: 2026/05/20 ZZ

#include "mex.h"
#include "blas.h"
#include <math.h>
#include <stddef.h>

static void require_double_matrix(const mxArray *a, const char *name) {
  if (!mxIsDouble(a) || mxGetNumberOfDimensions(a) != 2) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:Type", "%s must be a double 2D matrix.", name);
  }
}

static mxComplexDouble *to_complex_buffer(const mxArray *a, mwSize n, int *must_free) {
  mwSize i;
  if (mxIsComplex(a)) {
    *must_free = 0;
    return mxGetComplexDoubles(a);
  }
  *must_free = 1;
  {
    const double *r = mxGetDoubles(a);
    mxComplexDouble *buf = (mxComplexDouble *)mxCalloc(n, sizeof(mxComplexDouble));
    for (i = 0; i < n; ++i) {
      buf[i].real = r[i];
      buf[i].imag = 0.0;
    }
    return buf;
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  const mxArray *invL_in, *oldN_in, *psi_new_in, *phi_new_in, *psi_all_in, *phi_all_in;
  const double *invL = NULL;
  const mxComplexDouble *psi_new = NULL, *phi_new = NULL, *psi_all = NULL, *phi_all = NULL;
  mwSize ld, nmax, oldN, nb1, nb2, i, j;
  ptrdiff_t m_, n_, k_, lda, ldb, ldc;
  char trans_n = 'N', trans_c = 'C';
  double alpha_c[2] = {1.0, 0.0};
  double beta_c[2] = {0.0, 0.0};
  mxComplexDouble *m1 = NULL, *m2 = NULL;
  double *l21, *inv_new_old;
  double c2c2h, s, l22, invl22;
  double *c = NULL;
  int free_psi_new = 0, free_phi_new = 0, free_psi_all = 0, free_phi_all = 0;

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

  if (!mxIsDouble(invL_in) || mxIsComplex(invL_in) || mxGetNumberOfDimensions(invL_in) != 2) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:Type", "invL_CCH must be a real double 2D matrix.");
  }
  if (mxGetM(invL_in) != mxGetN(invL_in)) {
    mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:Size", "invL_CCH must be square.");
  }
  require_double_matrix(psi_new_in, "Psi_new");
  require_double_matrix(phi_new_in, "Phi_new");
  require_double_matrix(psi_all_in, "Psi_all");
  require_double_matrix(phi_all_in, "Phi_all");

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
  plhs[0] = mxCreateNumericMatrix(1, oldN, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
  plhs[2] = mxCreateNumericMatrix(1, oldN, mxDOUBLE_CLASS, mxREAL);
  l21 = mxGetDoubles(plhs[0]);
  inv_new_old = mxGetDoubles(plhs[2]);

  invL = mxGetDoubles(invL_in);
  psi_new = to_complex_buffer(psi_new_in, mxGetNumberOfElements(psi_new_in), &free_psi_new);
  phi_new = to_complex_buffer(phi_new_in, mxGetNumberOfElements(phi_new_in), &free_phi_new);
  psi_all = to_complex_buffer(psi_all_in, mxGetNumberOfElements(psi_all_in), &free_psi_all);
  phi_all = to_complex_buffer(phi_all_in, mxGetNumberOfElements(phi_all_in), &free_phi_all);

  {
    double npsi = 0.0;
    double nphi = 0.0;
    for (j = 0; j < nb1; ++j) {
      double xr = psi_new[j].real, xi = psi_new[j].imag;
      npsi += (double)xr * (double)xr + (double)xi * (double)xi;
    }
    for (j = 0; j < nb2; ++j) {
      double xr = phi_new[j].real, xi = phi_new[j].imag;
      nphi += (double)xr * (double)xr + (double)xi * (double)xi;
    }
    c2c2h = (double)(npsi * nphi);
  }

  if (oldN == 0) {
    s = c2c2h;
    if (s <= 0.0) {
      if (free_psi_new) mxFree((void *)psi_new);
      if (free_phi_new) mxFree((void *)phi_new);
      if (free_psi_all) mxFree((void *)psi_all);
      if (free_phi_all) mxFree((void *)phi_all);
      mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:SPD", "Non-positive Schur complement (oldN=0).");
    }
    mxGetDoubles(plhs[1])[0] = sqrt(s);
    if (free_psi_new) mxFree((void *)psi_new);
    if (free_phi_new) mxFree((void *)phi_new);
    if (free_psi_all) mxFree((void *)psi_all);
    if (free_phi_all) mxFree((void *)phi_all);
    return;
  }

  m1 = (mxComplexDouble *)mxCalloc((oldN > 0 ? oldN : 1), sizeof(mxComplexDouble));
  m2 = (mxComplexDouble *)mxCalloc((oldN > 0 ? oldN : 1), sizeof(mxComplexDouble));
  c = (double *)mxCalloc((oldN > 0 ? oldN : 1), sizeof(double));
  m_ = 1;
  n_ = (ptrdiff_t)oldN;
  k_ = (ptrdiff_t)nb1;
  lda = 1;
  ldb = (ptrdiff_t)mxGetM(psi_all_in);
  ldc = 1;
  zgemm(&trans_n, &trans_c, &m_, &n_, &k_, alpha_c, (double *)psi_new, &lda, (double *)psi_all, &ldb, beta_c, (double *)m1, &ldc);

  k_ = (ptrdiff_t)nb2;
  ldb = (ptrdiff_t)mxGetM(phi_all_in);
  zgemm(&trans_n, &trans_c, &m_, &n_, &k_, alpha_c, (double *)phi_new, &lda, (double *)phi_all, &ldb, beta_c, (double *)m2, &ldc);

  for (i = 0; i < oldN; ++i) {
    double ar = m1[i].real;
    double ai = -m1[i].imag;
    double br = m2[i].real;
    double bi = m2[i].imag;
    double ci = ar * bi + ai * br;
    double cr = ar * br - ai * bi;
    double tol = 1e-5 * (fabs(cr) + 1.0);
    c[i] = cr;
    if (fabs(ci) > tol) {
      mxFree(c);
      mxFree(m1);
      mxFree(m2);
      if (free_psi_new) mxFree((void *)psi_new);
      if (free_phi_new) mxFree((void *)phi_new);
      if (free_psi_all) mxFree((void *)psi_all);
      if (free_phi_all) mxFree((void *)phi_all);
      mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:Complex",
                        "C2C1H has significant imaginary part; fused rank-1 kernel expects near-real values.");
    }
  }

  for (i = 0; i < oldN; ++i) {
    double acc = 0.0;
    for (j = 0; j < oldN; ++j) {
      acc += (double)c[j] * (double)invL[i + j * ld];
    }
    l21[i] = (double)acc;
  }

  s = c2c2h;
  for (i = 0; i < oldN; ++i) {
    s -= l21[i] * l21[i];
  }
  if (s <= 0.0) {
    mxFree(c);
    mxFree(m1);
    mxFree(m2);
    if (free_psi_new) mxFree((void *)psi_new);
    if (free_phi_new) mxFree((void *)phi_new);
    if (free_psi_all) mxFree((void *)psi_all);
    if (free_phi_all) mxFree((void *)phi_all);
    mexErrMsgIdAndTxt("isdf_schur_rank1_prod_mex:SPD", "Non-positive Schur complement in rank-1 update.");
  }
  l22 = sqrt(s);
  mxGetDoubles(plhs[1])[0] = l22;
  invl22 = 1.0 / l22;

  for (j = 0; j < oldN; ++j) {
    double acc = 0.0;
    for (i = 0; i < oldN; ++i) {
      acc += (double)(invl22 * l21[i]) * (double)invL[i + j * ld];
    }
    inv_new_old[j] = (double)(-acc);
  }

  mxFree(c);
  mxFree(m1);
  mxFree(m2);
  if (free_psi_new) mxFree((void *)psi_new);
  if (free_phi_new) mxFree((void *)phi_new);
  if (free_psi_all) mxFree((void *)psi_all);
  if (free_phi_all) mxFree((void *)phi_all);
}
